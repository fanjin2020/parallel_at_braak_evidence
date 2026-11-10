"""Rebuild manuscript analyses from local authorized ADNI tables and public sources."""
import argparse
from collections import deque
import csv
from datetime import datetime
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time
import yaml
from source_io import sha256, local_cel_metadata
from linux_runtime import require_analysis_python, resolve_rscript
from offline_sources import load_cached_public

HERE = Path(__file__).resolve().parent
PACKAGE = HERE.parent
ADNI_FILES = ['ADSL.csv','APOERES.csv','CDR.csv','DATA_DOWNLOADED_DATE.csv','DXSUM.csv',
              'PTDEMOG.csv','UCBERKELEY_AMY_6MM.csv','UCBERKELEY_TAU_6MM.csv',
              'UCBERKELEY_TAUPVC_6MM.csv','UCBERKELEYFDG_8mm.csv','UCSFFSX7.csv']
CONFOUNDER_FILES = ['GDSCALE.csv','BACKMEDS.csv','MEDHIST.csv','MODHACH.csv','RECCMEDS.csv','INITHEALTH.csv']

def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--adni-raw-dir',type=Path,required=True)
    g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--adnimerge2-dir',type=Path,help='Already downloaded, unpacked ADNIMERGE2 package')
    g.add_argument('--confounder-raw-dir',type=Path,help='Existing raw CSV export, not dated_confounders.csv')
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--cache-dir',type=Path,default=PACKAGE/'cache/public')
    p.add_argument('--rscript',default='Rscript')
    p.add_argument('--go-gaf',type=Path)
    p.add_argument('--go-ontology',type=Path)
    p.add_argument('--gene-info',type=Path)
    p.add_argument('--bootstrap',type=int,default=1000)
    p.add_argument('--enhancement-bootstrap',type=int,default=500,
                   help='PET-day rematching and threshold-propagation bootstrap, matching the manuscript')
    p.add_argument('--cv-bootstrap',type=int,default=200)
    p.add_argument('--cv-repeats',type=int,default=20)
    p.add_argument('--imputations',type=int,default=100)
    p.add_argument('--threads',type=int,default=1)
    p.add_argument('--check-only',action='store_true')
    p.add_argument('--resume',action='store_true')
    p.add_argument('--keep-work',action='store_true',help='Retain temporary fit tables for debugging')
    return p

def preflight(a):
    missing=[str(a.adni_raw_dir/n) for n in ADNI_FILES if not (a.adni_raw_dir/n).is_file()]
    if a.adnimerge2_dir:
        missing += [str(a.adnimerge2_dir/'data'/Path(n).with_suffix('.rda')) for n in CONFOUNDER_FILES
                    if not (a.adnimerge2_dir/'data'/Path(n).with_suffix('.rda')).is_file()]
        if not (a.adnimerge2_dir/'DESCRIPTION').is_file():missing.append('ADNIMERGE2 DESCRIPTION')
    else:
        missing += [str(a.confounder_raw_dir/n) for n in CONFOUNDER_FILES if not (a.confounder_raw_dir/n).is_file()]
    if missing:raise FileNotFoundError('Missing local inputs:\n'+'\n'.join(missing))
    if any([a.go_gaf,a.go_ontology,a.gene_info]):
        raise ValueError('Offline mode uses the snapshots recorded in --cache-dir, not unrecorded overrides.')
    if a.enhancement_bootstrap<200:raise ValueError('enhancement-bootstrap must be at least 200')
    if a.bootstrap<200 or a.cv_bootstrap<100 or a.cv_repeats<2:raise ValueError('Insufficient resampling settings')
    if a.imputations<8 or a.imputations%4:raise ValueError('MI requires at least 8 draws, divisible by four')
    if a.threads<1:raise ValueError('threads must be positive')
    for name in ['numpy','pandas','scipy','statsmodels','patsy','yaml','matplotlib','openpyxl']:
        __import__(name)
    for p in [a.go_gaf,a.go_ontology,a.gene_info]:
        if p is not None and not p.is_file():raise FileNotFoundError(p)
    out=a.output_dir.resolve()
    for src in [a.adni_raw_dir,a.adnimerge2_dir,a.confounder_raw_dir]:
        if src and (out==src.resolve() or src.resolve() in out.parents):
            raise ValueError('Output must be outside the input directory')
    if out.exists() and any(out.iterdir()) and not a.resume:raise ValueError('Use a new output directory or --resume')

class Workflow:
    def __init__(self,a):
        from pipeline_runtime import Runner
        self.a=a;self.top=a.output_dir.resolve();self.root=self.top/'.work'
        self.root.mkdir(parents=True,exist_ok=True)
        self.env=os.environ.copy()
        self.env.update({k:str(a.threads) for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS']})
        self.env.update(PYTHONIOENCODING='utf-8',JAD_OFFLINE='1',BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS='FALSE')
        self.runner=Runner(self.top/'audit',self.root,self.env)
    def run(self,name,cmd):
        # A failed non-checkpointed primary job is restarted in its own output folder.
        state=self.top/'audit'/(name+'.json')
        if state.is_file() and json.loads(state.read_text()).get('status')!='complete':
            args=list(map(str,cmd))
            if '--output-dir' in args:
                target=Path(args[args.index('--output-dir')+1]).resolve()
                protected=self.root.resolve()
                if target==protected or protected not in target.parents or target.is_symlink():
                    raise ValueError('Unsafe partial-job directory')
                if target.exists():shutil.rmtree(target)
        self.runner.run(name,cmd,self.root,changed_only=True)
    def py(self,name,script,*args):self.run(name,[sys.executable,HERE/script,*args])
    def r(self,name,script,*args):self.run(name,[sys.executable,HERE/'run_r.py','--rscript',self.a.rscript,'--script',script,'--',*args])

def local_stage(source,dest,files):
    dest.mkdir(parents=True,exist_ok=True)
    for name in files:
        src=source/name;target=dest/name
        if target.exists():
            if sha256(src)!=sha256(target):raise ValueError('Staged input changed: '+name)
            continue
        # Hard links avoid copying large local tables; no input file is modified.
        try:os.link(src,target)
        except OSError:shutil.copy2(src,target)

def configure(root):
    c=yaml.safe_load((PACKAGE/'config/project.example.yaml').read_text(encoding='utf-8'))
    c['paths']['gse_subject_info']='data/public/gse131617/raw/published_supplement.xls'
    c['paths']['gse_subject_sheet']='Supple_Table_S8'
    (root/'config').mkdir(exist_ok=True)
    file=root/'config/project.yaml';file.write_text(yaml.safe_dump(c,sort_keys=False),encoding='utf-8')
    return file

def execute(a):
    require_analysis_python()
    preflight(a)
    a.rscript=resolve_rscript(a.rscript)
    sources=load_cached_public(a.cache_dir)
    if a.check_only:
        print('Local ADNI and public-source cache verified. No download or analysis has run.')
        subprocess.run([sys.executable,str(HERE/'run_r.py'),'--rscript',a.rscript,'--script','check_primary_environment.R'],check=True)
        print('Offline preflight passed. Full fits remain untested.')
        return
    w=Workflow(a);root=w.root
    settings={k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()}
    (root/'run_settings.json').write_text(json.dumps(settings,indent=2),encoding='utf-8')
    (root/'code_identity.json').write_text(json.dumps({p.name:sha256(p) for p in HERE.iterdir()
        if p.suffix in ['.py','.R']},indent=2),encoding='utf-8')
    cache=a.cache_dir.resolve();cache.mkdir(parents=True,exist_ok=True)
    w.r('dependencies','check_primary_environment.R')
    raw=root/'data/private/adni/raw'
    local_stage(a.adni_raw_dir.resolve(),raw,ADNI_FILES)
    confraw=root/'data/private/adni/confounder_raw'
    if a.adnimerge2_dir:
        w.r('export_local_ADNI','export_local_adni_tables.R','--package-dir',a.adnimerge2_dir.resolve(),'--output-dir',confraw)
    else:local_stage(a.confounder_raw_dir.resolve(),confraw,CONFOUNDER_FILES)
    private_inventory=[]
    for folder in [raw,confraw]:
        for f in sorted(folder.glob('*.csv')):
            private_inventory.append(dict(relative_path=f.relative_to(root).as_posix(),bytes=f.stat().st_size,sha256=sha256(f)))
    (root/'local_input_inventory.json').write_text(json.dumps(private_inventory,indent=2),encoding='utf-8')
    config=configure(root)
    (root/'public_sources_manifest.json').write_text(json.dumps(sources,indent=2),encoding='utf-8')
    upstream={key:Path(value) for key,value in sources['paths'].items()}
    public=root/'data/public/gse131617'
    (public/'raw').mkdir(parents=True,exist_ok=True)
    matrix=upstream['discovery_matrix'];donors=upstream['published_donor_supplement']
    shutil.copy2(matrix,public/'raw'/matrix.name)
    shutil.copy2(donors,public/'raw'/'published_supplement.xls')
    w.py('prepare_ADNI','prepare_adni.py','--config',config)
    w.py('original_ADNI','analyze_adni.py','--config',config)
    w.r('prepare_discovery','prepare_gse131617.R','--config',config)
    w.r('discovery_full_background','analyze_gse131617.R','--config',config)
    w.r('derive_targets','derive_discovery_targets.R','--project-root',root)
    w.py('original_clinical_extensions','clinical_sensitivity_models.py','--project-root',root,'--output-dir',root/'results/extensions')
    rev=root/'revision';enh=root/'enhancement';target=root/'targeted'
    rev.mkdir(exist_ok=True);enh.mkdir(exist_ok=True);target.mkdir(exist_ok=True)
    cel=upstream['cel_dir']
    import pandas as pd
    expected_samples=set(pd.read_csv(public/'sample_manifest.csv').sample_id.astype(str))
    downloaded_samples={r['sample_id'] for r in sources['cel_files']}
    if expected_samples!=downloaded_samples:raise ValueError('Downloaded CEL samples and rebuilt analysis manifest differ')
    local_cel_metadata(cel,root/'data/public/cel_metadata',sources['cel_files'])
    meta=root/'data/public/cel_metadata/sample_technical_metadata.csv'
    conf=root/'data/private/adni/confounders'
    w.py('prepare_dated_confounders','prepare_confounders.py','--raw-dir',confraw,'--output-dir',conf)
    common=['--project-root',root]
    w.py('incremental','analyze_tau_incremental.py',*common,'--output-dir',rev/'incremental','--bootstrap',a.bootstrap,'--repeats',a.cv_repeats)
    w.py('clinical','pet_staging_covariates.py',*common,'--output-dir',rev/'clinical','--data-derived-staging','--stage-bootstrap',a.bootstrap)
    w.py('longitudinal','longitudinal_model_comparison.py',*common,'--output-dir',rev/'longitudinal','--bootstrap',a.bootstrap,'--cv-bootstrap',a.cv_bootstrap)
    w.py('visit_process','visit_process_models.py',*common,'--output-dir',rev/'visit_process','--bootstrap',a.bootstrap)
    w.r('transcriptomics','tissue_models.R',*common,'--output-dir',rev/'transcriptomics','--sample-metadata',meta)
    ref=root/'data/public/cell_reference'
    w.py('cell_reference','prepare_cell_reference.py','--output-dir',ref,
        '--workbook',upstream['cell_reference_workbook'],'--gene-info',upstream['gene_info'])
    w.r('cell_composition','cell_composition_models.R',*common,'--reference-dir',ref,'--output-dir',rev/'cell_composition')
    w.py('additional_confounders','clinical_confounder_models.py',*common,'--output-dir',rev/'additional_confounders','--confounder-file',conf/'dated_confounders.csv')
    go=rev/'go_annotation'
    goargs=[]
    for key,value in [('--annotation-file',upstream['go_gaf']),('--ontology-file',upstream['go_ontology']),('--gene-info-file',upstream['gene_info'])]:
        goargs.extend([key,value])
    w.py('GO_annotation','prepare_go.py','--gene-ranks',rev/'transcriptomics/gene_level_ranks.csv','--output-dir',go,*goargs)
    w.r('pathways','pathway_models.R',*common,'--transcriptomic-results',rev/'transcriptomics','--go-dir',go,'--cell-results',rev/'cell_composition','--output-dir',rev/'pathways')
    for name in ['pathway_score_region_effects.csv','pathway_score_region_interactions.csv']:
        if not (rev/'pathways'/name).is_file():raise RuntimeError('Required pathway analysis did not run: '+name)
    w.r('pathway_LORO','pathway_leave_region_out.R',*common,'--transcriptomic-results',rev/'transcriptomics','--go-dir',go,'--pathway-results',rev/'pathways','--output-dir',rev/'pathway_loro')
    ext=root/'data/public/GSE44770'
    w.py('prepare_GSE44770','prepare_external_gse44770.py','--output-dir',ext,
        '--matrix-file',upstream['external_AD_matrix'],'--platform-file',upstream['external_AD_platform'])
    targets=rev/'transcriptomics/external_validation_targets.csv'
    w.r('external_GSE44770','external_cohort_models.R','--output-dir',rev/'external','--targets',targets,
        '--manifest',ext/'manifest.csv','--expression',ext/'expression_log2.csv','--covariates','age_at_death,sex_male,pmi_hours,tissue_ph,rin,batch',
        '--validation-kind','related_AD_vs_control','--source','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44770',
        '--independence-evidence','Independent HBTRC prefrontal AD/control cohort; not a matched Braak replication',
        '--go-dir',go,'--discovery-pathways',rev/'pathways/current_GO_all_tests.csv')
    w.r('enhanced_tissue','tissue_technical_models.R',*common,'--output-dir',enh/'tissue','--go-dir',go,'--sample-metadata',meta,'--cell-results',rev/'cell_composition')
    w.py('rematched_clinical','pet_day_clinical_models.py',*common,'--output-dir',enh/'clinical','--bootstrap',a.enhancement_bootstrap,'--repeats',a.cv_repeats)
    w.py('missingness','clinical_missingness.py',*common,'--output-dir',enh/'missingness','--imputations',20,'--bootstrap',a.bootstrap)
    w.r('CEL_QC','cel_qc_models.R',*common,'--cel-dir',cel,'--output-dir',root/'cel_qc','--go-dir',go,'--allow-qc-exclusion','--max-excluded-fraction','0.05')
    import pandas as pd
    qc=root/'cel_qc'
    # Dimensions are computed from this run, not enforced by selecting a favorable subset.
    n=pd.read_csv(qc/'retained_sample_manifest.csv').shape[0]
    nf=pd.read_csv(public/'feature_annotation.csv').shape[0]
    ng=pd.read_csv(qc/'fixed_gene_representatives.csv').shape[0]
    cov=pd.read_csv(qc/'GO_coverage.csv');ngo=int(cov.tested.astype(str).str.lower().isin(['true','1']).sum())
    w.r('QC_joint_models','qc_joint_models.R',*common,'--qc-results',qc,'--sample-metadata',meta,'--reference-dir',ref,'--go-dir',go,
        '--output-dir',target/'qc_joint','--expected-samples',n,'--expected-features',nf,'--expected-genes',ng,'--expected-go',ngo)
    w.py('matched_MI','mi_matched_estimands.py',*common,'--output-dir',target/'mi_matched','--imputations',a.imputations,'--audit-current-run')
    braak=root/'data/public/GSE106241'
    w.py('prepare_GSE106241','prepare_braak_external.py','--cache-dir',cache/'GSE106241','--output-dir',braak,
        '--soft-file',upstream['external_Braak_soft'],'--expression-file',upstream['external_Braak_matrix'],'--gene-info',upstream['gene_info'])
    w.r('external_GSE106241','validate_braak_external.R','--input-dir',braak,'--targets',targets,'--go-dir',go,'--output-dir',enh/'external')
    w.r('published_RIN','extract_published_rin.R','--source-xls',donors,'--output-dir',root/'published_rin')
    prov_args=['--package-dir',a.adnimerge2_dir.resolve()] if a.adnimerge2_dir else []
    w.py('reporting_tables','prepare_reporting_tables.py','--run-dir',root,*prov_args)
    env=dict(python=platform.python_version(),platform=platform.platform(),threads=a.threads,
             packages={p:importlib.metadata.version(p) for p in ['numpy','pandas','scipy','statsmodels','patsy','PyYAML','matplotlib','openpyxl']},
             settings={k:str(v) for k,v in vars(a).items()})
    (root/'run_environment.json').write_text(json.dumps(env,indent=2),encoding='utf-8')
    print('Primary models complete:',root,flush=True)

if __name__=='__main__':execute(parser().parse_args())
