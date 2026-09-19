"""Re-run the original and additional paper analyses into a new local directory."""
from pathlib import Path
import argparse
import copy
import importlib
import shutil
import subprocess
import sys

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--input-root',type=Path,help='Project containing the input data; defaults to this project')
    p.add_argument('--config',type=Path,default=Path('config/project.yaml'))
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--rscript',default=None,help='Rscript executable, e.g. conda run environment executable path')
    p.add_argument('--bootstrap',type=int,default=2000)
    p.add_argument('--selection-bootstrap',type=int,default=1000)
    p.add_argument('--seed',type=int,default=20260915)
    p.add_argument('--dry-run',action='store_true',help='Check inputs and print commands without writing results')
    a=p.parse_args()
    root=a.project_root.resolve()
    input_root=a.input_root.resolve() if a.input_root else root
    cfg=a.config if a.config.is_absolute() else root/a.config
    if not cfg.is_file():
        cfg=root/'config/project.example.yaml'
        print('Local config absent; using config/project.example.yaml',flush=True)
    try:
        import yaml
    except ImportError:
        p.error('PyYAML is missing in this Python environment. Activate the analysis environment or install PyYAML>=6,<7.')
    conf=yaml.safe_load(cfg.read_text(encoding='utf-8-sig'))
    paths=conf['paths']
    def source(key):
        v=Path(paths[key]); return v if v.is_absolute() else (root if key in ['gse_gene_symbols','go_bp_term_dictionary'] else input_root)/v
    if not source('gse_subject_info').is_file():
        candidates=list((input_root/'data/public/gse131617/raw').glob('GSE131617*Subject*info*.xlsx'))
        if len(candidates)==1:
            paths['gse_subject_info']=str(candidates[0])
            print('Using donor metadata: '+str(candidates[0]),flush=True)
    out=(a.output_dir if a.output_dir.is_absolute() else root/a.output_dir).resolve()
    if out==root or root/'data'==out or root/'data' in out.parents or out==input_root or input_root/'data'==out or input_root/'data' in out.parents:
        p.error('Choose a new output directory outside data/.')
    if out.exists() and any(out.iterdir()): p.error('Output directory must be empty; use a new run name.')
    if a.bootstrap<100 or a.selection_bootstrap<200: p.error('Bootstrap counts must be at least 100 and 200.')
    raw=source('adni_raw_dir')
    raw_names=['UCBERKELEY_AMY_6MM.csv','UCBERKELEY_TAU_6MM.csv','UCSFFSX7.csv','CDR.csv','APOERES.csv','PTDEMOG.csv','UCBERKELEYFDG_8mm.csv','DXSUM.csv','ADSL.csv','DATA_DOWNLOADED_DATE.csv']
    required=[raw/n for n in raw_names]
    required += [source(k) for k in ['gse_series_matrix','gse_subject_info','gse_gene_symbols','go_bp_term_dictionary']]
    missing=[str(f) for f in required if not f.is_file()]
    if missing: p.error('Missing inputs:\n'+'\n'.join(missing))
    for name in ['numpy','pandas','scipy','statsmodels','patsy']: importlib.import_module(name)
    r=a.rscript or conf.get('runtime',{}).get('rscript','Rscript')
    if not shutil.which(r): p.error('Rscript not found. Activate the R environment or pass --rscript /absolute/path/Rscript')
    runconf=copy.deepcopy(conf)

    runconf['paths']={k:str(source(k)) for k in paths}
    outputs={
      'adni_master':'data/private/adni/adni_pet_aligned_master.csv',
      'gse_feature_expression':'data/public/gse131617/expression_feature.csv.gz',
      'gse_feature_annotation':'data/public/gse131617/feature_annotation.csv',
      'gse_manifest':'data/public/gse131617/sample_manifest.csv',
      'gse_pathway_expression':'data/public/gse131617/pathway_expression.csv.gz',
      'output_dir':'results', 'adni_model_counts':'results/adni/adni_model_counts.csv',
      'adni_results':'results/adni/adni_results.csv','gse_results_dir':'results/gse131617',
      'pathway_results_dir':'results/pathways','gse_sensitivity_dir':'results/gse_sensitivity',
      'gse_robust_genes':'results/gse131617/braak_robust_genes.csv',
      'stable_pathways':'results/pathways/stable_pathways.csv',
      'manuscript_tables_dir':'results/manuscript_tables'}
    runconf['paths'].update({k:str(out/v) for k,v in outputs.items()})
    runconfig=out/'config/project.yaml'
    s=root/'scripts'
    jobs=[]
    for name in ['prepare_adni.py','analyze_adni.py','prepare_gse131617.R','analyze_gse131617.R','prepare_pathway_expression.R','analyze_pathways.R','summarize_uty_sensitivity.R','compile_results_tables.py']:
        jobs.append([r if name.endswith('.R') else sys.executable,str(s/name),'--config',str(runconfig)])
    jobs.append([sys.executable,str(s/'run_additional.py'),'--project-root',str(out),'--output-dir','results/extensions','--bootstrap',str(a.bootstrap),'--selection-bootstrap',str(a.selection_bootstrap),'--seed',str(a.seed)])
    jobs.append([sys.executable,str(s/'collect_paper_results.py'),'--run-dir',str(out)])
    for cmd in jobs:
        if not Path(cmd[1]).is_file(): p.error('Missing script: '+cmd[1])
    if a.dry_run:
        for cmd in jobs: print(subprocess.list2cmdline(cmd))
        print('Input and Python import checks complete. No analyses executed; R packages not yet checked.')
        return
    subprocess.run([r,'-e','pkgs <- c("yaml","readxl","DBI","RSQLite","limma","AnnotationDbi","org.Hs.eg.db","GO.db","huex10sttranscriptcluster.db"); missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]; if(length(missing)) stop(paste("Missing R packages:", paste(missing,collapse=", ")))'],check=True)
    runconfig.parent.mkdir(parents=True,exist_ok=True)
    runconfig.write_text(yaml.safe_dump(runconf,allow_unicode=True,sort_keys=False),encoding='utf-8')

    private=out/'data/private/adni/raw'; private.mkdir(parents=True)
    for name in ['CDR.csv','DXSUM.csv','ADSL.csv','DATA_DOWNLOADED_DATE.csv']:
        shutil.copy2(raw/name,private/name)
    (out/'.gitignore').write_text('*\n',encoding='utf-8')
    for i,cmd in enumerate(jobs,1):
        print(f'\n[{i}/{len(jobs)}] '+subprocess.list2cmdline(cmd),flush=True)
        subprocess.run(cmd,cwd=root,check=True)
    print('\nCompleted. Review '+str(out/'paper_results'))
    print('The run directory contains private inputs. Share only the reviewed paper_results folder.')

if __name__=='__main__': main()
