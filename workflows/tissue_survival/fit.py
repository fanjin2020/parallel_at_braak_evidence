"""Use local ADNI R-package tables and public GSE174367 counts for two follow-ups."""
import os
for name in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[name]='1'
import argparse
from datetime import datetime
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tarfile
import time

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
sys.path[:0]=[str(HERE),str(REPO/'scripts')]
EXECUTED_R_SCRIPTS=[]


def worker(command, log):
    for value in command:
        p=Path(str(value))
        if p.suffix=='.R' and p.is_file():
            EXECUTED_R_SCRIPTS.append(dict(file=p.relative_to(HERE.parents[1]).as_posix(),
                sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    with log.open('w',encoding='utf-8') as stream:
        proc=subprocess.Popen(list(map(str,command)),stdout=stream,stderr=subprocess.STDOUT)
        started=time.monotonic()
        while proc.poll() is None:
            try:
                proc.wait(timeout=20)
            except subprocess.TimeoutExpired:
                print(f'{log.stem}: running {int(time.monotonic()-started)} s; log={log}',flush=True)
            except KeyboardInterrupt:
                proc.terminate();proc.wait();raise
    print(log.read_text(encoding='utf-8',errors='replace')[-6000:],flush=True)
    if proc.returncode:
        raise RuntimeError('Step failed: '+str(log))


def unpack_private(archive, destination):
    allowed={'STUDYSUM.csv','REGISTRY.csv','ADVERSE.csv','RECADV.csv','TREATDIS.csv','OTDISC.csv',
             'DATADIC.csv','DATA_DOWNLOADED_DATE.csv','inventory.csv','R_sessionInfo.txt','SHA256.json'}
    with tarfile.open(archive,'r:gz') as tar:
        members={}
        for item in tar.getmembers():
            if item.isdir():
                continue
            parts=Path(item.name).parts
            if not item.isfile() or len(parts)!=2 or parts[0]!='ADNI_vital_local_20261101_PRIVATE' or parts[1] not in allowed:
                raise ValueError('Unexpected private archive member')
            if parts[1] in members:
                raise ValueError('Duplicate private archive member')
            members[parts[1]]=item
        if 'SHA256.json' not in members:
            raise ValueError('Private archive needs a checksum manifest')
        checks=json.load(tar.extractfile(members['SHA256.json']))
        if set(checks)!=set(members)-{'SHA256.json'}:
            raise ValueError('Private archive manifest does not cover all files')
        destination.mkdir(parents=True,mode=0o700,exist_ok=True)
        for name,item in members.items():
            target=destination/name
            if target.exists():
                if not target.is_file() or name not in checks or sha256_file(target)!=checks[name]:
                    if name=='SHA256.json' and target.is_file() and json.loads(target.read_text())==checks:
                        continue
                    raise ValueError('Existing extracted private input differs: '+name)
                continue
            temporary=target.with_suffix(target.suffix+'.partial')
            with tar.extractfile(item) as source,temporary.open('wb') as output:
                shutil.copyfileobj(source,output)
            temporary.chmod(0o600)
            if name in checks and sha256_file(temporary)!=checks[name]:
                raise ValueError('Private input checksum mismatch: '+name)
            temporary.replace(target)
    return destination


def sha256_file(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def main(a):
    if platform.system()!='Linux':raise RuntimeError('Run the analysis entry point on Linux; unit tests are platform independent')
    if sys.version_info[:2]!=(3,12):raise RuntimeError('Use Python 3.12 in paper2-r442')
    from pipeline_runtime import environment_versions
    from common import write_csv, write_json, sha256, NotEstimable
    packages=environment_versions(a.rscript)
    out=(a.output_dir or HERE/'runs'/('tissue_survival_'+datetime.now().strftime('%Y%m%d_%H%M%S'))).resolve()
    if out.exists() and any(out.iterdir()) and not a.resume:raise ValueError('Use --resume for an existing checked run')
    base=a.base_run.resolve()
    if not (base/'data/private/adni/adni_pet_aligned_master.csv').is_file() and (base/'.work').is_dir():base=base/'.work'
    if out==base or out in base.parents or base in out.parents:raise ValueError('Input and output directories must be separate')
    out.mkdir(parents=True,exist_ok=True);work=out/'_private';work.mkdir(mode=0o700,exist_ok=True);final=out/'results';final.mkdir(exist_ok=True)
    old_plan=json.loads((REPO/'config/survival_plan.json').read_text())
    plan=dict(version='2026-11-01.1',mnar=old_plan['mnar'],seed=old_plan['seed'],
              snrna=dict(cohort='GSE174367',contrast='Braak V-VI versus I-II',primary_cell_type='MG',
                         primary_covariates=['age','sex','batch','RIN'],sensitivity='PMI-complete donors with PMI adjustment',
                         exact_Braak0='reported as not estimable if no stage-0 donors',
                         cell_types=['MG','ASC','ODC','OPC','EX','INH','PER.END'],minimum_cells=20,minimum_group_donors=6,
                         target_multiplicity='Primary MG: 16 fixed hypotheses. All cell/model target family: 224 hypotheses, untested P=1.',
                         GO_multiplicity='4826 fixed sets per cell/model; untestable sets P=1. Additional all-cell/model BH reported.'))
    write_json(final/'Analysis_plan.json',plan);states=[];errors=[]
    def state(module,status,detail):
        states.append(dict(module=module,status=status,detail=detail));write_json(out/'WORKFLOW_STATUS.json',states)
        print(module+': '+status,flush=True)
    if a.only!='snrna':
        try:
            tables=a.adni_tables
            if tables is None:
                bundle=a.private_bundle
                default=HERE.parents[2]/'ADNI_vital_local_20261101_PRIVATE.tar.gz'
                if bundle is None and default.is_file():bundle=default
                if bundle:
                    tables=work/'adni_source_tables'
                    tables=unpack_private(bundle,tables)
                else:raise FileNotFoundError('Local ADNI event tables or private archive are required.')
            from vital_followup import run
            result=run(base,tables,work/'vital_followup',final,plan['mnar'],plan['seed'],a.diagnostic)
            ok=result['successful_scenarios']==result['planned_scenarios']
            state('ADNI_survival_MNAR','diagnostic_only' if a.diagnostic else 'completed' if ok else 'fit_failures',result)
            if not ok:errors.append('MNAR scenario fit failures')
        except NotEstimable as e:state('ADNI_survival_MNAR','not_estimable',str(e))
        except Exception as e:state('ADNI_survival_MNAR','failed',str(e));errors.append(str(e))
    else:state('ADNI_survival_MNAR','not_requested','')
    if a.only!='adni':
        try:
            from data_sources import obtain
            from morabito_pseudobulk import SPEC, prepare
            import pandas as pd
            go=a.go_dir or base/'revision/go_annotation'
            for name in ['go_gene_membership.csv','human_gene_dictionary.csv']:
                if not (go/name).is_file():raise FileNotFoundError('Missing retained annotation '+str(go/name))
            paths,provenance=obtain(SPEC,a.public_cache)
            report=prepare(paths['counts'],paths['metadata'],work/'GSE174367')
            write_json(final/'GSE174367_source_provenance.json',dict(files=provenance,aggregation=report,
                GO=[dict(file=name,sha256=sha256(go/name)) for name in ['go_gene_membership.csv','human_gene_dictionary.csv']]))
            worker([a.rscript,'--vanilla',HERE/'morabito_models.R',work/'GSE174367',go,final],out/'GSE174367.log')
            write_json(final/'Interaction_plan.json',dict(
                pair='OPC minus INH',stage_contrast='Braak V-VI minus I-II',
                interpretation='Exploratory follow-up selected after inspecting the preceding cell-type results; not a prespecified independent replication.',
                population='Same donors with at least 20 nuclei in each cell type; at least six low- and six late-stage donors.',
                design='Cell-type-specific stage and covariate effects; donor-blocked limma-voom with genome-wide duplicateCorrelation.',
                covariates='age, sex, batch and RIN; additional PMI-complete sensitivity',
                normalization='Within-cell-type TMM; expression required in both cell types, with full eligible expression background fitted.',
                multiplicity='16 targets and 4826 GO sets per model; global BH over two models, untested hypotheses P=1.',
                method_reference='https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf'))
            worker([a.rscript,'--vanilla',HERE/'celltype_interaction.R',work/'GSE174367',go,final],out/'interaction.log')
            statuses=pd.read_csv(final/'GSE174367_model_status.csv')
            for name in ['donor_cell_counts.csv','donors.csv']:
                shutil.copyfile(work/'GSE174367'/name,final/('GSE174367_'+name))
            state('GSE174367','completed_with_model_statuses' if statuses.status.eq('estimated').any() else 'not_estimable',
                  dict(models_estimated=int(statuses.status.eq('estimated').sum()),models_planned=len(statuses)))
        except Exception as e:state('GSE174367','failed',str(e));errors.append(str(e))
    else:state('GSE174367','not_requested','')
    write_json(final/'Environment.json',dict(python=sys.version,packages=packages,diagnostic_only=a.diagnostic))
    write_json(out/'R_invocations.json',EXECUTED_R_SCRIPTS)
    write_json(out/'COMPLETED.json',dict(workflow_finished=True,errors=errors))
    print('Review aggregate results:',final,flush=True)
    print('Keep _private and the ADNI source bundle off GitHub.',flush=True)
    if errors:raise SystemExit(1)


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base-run',type=Path,required=True)
    source=p.add_mutually_exclusive_group()
    source.add_argument('--adni-tables',type=Path)
    source.add_argument('--private-bundle',type=Path)
    p.add_argument('--go-dir',type=Path)
    p.add_argument('--public-cache',type=Path,default=HERE/'data/public/GSE174367')
    p.set_defaults(only='all',diagnostic=False)
    p.add_argument('--rscript',default='Rscript')
    p.add_argument('--output-dir',type=Path)
    p.add_argument('--resume',action='store_true')
    return p


if __name__=='__main__':main(parser().parse_args())
