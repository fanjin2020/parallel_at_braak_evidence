"""Recompute the manuscript from local inputs, then export one publication set."""
import argparse
import csv
from datetime import datetime
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
S = ROOT/'scripts'
sys.path.insert(0, str(S))
from pipeline_runtime import Runner, inventory, sha, signature, write_json, publish, recover_publication


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--check', action='store_true')
    p.add_argument('--resume', action='store_true')
    p.add_argument('--keep-work', action='store_true', help='Keep private fit checkpoints locally; never upload them')
    p.add_argument('--export-only', action='store_true', help='Regenerate publication summaries, figures and archives; no model fitting')
    p.add_argument('--from-run', type=Path, help='Completed V4 run, or code directory whose runs contain that run')
    p.add_argument('--server-root', type=Path, default=ROOT.parent)
    p.add_argument('--data-root', type=Path, help='Local input project, not an old statistical-result directory')
    for key in ['adni-raw-dir','cache-dir','public-cache','output-dir','prepared-inputs']:
        p.add_argument('--'+key, type=Path)
    for names in [('adnimerge2-dir','confounder-raw-dir'),('adni-tables','private-bundle')]:
        g = p.add_mutually_exclusive_group()
        for key in names: g.add_argument('--'+key, type=Path)
    return p


def choose(explicit, candidates, marker, label):
    choices = [explicit] if explicit else candidates
    for p in choices:
        if p is not None and (p/marker).is_file(): return p.resolve()
    raise FileNotFoundError(label+' missing. Specify its directory; no download will be attempted.')


def verify_code():
    with (ROOT/'CODE_SHA256SUMS.csv').open(encoding='utf-8-sig') as f:
        rows = list(csv.DictReader(f))
    for row in rows:
        path = ROOT/row['file']
        if not path.is_file() or sha(path) != row['sha256']:
            raise ValueError('Code checksum differs: '+row['file'])


def locate_inputs(a):
    home = a.server_root.resolve()
    data = a.data_root or home/'parallel_at_braak_evidence_linux_20261027_r2'
    historic = home/'parallel_at_braak_evidence_20261101_datafill_r1'
    sn = choose(a.public_cache, [data/'cache/public/GSE174367', data/'data/public/GSE174367',
        historic/'addons/remaining_20261101/data/public/GSE174367'],
        'GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5', 'GSE174367 cache (--public-cache)')
    event_candidates = [data/'private_inputs/adni_event_tables',
        historic/'addons/remaining_20261101/runs/datafill_20261101_192146/_private/adni_source_tables']
    tables = a.adni_tables or next((p for p in event_candidates if (p/'STUDYSUM.csv').is_file()), None)
    if a.private_bundle: tables = None
    bundle = a.private_bundle
    if tables:
        missing = [n for n in ['STUDYSUM','REGISTRY','ADVERSE','RECADV','TREATDIS','DATADIC','DATA_DOWNLOADED_DATE']
                   if not (tables/(n+'.csv')).is_file()]
        if missing: raise FileNotFoundError('Missing event tables: '+', '.join(missing))
        event_args = ['--adni-tables', tables.resolve()]
    else:
        if bundle is None:
            bundle = next((p/'ADNI_vital_local_20261101_PRIVATE.tar.gz' for p in [home,data,historic]
                           if (p/'ADNI_vital_local_20261101_PRIVATE.tar.gz').is_file()), None)
        if bundle is None or not bundle.is_file(): raise FileNotFoundError('Specify --adni-tables or --private-bundle.')
        event_args = ['--private-bundle', bundle.resolve()]
    result = dict(sn=sn, events=event_args, protected=[sn, event_args[1]])
    if a.prepared_inputs:
        if any([a.adni_raw_dir,a.adnimerge2_dir,a.confounder_raw_dir,a.cache_dir]):
            raise ValueError('--prepared-inputs cannot be combined with raw input options.')
        base = choose(None,[a.prepared_inputs,a.prepared_inputs/'.work'],
                      'data/private/adni/adni_pet_aligned_master.csv','Prepared input workspace')
        result.update(base=base, primary_args=[])
        result['protected'].append(base)
    else:
        raw = choose(a.adni_raw_dir,[data/'private_inputs/adni_local_inputs/primary_raw',data/'raw/raw'],
                     'CDR.csv','ADNI CSV (--adni-raw-dir)')
        cache = choose(a.cache_dir,[data/'cache/public'],'public_sources_manifest.json','Public cache (--cache-dir)')
        conf = a.confounder_raw_dir
        package = None if conf else choose(a.adnimerge2_dir,[data/'private_inputs/adni_local_inputs/ADNIMERGE2_local'],
                                           'DESCRIPTION','ADNIMERGE2 package (--adnimerge2-dir)')
        result['primary_args'] = ['--adni-raw-dir',raw,'--cache-dir',cache,
            '--confounder-raw-dir' if conf else '--adnimerge2-dir', (conf or package).resolve()]
        result['protected'] += [raw,cache,(conf or package).resolve()]
        from offline_sources import resolved_inputs
        # A cache manifest may refer to snapshots outside the cache directory.
        # Include those files in the resume signature, not only the manifest.
        result['protected'] += [p for p in resolved_inputs(cache)
                                if cache not in p.parents and p != cache]
    return result


def main(a):
    os.umask(0o077)
    for name in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS',
                 'NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
        os.environ[name]='1'
    os.environ.update(JAD_OFFLINE='1',PIP_NO_INDEX='1',BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS='FALSE')
    if platform.system() != 'Linux' or sys.version_info[:2] != (3,12):
        raise SystemExit('Use Python 3.12 and R 4.4.2 on Linux (paper2-r442).')
    verify_code()
    if a.export_only:
        from export_existing import run
        run(a,ROOT)
        return
    if a.from_run is not None:
        raise ValueError('--from-run is used only with --export-only')
    from prerequisites import check
    packages = check()
    rversion = subprocess.check_output(['Rscript','--vanilla','-e',
        'cat(R.version.string,"\\n");for(p in c("limma","edgeR","statmod","oligo","dtangle")) if(requireNamespace(p,quietly=TRUE))cat(p,as.character(packageVersion(p)),"\\n")'], text=True)
    inputs = locate_inputs(a)
    sys.path.insert(0, str(ROOT/'workflows/tissue_survival'))
    from local_counts import obtain
    from morabito_pseudobulk import SPEC
    obtain(SPEC, inputs['sn'])
    if inputs['events'][0] == '--private-bundle':
        from private_bundle import check_bundle
        check_bundle(inputs['events'][1])
    out = (a.output_dir or ROOT/'runs'/('v4_'+datetime.now().strftime('%Y%m%d_%H%M%S'))).resolve()
    if a.resume and not a.output_dir: raise ValueError('--resume requires --output-dir.')
    for p in inputs['protected']:
        if p==out or p in out.parents or out in p.parents: raise ValueError('Output overlaps an input directory.')
    if ROOT==out or out in ROOT.parents: raise ValueError('Output cannot replace the code directory.')
    work = out/'.work'
    raw_args = [*inputs['primary_args'], '--output-dir',work/'primary']
    if a.check:
        if inputs['primary_args']:
            subprocess.run([sys.executable,str(S/'run_from_sources.py'),*map(str,raw_args),'--check-only'],check=True)
        else:
            from manuscript_bindings import retained_counts
            retained_counts(inputs['base'])
        print('OFFLINE_CHECK_OK. No model fits or downloads performed.')
        return
    if out.exists() and not a.resume: raise ValueError('Existing output: use --resume or select a new directory.')
    out.mkdir(parents=True,exist_ok=True)
    import fcntl
    lock = (out/'RUN.lock').open('a')
    fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
    print('Verifying input and code fingerprints...',flush=True)
    fingerprints = {str(i): inventory(p) for i,p in enumerate(inputs['protected'])}
    env = dict(python=sys.version, packages=packages, R=rversion, threads=1)
    sig = dict(inputs=fingerprints,locations=list(map(str,inputs['protected'])),
               code_manifest=sha(ROOT/'CODE_SHA256SUMS.csv'),environment=env,
               mode='prepared_inputs' if a.prepared_inputs else 'local_raw_sources')
    signature(out/'audit/private_run_signature.json',sig,a.resume)
    write_json(out/'audit/environment.json',env)
    shutil.copyfile(ROOT/'CODE_SHA256SUMS.csv',out/'audit/CODE_SHA256SUMS.csv')
    completed=recover_publication(out)
    if completed:
        if not a.keep_work and completed.get('manuscript_match') and work.exists():
            if work.is_symlink() or work.resolve().parent!=out: raise ValueError('Unsafe work path')
            shutil.rmtree(work)
        print('Run already complete:',out/completed['output_directory']);return
    work.mkdir(exist_ok=True)
    runner=Runner(out/'audit/steps',ROOT)
    if inputs['primary_args']:
        cmd=[sys.executable,'-u',S/'run_from_sources.py',*raw_args,'--resume']
        runner.run('primary',cmd,work/'primary')
        base=work/'primary/.work'
    else: base=inputs['base']
    models=work/'models';models.mkdir(exist_ok=True)
    for name in ['external','clinical','mnar','prediction']:
        runner.run(name,[sys.executable,'-u',S/'analysis_tasks.py',name,'--primary',base,'--output',models/name],models/name)
    for name,script in [('tissue','tissue_influence.R'),('sex_interaction','sex_interaction.R')]:
        runner.run(name,['Rscript','--vanilla',S/script,'--project-root',base,'--output-dir',models/name],models/name)
    runner.run('observation',[sys.executable,'-u',S/'analysis_tasks.py','observation','--primary',base,
                              '--output',work/'observation'],work/'observation')
    runner.run('survival_single_nucleus',[sys.executable,'-u',ROOT/'workflows/tissue_survival/fit.py',
        '--base-run',base,'--go-dir',base/'revision/go_annotation','--public-cache',inputs['sn'],
        *inputs['events'],'--output-dir',work/'survival','--resume'],work/'survival')
    runner.run('presentation',[sys.executable,'-u',S/'analysis_tasks.py','presentation','--primary',base,
        '--models',models,'--output',work/'presentation'],work/'presentation')
    stage=out/'.submission_build'
    if stage.exists():
        if stage.is_symlink() or stage.resolve().parent!=out: raise ValueError('Unsafe staging directory')
        shutil.rmtree(stage)
    from publication_export import export
    report=export(base,models,work/'presentation',work/'observation',work/'survival/results',stage,out/'audit')
    name='submission' if report['manuscript_match'] else 'results_to_review'
    nested=work/'primary/audit'
    if nested.is_dir():shutil.copytree(nested,out/'audit/primary_steps',dirs_exist_ok=True)
    # Commit the completion record before cleanup so an interrupted cleanup is resumable.
    publish(stage,out,name,report)
    if not a.keep_work and report['manuscript_match']:
        if work.is_symlink() or work.resolve().parent!=out: raise ValueError('Unsafe work path')
        shutil.rmtree(work)
    print('Outputs:',out/name)
    print('Local audit records:',out/'audit','(do not upload the whole run folder)')


if __name__=='__main__': main(parser().parse_args())
