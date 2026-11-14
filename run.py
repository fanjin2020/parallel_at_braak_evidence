"""Run the V5 analyses from locally available source data."""
import argparse
from datetime import datetime
import json
import os
from pathlib import Path
import shutil
import sys

ROOT=Path(__file__).resolve().parent
S=ROOT/'scripts'
sys.path.insert(0,str(S))
from pipeline_runtime import Runner, inventory, signature, write_json, publish, recover_publication
from pipeline_runtime import code_identity, require_analysis_python, environment_versions, sha
from data_sources import COUNTS, main_cache_manifest


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-root',type=Path,help='Existing source-data project (for example V4); no data are downloaded')
    p.add_argument('--output-dir',type=Path)
    p.add_argument('--resume',action='store_true')
    p.add_argument('--export-only',action='store_true',help='Export saved fitted results without refitting')
    p.add_argument('--from-run',type=Path,help='Completed run used only with --export-only')
    p.add_argument('--discard-work',action='store_true',help='Delete private model checkpoints after successful export')
    p.add_argument('--adni-raw-dir',type=Path)
    c=p.add_mutually_exclusive_group()
    c.add_argument('--adnimerge2-dir',type=Path)
    c.add_argument('--confounder-raw-dir',type=Path)
    p.add_argument('--cache-dir',type=Path)
    p.add_argument('--public-cache',type=Path)
    e=p.add_mutually_exclusive_group()
    e.add_argument('--adni-tables',type=Path)
    e.add_argument('--private-bundle',type=Path)
    return p


def choose(explicit,candidates,required,label,main_cache=False,identity_optional=()):
    choices=[explicit] if explicit is not None else candidates
    found=[];identities={};reasons=[]
    for value in choices:
        if value is None:continue
        p=Path(value).resolve()
        if p in found:continue
        if not p.is_dir():continue
        missing=[name for name in required if not (p/name).is_file()]
        if missing:
            if explicit is not None:reasons.append(str(p)+': missing '+', '.join(missing))
            continue
        if main_cache:
            try:
                manifest=main_cache_manifest(p)
            except (ValueError,FileNotFoundError) as error:
                reasons.append(str(error));continue
            identities[p]=json.dumps(dict(
                sources=sorted((r['key'],r['sha256'],int(r['bytes'])) for r in manifest['sources']),
                cel_files=sorted((r['sample_id'],r['sha256'],int(r['bytes'])) for r in manifest['cel_files'])),sort_keys=True)
        found.append(p)
    if not found:
        detail=('\n'+'\n'.join(reasons)) if reasons else ''
        raise FileNotFoundError(label+' not found in the selected data project. Supply its explicit path. No download was attempted.'+detail)
    if len(found)>1:
        if not main_cache:
            names=[*required,*identity_optional]
            identities={p:tuple((name,sha(p/name) if (p/name).is_file() else None)
                               for name in names) for p in found}
        if len(set(identities.values()))!=1:
            raise ValueError('Different input snapshots for '+label+'; select one explicitly:\n'+'\n'.join(map(str,found)))
    return found[0]


def locate_inputs(a):
    if a.data_root is not None and not a.data_root.is_dir():
        raise FileNotFoundError('Source-data project not found: '+str(a.data_root))
    roots=[(a.data_root or ROOT).resolve()]
    saved=[];saved_by_root={}
    for root in roots:
        saved_by_root[root]=[]
        runs=root/'runs'
        if not runs.is_dir():continue
        signatures=sorted(runs.glob('*/audit/private_run_signature.json'),
                          key=lambda p:p.stat().st_mtime_ns,reverse=True)
        for path in signatures:
            try:
                values=[Path(p) for p in json.loads(path.read_text(encoding='utf-8')).get('locations',[])]
            except (ValueError,TypeError) as error:
                raise ValueError('Invalid saved input locations: '+str(path)) from error
            saved.extend(values);saved_by_root[root].extend(values)
    def choices(*suffixes):
        return list(dict.fromkeys(p for root in roots
                    for p in [*[root/part for part in suffixes],*saved_by_root[root]]))
    from run_from_sources import ADNI_FILES, CONFOUNDER_FILES
    raw=choose(a.adni_raw_dir,choices('private_inputs/adni_local_inputs/primary_raw','raw/raw',
                                    'data/private/adni/raw','.'),ADNI_FILES,'--adni-raw-dir')
    cache=choose(a.cache_dir,choices('cache/public','.'),['public_sources_manifest.json'],'--cache-dir',main_cache=True)
    sn=choose(a.public_cache,choices('cache/public/GSE174367','data/public/GSE174367','.')+
              [cache/'GSE174367'],list(COUNTS),'--public-cache')
    conf=a.confounder_raw_dir
    package=None
    if conf is not None:
        conf=choose(conf,[],CONFOUNDER_FILES,'--confounder-raw-dir')
    else:
        try:
            package=choose(a.adnimerge2_dir,choices('private_inputs/adni_local_inputs/ADNIMERGE2_local',
                           'ADNIMERGE2_local','ADNIMERGE2','.'),
                           ['DESCRIPTION',*[str(Path('data')/Path(n).with_suffix('.rda')) for n in CONFOUNDER_FILES]],'--adnimerge2-dir')
        except FileNotFoundError:
            if a.adnimerge2_dir is not None:raise
            conf=choose(None,choices('private_inputs/adni_local_inputs/confounder_raw',
                        'data/private/adni/confounder_raw','confounder_raw','.'),
                        CONFOUNDER_FILES,'--confounder-raw-dir (or --adnimerge2-dir)')
    event_args=None
    if a.private_bundle is not None:
        if not a.private_bundle.is_file():raise FileNotFoundError(a.private_bundle)
        event_args=['--private-bundle',a.private_bundle.resolve()]
    else:
        try:
            tables=choose(a.adni_tables,choices('private_inputs/adni_event_tables',
                'private_inputs/adni_local_inputs/event_tables','adni_event_tables','.'),
                [n+'.csv' for n in ['STUDYSUM','REGISTRY','ADVERSE','RECADV','TREATDIS','DATADIC','DATA_DOWNLOADED_DATE']],
                '--adni-tables',identity_optional=['OTDISC.csv'])
            event_args=['--adni-tables',tables]
        except FileNotFoundError:
            if a.adni_tables is not None:raise
            bundles=[p for p in saved if p.is_file() and p.name.endswith('_PRIVATE.tar.gz')]
            bundles.extend(root/'ADNI_vital_local_20261101_PRIVATE.tar.gz' for root in roots)
            bundles=list(dict.fromkeys(p.resolve() for p in bundles if p.is_file()))
            if not bundles:raise FileNotFoundError('Supply --adni-tables or --private-bundle; ADNI is never downloaded.')
            if len(bundles)>1 and len({sha(p) for p in bundles})>1:
                raise ValueError('Different private input bundles; select --private-bundle explicitly:\n'+'\n'.join(map(str,bundles)))
            event_args=['--private-bundle',bundles[0]]
    from data_sources import resolved_inputs
    protected=[raw,conf or package,cache,sn,event_args[1]]
    protected.extend(p for p in resolved_inputs(cache) if cache not in p.parents and p!=cache)
    primary=['--adni-raw-dir',raw,'--cache-dir',cache,
             '--confounder-raw-dir' if conf else '--adnimerge2-dir',conf or package]
    for label,value in [('ADNI',raw),('Public cache',cache),('Single-nucleus cache',sn),
                        ('ADNI supplementary inputs',conf or package),('ADNI event inputs',event_args[1])]:
        print(label+': '+str(value),flush=True)
    return dict(primary_args=primary,sn=sn,events=event_args,protected=list(dict.fromkeys(protected)))


def main(a):
    os.umask(0o077)
    for name in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS',
                 'NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
        os.environ[name]='1'
    os.environ.update(JAD_OFFLINE='1',PIP_NO_INDEX='1',BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS='FALSE')
    require_analysis_python()
    if a.export_only:
        from export_existing import run
        run(a,ROOT)
        return
    if a.from_run is not None:raise ValueError('--from-run requires --export-only')
    if a.resume and a.output_dir is None:raise ValueError('--resume requires --output-dir')
    env=environment_versions()
    inputs=locate_inputs(a)
    print('Reading local input identities; no downloads or model fitting at this stage.',flush=True)
    for name,digest in COUNTS.items():
        if sha(inputs['sn']/name)!=digest:raise ValueError('Cached GSE174367 file changed: '+name)
    out=(a.output_dir or ROOT/'runs'/('v5_'+datetime.now().strftime('%Y%m%d_%H%M%S'))).resolve()
    for p in inputs['protected']:
        if p==out or p in out.parents or out in p.parents:raise ValueError('Output overlaps a source-data path')
    if out==ROOT or out in ROOT.parents:raise ValueError('Output cannot replace the code directory')
    if out.exists() and not a.resume:raise ValueError('Use a new output directory or --resume')
    out.mkdir(parents=True,exist_ok=True)
    import fcntl
    lock=(out/'RUN.lock').open('a')
    fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
    sig=dict(inputs={str(i):inventory(p) for i,p in enumerate(inputs['protected'])},
             locations=list(map(str,inputs['protected'])),code=code_identity(ROOT),environment=env)
    signature(out/'audit/private_run_signature.json',sig,a.resume)
    write_json(out/'audit/environment.json',env)
    completed=recover_publication(out)
    if completed:
        print('Run already complete:',out/completed['output_directory'])
        return
    work=out/'.work';work.mkdir(exist_ok=True)
    runner=Runner(out/'audit/steps',ROOT)
    primary_args=[*inputs['primary_args'],'--output-dir',work/'primary','--resume']
    runner.run('primary',[sys.executable,'-u',S/'run_from_sources.py',*primary_args],work/'primary')
    base=work/'primary/.work'
    models=work/'models';models.mkdir(exist_ok=True)
    for name in ['external','clinical','mnar','prediction']:
        runner.run(name,[sys.executable,'-u',S/'analysis_tasks.py',name,
                        '--primary',base,'--output',models/name],models/name)
    for name,script in [('tissue','tissue_influence.R'),('sex_interaction','sex_interaction.R')]:
        runner.run(name,['Rscript','--vanilla',S/script,'--project-root',base,
                        '--output-dir',models/name],models/name)
    runner.run('observation',[sys.executable,'-u',S/'analysis_tasks.py','observation',
                            '--primary',base,'--output',work/'observation'],work/'observation')
    runner.run('survival_single_nucleus',[sys.executable,'-u',ROOT/'workflows/tissue_survival/fit.py',
        '--base-run',base,'--go-dir',base/'revision/go_annotation','--public-cache',inputs['sn'],
        *inputs['events'],'--output-dir',work/'survival','--resume'],work/'survival')
    runner.run('presentation',[sys.executable,'-u',S/'analysis_tasks.py','presentation','--primary',base,
                              '--models',models,'--output',work/'presentation'],work/'presentation')
    stage=out/'.submission_build'
    if stage.exists():
        if stage.is_symlink() or stage.resolve().parent!=out:raise ValueError('Unsafe staging path')
        shutil.rmtree(stage)
    from publication_export import export
    report=export(base,models,work/'presentation',work/'observation',work/'survival/results',stage,out/'audit')
    publish(stage,out,'submission',report)
    if a.discard_work:
        if work.is_symlink() or work.resolve().parent!=out:raise ValueError('Unsafe work path')
        shutil.rmtree(work)
    print('Final publication files:',out/'submission')
    print('Private checkpoints and execution logs are local only:',out)


if __name__=='__main__':
    main(parser().parse_args())
