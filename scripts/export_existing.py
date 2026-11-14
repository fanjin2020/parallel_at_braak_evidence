"""Export a completed model workspace without running analyses or downloads."""
from datetime import datetime
import json
from pathlib import Path
import shutil
from pipeline_runtime import write_json,publish,code_identity


def fitted_workspace(path):
    path=Path(path)
    return (path/'COMPLETED.json').is_file() and all((path/'.work'/n).is_dir()
        for n in ['models','presentation','observation','survival/results'])


def locate_run(path):
    path=Path(path).resolve()
    if fitted_workspace(path):return path
    runs=path/'runs'
    choices=[p for p in runs.iterdir() if p.is_dir() and fitted_workspace(p)] if runs.is_dir() else []
    if choices:return max(choices,key=lambda p:(p/'COMPLETED.json').stat().st_mtime_ns)
    raise FileNotFoundError('No completed model workspace at '+str(path)+
        '. Supply a run containing .work and COMPLETED.json, not an exported-data folder.')


def existing_sources(run):
    run=Path(run);work=run/'.work'
    record=json.loads((run/'COMPLETED.json').read_text(encoding='utf-8'))
    if not record.get('completed'):raise ValueError('The source run is incomplete')
    base=work/'primary/.work'
    if not (base/'targeted/qc_joint/features_baseline.csv').is_file():
        path=run/'audit/private_run_signature.json'
        signature=json.loads(path.read_text(encoding='utf-8')) if path.is_file() else {}
        choices=[Path(p) for p in signature.get('locations',[])
                 if (Path(p)/'targeted/qc_joint/features_baseline.csv').is_file()]
        if len(choices)!=1:raise FileNotFoundError('Cannot locate the fitted primary workspace')
        base=choices[0]
    return dict(base=base,models=work/'models',presentation=work/'presentation',
                observation=work/'observation',survival=work/'survival/results')


def run(args,code_root):
    if args.resume:raise ValueError('Export-only creates a new output; --resume is not supported')
    names=['adni_raw_dir','adnimerge2_dir','confounder_raw_dir','cache_dir','public_cache',
           'adni_tables','private_bundle','data_root']
    if any(getattr(args,k,None) is not None for k in names):
        raise ValueError('Export-only reads a completed run; do not supply raw-input arguments')
    code_root=Path(code_root)
    previous=locate_run(args.from_run or code_root)
    sources=existing_sources(previous)
    out=(args.output_dir or code_root/'runs'/('v5_export_'+datetime.now().strftime('%Y%m%d_%H%M%S'))).resolve()
    if out.exists():raise ValueError('Choose a new export output directory')
    for protected in [previous,code_root,*sources.values()]:
        protected=protected.resolve()
        if out==protected or out in protected.parents or (protected!=code_root and protected in out.parents):
            raise ValueError('Export output overlaps its source')
    out.mkdir(parents=True)
    audit=out/'audit';audit.mkdir()
    presentation=out/'.work/presentation'
    shutil.copytree(sources['presentation'],presentation)
    from manuscript_figures import prepare
    prepare(sources['base'],sources['models'],presentation/'publication_tables')
    write_json(audit/'EXPORT_INPUTS.json',dict(source_run=str(previous),model_fits_started=0,
               downloads_started=0,code=code_identity(code_root)))
    from publication_export import export
    stage=out/'.submission_build'
    report=export(sources['base'],sources['models'],presentation,sources['observation'],
                  sources['survival'],stage,audit)
    report.update(mode='export_only',model_fits_started=0,downloads_started=0)
    publish(stage,out,'submission',report)
    if args.discard_work:
        work=out/'.work'
        if work.is_symlink() or work.resolve().parent!=out:raise ValueError('Unsafe work path')
        shutil.rmtree(work)
    print('Final publication files:',out/'submission',flush=True)
