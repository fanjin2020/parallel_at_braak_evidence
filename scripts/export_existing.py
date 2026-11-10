"""Rebuild publication files from a completed V4 run without fitting models."""
import json
from datetime import datetime
from pathlib import Path
import shutil
import zipfile
import io
import pandas as pd
from pipeline_runtime import inventory, sha, write_json, publish
from source_schema import read_csv_exact


def fitted_workspace(path):
    path = Path(path)
    return ((path/'COMPLETED.json').is_file() and
            all((path/'.work'/name).is_dir() for name in
                ['models','presentation','observation','survival/results']))


def locate_run(path):
    path = Path(path).resolve()
    if fitted_workspace(path):
        return path
    runs = path/'runs'
    candidates = [p for p in runs.iterdir() if p.is_dir() and fitted_workspace(p)] if runs.is_dir() else []
    if candidates:
        return max(candidates,key=lambda p:(p/'COMPLETED.json').stat().st_mtime_ns)
    raise FileNotFoundError('No completed V4 model workspace under '+str(path)+
                           '. Supply the original run directory containing .work, audit and COMPLETED.json. '
                           'A publication-only results_to_review folder is insufficient; no models will be rerun.')


def resolve_source(source,base,models,presentation):
    if source.startswith('review/'):
        return models/source.removeprefix('review/')
    if source.startswith('publication/'):
        return presentation/'publication_tables'/source.removeprefix('publication/')
    if source.startswith('closeout/sex_interaction/'):
        return models/'sex_interaction'/Path(source).name
    if source.startswith('closeout/figure_sources/'):
        return presentation/'figure_sources'/Path(source).name
    return base/source


def existing_sources(run):
    run = Path(run); work = run/'.work'
    record = json.loads((run/'COMPLETED.json').read_text(encoding='utf-8'))
    name = record['output_directory']
    if name not in {'submission','results_to_review'} or not record.get('completed'):
        raise ValueError('The original run has no valid completion record')
    public = run/name
    print('Verifying the completed publication files and stored model sources...',flush=True)
    if inventory(public) != record['outputs']:
        raise ValueError('Original publication files changed; refusing to combine unverified files')
    base = work/'primary/.work'
    if not (base/'targeted/qc_joint/features_baseline.csv').is_file():
        sig = run/'audit/private_run_signature.json'
        value = json.loads(sig.read_text(encoding='utf-8')) if sig.is_file() else {}
        candidates = [Path(p) for p in value.get('locations',[])
                      if (Path(p)/'targeted/qc_joint/features_baseline.csv').is_file()]
        if len(candidates) != 1:
            raise FileNotFoundError('Cannot locate the completed primary fit workspace; no raw-data refit was started')
        base = candidates[0]
    sources = dict(base=base,models=work/'models',presentation=work/'presentation',
                   observation=work/'observation',survival=work/'survival/results')
    index = read_csv_exact(public/'SOURCE_INDEX.csv')
    checked = []
    with zipfile.ZipFile(public/'Supplementary_Data_3.zip') as z:
        for row in index.to_dict('records'):
            if row['record_id'].startswith('D3_'):
                name = row['file'].split('.zip/',1)[1]
                parent = sources['observation'] if name.startswith('Observation_') else sources['survival']
                actual = read_csv_exact(parent/name)
                exported = pd.read_csv(io.BytesIO(z.read(name)),dtype=str,keep_default_na=False)
                actual = actual[list(exported)]
                if not actual.reset_index(drop=True).equals(exported):
                    raise ValueError('Stored model output differs from its completed export: '+name)
            else:
                path = resolve_source(row['source'],base,sources['models'],sources['presentation'])
                if not path.is_file() or sha(path) != row['source_sha256']:
                    raise ValueError('Stored source changed or is missing: '+row['source'])
            checked.append(row['record_id'])
            if len(checked)%25==0 or len(checked)==len(index):
                print(f'  Stored sources verified: {len(checked)}/{len(index)}',flush=True)
    a = read_csv_exact(base/'targeted/qc_joint/features_baseline.csv')
    b = read_csv_exact(public/'Supplementary_Data_2.csv')
    if not a.equals(b):
        raise ValueError('Stored QC baseline differs from completed Data 2')
    # Full processing outputs were local-only. Verify their saved step hashes when available.
    primary = run/'audit/steps/primary.json'
    if primary.is_file():
        saved = json.loads(primary.read_text(encoding='utf-8'))
        if saved.get('status') != 'complete':
            raise ValueError('Primary model workflow was not completed')
        for name,digest in saved['outputs'].items():
            if '/cel_qc/features_' in name or name.endswith('/cel_qc/analysis_status.csv'):
                path = work/'primary'/name
                if not path.is_file() or sha(path) != digest:
                    raise ValueError('Full-feature processing result changed: '+name)
    return sources, checked


def run(args,code_root):
    if args.resume or args.check or args.prepared_inputs:
        raise ValueError('--export-only uses a new output directory and cannot be combined with --resume, --check or --prepared-inputs')
    input_options = ['adni_raw_dir','adnimerge2_dir','confounder_raw_dir','cache_dir','public_cache','adni_tables','private_bundle','data_root']
    if any(getattr(args,key,None) is not None for key in input_options):
        raise ValueError('--export-only reads its completed run; do not supply raw-data options')
    code_root = Path(code_root)
    previous = locate_run(args.from_run or code_root)
    print('Existing completed run:',previous,flush=True)
    print('Mode: publication export only. No downloads, installations or model fits.',flush=True)
    sources,checked = existing_sources(previous)
    print('Stored results verified. Rebuilding publication files only.',flush=True)
    out = (args.output_dir or code_root/'runs'/('v4_export_'+datetime.now().strftime('%Y%m%d_%H%M%S'))).resolve()
    if out.exists():
        raise ValueError('Export output already exists; choose a new directory')
    for protected in [previous,code_root,*sources.values()]:
        if out == protected or out in protected.parents:
            raise ValueError('Export output would replace an input/code directory')
    if previous in out.parents or any(p in out.parents for p in sources.values()):
        raise ValueError('Export output cannot be inside the original model workspace')
    out.mkdir(parents=True)
    audit=out/'audit';audit.mkdir()
    presentation=out/'.work/presentation'
    shutil.copytree(sources['presentation'],presentation)
    write_json(audit/'REUSED_MODEL_RUN.json',dict(mode='export_only',source_run=str(previous),
        verified_export_records=checked,model_fits_started=0,downloads_started=0,
        source_completion_sha256=sha(previous/'COMPLETED.json'),code_manifest_sha256=sha(code_root/'CODE_SHA256SUMS.csv')))
    shutil.copyfile(code_root/'CODE_SHA256SUMS.csv',audit/'CODE_SHA256SUMS.csv')
    from publication_export import export
    stage=out/'.submission_build'
    report=export(sources['base'],sources['models'],presentation,sources['observation'],sources['survival'],stage,audit)
    report.update(mode='export_only',model_fits_started=0,downloads_started=0)
    name='submission' if report['manuscript_match'] else 'results_to_review'
    publish(stage,out,name,report)
    if not args.keep_work:
        work=out/'.work'
        if work.is_symlink() or work.resolve().parent != out:
            raise ValueError('Unsafe export staging path')
        shutil.rmtree(work)
    print('Outputs:',out/name,flush=True)
    print('Original model run preserved:',previous,flush=True)
    print('Reference match:',report['manuscript_match'],'; Word documents were not changed.',flush=True)
