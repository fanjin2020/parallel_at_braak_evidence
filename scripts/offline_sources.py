"""Verify an existing public-source cache without making network requests."""
import csv
import json
from pathlib import Path
from public_sources import SOURCES, matrix_sample_ids, cel_selection, gaf_release
from source_io import sha256


def checked(path, record):
    if not path.is_file():
        raise FileNotFoundError(f'Missing cached input: {path}. Offline mode will not download it.')
    if path.stat().st_size != int(record['bytes']) or sha256(path) != record['sha256']:
        raise ValueError(f'Cached input changed: {path}. Restore the recorded snapshot before analysis.')
    return path.resolve()


def resolve_path(cache, saved, relative):
    local = cache / relative
    if local.is_file():
        return local
    original = Path(saved)
    if original.is_file():
        return original
    raise FileNotFoundError(f'Missing cached input: {local}. No download attempted.')


def resolved_inputs(cache):
    """Resolve the exact files used by the cache, including relocated snapshots."""
    cache = Path(cache).resolve()
    saved = json.loads((cache/'public_sources_manifest.json').read_text(encoding='utf-8'))
    paths = [resolve_path(cache, saved['paths'][key], relative).resolve()
             for key, relative, _ in SOURCES]
    paths.append(resolve_path(cache, saved['paths']['go_ontology'],
                              'annotations/go-basic.obo').resolve())
    cel = cache/'GSE131617/CEL'
    if not cel.is_dir():
        cel = Path(saved['paths']['cel_dir'])
    if not cel.is_dir():
        raise FileNotFoundError('Cached CEL directory missing: '+str(cel))
    paths.append(cel.resolve())
    return list(dict.fromkeys(paths))


def load_cached_public(cache, progress=True):
    cache = Path(cache).resolve()
    manifest = cache / 'public_sources_manifest.json'
    if not manifest.is_file():
        raise FileNotFoundError(f'Existing public cache manifest not found: {manifest}. '
                                'Point --cache-dir to the cache used by the previous complete run.')
    saved = json.loads(manifest.read_text(encoding='utf-8'))
    records = {r['key']: r for r in saved['sources']}
    if len(records) != len(saved['sources']):
        raise ValueError('Duplicate source keys in the cache manifest')
    paths = {}
    for key, relative, _ in SOURCES:
        if key not in records:
            raise ValueError(f'No recorded checksum for {key}')
        paths[key] = checked(resolve_path(cache, saved['paths'][key], relative), records[key])
        if progress:
            print(f'Local public input verified: {key}', flush=True)
    paths['go_ontology'] = checked(resolve_path(cache, saved['paths']['go_ontology'],
                                               'annotations/go-basic.obo'), records['go_ontology'])
    release = gaf_release(paths['go_gaf'])
    versions = [x.strip() for x in paths['go_ontology'].read_text(encoding='utf-8').splitlines()
                if x.startswith('data-version:')]
    if len(versions) != 1 or release not in versions[0]:
        raise ValueError('Cached GAF and ontology release do not match')
    ids = matrix_sample_ids(paths['discovery_matrix'])
    expected = cel_selection(ids, paths['cel_file_list'])
    cel_dir = cache / 'GSE131617/CEL'
    if not cel_dir.is_dir():
        cel_dir = Path(saved['paths']['cel_dir'])
    if not cel_dir.is_dir():
        raise FileNotFoundError(f'Cached CEL directory missing: {cel_dir}')
    cel_records = {r['sample_id']: r for r in saved.get('cel_files', [])}
    if set(cel_records) != set(ids) or len(cel_records) != len(saved.get('cel_files', [])):
        raise ValueError('Cache manifest must record one CEL file for every discovery sample')
    verified = []
    for i, r in enumerate(expected, 1):
        sample = r['Name'].split('_')[0]
        entry = cel_records[sample]
        if Path(entry['filename']).name != r['Name']:
            raise ValueError('CEL manifest filename differs from the official list: '+sample)
        path = checked(cel_dir / r['Name'], entry)
        if path.stat().st_size != int(r['Size']):
            raise ValueError(f'CEL size differs from the official list: {sample}')
        verified.append({**entry, 'filename': str(path)})
        if progress and (i % 20 == 0 or i == len(expected)):
            print(f'Local CEL checksum verification: {i}/{len(expected)}', flush=True)
    paths['cel_dir'] = cel_dir.resolve()
    return {**saved, 'paths': {k: str(v) for k, v in paths.items()},
            'cel_files': verified, 'offline_verified': True, 'network_requests': 0}


def write_cel_index(record, output):
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    rows = record['cel_files']
    with output.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
