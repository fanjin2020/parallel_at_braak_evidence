"""Local public-source loading and source identifiers; no network calls."""
import csv, gzip, re
from pathlib import Path
SOURCES = [
    ('discovery_matrix','GSE131617/GSE131617-GPL5175_series_matrix.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131617/matrix/GSE131617-GPL5175_series_matrix.txt.gz'),
    ('published_donor_supplement','GSE131617/published_supplement.xls','https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Ftp.2014.35/MediaObjects/41398_2014_BFtp201435_MOESM249_ESM.xls'),
    ('cel_file_list','GSE131617/filelist.txt','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131617/suppl/filelist.txt'),
    ('external_AD_matrix','GSE44770/GSE44770_series_matrix.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE44nnn/GSE44770/matrix/GSE44770_series_matrix.txt.gz'),
    ('external_AD_platform','GSE44770/GPL4372_platform_table.soft.gz','https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL4nnn/GPL4372/soft/GPL4372_family.soft.gz'),
    ('external_Braak_soft','GSE106241/GSE106241.soft.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE106nnn/GSE106241/soft/GSE106241_family.soft.gz'),
    ('external_Braak_matrix','GSE106241/GSE106241_Data_matrix_full.csv.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE106nnn/GSE106241/suppl/GSE106241_Data_matrix_full.csv.gz'),
    ('cell_reference_workbook','cell_reference/Sutton2022_Data5.xlsx','https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-022-28655-4/MediaObjects/41467_2022_28655_MOESM8_ESM.xlsx'),
    ('gene_info','annotations/Homo_sapiens.gene_info.gz','https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz'),
    ('go_gaf','annotations/goa_human.gaf.gz','https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz'),
]

def matrix_sample_ids(path):
    with gzip.open(path,'rt',encoding='utf-8') as f:
        for line in f:
            if line.startswith('!Sample_geo_accession\t'):
                ids=next(csv.reader([line],delimiter='\t'))[1:]
                if not ids or len(set(ids))!=len(ids) or not all(re.fullmatch(r'GSM\d+',s) for s in ids):
                    raise ValueError('Invalid or duplicate discovery GEO sample IDs')
                return ids
    raise ValueError('No discovery sample-accession row in the upstream matrix')

def cel_selection(sample_ids,file_list):
    ids=set(sample_ids)
    with Path(file_list).open(encoding='utf-8-sig') as f:rows=list(csv.DictReader(f,delimiter='\t'))
    selected=[r for r in rows if r['Type'].strip()=='CEL' and r['Name'].split('_')[0] in ids]
    if len(selected)!=len(ids) or {r['Name'].split('_')[0] for r in selected}!=ids:
        raise ValueError('Samples do not map one-to-one to the official CEL file list')
    for r in selected:
        if Path(r['Name']).name!=r['Name'] or int(r['Size'])<=0:raise ValueError('Unsafe CEL file list entry')
    return sorted(selected,key=lambda r:r['Name'])

def gaf_release(gaf):
    with gzip.open(gaf,'rt',encoding='utf-8') as f:
        for line in f:
            if not line.startswith('!'):break
            if line.lower().startswith('!go-version:'):
                match=re.search(r'\d{4}-\d{2}-\d{2}',line)
                if match:return match[0]
    raise ValueError('The GAF does not declare a GO release; supply documented matching GAF/OBO sources')



import csv
import json
from pathlib import Path
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


def main_cache_manifest(cache):
    path=Path(cache)/'public_sources_manifest.json'
    if not path.is_file():raise FileNotFoundError('Main public cache manifest missing: '+str(path))
    try:
        saved=json.loads(path.read_text(encoding='utf-8'))
        needed={key for key,_,_ in SOURCES}|{'go_ontology'}
        records=saved['sources'];paths=saved['paths'];cells=saved['cel_files']
        keys=[r['key'] for r in records]
        if len(keys)!=len(set(keys)) or not needed.issubset(keys):raise ValueError('source keys')
        if not (needed|{'cel_dir'}).issubset(paths):raise ValueError('source paths')
        if not cells:raise ValueError('CEL records')
        for row in [*records,*cells]:
            if not re.fullmatch(r'[0-9a-fA-F]{64}',row['sha256']) or int(row['bytes'])<=0:
                raise ValueError('source identity')
        if any(not r.get('sample_id') for r in cells):raise ValueError('CEL sample IDs')
    except (ValueError,KeyError,TypeError) as error:
        raise ValueError('Not a complete main public-data manifest: '+str(path)+
                         '. A single-nucleus manifest cannot be used as --cache-dir.') from error
    return saved


def resolved_inputs(cache):
    """Resolve the exact files used by the cache, including relocated snapshots."""
    cache = Path(cache).resolve()
    saved = main_cache_manifest(cache)
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
    saved = main_cache_manifest(cache)
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

from pathlib import Path
COUNTS = {'GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5': '6ba98a1af8772af08c8cdd5e5e63eebb0890daed48cac2e1b8961dcc67069b77',
          'GSE174367_snRNA-seq_cell_meta.csv.gz': 'ab1a029deb43196e2bb1fea5907d838885750cfff3850c600c07588c7c7cdb2b'}
from source_io import sha256

def obtain(spec,cache):
    paths={};provenance={}
    for item in spec['files']:
        name=item['filename'];file=Path(cache)/name
        if not file.is_file():raise FileNotFoundError(f'Missing cached public file: {file}. Nothing was downloaded.')
        digest=sha256(file)
        if digest!=COUNTS[name]:raise ValueError('GSE174367 checksum mismatch: '+name)
        paths[item['id']]=file
        provenance[item['id']]=dict(filename=name,url=item['url'],sha256=digest,bytes=file.stat().st_size)
    return paths,provenance

