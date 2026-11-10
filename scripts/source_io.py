import csv, gzip, hashlib, json
from pathlib import Path
def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def local_cel_metadata(cel_dir, out, verified_records=None):
    from cel_header import parse_header
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    records = []
    if verified_records is None:
        with (Path(cel_dir) / 'downloaded_arrays.csv').open(encoding='utf-8-sig') as f:
            verified_records = list(csv.DictReader(f))
    for row in verified_records:
        p = Path(cel_dir) / Path(row['filename']).name
        if not p.is_file() or sha256(p) != row['sha256']:
            raise ValueError(f'CEL checksum changed or file missing: {p.name}')
        with gzip.open(p, 'rb') as g:
            header = parse_header(g.read(100024))
        records.append(dict(sample_id=row['sample_id'], array_scan_date=header['array_scan_date'],
                            scanner_id=header['scanner_id']))
    if len({r['sample_id'] for r in records}) != len(records):
        raise ValueError('Duplicate CEL metadata')
    with (out / 'sample_technical_metadata.csv').open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['sample_id','array_scan_date','scanner_id'])
        w.writeheader()
        w.writerows(records)
    (out/'metadata_scope.json').write_text(json.dumps(dict(source='Raw CEL DAT headers',
        scan_date_is_not='RNA quality, extraction batch, hybridization batch, or cause of death'),indent=2),encoding='utf-8')

