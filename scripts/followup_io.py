import gzip, json, os
from pathlib import Path
from source_io import sha256

def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(path) + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False), encoding='utf-8')
    os.replace(temporary, path)


def soft_samples(path):
    """Parse GEO SOFT records; retain all characteristics without guessing field meanings."""
    opener = gzip.open if str(path).endswith('.gz') else open
    samples, current = [], None
    with opener(path, 'rt', encoding='utf-8', errors='strict') as f:
        for line in f:
            if line.startswith('^SAMPLE = '):
                if current is not None:
                    samples.append(current)
                current = {'sample_id': line.strip().split(' = ', 1)[1]}
            elif current is not None and line.startswith('!Sample_') and ' = ' in line:
                key, value = line.rstrip('\r\n').split(' = ', 1)
                key = key.removeprefix('!Sample_')
                if key.startswith('characteristics_ch') and ': ' in value:
                    key, value = value.split(': ', 1)
                    key = key.strip().lower()
                if key in current and current[key] != value:
                    current[key] = str(current[key]) + ' | ' + value
                else:
                    current[key] = value
    if current is not None:
        samples.append(current)
    return samples


def series_url(accession, area, filename):
    return f'https://ftp.ncbi.nlm.nih.gov/geo/series/{accession[:-3]}nnn/{accession}/{area}/{filename}'

