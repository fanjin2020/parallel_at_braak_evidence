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
