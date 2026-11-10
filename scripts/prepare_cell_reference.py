"""Prepare published adult-brain single-nucleus reference profiles."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import pandas as pd

REFERENCE = 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-022-28655-4/MediaObjects/41467_2022_28655_MOESM8_ESM.xlsx'
GENES = 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz'
CELL_TYPES = ['Astrocytes','Endothelia','Microglia','Oligodendrocytes','OPCs','Neurons']


def run(out, workbook=None, gene_info=None):
    out.mkdir(parents=True,exist_ok=True)
    if (out/'reference_provenance.json').exists():raise ValueError('Choose an empty reference output directory')
    book=Path(workbook) if workbook else out/'Sutton2022_Data5.xlsx'
    info=Path(gene_info) if gene_info else out/'human_gene_info.gz'
    if not workbook or not gene_info:raise ValueError('Supply local reference workbook and gene annotation.')
    mapping={}
    with gzip.open(info,'rt',encoding='utf-8') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['#tax_id']!='9606':continue
            for ref in r['dbXrefs'].split('|'):
                if ref.startswith('Ensembl:ENSG'):
                    mapping.setdefault(ref.split(':')[1],set()).add(r['GeneID'])
    records=[]
    for sheet,source in [('CA','Hodge et al. 2019; adult human middle temporal gyrus; SmartSeq2 snRNA-seq'),('NG','Nagy et al. 2020; post-mortem prefrontal cortex; 10X snRNA-seq')]:
        frame=pd.read_excel(book,sheet_name=sheet).rename(columns=lambda x:'ensembl_id' if str(x).startswith('Unnamed:') else x)
        if not set(CELL_TYPES).issubset(frame):raise ValueError('Reference cell-type headers changed')
        frame['entrez_id']=frame.ensembl_id.map(lambda g: next(iter(mapping[g])) if g in mapping and len(mapping[g])==1 else None)
        frame=frame.dropna(subset=['entrez_id'])
        # Multiple Ensembl features are combined without looking at Braak labels or effects.
        ref=frame.groupby('entrez_id')[CELL_TYPES].median().sort_index()
        if ref.isna().any().any() or (ref<0).any().any():raise ValueError('Invalid linear reference expression')
        ref.to_csv(out/f'reference_{sheet}.csv')
        records.append(dict(reference=sheet,source=source,genes=len(ref),cell_types=CELL_TYPES,
                            expression='Published normalized linear expression; not cell-level counts',
                            limitation='Cross-platform and cross-region sensitivity reference, not a measured cell count or a matched validation cohort'))
    record=dict(publication='https://doi.org/10.1038/s41467-022-28655-4',url=REFERENCE,
                workbook_sha256=hashlib.sha256(book.read_bytes()).hexdigest(),
                dictionary_url=GENES,dictionary_sha256=hashlib.sha256(info.read_bytes()).hexdigest(),references=records,
                mapping='Unique human Ensembl to Entrez from NCBI; median if multiple Ensembl identifiers map to one Entrez; no alias guessing')
    (out/'reference_provenance.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    print([(r['reference'],r['genes']) for r in records])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--workbook',type=Path);p.add_argument('--gene-info',type=Path)
    a=p.parse_args();run(a.output_dir,a.workbook,a.gene_info)
