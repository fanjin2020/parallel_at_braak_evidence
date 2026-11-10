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

