"""Source dimensions and manuscript table bindings."""
import csv
import gzip
from pathlib import Path

COVERAGE='targeted/qc_joint/scan_date_braak_region_counts.csv'
FEATURES='targeted/qc_joint/features_baseline.csv'
GENES='data/public/GSE106241/expression_gene.csv.gz'

def retained_counts(work):
    with (work/COVERAGE).open(encoding='utf-8-sig',newline='') as stream:
        reader=csv.DictReader(stream)
        if reader.fieldnames!=['Var1','Var2','Var3','Freq']:raise ValueError('Unknown QC scan coverage schema')
        keys=set();counts=0
        for r in reader:
            k=tuple(r[n] for n in ('Var1','Var2','Var3'));v=int(r['Freq'])
            if k in keys or not all(k) or v<0:raise ValueError('Invalid QC scan coverage cell')
            keys.add(k);counts+=v
        if len(keys)!=216 or [len({k[i] for k in keys}) for i in range(3)]!=[18,4,3]:
            raise ValueError('Expected 18 scan dates x 4 Braak stages x 3 regions')
    with gzip.open(work/GENES,'rt',encoding='utf-8-sig',newline='') as stream:
        rows=csv.reader(stream);header=next(rows)
        if header[0]!='entrez_id' or len(header)!=61 or len(set(header))!=61:raise ValueError('Expected 60 GSE106241 donors')
        ids=set()
        for r in rows:
            if len(r)!=len(header) or not r[0] or r[0] in ids:raise ValueError('Invalid GSE106241 expression matrix')
            ids.add(r[0])
    return [('S19/QC-retained scan-date coverage',212,counts),('S21B/mapped genes',16597,len(ids))]

def table_bindings(template, registered):
    with template.open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
    if len(rows)!=106 or len({(r['table'],r['panel'],r['record_id']) for r in rows})!=len(rows):
        raise ValueError('Incomplete or duplicated November 6 table bindings')
    for r in rows:
        key,src=r['record_id'],r['source']
        if key.startswith('D3_'):
            if src!='Supplementary_Data_3/'+key[3:]+'.csv':raise ValueError(src)
        elif key=='Supplementary_Data_2':
            if src!=FEATURES:raise ValueError(src)
        elif registered.get(src)!=key:raise ValueError('Unresolved manuscript table binding: '+str(r))
    return rows

