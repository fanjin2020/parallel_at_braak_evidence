"""Aggregate public GSE174367 counts by donor and author-assigned cell type."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

CELL_TYPES = ['MG', 'ASC', 'ODC', 'OPC', 'EX', 'INH', 'PER.END']
TARGETS = ['MTHFD2','OXTR','HOPX','CXCL10','IRF1','CD14','TAP1','PTPN6',
           'TNFAIP2','ICAM1','HAMP','ITPKC','JAK3','GMFG','HCK','CEBPB']
PREFIX = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174367/suppl/'
SPEC = {'files': [
    dict(id='metadata', filename='GSE174367_snRNA-seq_cell_meta.csv.gz',
         url=PREFIX+'GSE174367_snRNA-seq_cell_meta.csv.gz'),
    dict(id='counts', filename='GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5',
         url=PREFIX+'GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5')]}


def decode(values):
    return np.array([x.decode('utf-8') if isinstance(x, bytes) else str(x) for x in values])


def prepare(matrix_path, metadata_path, out):
    import h5py
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    cells = pd.read_csv(metadata_path)
    required = ['Barcode','SampleID','Cell.Type','Age','Sex','Batch','RIN','PMI','Tangle.Stage']
    if not set(required).issubset(cells):
        raise ValueError('Missing public cell metadata fields')
    if cells.Barcode.isna().any() or cells.Barcode.duplicated().any():
        raise ValueError('Missing or duplicate cell barcodes')
    fields = ['SampleID','Age','Sex','Batch','RIN','PMI','Tangle.Stage']
    donors = cells[fields].drop_duplicates().sort_values('SampleID').reset_index(drop=True)
    if donors.SampleID.duplicated().any():
        raise ValueError('Inconsistent donor-level covariates')
    donors['braak'] = pd.to_numeric(donors['Tangle.Stage'].str.extract(r'^Stage ([0-6])$')[0], errors='raise')
    if donors[['SampleID','Age','Sex','Batch','RIN','braak']].isna().any().any():
        raise ValueError('Incomplete primary donor covariates')
    if not donors.Sex.isin(['F','M']).all():
        raise ValueError('Unrecognized sex labels')
    donors['sex_male'] = donors.Sex.eq('M').astype(int)
    donors = donors.rename(columns={'SampleID':'donor_id','Age':'age','Batch':'batch','RIN':'rin','PMI':'pmi'})
    sample_index = {x:i for i,x in enumerate(donors.donor_id)}
    type_index = {x:i for i,x in enumerate(CELL_TYPES)}
    nc, nd = len(CELL_TYPES), len(donors)
    with h5py.File(matrix_path, 'r') as handle:
        matrix = handle['matrix']
        shape = tuple(matrix['shape'][:])
        barcodes = decode(matrix['barcodes'][:])
        ids, symbols = decode(matrix['features/id'][:]), decode(matrix['features/name'][:])
        if len(shape)!=2 or shape != (len(ids), len(barcodes)):
            raise ValueError('Unexpected count matrix shape')
        if len(set(barcodes)) != len(barcodes) or len(set(ids)) != len(ids):
            raise ValueError('Duplicate matrix feature IDs or barcodes')
        missing = set(cells.Barcode)-set(barcodes)
        if missing:
            raise ValueError(f'{len(missing)} metadata cells are absent from the count matrix')
        aligned = cells.set_index('Barcode').reindex(barcodes)
        group = np.full(len(barcodes), -1, dtype=int)
        for i,(sample,kind) in enumerate(zip(aligned.SampleID, aligned['Cell.Type'])):
            if sample in sample_index and kind in type_index:
                group[i] = type_index[kind]*nd+sample_index[sample]
        totals = np.zeros((len(ids), nc*nd), dtype=np.int64)
        pointers = matrix['indptr'][:]
        if len(pointers)!=len(barcodes)+1 or pointers[0]!=0 or np.any(np.diff(pointers)<0):
            raise ValueError('Invalid CSC pointers')
        if int(pointers[-1]) != len(matrix['data']) or len(matrix['indices']) != len(matrix['data']):
            raise ValueError('CSC array lengths differ')
        for start in range(0,len(barcodes),2000):
            end = min(start+2000,len(barcodes)); left,right = int(pointers[start]),int(pointers[end])
            values = matrix['data'][left:right]
            indices = matrix['indices'][left:right].astype(np.int64)
            if not np.isfinite(values).all() or np.any(values<0) or np.any(values!=np.floor(values)):
                raise ValueError('Expected nonnegative integer raw counts, not normalized expression')
            if np.any(indices<0) or np.any(indices>=len(ids)):
                raise ValueError('Invalid CSC row indices')
            groups = np.repeat(group[start:end], np.diff(pointers[start:end+1]))
            keep = groups>=0
            np.add.at(totals, (indices[keep],groups[keep]), values[keep].astype(np.int64))
            print(f'GSE174367 pseudobulk: {end}/{len(barcodes)} matrix columns',flush=True)
    count_table = pd.crosstab(cells.SampleID,cells['Cell.Type']).reindex(index=donors.donor_id,columns=CELL_TYPES,fill_value=0)
    count_table.index.name='donor_id'
    count_table.to_csv(out/'donor_cell_counts.csv')
    donors.to_csv(out/'donors.csv',index=False)
    pd.DataFrame({'gene_id':ids,'gene_symbol':symbols}).to_csv(out/'genes.csv',index=False)
    for i,kind in enumerate(CELL_TYPES):
        matrix = pd.DataFrame(totals[:,i*nd:(i+1)*nd],index=ids,columns=donors.donor_id)
        matrix.index.name='gene_id'; matrix.to_csv(out/('counts_'+kind+'.csv.gz'))
    report = dict(donors=len(donors),metadata_nuclei=len(cells),matrix_nuclei=len(barcodes),
                  matrix_cells_not_in_metadata=int(aligned.SampleID.isna().sum()),
                  cells_aggregated=int((group>=0).sum()),features=len(ids),
                  stage_counts=donors.braak.value_counts().sort_index().to_dict(),
                  cell_types=CELL_TYPES,primary_cell_type='MG',minimum_cells=20,
                  scope='Independent prefrontal single-nucleus context; V-VI versus I-II, not exact Braak-0 replication')
    (out/'aggregation.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    return report
