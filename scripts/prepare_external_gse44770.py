"""Prepare public HBTRC prefrontal AD/control data for a related-cohort check."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import numpy as np
import pandas as pd

MATRIX='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE44nnn/GSE44770/matrix/GSE44770_series_matrix.txt.gz'
PLATFORM='https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL4nnn/GPL4372/soft/GPL4372_family.soft.gz'


def series_matrix(path):
    metadata={};characteristics={}
    with gzip.open(path,'rt',encoding='utf-8') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'):
                table=pd.read_csv(f,sep='\t',comment='!',dtype={'ID_REF':str},na_values=['null','NULL'])
                break
            values=next(csv.reader([line],delimiter='\t'))
            if not values:continue
            if values[0]=='!Sample_characteristics_ch2':
                pairs=[x.split(': ',1) for x in values[1:]]
                if any(len(x)!=2 for x in pairs) or len({x[0] for x in pairs})!=1:raise ValueError('Inconsistent sample metadata keys')
                key=pairs[0][0]
                if key in characteristics:raise ValueError('Repeated sample metadata key')
                characteristics[key]=[x[1] for x in pairs]
            elif values[0].startswith('!Sample_'):metadata[values[0][8:]]=values[1:]
        else:raise ValueError('GEO expression table not found')
    return table,metadata,characteristics


def platform_table(path):
    with gzip.open(path,'rt',encoding='utf-8') as f:
        for line in f:
            if line.startswith('!platform_table_begin'):
                header=next(csv.reader([next(f)],delimiter='\t'))
                rows=[]
                for line in f:
                    if line.startswith('!platform_table_end'):return pd.DataFrame(rows,columns=header)
                    values=next(csv.reader([line],delimiter='\t'))
                    if len(values)!=len(header):raise ValueError('Malformed platform row')
                    rows.append(values)
    raise ValueError('Platform table missing')


def run(out,matrix_file=None,platform_file=None):
    out.mkdir(parents=True,exist_ok=True)
    if (out/'manifest.csv').exists():raise ValueError('Prepared output exists; choose a new directory')
    matrix=out/'GSE44770_series_matrix.txt.gz'
    if not matrix_file or not platform_file:raise ValueError('Supply cached GSE44770 matrix and GPL4372 annotation.')
    shutil.copy2(matrix_file,matrix)
    platform=out/'GPL4372_platform_table.soft.gz';shutil.copy2(platform_file,platform)
    table,meta,chars=series_matrix(matrix)
    required={'disease','age','gender','pmi','ph','rin','batch'}
    if not required.issubset(chars):raise ValueError('Required external confounders missing: '+str(required-set(chars)))
    if set(meta.get('platform_id',[]))!={'GPL4372'}:raise ValueError('Unexpected expression platform')
    if not set(meta.get('label_ch2',[]))=={'Cy5'} or not set(meta.get('label_ch1',[]))=={'Cy3'}:raise ValueError('Unexpected channel orientation')
    m=pd.DataFrame(chars)
    m['sample_id']=meta['geo_accession'];m['donor_key']=meta['title']
    if m.sample_id.duplicated().any() or m.donor_key.duplicated().any():raise ValueError('PFC donor/sample keys are not unique')
    m['group']=m.disease.map({'N':'reference','A':'case'})
    selected=m.group.notna();excluded_other_disease=int((~selected).sum());m=m.loc[selected].copy()
    m['age_at_death']=pd.to_numeric(m.age,errors='coerce');m['sex_male']=m.gender.map({'M':1,'F':0})
    for new,old in [('pmi_hours','pmi'),('tissue_ph','ph'),('rin','rin')]:m[new]=pd.to_numeric(m[old],errors='coerce')
    m['batch']=m.batch.replace({'NA':np.nan,'':np.nan}).map(lambda x:'batch_'+str(x) if pd.notna(x) else np.nan)
    columns=['sample_id','donor_key','group','age_at_death','sex_male','pmi_hours','tissue_ph','rin','batch']
    m=m[columns]
    if not set(m.sample_id).issubset(table):raise ValueError('External matrix/sample mismatch')
    annotation=platform_table(platform)
    if not {'ID','EntrezGeneID'}.issubset(annotation):raise ValueError('Unexpected GPL4372 columns')
    annotation=annotation.loc[annotation.EntrezGeneID.fillna('').str.fullmatch(r'\d+')].copy()
    if annotation.ID.duplicated().any():raise ValueError('Duplicate platform probe mapping')
    if table.ID_REF.duplicated().any():raise ValueError('Duplicate GEO matrix probe')
    expression=table.set_index('ID_REF')[list(m.sample_id)].apply(pd.to_numeric,errors='coerce')
    mapping=annotation.set_index('ID').EntrezGeneID.reindex(expression.index)
    good=mapping.notna()&np.isfinite(expression).all(axis=1)
    # GEO sample VALUE: normalized log10(Cy5/Cy3), test/reference. Do not log again.
    expression=expression.loc[good]*np.log2(10)
    expression['entrez_id']=mapping.loc[good]
    expression=expression.groupby('entrez_id',sort=True).median()
    if len(expression)<1000:raise ValueError('Too few uniquely mapped complete genes')
    expression.to_csv(out/'expression_log2.csv',index=True)
    m.to_csv(out/'manifest.csv',index=False)
    source='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44770'
    independence='GSE44770 is the HBTRC PFC donor collection; discovery GSE131617 comprises Japanese brain donors. Separate accession/study donor collections as described in GEO; no participant identity re-identification attempted.'
    config={'manifest':'manifest.csv','expression':'expression_log2.csv','covariates':'age_at_death,sex_male,pmi_hours,tissue_ph,rin,batch','source':source,'independence-evidence':independence,'validation-kind':'related_AD_vs_control'}
    (out/'external_config.json').write_text(json.dumps(config,indent=2),encoding='utf-8')
    hashes={}
    for path in [matrix,platform]:
        with path.open('rb') as f:hashes[path.name]=hashlib.file_digest(f,'sha256').hexdigest()
    info=dict(matrix_url=MATRIX,platform_url=PLATFORM,platform_download='Stream only through platform_table_end; exclude unrelated sample records in the large platform family archive',sample_value_definition_source='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1090501',
              scale='Submitted normalized log10(Cy5/Cy3) test/reference multiplied by log2(10); no re-normalization',
              probe_mapping='Only single Entrez mappings; median across all retained probes per gene, independent of group effects',
              excluded_probe_rows=int((~good).sum()),retained_genes=len(expression),samples=len(m),excluded_other_disease=excluded_other_disease,
              group_counts=m.group.value_counts().to_dict(),missing_covariates=m.isna().sum().to_dict(),sha256=hashes,
              comparison='AD versus control in PFC; not Braak V-VI versus 0, not cell-resolved validation',
              data_processing=meta.get('data_processing',[''])[0])
    (out/'preparation_specification.json').write_text(json.dumps(info,indent=2),encoding='utf-8')
    print('Prepared external related-cohort inputs:',len(m),'samples;',len(expression),'genes',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--matrix-file',type=Path);p.add_argument('--platform-file',type=Path)
    a=p.parse_args();run(a.output_dir.resolve(),a.matrix_file,a.platform_file)
