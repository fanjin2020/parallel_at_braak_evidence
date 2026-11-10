"""Prepare an independently collected temporal-cortex Braak series from GEO."""
import argparse
import gzip
import io
from pathlib import Path
import numpy as np
import pandas as pd
from followup_io import soft_samples, series_url, write_json, sha256

ACCESSION = 'GSE106241'
GENE_INFO = 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz'


def platform_table(path):
    lines, active = [], False
    with gzip.open(path, 'rt', encoding='utf-8') as f:
        for line in f:
            if line.startswith('!platform_table_begin'):
                active = True
            elif line.startswith('!platform_table_end'):
                break
            elif active:
                lines.append(line)
    return pd.read_csv(io.StringIO(''.join(lines)), sep='\t', dtype=str, keep_default_na=False)


def prepare(cache, out, soft_file=None, expression_file=None, gene_info=None):
    out.mkdir(parents=True, exist_ok=True)
    if not all((soft_file,expression_file,gene_info)):raise ValueError('Local SOFT, matrix and gene annotation files are required.')
    soft,matrix,info=map(Path,(soft_file,expression_file,gene_info))
    raw = pd.DataFrame(soft_samples(soft))
    raw.to_csv(out/'GEO_sample_metadata.csv', index=False)
    apoe = 'apoe (additive) (0=0 copies of apoee4, 1=1 copy of apoee4, 2=two copies of apoee4)'
    fields = ['sample_id','title','braak stage','brain region','age','sex (1= female, 2= male)',apoe,'rna quality']
    if not set(fields).issubset(raw):
        raise ValueError('GEO metadata schema changed; review field mapping')
    m = pd.DataFrame(dict(sample_id=raw.sample_id, donor_key=raw.title,
        brain_region=raw['brain region'], braak_stage=pd.to_numeric(raw['braak stage'], errors='raise'),
        age_at_death=pd.to_numeric(raw.age, errors='raise'),
        sex_male=raw['sex (1= female, 2= male)'].map({'1':0,'2':1}),
        apoe4_count=pd.to_numeric(raw[apoe], errors='coerce'),
        rna_quality_category=raw['rna quality'].replace({'NA':np.nan,'':np.nan})))
    if m.sample_id.duplicated().any() or m.donor_key.duplicated().any() or m.sex_male.isna().any():
        raise ValueError('Donor/sample/sex metadata failed validation')
    if not m.braak_stage.isin(range(7)).all() or m.brain_region.nunique()!=1:
        raise ValueError('Unexpected Braak stages or tissue regions')
    m['braak_group'] = m.braak_stage.map({0:'0',1:'I-II',2:'I-II',3:'III-IV',4:'III-IV',5:'V-VI',6:'V-VI'})
    m.to_csv(out/'manifest.csv', index=False)
    values = pd.read_csv(matrix, index_col=0)
    if values.index.duplicated().any() or set(values.columns)!=set(m.donor_key):
        raise ValueError('Expression matrix does not match donor titles one-to-one')
    values=values.loc[:,m.donor_key]
    if not np.isfinite(values.to_numpy()).all() or not (values.to_numpy()>0).all():
        raise ValueError('Expected finite positive linear-scale normalized intensities')
    annotation = platform_table(soft)
    if not {'ID','GeneSymbol'}.issubset(annotation) or annotation.ID.duplicated().any():
        raise ValueError('Unexpected GPL24170 annotation')
    genes = pd.read_csv(info, sep='\t', usecols=['GeneID','Symbol'], dtype=str)
    unique_symbols = genes.groupby('Symbol').GeneID.nunique()
    genes=genes[genes.Symbol.isin(unique_symbols[unique_symbols.eq(1)].index)].drop_duplicates('Symbol')
    mapping=annotation[['ID','GeneSymbol']].merge(genes, left_on='GeneSymbol', right_on='Symbol', how='left', validate='many_to_one')
    mapping['available_in_expression']=mapping.ID.isin(values.index)
    mapping['mapping_status']=np.where(mapping.GeneID.notna(),'exact_unique_official_symbol','unmapped_or_ambiguous')
    mapping.to_csv(out/'probe_mapping_audit.csv',index=False)
    mapped=mapping.loc[mapping.GeneID.notna() & mapping.available_in_expression].set_index('ID')
    # Same mapping and outcome-blind aggregation for every donor and every tested gene.
    log=np.log2(values.loc[mapped.index])
    log['entrez_id']=mapped.GeneID
    gene=log.groupby('entrez_id',sort=True).median()
    gene.columns=m.sample_id
    gene.reset_index().to_csv(out/'expression_gene.csv.gz',index=False)
    m.groupby('braak_group',sort=False).agg(n_donors=('donor_key','nunique'),age_min=('age_at_death','min'),
        age_max=('age_at_death','max'),n_male=('sex_male','sum')).reset_index().to_csv(out/'cohort_coverage.csv',index=False)
    write_json(out/'provenance.json',dict(accession=ACCESSION,source='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106241',
        tissue='Inferior temporal cortex; related to discovery TC, not identical anatomical sampling established',
        independence='GEO describes University of Eastern Finland temporal-cortex cohort (PMID 30557660); different source from GSE131617. Cross-study donor identifiers cannot be globally linked.',
        contrast='Braak V-VI versus 0, with intermediate stages retained in model',
        expression='Published background-corrected, quantile-normalized positive probe intensities; log2 then median over uniquely mapped probes per gene. No case-control-based probe selection.',
        rna_quality='Published categorical field, NOT relabeled as numeric RIN',
        missing_covariates=['PMI','verified experimental batch','cause_of_death'],
        n_samples=len(m),n_genes=len(gene),input_hashes={p.name:sha256(p) for p in [soft,matrix,info]}))
    print(f'GSE106241 prepared: {len(m)} donors, {len(gene)} mapped genes',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cache-dir',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--soft-file',type=Path);p.add_argument('--expression-file',type=Path);p.add_argument('--gene-info',type=Path)
    a=p.parse_args();prepare(a.cache_dir,a.output_dir,a.soft_file,a.expression_file,a.gene_info)
