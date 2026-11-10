import pandas as pd
IDENTIFIERS={'rid','ptid','subject','subject_id','participant_id','patient_id','sample_id','donor_key','loniuid'}
PUBLIC_GEO_TABLES={'cel_qc/array_quality.csv','cel_qc/technical_exclusions.csv'}
PATH_PATTERN=r'(?<![A-Za-z0-9/:])(?:[A-Za-z]:[\\/]|/(?:home|Users|tmp|mnt)/)'

def check_public(df,source,public_manifest=None):
    """Permit GEO identifiers only in named public QC tables with verified linkage."""
    blocked=IDENTIFIERS.intersection(str(c).strip().lower() for c in df.columns)
    if source in PUBLIC_GEO_TABLES and public_manifest is not None:
        ids=df['sample_id'].astype(str)
        if not ids.str.fullmatch(r'GSM[0-9]+').all() or not ids.isin(public_manifest.sample_id).all():
            raise ValueError('Unverified public GEO identifier in '+source)
        if 'donor_key' in df:
            lookup=public_manifest.set_index('sample_id').donor_key
            if not ids.map(lookup).astype(str).eq(df.donor_key.astype(str)).all():
                raise ValueError('Public donor/GSM linkage differs in '+source)
        blocked-= {'sample_id','donor_key'}
    if blocked:raise ValueError(f'Individual-level key in export {source}: {sorted(blocked)}')
    for col in df:
        if df[col].dropna().astype(str).str.contains(PATH_PATTERN,regex=True).any():
            raise ValueError(f'Local path in proposed public source: {source}, {col}')


def read_csv_exact(path):
    # Strings preserve every stored decimal and missing-value token during grouping.
    return pd.read_csv(path,dtype=str,keep_default_na=False)


def csv_bytes(df):
    return df.to_csv(index=False,lineterminator='\n').encode('utf-8')

