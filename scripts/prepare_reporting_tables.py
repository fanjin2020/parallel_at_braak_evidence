"""Build descriptive and display tables from this run without changing fitted results."""
import argparse
import json
from pathlib import Path
import pandas as pd
from source_io import sha256

def run(root,package_dir=None):
    out=root/'results/descriptive_tables';out.mkdir(parents=True,exist_ok=True)
    m=pd.read_csv(root/'data/public/gse131617/sample_manifest.csv').drop_duplicates('donor_key')
    rows=[]
    for stage,g in m.groupby('braak_stage',sort=False):
        rows.append(dict(braak_stage=stage,donors=len(g),male=int(g.sex_male.eq(1).sum()),
            female=int(g.sex_male.eq(0).sum()),apoe4_carriers=int(g.apoe4_carrier.eq(1).sum()),
            age_median=g.age_at_death.median(),pmi_median=g.pmi_hours.median()))
    pd.DataFrame(rows).to_csv(out/'S1_donor_characteristics.csv',index=False)
    a=pd.read_csv(root/'data/private/adni/adni_pet_aligned_master.csv')
    rows=[]
    for v in ['age_at_tau','sex_male','education_years','apoe4_carrier','tau_temporal_suvr','hippocampus_icv','cdrsb',
              'hippocampus_icv_days_from_tau','cdrsb_days_from_tau']:
        x=a.loc[a.at_stage.eq('A-T+'),v].dropna()
        rows.append(dict(variable=v,n=len(x),mean=x.mean(),sd=x.std(),median=x.median(),minimum=x.min(),maximum=x.max()))
    pd.DataFrame(rows).to_csv(out/'S3_ATminus_Tplus_description.csv',index=False)
    go=pd.read_csv(root/'revision/pathways/current_GO_significant_grouped.csv')
    go=go.loc[go.display_representative.astype(str).str.lower().eq('true')].sort_values(['q','p','go_id']).head(12)
    go.to_csv(root/'revision/pathways/GO_display_12.csv',index=False)
    qc_go_display(root)
    adni_provenance(root,package_dir)
    definitions(root)

def qc_go_display(root):
    out=root/'reporting/QC_GO_display';out.mkdir(parents=True,exist_ok=True)
    tests=pd.read_csv(root/'targeted/qc_joint/GO_all_models.csv')
    terms=pd.read_csv(root/'revision/go_annotation/go_terms.csv')
    members=pd.read_csv(root/'revision/go_annotation/go_gene_membership.csv',dtype=str)
    sets={key:set(g.entrez_id) for key,g in members.groupby('go_id',sort=False)}
    d=tests.loc[tests.model.eq('baseline')].copy()
    if d.go_id.duplicated().any() or set(d.go_id)!=set(terms.go_id):
        raise ValueError('QC baseline GO family differs from the recorded eligible family')
    d=d.merge(terms[['go_id','name']],on='go_id',validate='one_to_one')
    d=d.loc[d.q.lt(.05)].sort_values(['q','p','go_id'],kind='stable').reset_index(drop=True)
    representatives=[];groups=[];flags=[]
    for row in d.itertuples():
        group=None
        for i,rep in enumerate(representatives,1):
            a,b=sets[row.go_id],sets[rep]
            if len(a&b)/len(a|b)>=.5:
                group=i;break
        is_rep=group is None
        if is_rep:
            representatives.append(row.go_id);group=len(representatives)
        groups.append(group);flags.append(is_rep)
    d['overlap_group']=groups;d['display_representative']=flags
    d.to_csv(out/'QC_baseline_GO_grouped.csv',index=False)
    d.loc[d.display_representative].head(12).to_csv(out/'QC_baseline_GO_display_12.csv',index=False)

def adni_provenance(root,package_dir=None):
    raw=root/'data/private/adni/raw'
    info={}
    if package_dir:
        description=Path(package_dir)/'DESCRIPTION'
        for line in description.read_text(encoding='utf-8').splitlines():
            if ':' in line and not line.startswith((' ','\t')):
                key,value=line.split(':',1)
                if key in {'Package','Version','Date','Date/Publication','Built','RemoteSha','GithubSHA1'}:
                    info[key]=value.strip()
        info['DESCRIPTION_sha256']=sha256(description)
    date=pd.read_csv(raw/'DATA_DOWNLOADED_DATE.csv',dtype=str)
    pet=[]
    for name in ['UCBERKELEY_AMY_6MM.csv','UCBERKELEY_TAU_6MM.csv','UCBERKELEY_TAUPVC_6MM.csv']:
        cols=pd.read_csv(raw/name,nrows=0).columns
        selected=[c for c in ['update_stamp','PROCESSDATE','IMAGE_RESOLUTION','TRACER'] if c in cols]
        d=pd.read_csv(raw/name,usecols=selected,dtype=str)
        entry=dict(table=name,sha256=sha256(raw/name),rows=len(d))
        for col in selected:
            values=sorted(d[col].dropna().unique())
            entry[col]=values if col in {'TRACER','IMAGE_RESOLUTION'} else dict(min=min(values,default=None),max=max(values,default=None))
        pet.append(entry)
    config=json.loads((root/'run_settings.json').read_text(encoding='utf-8'))
    payload=dict(record_id='SD242',package=info or {'status':'Raw CSV input; no package DESCRIPTION supplied'},
                 download_date_records=date.to_dict('records'),PET_tables=pet,
                 bootstrap=config.get('bootstrap'),enhancement_bootstrap=config.get('enhancement_bootstrap'),
                 scope='Local data snapshot identifiers; update_stamp is not a PET processing-release identifier.')
    dest=root/'reporting/ADNI_snapshot_provenance.json';dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(payload,indent=2),encoding='utf-8')

def source_definitions():
    return pd.DataFrame([
        dict(source='ADNI',access='Authorized local data only',role='Clinical PET-aligned analyses',version_source='local_input_inventory.json; supplied package DESCRIPTION'),
        dict(source='GSE131617',access='Public GEO and original published supplement',role='Discovery cortical expression and CEL QC',version_source='download receipts and annotation package versions'),
        dict(source='GSE44770',access='Public GEO',role='Related AD/control external comparison',version_source='preparation provenance'),
        dict(source='GSE106241',access='Public GEO',role='Independent temporal-cortex Braak comparison',version_source='preparation provenance'),
        dict(source='GSE174367',access='Public GEO',role='Donor-level single-nucleus pseudobulk and OPC-INH interaction',version_source='GSE174367_source_provenance.json; model status and fixed testing families'),
        dict(source='GO',access='Public GOA, GO and NCBI',role='Genome-wide pathway testing',version_source='revision/go_annotation/annotation_versions.json'),
        dict(source='CA and NG references',access='Published public reference workbook',role='Estimated cell-composition sensitivity',version_source='reference_provenance.json'),
    ])

def definitions(root):
    defs=root/'definitions';defs.mkdir(exist_ok=True)
    source_definitions().to_csv(defs/'data_sources.csv',index=False)
    pd.DataFrame([
        dict(analysis='Confirmatory A/T',definition='A+T- and A+T+ versus A-T-; two outcomes; four-test BH family'),
        dict(analysis='Exploratory clinical',definition='Two direct contrasts and two within-A+ slopes; separate four-test BH family'),
        dict(analysis='Sensitivity analyses',definition='Covariate, missingness, technical and external-cohort analyses'),
        dict(analysis='CEL QC',definition='Technical exclusions recorded before refitting; provided-matrix and reprocessed RMA comparisons'),
        dict(analysis='Missing data limits',definition='No fabricated individual RIN, experimental batch, cause of death, or unverified survival data'),
    ]).to_csv(defs/'analysis_roles.csv',index=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run-dir',type=Path,required=True);p.add_argument('--package-dir',type=Path)
    a=p.parse_args();run(a.run_dir,a.package_dir)
