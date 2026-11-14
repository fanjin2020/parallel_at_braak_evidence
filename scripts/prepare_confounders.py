"""Prepare dated ADNI covariates from exported ADNIMERGE2 tables."""
import argparse, json, re
from pathlib import Path
import numpy as np
import pandas as pd

HISTORY=['MHPSYCH','MH2NEURL','MH4CARD','MH5RESP','MH6HEPAT','MH7DERM','MH8MUSCL','MH9ENDO','MH10GAST','MH11HEMA','MH12RENA','MH17MALI']
KEYMED={'Aricept','Cognex','Exelon','Namenda','Razadyne','Anti-depressant medication','Other behavioral medication','None of the above'}
DRUGS={
 'antidementia':r'\b(?:aricept|donepezil|donepezi|exelon|rivastigmine|razadyne|reminyl|galantamine|galanthamine|namenda|memantine|namzaric|cognex|tacrine)\b',
 'antidepressant':r'\b(?:sertraline|zoloft|fluoxetine|prozac|citalopram|celexa|escitalopram|lexapro|paroxetine|paxil|fluvoxamine|luvox|venlafaxine|effexor|desvenlafaxine|pristiq|duloxetine|cymbalta|bupropion|wellbutrin|mirtazapine|remeron|trazodone|desyrel|amitriptyline|elavil|nortriptyline|pamelor|imipramine|tofranil|clomipramine|anafranil|doxepin|sinequan|desipramine|norpramin|phenelzine|nardil|tranylcypromine|parnate|vortioxetine|trintellix|vilazodone|viibryd)\b'}

def yesno(s):
    return s.astype('string').str.strip().str.lower().map({'yes':1.,'no':0.,'1':1.,'0':0.,'1.0':1.,'0.0':0.})

def keymed(value):
    if pd.isna(value):return np.nan,np.nan
    parts=set(str(value).split('|'))
    if not parts.issubset(KEYMED) or ('None of the above' in parts and len(parts)>1):return np.nan,np.nan
    return float(bool(parts&{'Aricept','Cognex','Exelon','Namenda','Razadyne'})),float('Anti-depressant medication' in parts)

def load(raw,name,checks):
    p=raw/(name+'.csv');d=pd.read_csv(p,low_memory=False)
    d['rid']=pd.to_numeric(d.RID,errors='raise').astype('Int64').astype(str)
    d['assessment_date']=pd.to_datetime(d.VISDATE,errors='coerce')
    bad=d.assessment_date.isna()
    if 'HAS_QC_ERROR' in d:
        bad|=d.HAS_QC_ERROR.astype(str).str.lower().isin(['1','1.0','has qc error'])
    checks.append(dict(table=name,total_rows=len(d),excluded_missing_date_or_qc=int(bad.sum())))
    return d.loc[~bad].copy()

def run(raw,out,health_file=None):
    out.mkdir(parents=True,exist_ok=False);checks=[];pieces=[]
    gd=load(raw,'GDSCALE',checks);gd['depression_score']=pd.to_numeric(gd.GDTOTAL,errors='coerce').where(lambda s:s.between(0,15))
    pieces.append(gd[['rid','assessment_date','depression_score']])
    b=load(raw,'BACKMEDS',checks);vals=b.KEYMED.map(keymed)
    b['antidementia_medication']=[v[0] for v in vals];b['antidepressant_medication']=[v[1] for v in vals]
    pieces.append(b[['rid','assessment_date','antidementia_medication','antidepressant_medication']])
    mh=load(raw,'MEDHIST',checks)
    for c in HISTORY:mh[c]=yesno(mh[c])
    mh['vascular_history']=mh.MH4CARD
    mh['comorbidity_count']=mh[HISTORY].sum(axis=1,min_count=len(HISTORY))
    mh['prior_cardiovascular_history']=mh.MH4CARD
    mh['prior_medical_domain_count']=mh.comorbidity_count
    mh['prior_other_neurologic_history']=mh.MH2NEURL
    pieces.append(mh[['rid','assessment_date','vascular_history','comorbidity_count','prior_cardiovascular_history','prior_medical_domain_count','prior_other_neurologic_history']])
    h=load(raw,'MODHACH',checks);h['ischemic_score']=pd.to_numeric(h.HMSCORE,errors='coerce').where(lambda s:s.between(0,12))
    h['stroke_history']=h.HMSTROKE.astype('string').str.lower().map({'absent':0.,'present - 2 points':1.,'0':0.,'2':1.})
    pieces.append(h[['rid','assessment_date','ischemic_score','stroke_history']])
    health=health_file or raw/'INITHEALTH.csv'
    if health.exists():
        ih=load(health.parent,health.stem,checks)
        ih['present']=yesno(ih.IHPRESENT)
        categories={'Psychiatric','Neurologic (other than Cognitive Disorder)','Cardiovascular','Respiratory','Hepatic','Dermatologic-Connective Tissue','Musculoskeletal','Endocrine-Metabolic','Gastrointestinal','Hematopoietic-Lymphatic','Renal-Genitourinary','Malignancy'}
        records=[]
        for (rid,day),g in ih.groupby(['rid','assessment_date']):
            row=dict(rid=rid,assessment_date=day)
            if g.present.isna().any() or g.IHSYMPTOM.isna().any():
                row.update(initial_recorded_cv=np.nan,initial_recorded_neuro=np.nan,initial_recorded_domain_count=np.nan)
            else:
                positive=set(g.loc[g.present.eq(1),'IHSYMPTOM'])
                row.update(initial_recorded_cv=float('Cardiovascular' in positive),initial_recorded_neuro=float('Neurologic (other than Cognitive Disorder)' in positive),initial_recorded_domain_count=len(positive&categories))
            records.append(row)
        pieces.append(pd.DataFrame(records))
    table=pd.concat(pieces,ignore_index=True)
    table.to_csv(out/'dated_confounders.csv',index=False)
    # Record exact drug-string matches separately; absence of a matching record is not no treatment.
    meds=load(raw,'RECCMEDS',checks);names=(meds.CMMED.fillna('')+' '+meds.CMMEDO.fillna('')).str.strip()
    pd.DataFrame(checks).to_csv(out/'preparation_counts.csv',index=False)
    mapping=[]
    for group,pattern in DRUGS.items():
        for name,count in names[names.str.contains(pattern,case=False,regex=True)].value_counts().items():mapping.append(dict(drug_class=group,recorded_name=name,n_records=int(count)))
    pd.DataFrame(mapping).to_csv(out/'recorded_drug_name_matches.csv',index=False)
    meds['mapped_antidementia']=names.str.contains(DRUGS['antidementia'],case=False,regex=True)
    meds['mapped_antidepressant']=names.str.contains(DRUGS['antidepressant'],case=False,regex=True)
    meds.loc[meds.mapped_antidementia|meds.mapped_antidepressant,['rid','assessment_date','CMBGN','CMEND','CMCONT','mapped_antidementia','mapped_antidepressant']].to_csv(out/'restricted_medication_crosscheck.csv',index=False)
    spec=dict(sources='Existing authorized ADNIMERGE2 export; see source_inventory.csv for source version and files',
        assessment_dates='ADNIMERGE2 VISDATE: assessment-specific date where supplied, otherwise REGISTRY EXAMDATE; never USERDATE',
        depression='GDS total 0-15, latest pre-PET observation within 180 days',
        medication='BACKMEDS explicit current-use checklist at latest pre-PET visit within 180 days; not verified adherence on PET day. RECCMEDS name matches are corroborative records, not used to turn missing data into non-use.',
        cardiovascular='MEDHIST MH4CARD is cardiovascular history, not a cerebrovascular imaging measure',
        comorbidity='Count of 12 recorded medical-history domains, not Charlson or count of distinct diseases',
        counted_domains=HISTORY,prior_history='Variables prefixed prior_ use the most recent available assessment before PET without 180-day restriction; historical status only, cannot rule out subsequent disease',
        initial_health='INITHEALTH documented conditions present in the three months before screening. The 12-domain count and flags reflect recorded conditions, not an adjudicated comorbidity index; zero means no qualifying condition recorded, not proof of disease absence. Requires nonmissing category and present/absent status for all records at that assessment.',
        ischemic_score='MODHACH total 0-12; stroke is the explicit history-of-stroke item, not MRI vascular burden',
        restrictions='Participant-level files must remain local and must not be included in the public code archive')
    (out/'coding_specification.json').write_text(json.dumps(spec,indent=2),encoding='utf-8')
    print('Prepared',len(table),'dated covariate records',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--raw-dir',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--health-file',type=Path)
    a=p.parse_args();run(a.raw_dir,a.output_dir,a.health_file)
