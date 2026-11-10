"""Date-aligned clinical-confounder sensitivity with explicit missing-input reporting."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import pet_staging_covariates as adni
import cross_sectional as cross
import longitudinal_model_comparison as longitudinal
import longitudinal_sensitivity as sensitivity
import analyze_tau_incremental as inc
import longitudinal as core

FIELDS=['depression_score','vascular_history','comorbidity_count','antidementia_medication','antidepressant_medication','ischemic_score','stroke_history','prior_cardiovascular_history','prior_medical_domain_count','prior_other_neurologic_history','initial_recorded_cv','initial_recorded_neuro','initial_recorded_domain_count']
REQUESTS=[('depression_score','GDSCALE.csv','RID, EXAMDATE, GDTOTAL','Nearest dated GDS before tau PET, within 180 days; valid 0-15'),
          ('vascular_history','MEDHIST.csv and dated diagnosis/history updates','RID, visit identifier, clinical assessment date, cardiovascular/cerebrovascular items','Verify coded definitions and distinguish cardiovascular from cerebrovascular disease; do not substitute data-entry USERDATE'),
          ('comorbidity_count','MEDHIST.csv and dated diagnosis/history updates','RID, clinical assessment date, clinically defined condition indicators','Predefine counted conditions; an absent file or unanswered item is not no disease'),
          ('antidementia_medication','RECCMEDS.csv','RID, CMMED, CMBGN, CMEND and ongoing-use fields','Determine active treatment at tau PET using validated drug mapping and dates; unknown dates remain missing'),
          ('antidepressant_medication','RECCMEDS.csv','RID, CMMED, CMBGN, CMEND and ongoing-use fields','Same temporal and drug-coding checks; do not identify medication solely from a vague free-text indication')]


def align(table,master):
    required={'rid','assessment_date'}
    if not required.issubset(table):raise ValueError('Confounder file needs rid and assessment_date')
    available=[v for v in FIELDS if v in table]
    if not available:raise ValueError('No supported confounder variables supplied')
    t=table.copy();t['rid']=t.rid.map(core.rid);t['assessment_date']=inc.dates(t.assessment_date)
    if t.assessment_date.isna().any():raise ValueError('Missing/unparseable assessment dates; do not use USERDATE as assessment date')
    for v in available:
        raw=t[v];t[v]=pd.to_numeric(raw,errors='coerce')
        if (raw.notna() & t[v].isna()).any():raise ValueError('Unmapped nonnumeric code in '+v)
        if (t[v].dropna()<0).any():raise ValueError('Negative/unknown codes must be missing: '+v)
        if v.endswith('medication') or v.endswith('history') or v in ['initial_recorded_cv','initial_recorded_neuro']:
            if not set(t[v].dropna()).issubset({0,1}):raise ValueError('Binary 0/1 required: '+v)
        if v=='depression_score' and t[v].gt(15).any():raise ValueError('GDS score must be 0-15')
    rows=[]
    for p in master:
        g=t.loc[t.rid.eq(p['rid'])].copy();g['gap']=(pd.Timestamp(p['anchor'])-g.assessment_date).dt.days
        g=g.loc[g.gap.ge(0)]
        row={'rid':p['rid']}
        for v in available:
            h=g.loc[g[v].notna() & (g.gap.ge(0) if v.startswith('prior_') else g.gap.le(180))].sort_values(['gap','assessment_date'])
            if len(h):
                same=h.loc[h.gap.eq(h.iloc[0].gap),v].unique()
                row[v]=same[0] if len(same)==1 else np.nan
                row[v+'_lag_days']=float(h.iloc[0].gap) if len(same)==1 else np.nan
            else:row[v]=np.nan
        rows.append(row)
    return pd.DataFrame(rows),available


def run(root,out,file=None,requested_groups=None):
    out.mkdir(parents=True,exist_ok=False);raw=root/'data/private/adni/raw'
    pd.DataFrame(REQUESTS,columns=['variable','source_files','required_fields','coding_rule']).to_csv(out/'required_source_data.csv',index=False,encoding='utf-8-sig')
    if file is None:
        gds=raw/'GDSCALE.csv'
        if not gds.exists():
            pd.DataFrame([dict(component=v,status='requires_restricted_source_data',reason='No verified dated confounder data supplied; not coded as absent') for v in FIELDS]).to_csv(out/'analysis_status.csv',index=False)
            return
        table=pd.read_csv(gds,low_memory=False)
        if not {'RID','EXAMDATE','GDTOTAL'}.issubset(table):raise ValueError('GDSCALE schema differs; review the data dictionary')
        table=table[['RID','EXAMDATE','GDTOTAL']].rename(columns={'RID':'rid','EXAMDATE':'assessment_date','GDTOTAL':'depression_score'})
        table.loc[~pd.to_numeric(table.depression_score,errors='coerce').between(0,15),'depression_score']=np.nan
    else:table=pd.read_csv(file,dtype={'rid':str})
    data,master,_=adni.load_data(root);aligned,fields=align(table,master);data=data.merge(aligned,on='rid',validate='one_to_one')
    status=[];results=[];coverage=[];diagnostics=[]
    for v in FIELDS:
        n=int(data[v].notna().sum()) if v in fields else 0
        lag=data.get(v+'_lag_days',pd.Series(dtype=float)).dropna()
        coverage.append(dict(variable=v,n_available=n,denominator=len(data),window='any prior visit (historical status)' if v.startswith('prior_') else 'prePET 0-180 days',median_lag_days=lag.median(),max_lag_days=lag.max()))
        if v not in fields:status.append(dict(component=v,status='requires_source_data',reason='Variable not supplied'))
    people,long,cdr=longitudinal.prepare(root);long=long.merge(aligned,on='rid',validate='many_to_one')
    cov=inc.DEMOGRAPHIC+' + centiloid10 + C(diagnosis) + hip1000'
    groups={v:[v] for v in fields}
    groups.update(depression_and_medications=['depression_score','antidementia_medication','antidepressant_medication'],
                  depression_medications_ischemic=['depression_score','antidementia_medication','antidepressant_medication','ischemic_score'],
                  historical_multimorbidity=['prior_cardiovascular_history','prior_medical_domain_count','prior_other_neurologic_history'],
                  initial_health_documentation=['initial_recorded_cv','initial_recorded_neuro','initial_recorded_domain_count'])
    if requested_groups is not None and not set(requested_groups).issubset(groups):raise ValueError('Unknown or unavailable requested confounder group')
    for group,used in groups.items():
        if requested_groups is not None and group not in requested_groups:continue
        if not set(used).issubset(fields):continue
        print('Confounder group:',group,flush=True)
        for y in cross.OUTCOMES:
            d=cross.eligible(data,y);d=d.loc[d.amyloid_positive.eq(1)].dropna(subset=['centiloids']+used)
            varying=[v for v in used if d[v].nunique()>1]
            if not varying or len(d)<30:
                status.append(dict(component=group,outcome=y,status='not_estimable',reason='Insufficient sample or no covariate variation',n=len(d)));continue
            for label,extra in [('same_sample_reference',''),('additional_confounders',' + '+' + '.join(varying))]:
                try:
                    result=adni.fit_effect(d,y,'tau_per_01 + centiloid10 + '+adni.COV+extra,'tau_per_01')
                    results.append(dict(group=group,analysis=label,outcome=y,adjusted_fields=';'.join(varying),**result))
                    status.append(dict(component=group+'__'+label,outcome=y,status='completed',n=len(d)))
                except (ValueError,np.linalg.LinAlgError) as e:status.append(dict(component=group+'__'+label,outcome=y,status='not_estimable',reason=str(e)))
        d=long.dropna(subset=used+['centiloids','hippocampus_icv']);d=d.loc[d.diagnosis.ne('Unavailable')&d.mri_match_status.eq('matched')]
        varying=[v for v in used if d[v].nunique()>1]
        if not varying or d.rid.nunique()<40:
            status.append(dict(component=group,outcome='longitudinal_CDRSB',status='not_estimable',reason='Fewer than 40 participants or no covariate variation',n=d.rid.nunique()));continue
        for label,extra in [('same_sample_reference',''),('additional_confounders',' + '+' + '.join(varying))]:
            try:
                fit,X,diagnostic=sensitivity.fit_gee(d,'time_years * (tau01 + '+cov+extra+')',True)
                diagnostics.append(dict(group=group,analysis=label,adjusted_fields=';'.join(varying),**diagnostic))
                status.append(dict(component=group+'__'+label,outcome='longitudinal_CDRSB',status='completed',n=d.rid.nunique()))
                results.extend(dict(group=group,analysis=label,outcome='longitudinal_CDRSB',adjusted_fields=';'.join(varying),**r) for r in sensitivity.standardized_change(fit,X.design_info,d,True,label))
            except (ValueError,np.linalg.LinAlgError) as e:status.append(dict(component=group+'__'+label,outcome='longitudinal_CDRSB',status='not_estimable',reason=str(e)))
    pd.DataFrame(results).to_csv(out/'additional_confounder_results.csv',index=False)
    pd.DataFrame(coverage).to_csv(out/'confounder_coverage.csv',index=False)
    pd.DataFrame(diagnostics).to_csv(out/'longitudinal_diagnostics.csv',index=False)
    pd.DataFrame(status or [dict(component='supplied_confounders',status='completed')]).to_csv(out/'analysis_status.csv',index=False)
    (out/'analysis_specification.json').write_text(json.dumps(dict(alignment='Latest nonmissing unambiguous pre-PET value within 180 days; prior_ variables explicitly use historical assessments without 180-day restriction',
                   covariates=fields,selection='Same complete sample before/after adjustment; no absent-as-zero imputation',
                   inference='Post hoc sensitivity; additional variables and medication indications may be consequences of disease; nominal conditional intervals, not causal effects'),indent=2),encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--project-root',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--confounder-file',type=Path)
    p.add_argument('--groups',nargs='+',help='Optional subset of named confounder groups')
    a=p.parse_args();run(a.project_root,a.output_dir,a.confounder_file,a.groups)
