"""History-dependent visit-observation sensitivity on fixed annual windows."""
import argparse
from datetime import timedelta
import json
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm
import longitudinal as core
import longitudinal_sensitivity as ls
import longitudinal_model_comparison as rev
import analyze_tau_incremental as inc

BASE=inc.DEMOGRAPHIC+' + baseline_cdrsb + C(diagnosis) + tau01 + centiloid10'
HISTORY='last_cdr + prior_visit_count + years_since_last + previous_window_observed'


def person_period(people,cdr,cutoff):
    rows=[];flow=[]
    for p in people:
        visits=sorted((d,next(iter(v))) for d,v in cdr.get(p['rid'],{}).items() if len(v)==1 and d<=cutoff)
        previous_observed=0
        for year in (1,2,3):
            start=p['anchor']+timedelta(days=(year-.5)*365.25)
            end=p['anchor']+timedelta(days=(year+.5)*365.25)
            if end>cutoff:
                flow.append(dict(year=year,reason='window_not_complete_at_data_cutoff'));continue
            history=[(d,v) for d,v in visits if d<=start]
            if not history:
                flow.append(dict(year=year,reason='no_prior_valid_cdr'));continue
            available=[(d,v) for d,v in visits if start<d<=end and d>p['baseline_date']]
            ordered=sorted(available,key=lambda x:(abs((x[0]-p['anchor']).days/365.25-year),x[0]))
            chosen=ordered[0] if ordered else None
            last_date,last_value=history[-1]
            rows.append(dict(**p,year=year,observed=int(chosen is not None),last_cdr=last_value,
                             prior_visit_count=sum(d>p['baseline_date'] for d,v in history),
                             years_since_last=(start-last_date).days/365.25,previous_window_observed=previous_observed,
                             outcome_cdrsb=chosen[1] if chosen else np.nan,
                             visit_offset=(chosen[0]-p['anchor']).days/365.25-year if chosen else np.nan))
            previous_observed=int(chosen is not None)
    return inc.scaled(pd.DataFrame(rows)),pd.DataFrame(flow)


def fit_probability(data,rhs):
    X=patsy.dmatrix(rhs,data,return_type='dataframe',NA_action='raise')
    if np.linalg.matrix_rank(X)!=X.shape[1] or data.observed.nunique()!=2:raise ValueError('Observation design not estimable')
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter('always');f=sm.GLM(data.observed,X,family=sm.families.Binomial()).fit(maxiter=200)
    prob=np.asarray(f.predict(X))
    if not f.converged or not np.isfinite(prob).all() or any('Separation' in type(w.message).__name__ for w in ws):raise ValueError('Observation fit failed')
    if (prob<.001).any() or (prob>.999).any():raise ValueError('Observation probabilities outside [0.001,0.999]; positivity failure')
    return prob,f


def estimate(data,weighted=True):
    rhs='C(year) + '+BASE+' + '+HISTORY
    prob,selection=fit_probability(data,rhs)
    observed=data.observed.eq(1);d=data.loc[observed].copy()
    marginal=data.groupby('year')['observed'].mean()
    raw=d.year.map(marginal).to_numpy()/prob[observed]
    cut=np.quantile(raw,[.01,.99]);weights=np.clip(raw,*cut) if weighted else np.ones(len(d))
    model_rhs='C(year) * tau01 + '+inc.DEMOGRAPHIC+' + baseline_cdrsb + C(diagnosis) + centiloid10 + visit_offset'
    X=patsy.dmatrix(model_rhs,d,return_type='dataframe',NA_action='raise')
    if np.linalg.matrix_rank(X)!=X.shape[1] or d.rid.nunique()<40:raise ValueError('Outcome design not estimable')
    f=sm.GLM(d.outcome_cdrsb/18,X,family=sm.families.Binomial(),freq_weights=weights).fit(maxiter=200)
    if not f.converged or not np.isfinite(f.params).all():raise ValueError('Outcome fit failed')
    effects=[]
    for year in (1,2,3):
        target=data.loc[data.year.eq(year)].drop_duplicates('rid').copy()
        target['visit_offset']=0.
        v=inc.tau_contrast(f,X.design_info,target,True)
        effects.append(dict(year=year,estimate=v['difference_per_01'],n_target=len(target),n_observed=int((d.year==year).sum())))
    diagnostic=dict(n_people=data.rid.nunique(),n_person_periods=len(data),n_observed=len(d),p_min=prob.min(),p_max=prob.max(),
                    weight_max_raw=raw.max(),weight_max_used=weights.max(),trim_low=cut[0],trim_high=cut[1],effective_observations=weights.sum()**2/(weights@weights))
    coefficients=pd.DataFrame(dict(term=selection.params.index,beta=selection.params.values,p=selection.pvalues.values))
    return effects,diagnostic,coefficients


def run(root,out,bootstrap,seed):
    out.mkdir(parents=True,exist_ok=False)
    people,cdr=rev.prepare_followup_descriptions(root)
    master,_=core.load_master(root/'data/private/adni/adni_pet_aligned_master.csv')
    cl,_=inc.match_centiloids(master,root/'data/private/adni/raw/UCBERKELEY_AMY_6MM.csv')
    for p in people:p['centiloids']=cl[p['rid']]['centiloids']
    people=[p for p in people if np.isfinite(p['centiloids']) and p['diagnosis']!='Unavailable']
    cutoff,source=ls.cutoff_date(root/'data/private/adni/raw')
    data,flow=person_period(people,cdr,cutoff)
    # A single target population must have administrative opportunity for all three windows.
    complete_ids=data.groupby('rid').year.nunique();data=data.loc[data.rid.isin(complete_ids[complete_ids==3].index)].copy()
    support=data.groupby('year').agg(n_target=('rid','nunique'),n_observed=('observed','sum')).reset_index()
    support.to_csv(out/'annual_window_support.csv',index=False)
    if len(flow):flow.value_counts().rename('n').reset_index().to_csv(out/'calendar_exclusions.csv',index=False)
    results=[];diagnostics=[];draws=[];statuses=[]
    for label,weighted in [('unweighted_windows',False),('history_weighted_windows',True)]:
        try:
            r,d,c=estimate(data,weighted);diagnostics.append(dict(model=label,**d));c.to_csv(out/f'{label}_visit_coefficients.csv',index=False)
            results.extend(dict(model=label,**x) for x in r)
        except (ValueError,np.linalg.LinAlgError) as e:statuses.append(dict(component=label,status='not_estimable',reason=str(e)))
    if results:
        rng=np.random.default_rng(seed);ids=data.rid.unique();groups={r:g for r,g in data.groupby('rid')}
        for b in range(bootstrap):
            sample=[]
            for copy,rid in enumerate(rng.choice(ids,len(ids),replace=True)):
                g=groups[rid].copy();g['rid']=str(copy);sample.append(g)
            sample=pd.concat(sample,ignore_index=True)
            for label,weighted in [('unweighted_windows',False),('history_weighted_windows',True)]:
                try:
                    es,_,_=estimate(sample,weighted);draws.extend(dict(replicate=b+1,model=label,status='ok',**x) for x in es)
                except (ValueError,np.linalg.LinAlgError) as e:draws.append(dict(replicate=b+1,model=label,status='failed',reason=str(e)))
            if (b+1)%25==0:print('Visit-process bootstrap',b+1,'/',bootstrap,flush=True)
        boot=pd.DataFrame(draws);boot.to_csv(out/'bootstrap_replicates.csv',index=False)
        for r in results:
            v=boot.loc[boot.model.eq(r['model']) & boot.status.eq('ok') & boot.year.eq(r['year']),'estimate']
            r.update(bootstrap_valid=len(v),bootstrap_requested=bootstrap,status='ok' if len(v)>=.9*bootstrap else 'bootstrap_unstable')
            if len(v)>=.9*bootstrap:r.update(ci_low=v.quantile(.025),ci_high=v.quantile(.975))
    pd.DataFrame(results).to_csv(out/'visit_process_results.csv',index=False)
    pd.DataFrame(diagnostics).to_csv(out/'weight_diagnostics.csv',index=False)
    pd.DataFrame(statuses or [dict(component='visit_process',status='completed')]).to_csv(out/'analysis_status.csv',index=False)
    spec=dict(cutoff=str(cutoff),cutoff_source=source,windows='Non-overlapping (year-0.5,year+0.5] relative to tau PET; nearest observed CDR in each window',
              target='A+ complete baseline/CL predictors with all three windows administratively observable; 3-year window extends to 3.5 years, unlike original 3-year follow-up analysis',
              selection='P(observed in window | year, baseline predictors, last CDR strictly before window, previous-window observed, count of prior visits, time since last)',
              weights='Period-specific marginal observation probability divided by modeled conditional probability, clipped at observed 1st/99th percentiles; not a cumulative survival weight',
              inference='Participant bootstrap refits both models and trimming. No person IDs exported.',
              estimand='Conditional CDR-SB level contrast per 0.1 SUVR normalized IQR, baseline severity held fixed; not an annual slope',
              limitations='History-dependent observation sensitivity, not correction for unobserved worsening, death, outcome-dependent within-window scheduling or non-ignorable missingness; last CDR is observed history, not imputed outcome')
    (out/'analysis_specification.json').write_text(json.dumps(spec,indent=2),encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--project-root',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--bootstrap',type=int,default=1000);p.add_argument('--seed',type=int,default=20261017)
    a=p.parse_args()
    if a.bootstrap<100:p.error('Use at least 100 replicates for a diagnostic run; 1000 for final run')
    run(a.project_root,a.output_dir,a.bootstrap,a.seed)
