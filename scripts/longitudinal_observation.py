"""Common-sample observation-process and stage-support sensitivity analyses."""
from datetime import timedelta
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm
import longitudinal_model_comparison as rev
import longitudinal_sensitivity as ls
import analyze_tau_incremental as inc
from analysis_utils import csv,dump,checkpoint,NotEstimable,bootstrap_people

COV=inc.DEMOGRAPHIC+' + centiloid10 + C(diagnosis) + hip1000'
RHS='time_years * (tau01 + '+COV+')'

class ObservationNotEstimable(NotEstimable):
    def __init__(self,reason,diagnostics):
        super().__init__(reason)
        self.diagnostics=dict(diagnostics,status='not_estimable',reason=reason)

def probability_diagnostics(pred,observed):
    pred=np.asarray(pred,dtype=float);observed=np.asarray(observed,dtype=int)
    if pred.shape!=observed.shape or not np.isfinite(pred).all():
        raise ValueError('Invalid observation probability diagnostics')
    low=pred<.001;high=pred>.999
    return dict(p_min=float(pred.min()),p_max=float(pred.max()),positivity_lower=.001,
        positivity_upper=.999,n_probability_below=int(low.sum()),n_probability_above=int(high.sum()),
        below_observed=int((low&(observed==1)).sum()),below_unobserved=int((low&(observed==0)).sum()),
        above_observed=int((high&(observed==1)).sum()),above_unobserved=int((high&(observed==0)).sum()))

def fit_gee(d,rhs=RHS,weights=None):
    X=patsy.dmatrix(rhs,d,return_type='dataframe',NA_action='raise')
    if d.rid.nunique()<40 or np.linalg.matrix_rank(X)!=X.shape[1]:raise NotEstimable('Trajectory model rank/support failure')
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter('always')
        f=sm.GEE(d.cdr.to_numpy()/18,X,groups=d.rid.to_numpy(),family=sm.families.Binomial(),cov_struct=sm.cov_struct.Independence(),weights=weights).fit(maxiter=200)
    if not f.converged or not np.isfinite(f.params).all() or not np.isfinite(f.cov_params()).all().all():raise NotEstimable('Trajectory GEE failed')
    return f,X.design_info

def common_data(root):
    people,long,cdr=rev.prepare(root)
    common=long.dropna(subset=['centiloids','hippocampus_icv','baseline_cdrsb']).copy()
    common=common.loc[common.diagnosis.ne('Unavailable')&common.mri_match_status.eq('matched')].reset_index(drop=True)
    if common.rid.nunique()<40:raise NotEstimable('Insufficient common-sample participants')
    return people,common,cdr

def change(f,info,target,label):
    return next(r for r in ls.standardized_change(f,info,target,True,label) if r['year']==2)

def observation_periods(common,cdr,cutoff):
    rows=[];mapping={}
    for key,p in common.drop_duplicates('rid').set_index('rid').iterrows():
        start0=max(p.anchor,p.baseline_date);end0=min(p.anchor+timedelta(days=3*365.25),cutoff)
        visits=sorted((day,next(iter(v))) for day,v in cdr.get(key,{}).items() if len(v)==1 and day<=cutoff)
        for index in range(7):
            left=max(start0,p.baseline_date+timedelta(days=index*.5*365.25))
            right=min(end0,p.baseline_date+timedelta(days=(index+1)*.5*365.25))
            if right<=left:continue
            history=[(day,value) for day,value in visits if day<=left]
            if not history:raise NotEstimable('No prior observed CDR for an eligible observation window')
            d=common.loc[common.rid.eq(key)&common.time_years.gt(0)]
            time=np.rint(d.time_years*365.25)
            selected=d.index[(time>(left-p.baseline_date).days)&(time<=(right-p.baseline_date).days)]
            for row in selected:mapping[row]=(key,index)
            rows.append(dict(**p.to_dict(),rid=key,period=index,observed=int(len(selected)>0),
                visits=len(selected),duration_years=(right-left).days/365.25,last_cdr=history[-1][1],
                prior_visits=sum(day>p.baseline_date for day,value in history),years_since_last=(left-history[-1][0]).days/365.25))
    periods=pd.DataFrame(rows)
    expected=set(common.index[common.time_years.gt(0)])
    if set(mapping)!=expected:raise ValueError('Observation windows failed to cover every primary follow-up row exactly once')
    long=common.copy();long['period']=[mapping[i][1] if i in mapping else -1 for i in long.index]
    return periods,long

def observation_weights(periods,long,trim):
    rhs='C(period) + '+COV+' + tau01 + baseline_cdrsb + last_cdr + prior_visits + years_since_last + duration_years'
    # Periods always observed carry weight one, not a fictitious fitted Bernoulli probability.
    varying=periods.groupby('period')['observed'].nunique();active=periods.period.isin(varying[varying>1].index)
    p=np.ones(len(periods));fit_data=periods.loc[active]
    by_period=periods.groupby('period')['observed'].agg(['min','max'])
    diag=dict(n_people=int(long.rid.nunique()),n_observations=len(long),n_periods=len(periods),
        n_observed_periods=int(periods.observed.sum()),n_variable_windows=len(fit_data),
        n_constant_zero_periods=int((by_period['max']==0).sum()),
        n_constant_one_periods=int((by_period['min']==1).sum()),
        n_baseline_observations=int(long.period.eq(-1).sum()),
        n_followup_observations=int(long.period.ne(-1).sum()),
        bootstrap_requested=0,bootstrap_attempted=0,bootstrap_valid=0,
        status='completed',reason='')
    def fail(reason):raise ObservationNotEstimable(reason,diag)
    if len(fit_data)<40 or fit_data.observed.nunique()!=2:fail('Insufficient variable observation windows')
    X=patsy.dmatrix(rhs,fit_data,return_type='dataframe',NA_action='raise')
    diag.update(design_rank=int(np.linalg.matrix_rank(X)),design_columns=X.shape[1])
    if diag['design_rank']!=X.shape[1]:fail('Observation model rank deficient; no covariates silently dropped')
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter('always');f=sm.GLM(fit_data.observed,X,family=sm.families.Binomial()).fit(maxiter=200)
    pred=np.asarray(f.predict(X))
    diag.update(converged=bool(f.converged),separation_warning=any('Separation' in type(w.message).__name__ for w in ws))
    if not f.converged or not np.isfinite(pred).all() or diag['separation_warning']:fail('Observation model separation/convergence failure')
    diag.update(probability_diagnostics(pred,fit_data.observed))
    if pred.min()<.001 or pred.max()>.999:fail(f'Observation positivity failure: probability range [{pred.min():.6g},{pred.max():.6g}]; {int(((pred<.001)|(pred>.999)).sum())} windows outside [0.001,0.999]. No probability clipping or result-based model simplification applied.')
    p[active]=pred;periods=periods.copy();periods['probability']=p
    marginal=periods.groupby('period')['observed'].mean();periods['raw_weight']=periods.period.map(marginal)/p
    lookup=periods.set_index(['rid','period']).raw_weight
    raw=np.array([1. if r.period==-1 else lookup.loc[(r.rid,r.period)] for r in long.itertuples()])
    follow=long.period.ne(-1).to_numpy();cuts=np.quantile(raw[follow],trim) if trim else [raw[follow].min(),raw[follow].max()]
    weights=raw.copy();weights[follow]=np.clip(raw[follow],*cuts)
    diag.update(trim_low=float(cuts[0]),trim_high=float(cuts[1]),
        raw_weight_max=float(raw.max()),used_weight_max=float(weights.max()),effective_observations=float(weights.sum()**2/(weights@weights)),
        n_people=int(long.rid.nunique()),n_observations=len(long),n_periods=len(periods),n_observed_periods=int(periods.observed.sum()))
    return weights,diag

def paired_estimates(periods,long):
    result=[];diags=[]
    f,info=fit_gee(long);result.append(change(f,info,long,'unweighted_primary_sample'))
    for label,trim in [('untrimmed',None),('trim_01_99',[.01,.99]),('trim_05_95',[.05,.95])]:
        weights,diag=observation_weights(periods,long,trim);f,info=fit_gee(long,weights=weights)
        result.append(change(f,info,long,label));diags.append(dict(model=label,**diag))
    return result,diags

def run(root,out,bootstrap=500,seed=20261029):
    root=Path(root);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    people,common,cdr=common_data(root);cutoff,_=ls.cutoff_date(root/'data/private/adni/raw')
    states=[];results=[];tests=[];diagnostics=[]
    periods,long=observation_periods(common,cdr,cutoff)
    csv(out/'observation_support.csv',periods.groupby('period').agg(people=('rid','nunique'),observed_people=('observed','sum'),visits=('visits','sum')).reset_index())
    try:
        points,diagnostics=paired_estimates(periods,long);draws=[]
        groups=dict(tuple(periods.groupby('rid')))
        for b in range(bootstrap):
            def one():
                rng=np.random.default_rng(seed+b);sample=bootstrap_people(long,rng);pp=[]
                for row in sample[['rid','original_rid']].drop_duplicates().itertuples():
                    z=groups[row.original_rid].copy();z['rid']=row.rid;pp.append(z)
                try:
                    es,_=paired_estimates(pd.concat(pp,ignore_index=True),sample)
                    return [dict(replicate=b+1,status='completed',**r) for r in es]
                except (NotEstimable,np.linalg.LinAlgError) as e:return [dict(replicate=b+1,status='not_estimable',reason=str(e))]
            draws.extend(checkpoint(out,'weight_boot_'+str(b+1),one))
            if (b+1)%10==0:print('Common-sample observation bootstrap',b+1,'/',bootstrap,flush=True)
        draw=pd.DataFrame(draws)
        for row in points:
            v=draw.loc[draw.status.eq('completed')&draw.model.eq(row['model']),'estimate'] if 'model' in draw else pd.Series(dtype=float)
            row.update(bootstrap_valid=len(v),bootstrap_requested=bootstrap,status='completed' if len(v)>=.9*bootstrap else 'bootstrap_unstable')
            row['ci_low']=float(v.quantile(.025)) if len(v)>=.9*bootstrap else np.nan;row['ci_high']=float(v.quantile(.975)) if len(v)>=.9*bootstrap else np.nan
        csv(out/'observation_results.csv',points);csv(out/'observation_bootstrap.csv',draw)
        for row in diagnostics:
            point=next(r for r in points if r['model']==row['model'])
            row.update(bootstrap_requested=bootstrap,bootstrap_attempted=bootstrap,bootstrap_valid=point['bootstrap_valid'],status=point['status'])
        states.append(dict(component='common_observation',status='completed' if all(r['status']=='completed' for r in points) else 'bootstrap_unstable',reason='Period-observation sensitivity within the selected common sample, not full selection correction'))
    except (NotEstimable,np.linalg.LinAlgError) as e:
        states.append(dict(component='common_observation',status='not_estimable',reason=str(e)))
        if isinstance(e,ObservationNotEstimable):
            diagnostics=[dict(e.diagnostics,model='observation_probability_model',bootstrap_requested=bootstrap)]
    csv(out/'weight_diagnostics.csv',diagnostics)
    support=[];people=common.drop_duplicates('rid')
    for stage,g in people.groupby('diagnosis'):
        support.append(dict(diagnosis=stage,n=len(g),tau_min=g.tau_temporal_suvr.min(),tau_q05=g.tau_temporal_suvr.quantile(.05),tau_q25=g.tau_temporal_suvr.quantile(.25),tau_q75=g.tau_temporal_suvr.quantile(.75),tau_q95=g.tau_temporal_suvr.quantile(.95),tau_max=g.tau_temporal_suvr.max()))
    s=pd.DataFrame(support);csv(out/'stage_tau_support.csv',s)
    try:
        if len(s)!=3 or s.n.min()<10:raise NotEstimable('Fewer than ten participants in a clinical stage')
        rhs=RHS+' + tau01:C(diagnosis) + time_years:tau01:C(diagnosis)'
        f,info=fit_gee(common,rhs)
        terms=[k for k in f.params.index if 'time_years' in k and 'tau01' in k and 'C(diagnosis)' in k]
        if len(terms)!=2:raise ValueError('Expected two time-by-tau-by-stage terms')
        test=ls.robust_joint_test(f,terms,'time_tau_stage_heterogeneity');tests.append(test)
        for stage,g in common.groupby('diagnosis'):results.append(dict(scope='stage_specific',diagnosis=stage,**change(f,info,g,'stage_interaction')))
        states.append(dict(component='stage_interaction',status='completed',reason='One joint exploratory interaction test; stage contrasts use their own observed IQRs'))
    except (NotEstimable,np.linalg.LinAlgError) as e:states.append(dict(component='stage_interaction',status='not_estimable',reason=str(e)))
    lower,upper=float(s.tau_q05.max()),float(s.tau_q95.min())
    trimmed=common.loc[common.tau_temporal_suvr.between(lower,upper)].copy()
    try:
        counts=trimmed.drop_duplicates('rid').diagnosis.value_counts()
        if lower>=upper or len(counts)!=3 or counts.min()<10:raise NotEstimable('Insufficient central-90%-range overlap across all three stages')
        f,info=fit_gee(trimmed);results.append(dict(scope='common_tau_support',diagnosis='all',**change(f,info,trimmed,'support_restricted')))
        states.append(dict(component='support_restricted',status='completed',reason='Restricted estimand, not a replacement primary analysis'))
    except (NotEstimable,np.linalg.LinAlgError) as e:states.append(dict(component='support_restricted',status='not_estimable',reason=str(e)))
    csv(out/'stage_contrasts.csv',results);csv(out/'stage_interaction_tests.csv',tests);csv(out/'analysis_status.csv',states)
    dump(out/'specification.json',dict(n_people=int(common.rid.nunique()),n_observations=len(common),bootstrap=bootstrap,
      weights='Half-year baseline-CDR windows clipped to PET follow-up horizon and administrative cutoff; baseline observations weight one. Each actual follow-up retains its original outcome and time. History precedes each window. Weights are period-specific, not cumulative.',
      limitations='Conditions on the original common sample. Does not correct initial sample selection, death, within-window visit frequency or unmeasured outcome-dependent attendance. Multiple visits in a window share a weight.',
      stages='Same primary covariate set, two-degree-of-freedom time-by-tau-by-diagnosis test. No inference of equal effects from a null test.',
      overlap_lower=lower,overlap_upper=upper,overlap_rule='Intersection of stage-specific 5th to 95th tau percentiles, selected without outcomes.'))

