"""Refitted calibration uncertainty and site-grouped internal validation."""
from pathlib import Path
import numpy as np
import pandas as pd
import patsy
import analyze_tau_incremental as inc
import pet_day_clinical_models as timing
import pet_staging_covariates as adni
from analysis_utils import csv,dump,checkpoint,NotEstimable,bootstrap_people

MODELS={k:v.replace(' + visit_offset','') for k,v in inc.MODEL_RHS.items() if k in ['clinical_amyloid_mri','clinical_amyloid_mri_tau']}

def folds_for(groups,rng,k=5):
    counts=pd.Series(groups).value_counts()
    if len(counts)<k:raise NotEstimable('Fewer than five independent validation groups')
    order=list(rng.permutation(counts.index));order.sort(key=lambda x:-counts[x])
    sizes=np.zeros(k,int);assigned={}
    for g in order:
        f=int(rng.choice(np.flatnonzero(sizes==sizes.min())));assigned[g]=f;sizes[f]+=int(counts[g])
    return np.array([assigned[g] for g in groups])

def validation(d,repeats,seed,mode):
    if len(d)<40:raise NotEstimable('Fewer than 40 observed endpoints')
    groups=d.site.astype(str) if mode=='site' else d.get('original_rid',d.rid).astype(str)
    if mode=='site' and (groups.isin(['Unavailable','nan','']).any() or groups.nunique()<5):raise NotEstimable('Site validation requires known sites and at least five sites')
    y=d.outcome_cdrsb.to_numpy();estimates=[];predictions=[]
    for rep in range(repeats):
        assignment=folds_for(groups,np.random.default_rng(seed+rep))
        for family in ['gaussian','fractional_logit']:
          for label,rhs in MODELS.items():
            pred=np.full(len(d),np.nan)
            for fold in range(5):
                train=assignment!=fold;test=~train
                if not test.any() or set(groups[train])&set(groups[test]):raise AssertionError('Validation split leakage')
                # Build category coding from training data only. Unseen held-out levels make this split non-estimable.
                X=patsy.dmatrix(rhs,d.loc[train],return_type='dataframe',NA_action='raise')
                try:T=patsy.build_design_matrices([X.design_info],d.loc[test],NA_action='raise')[0]
                except patsy.PatsyError as e:raise NotEstimable('Held-out categorical level absent from training data') from e
                if groups[train].nunique()<2 or d.loc[train].get('original_rid',d.loc[train].rid).nunique()<=X.shape[1]+5:raise NotEstimable('Insufficient independent training participants')
                try:f=inc.fit_endpoint(X,y[train],family=='fractional_logit')
                except ValueError as e:raise NotEstimable(str(e)) from e
                pred[test]=inc.predict(f,T,family=='fractional_logit')
            if not np.isfinite(pred).all() or np.ptp(pred)<1e-10:raise NotEstimable('Invalid or constant out-of-fold predictions')
            estimates.append(dict(repeat=rep+1,family=family,model=label,**inc.metrics(y,pred)))
            predictions.append(dict(family=family,model=label,pred=pred))
    result=pd.DataFrame(estimates).groupby(['family','model'],as_index=False).mean(numeric_only=True).drop(columns='repeat')
    paired=[]
    for family,g in result.groupby('family'):
        a=g.set_index('model');pair=dict(family=family,model='tau_minus_no_tau')
        for c in ['rmse','mae','r2','mean_error']:pair[c]=a.loc['clinical_amyloid_mri_tau',c]-a.loc['clinical_amyloid_mri',c]
        paired.append(pair)
    bins=[]
    for family in ['gaussian','fractional_logit']:
      for model in MODELS:
        pred=np.mean([r['pred'] for r in predictions if r['family']==family and r['model']==model],axis=0)
        # Do not export individual predictions or identifiers.
        order=np.argsort(pred,kind='stable')
        for index,indices in enumerate(np.array_split(order,5)):
            bins.append(dict(family=family,model=model,bin=index+1,n=len(indices),predicted=float(pred[indices].mean()),observed=float(y[indices].mean())))
    return pd.concat([result,pd.DataFrame(paired)],ignore_index=True).to_dict('records'),bins

def bootstrap_sites(d,rng):
    groups=dict(tuple(d.groupby('site')));sites=list(groups);blocks=[]
    for index,key in enumerate(rng.choice(sites,len(sites),replace=True)):
        z=groups[key].copy();z['original_rid']=z.rid;z['rid']=str(index)+'_'+z.rid.astype(str);blocks.append(z)
    return pd.concat(blocks,ignore_index=True)

def run(root,out,bootstrap=200,repeats=5,seed=20261029):
    root=Path(root);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    master,_,_=adni.load_data(root);sites=master.set_index('rid').site
    states=[];summary=[];binrows=[];draws=[];flows=[]
    for cutoff in [0,180]:
        people,cdr,flow,lags=timing.match_landmark(root,cutoff)
        d,flow2=inc.endpoint_cohort(people,cdr,2,180,cutoff);flows.extend(flow);flows.extend([dict(**r,landmark_days=cutoff) for r in flow2])
        if d.empty:
            states.append(dict(component=f'landmark_{cutoff}',status='not_estimable',reason='No eligible observed endpoints'));continue
        d=d.reset_index(drop=True);d['site']=d.rid.map(sites).fillna('Unavailable')
        for mode in ['participant','site']:
            label=f'landmark_{cutoff}_{mode}'
            try:point,bins=validation(d,repeats,seed,mode)
            except NotEstimable as e:
                states.append(dict(component=label,status='not_estimable',reason=str(e)));continue
            binrows.extend([dict(component=label,**r) for r in bins]);local=[]
            for b in range(bootstrap):
                def one():
                    rng=np.random.default_rng(seed+10000+cutoff+b)
                    sample=bootstrap_sites(d,rng) if mode=='site' else bootstrap_people(d,rng)
                    try:
                        metrics,_=validation(sample,repeats,seed+100000+b,mode)
                        return [dict(replicate=b+1,status='completed',**r) for r in metrics]
                    except (NotEstimable,np.linalg.LinAlgError) as e:return [dict(replicate=b+1,status='not_estimable',reason=str(e))]
                local.extend(checkpoint(out,label+'_boot'+str(b+1),one))
                if (b+1)%10==0:print(label,'calibration bootstrap',b+1,'/',bootstrap,flush=True)
            dr=pd.DataFrame(local);draws.extend([dict(component=label,**r) for r in local]);valid=int(dr.loc[dr.status.eq('completed'),'replicate'].nunique())
            for row in point:
                for metric,value in row.items():
                    if metric in ['family','model'] or not np.isfinite(value):continue
                    v=dr.loc[dr.status.eq('completed')&dr.family.eq(row['family'])&dr.model.eq(row['model']),metric] if 'family' in dr and metric in dr else pd.Series(dtype=float)
                    ok=len(v.dropna())>=.9*bootstrap
                    summary.append(dict(component=label,family=row['family'],model=row['model'],metric=metric,estimate=value,
                        ci_low=float(v.quantile(.025)) if ok else np.nan,ci_high=float(v.quantile(.975)) if ok else np.nan,
                        bootstrap_valid=len(v.dropna()),bootstrap_requested=bootstrap,n=len(d),n_sites=d.site.nunique(),status='completed' if ok else 'bootstrap_unstable'))
            states.append(dict(component=label,status='completed' if valid>=.9*bootstrap else 'bootstrap_unstable',reason='Internal validation; not independent external clinical validation',bootstrap_valid=valid,bootstrap_requested=bootstrap))
            csv(out/'analysis_status.csv',states)
    csv(out/'calibration_uncertainty.csv',summary);csv(out/'calibration_bins.csv',binrows);csv(out/'calibration_bootstrap.csv',draws);csv(out/'cohort_flow.csv',flows);csv(out/'analysis_status.csv',states)
    dump(out/'specification.json',dict(bootstrap=bootstrap,repeats=repeats,folds=5,
      outcome='Observed CDR-SB nearest two years within +/-180 days; a continuous bounded outcome, not a conversion event.',
      predictors='Measured by PET day or the fixed 180-day landmark. Outcome-visit offset is excluded.',
      uncertainty='Percentile intervals from refitting the entire paired validation procedure; participant bootstrap for participant CV and site-cluster bootstrap for site CV. Duplicate participants/sites remain in a single fold. At least 90% valid resamples required.',
      calibration='Continuous-outcome OLS intercept and slope fitted to held-out predictions, averaged across repeats; ideal intercept 0, slope 1. Bins use mean repeated out-of-fold prediction, no individual rows exported.',
      limitations='No independent external clinical cohort. Site grouping is conditional on data support; failure is reported. No calibration benefit inferred solely from a nonsignificant test.'))

