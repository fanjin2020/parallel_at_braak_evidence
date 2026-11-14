"""Common-sample trajectory contrasts and observed-endpoint sensitivities."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from types import SimpleNamespace
import warnings
import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm
import longitudinal as core
import longitudinal_sensitivity as diag
import analyze_tau_incremental as inc
import pet_staging_covariates as adni


def prepare_followup_descriptions(root):
    """Visit summaries do not require PET ROI or MRI refitting inputs."""
    from datetime import timedelta
    raw=root/'data/private/adni/raw'
    master,_=core.load_master(root/'data/private/adni/adni_pet_aligned_master.csv')
    ids={p['rid'] for p in master}
    cdr,_=core.load_visits(raw/'CDR.csv',ids,'cdr')
    dx,_=core.load_visits(raw/'DXSUM.csv',ids,'dx')
    people,_,_=core.make_cohort(master,cdr)
    for p in people:
        d,_=core.nearest(dx.get(p['rid'],{}),p['anchor'],180)
        p['diagnosis']={1:'CN',2:'MCI',3:'Dementia'}[d[1]] if d else 'Unavailable'
        p['baseline_date']=p['anchor']+timedelta(days=p['baseline_offset_days'])
    return people,cdr


def prepare(root):
    df,master,_=adni.load_data(root)
    raw=root/'data/private/adni/raw'; ids={p['rid'] for p in master}
    cdr,_=core.load_visits(raw/'CDR.csv',ids,'cdr')
    dx,_=core.load_visits(raw/'DXSUM.csv',ids,'dx')
    eligible,long,_=core.make_cohort(master,cdr)
    mri,_=inc.match_mri_dates(master,raw/'UCSFFSX7.csv')
    lookup=df.set_index('rid').to_dict('index')
    people=[]
    for p in eligible:
        d,s=core.nearest(dx.get(p['rid'],{}),p['anchor'],180)
        z={**p,**lookup[p['rid']],**mri[p['rid']]}
        z['diagnosis']={1:'CN',2:'MCI',3:'Dementia'}[d[1]] if d else 'Unavailable'
        z['diagnosis_date']=d[0] if d else None
        from datetime import timedelta
        z['baseline_date']=p['anchor']+timedelta(days=p['baseline_offset_days'])
        people.append(z)
    lookup={p['rid']:p for p in people}
    long=pd.DataFrame([{**lookup[r['rid']], 'cdr':r['cdr'],'time_years':r['time_years']} for r in long])
    long=inc.scaled(long)
    long['tau_years']=long.time_years+long.baseline_offset_days/365.25
    return people,long,cdr


def observation_bootstrap(people,endpoint,year,replicates,seed):
    table=pd.DataFrame(people)
    predictors=inc.PREDICTOR_COLUMNS
    table=table.dropna(subset=predictors)
    table=table.loc[table.diagnosis.ne('Unavailable') & table.mri_match_status.eq('matched')].copy()
    from datetime import timedelta
    valid=[]
    for p in table.to_dict('records'):
        landmark=p['anchor']+timedelta(days=180)
        dates=[core.iso_date(p['amy_pet_date']),p['baseline_date'],p['mri_date'],p['diagnosis_date']]
        valid.append(all(d is not None and d<=landmark for d in dates))
    table=inc.scaled(table.loc[valid]).reset_index(drop=True)
    observed=table.rid.isin(endpoint.rid)
    # Endpoint-observation weights target complete, landmark-eligible predictors only.
    sel_rhs=inc.DEMOGRAPHIC+' + baseline_cdrsb + C(diagnosis) + centiloid10 + hip1000 + tau01'
    merged=table.merge(endpoint[['rid','outcome_cdrsb','visit_offset']],on='rid',how='left',validate='one_to_one')
    def one(sample):
        observed=sample.outcome_cdrsb.notna()
        w,wd=adni.observation_weights(sample,observed,sel_rhs)
        d=sample.loc[observed]
        x=patsy.dmatrix(inc.MODEL_RHS['clinical_amyloid_mri_tau'],d,return_type='dataframe',NA_action='raise')
        if np.linalg.matrix_rank(x)!=x.shape[1]: raise ValueError('Endpoint design rank deficient')
        f=sm.GLM(d.outcome_cdrsb/18,x,family=sm.families.Binomial(),freq_weights=w).fit(maxiter=200)
        if not f.converged or not np.isfinite(f.params).all(): raise ValueError('Weighted endpoint did not converge')
        contrast=inc.tau_contrast(f,x.design_info,d,True)
        return contrast['difference_per_01'],wd
    estimate,wd=one(merged)
    rng=np.random.default_rng(seed); draws=[]; failed=0
    for b in range(replicates):
        try:
            value,_=one(merged.iloc[rng.integers(0,len(merged),len(merged))].reset_index(drop=True)); draws.append(value)
        except (ValueError,np.linalg.LinAlgError): failed+=1
    ok=len(draws)>=.9*replicates
    ci=np.quantile(draws,[.025,.975]) if ok else [np.nan,np.nan]
    return dict(year=year,estimate=estimate,ci_low=ci[0],ci_high=ci[1],n_observed=int(observed.sum()),
                n_target=len(table),bootstrap_requested=replicates,bootstrap_valid=len(draws),
                status='ok' if ok else 'too_many_failed_bootstraps',
                interval='Participant bootstrap refits observation and endpoint models and trimming cutoffs',
                estimand='Conditional observed-endpoint level contrast; baseline severity held fixed; per 0.1 SUVR normalized IQR',
                limitation='Measured-predictor MAR sensitivity within the complete-predictor target; bootstrap refits weights but does not address MNAR or full analytic selection',
                **{k:v for k,v in wd.items() if k not in ['n_target','n_observed','limitation']})


def refit_cv_uncertainty(df,reps,folds,seed):
    """Subject pairs bootstrap; duplicate copies always stay in the same fold."""
    rng=np.random.default_rng(seed); rows=[]
    for b in range(reps):
        sampled=df.iloc[rng.integers(0,len(df),len(df))].copy().reset_index(drop=True)
        unique=sampled.rid.unique(); fold={r:i%folds for i,r in enumerate(rng.permutation(unique))}
        k=sampled.rid.map(fold).to_numpy()
        for family in ['gaussian','fractional_logit']:
            preds=np.full((len(sampled),2),np.nan); status='ok'
            try:
                for f in range(folds):
                    train,test=sampled.loc[k!=f],sampled.loc[k==f]
                    if set(train.rid)&set(test.rid): raise AssertionError('Bootstrap ID leakage')
                    for j,label in enumerate(['clinical_amyloid_mri','clinical_amyloid_mri_tau']):
                        x=patsy.dmatrix(inc.MODEL_RHS[label],train,return_type='dataframe')
                        xt=patsy.build_design_matrices([x.design_info],test)[0]
                        fitted=inc.fit_endpoint(x,train.outcome_cdrsb,family=='fractional_logit')
                        preds[k==f,j]=inc.predict(fitted,xt,family=='fractional_logit')
                if not np.isfinite(preds).all(): raise ValueError('Nonfinite prediction')
                values=[inc.metrics(sampled.outcome_cdrsb,preds[:,j]) for j in range(2)]
                rows.append(dict(replicate=b+1,family=family,status=status,
                                 delta_rmse=values[1]['rmse']-values[0]['rmse'],delta_mae=values[1]['mae']-values[0]['mae']))
            except (ValueError,np.linalg.LinAlgError,patsy.PatsyError) as e:
                rows.append(dict(replicate=b+1,family=family,status='failed',reason=str(e).splitlines()[0]))
        if (b+1)%25==0: print('Refitted CV bootstrap',b+1,'/',reps,flush=True)
    data=pd.DataFrame(rows); summaries=[]
    for family,g in data.groupby('family'):
        valid=g.loc[g.status.eq('ok')]
        row=dict(family=family,requested=reps,valid=len(valid),
                 status='ok' if len(valid)>=.9*reps else 'too_many_failed_bootstraps',
                 interpretation='Exploratory bootstrap distribution of refitted grouped CV; distinct from conditional fixed-OOF intervals; not external validation')
        for metric in ['delta_rmse','delta_mae']:
            if len(valid)>=.9*reps:
                row.update({metric+'_q025':valid[metric].quantile(.025),metric+'_q975':valid[metric].quantile(.975)})
        summaries.append(row)
    return summaries,rows


def followup_tables(people,cdr):
    support=[];duration=[];near=[];periods=[]
    q=np.quantile([p['tau_temporal_suvr'] for p in people],[.25,.75])
    for grouping in ['all','diagnosis','tau_group']:
        for p in people:
            group='all' if grouping=='all' else p['diagnosis'] if grouping=='diagnosis' else ('low' if p['tau_temporal_suvr']<=q[0] else 'high' if p['tau_temporal_suvr']>=q[1] else 'middle')
            times=[(d-p['anchor']).days/365.25 for d,v in cdr.get(p['rid'],{}).items() if len(v)==1 and d>max(p['anchor'],p['baseline_date'])]
            for clock in ['tau_PET','baseline_CDR']:
                shift=0 if clock=='tau_PET' else p['baseline_offset_days']/365.25
                for scope,values in [('all_available',times),('analysis_3years',[t for t in times if t<=3])]:
                    if values:duration.append(dict(grouping=grouping,group=group,scope=scope,clock=clock,last_year=max(values)-shift))
                selected=[t-shift for t in times if t<=3]
                for month in range(36 if clock=='tau_PET' else 42):
                    here=[t for t in selected if month/12<t<=(month+1)/12]
                    support.append(dict(grouping=grouping,group=group,clock=clock,month=month+1,people=int(bool(here)),visits=len(here)))
                for label,lower,upper in [('0_to_1',0,1),('1_to_2',1,2),('over_2',2,np.inf)]:
                    here=[t for t in selected if lower<t<=upper]
                    periods.append(dict(grouping=grouping,group=group,clock=clock,period=label,people=int(bool(here)),visits=len(here)))
            for year in [1,2,3]:
                for window in [90,180]:
                    here=[t for t in times if abs(t-year)*365.25<=window]
                    near.append(dict(grouping=grouping,group=group,year=year,window_days=window,clock='tau_PET',people=int(bool(here)),visits=len(here),denominator=1))
    dist=[];d=pd.DataFrame(duration)
    for (grouping,group,scope,clock),g in d.groupby(['grouping','group','scope','clock']):
        dist.append(dict(grouping=grouping,group=group,scope=scope,clock=clock,n=len(g),median=g.last_year.median(),q25=g.last_year.quantile(.25),q75=g.last_year.quantile(.75),maximum=g.last_year.max()))
    support=pd.DataFrame(support).groupby(['grouping','group','clock','month'],as_index=False)[['people','visits']].sum()
    near=pd.DataFrame(near).groupby(['grouping','group','year','window_days','clock'],as_index=False)[['people','visits','denominator']].sum()
    periods=pd.DataFrame(periods).groupby(['grouping','group','clock','period'],as_index=False)[['people','visits']].sum()
    return {'followup_duration':dist,'monthly_visit_support':support,'observed_visit_windows':near,'observed_visit_periods':periods}


def run(root,out,bootstrap=200,cv_bootstrap=100,seed=20261017):
    people,long,cdr=prepare(root)
    common=long.dropna(subset=['centiloids','hippocampus_icv','baseline_cdrsb']).copy()
    common=common.loc[common.diagnosis.ne('Unavailable') & common.mri_match_status.eq('matched')]
    rows,coefficients,diagnostics,status=[],[],[],[]
    chain=[('demographic',inc.DEMOGRAPHIC),('amyloid',inc.DEMOGRAPHIC+' + centiloid10'),
           ('stage',inc.DEMOGRAPHIC+' + centiloid10 + C(diagnosis)'),
           ('MRI',inc.DEMOGRAPHIC+' + centiloid10 + C(diagnosis) + hip1000')]
    for family in ['gaussian','fractional_logit']:
        for label,cov in chain:
            rhs='time_years * (tau01 + '+cov+')'
            try:
                fit,x,d=diag.fit_gee(common,rhs,family=='fractional_logit')
                r=diag.standardized_change(fit,x.design_info,common,family=='fractional_logit',family+'_'+label)
                rows.extend(dict(**z,family=family,adjustment=label,clock='Years since selected baseline CDR',estimand='Marginal modeled change contrast; same complete sample and tau IQR') for z in r)
                diagnostics.append(dict(model=family+'_'+label,formula=rhs,**d))
                ci=fit.conf_int()
                coefficients.extend(dict(model=family+'_'+label,term=t,beta=fit.params[t],ci_low=ci.loc[t].iloc[0],ci_high=ci.loc[t].iloc[1],nominal_p=fit.pvalues[t]) for t in fit.params.index)
            except (ValueError,np.linalg.LinAlgError) as e: status.append(dict(component=family+'_'+label,status='failed',reason=str(e)))
    conditional=common.loc[common.time_years.gt(0)].copy()
    for family in ['gaussian','fractional_logit']:
        rhs='time_years * (tau01 + '+chain[-1][1]+' + baseline_cdrsb)'
        try:
            fit,x,d=diag.fit_gee(conditional,rhs,family=='fractional_logit')
            for year in [1,2]:
                p=conditional.drop_duplicates('rid').copy(); p['visit_offset']=0
                # tau_contrast keeps this fixed target time and baseline CDR constant.
                p['time_years']=year
                r=inc.tau_contrast(fit,x.design_info,p,family=='fractional_logit')
                rows.append(dict(model=family+'_conditional_followup_only',family=family,year=year,
                                 estimate=r['difference_per_01'],ci_low=r['ci_low_per_01'],ci_high=r['ci_high_per_01'],
                                 n_standardized_people=len(p),clock='Years since selected baseline CDR',
                                 estimand='Conditional level contrast holding baseline CDR fixed; not marginal change contrast'))
            diagnostics.append(dict(model=family+'_conditional_followup_only',formula=rhs,**d))
        except (ValueError,np.linalg.LinAlgError) as e: status.append(dict(component=family+'_conditional_followup',status='failed',reason=str(e)))
    try:
        rhs='time_years * (tau01 + '+chain[-1][1]+')'
        x=patsy.dmatrix(rhs,common,return_type='dataframe')
        z=np.column_stack([np.ones(len(common)),common.time_years])
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            fit=sm.MixedLM(common.cdr,x,groups=common.rid,exog_re=z).fit(reml=False,method='lbfgs',maxiter=500,disp=False)
        if not fit.converged or any('Hessian' in str(t.message) or 'singular' in str(t.message) for t in w):
            raise ValueError('Mixed model convergence/Hessian/singularity issue')
        proxy=SimpleNamespace(params=fit.fe_params,cov_params=lambda:fit.cov_params().loc[x.columns,x.columns])
        rows.extend(dict(**r,family='mixed_gaussian',clock='Years since selected baseline CDR',estimand='Fixed-effect change contrast; random intercept and slope sensitivity') for r in diag.standardized_change(proxy,x.design_info,common,False,'mixed_intercept_slope'))
    except (ValueError,np.linalg.LinAlgError) as e: status.append(dict(component='mixed_intercept_slope',status='not_estimable',reason=str(e)))
    endpoint_rows=[]; flow=[]
    for year in [1,2]:
        for window in [90,180]:
            endpoint,f=inc.endpoint_cohort(people,cdr,year,window)
            flow.extend(dict(**r,window_days=window,landmark_days=180) for r in f)
            m,c,t,fail=inc.endpoint_models(endpoint,year)
            endpoint_rows.extend(dict(**r,window_days=window,landmark_days=180,clock='Years since tau PET',estimand='Conditional observed-endpoint level contrast') for r in t)
            status.extend(fail)
            if window==180:
                try:
                    w=observation_bootstrap(people,endpoint,year,bootstrap,seed+year)
                    adni.write(out,'endpoint_observation_weighting_'+str(year),[w])
                except (ValueError,np.linalg.LinAlgError) as e: status.append(dict(component=f'endpoint_{year}_observation_weighting',status='not_estimable',reason=str(e)))
                if year==2 and cv_bootstrap>0:
                    summary,replicates=refit_cv_uncertainty(endpoint,cv_bootstrap,5,seed)
                    adni.write(out,'refitted_cv_bootstrap_summary',summary)
                    adni.write(out,'refitted_cv_bootstrap_replicates',replicates)
        pre,f=inc.endpoint_cohort(people,cdr,year,180,landmark_days=0)
        m,c,t,fail=inc.endpoint_models(pre,year)
        flow.extend(dict(**r,window_days=180,landmark_days=0) for r in f)
        endpoint_rows.extend(dict(**r,window_days=180,landmark_days=0,clock='Years since tau PET',estimand='PET-date information subset; no re-selection of original baseline scans') for r in t)
        status.extend(dict(**r,scope='PET_date_information') for r in fail)
    for name,table in followup_tables(people,cdr).items():adni.write(out,name,table)
    for name,table in [('common_sample_trajectory_contrasts',rows),('trajectory_coefficients',coefficients),('trajectory_diagnostics',diagnostics),
                       ('observed_endpoint_contrasts',endpoint_rows),('endpoint_flow',flow),
                       ('analysis_status',status or [dict(component='requested_estimable_models',status='ok')])]:adni.write(out,name,table)
    spec=dict(common_people=common.rid.nunique(),common_observations=len(common),
              primary_exploratory_trajectory='fractional_logit_MRI on common sample; other repeated models sensitivity',
              primary_new_prediction='Use step1 fractional-logit observed endpoint comparison; not same estimand as trajectories',
              clocks='Trajectory: selected baseline CDR date. Observed endpoints: tau PET date; landmark tau+180d.',
              limitations='All revision-added analyses exploratory; intervals do not encompass full analytic selection. Observation weights require measured-predictor MAR and positivity; not visit-intensity weighting.')
    (out/'analysis_specification.json').write_text(json.dumps(spec,indent=2,default=int),encoding='utf-8')
    print('Longitudinal revision complete;',len(status),'items need inspection',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--bootstrap',type=int,default=200);p.add_argument('--cv-bootstrap',type=int,default=100)
    p.add_argument('--seed',type=int,default=20261017)
    p.add_argument('--descriptions-only',action='store_true')
    a=p.parse_args()
    if a.bootstrap<100 or a.cv_bootstrap<0: p.error('bootstrap must be >=100; cv-bootstrap >=0')
    out=adni.output_dir(a.output_dir)
    if a.descriptions_only:
        people,cdr=prepare_followup_descriptions(a.project_root.resolve())
        for name,table in followup_tables(people,cdr).items():adni.write(out,name,table)
    else:run(a.project_root.resolve(),out,a.bootstrap,a.cv_bootstrap,a.seed)
