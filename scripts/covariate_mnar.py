"""Predefined pattern-mixture stress tests for missing covariates only."""
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.special import expit,logit
import mi_matched_estimands as mi
from analysis_utils import csv,dump,checkpoint,NotEstimable

BINARY={'sex_male','apoe4_carrier'}

def scenarios(missing):
    result=[dict(name='MAR',shifts={})]
    for variable in missing:
        deltas=[-np.log(2),np.log(2)] if variable in BINARY else [-1.,-.5,.5,1.]
        for delta in deltas:
            result.append(dict(name=variable+('_logOR_' if variable in BINARY else '_residualSD_')+str(round(delta,6)),shifts={variable:delta}))
    if len(missing)>1:
        for sign in [-1,1]:
            result.append(dict(name='joint_'+('lower' if sign<0 else 'higher'),shifts={k:sign*(np.log(2) if k in BINARY else .5) for k in missing}))
    return result

def sensitivity_parameters(completed,original,missing):
    result={}
    for k in missing:
        mask=original[k].notna();columns=[c for c in completed if c!=k and completed.loc[mask,c].nunique()>1]
        X=sm.add_constant(completed[columns],has_constant='add').astype(float)
        if np.linalg.matrix_rank(X.loc[mask])!=X.shape[1]:raise NotEstimable('Sensitivity covariate model rank deficient: '+k)
        if k in BINARY:
            with warnings.catch_warnings(record=True) as ws:
                warnings.simplefilter('always');f=sm.GLM(completed.loc[mask,k],X.loc[mask],family=sm.families.Binomial()).fit(maxiter=200)
            if not f.converged or any('Separation' in type(w.message).__name__ for w in ws):raise NotEstimable('Sensitivity binary model separation: '+k)
            p=np.asarray(f.predict(X));result[k]=np.clip(p,1e-6,1-1e-6)
        else:
            f=sm.OLS(completed.loc[mask,k],X.loc[mask]).fit()
            sd=float(np.sqrt(f.ssr/f.df_resid))
            if not np.isfinite(sd) or sd<=0:raise NotEstimable('Invalid residual SD: '+k)
            result[k]=sd
    return result

def shifted(completed,original,parameters,shifts,uniforms):
    z=completed.copy();clipped=0
    for k,delta in shifts.items():
        missing=original[k].isna().to_numpy();values=z[k].to_numpy().copy()
        if k in BINARY:
            p=parameters[k];changed=expit(logit(p)+delta)
            flip=(changed-p)/(1-p) if delta>0 else (p-changed)/p
            candidates=missing & (values==(0 if delta>0 else 1)) & (uniforms[k]<np.clip(flip,0,1))
            values[candidates]=1 if delta>0 else 0
        else:
            values[missing]+=delta*parameters[k]
            # Bounds are the observed target-sample covariate support, not outcome-based cutoffs.
            lower,upper=original[k].min(),original[k].max()
            clipped+=int(((values[missing]<lower)|(values[missing]>upper)).sum())
            values[missing]=np.clip(values[missing],lower,upper)
        z[k]=values
    mi.assert_observed_preserved(z,original)
    for k in BINARY:
        if not z[k].isin([0,1]).all():raise AssertionError('Binary support changed')
    return z,clipped

def run(root,out,imputations=100,seed=20261029,diagnostic=False):
    root=Path(root);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(root/'data/private/adni/adni_pet_aligned_master.csv',dtype={'rid':str});d.rid=d.rid.map(mi.rid);mi.validate_master(d)
    d,_=mi.add_centiloid(d,root/'data/private/adni/raw/UCBERKELEY_AMY_6MM.csv')
    chains=4;per=imputations//chains
    if imputations%4 or per<2:raise ValueError('Imputations must be divisible by four, at least eight')
    pooled=[];all_draws=[];states=[];scenario_specs=[];job=0
    for adjusted in [False,True]:
      for model in ['AT','within_Aplus']:
       for y in mi.OUTCOMES:
        job+=1;label=f'{model}_{y}_{int(adjusted)}';print('Covariate MNAR',job,'/8',label,flush=True)
        z=mi.eligible(d,y,model,adjusted)
        fixed=[y,'tau_per_01']+(['group_AplusTminus','group_AplusTplus'] if model=='AT' else [])+(['centiloid10'] if adjusted else [])
        original=z[mi.COV+fixed].astype(float);missing=[k for k in mi.COV if original[k].isna().any()]
        specs=scenarios(missing)
        for s in specs:scenario_specs.append(dict(component=label,scenario=s['name'],shifts=str(s['shifts']),n_target=len(z),missing_fields=';'.join(missing)))
        if not missing:
            estimates=mi.fit(z,y,model,adjusted)
            for c,v in estimates.items():pooled.append(dict(component=label,outcome=y,model=model,adjustment='Centiloid' if adjusted else 'original_covariates',scenario='MAR',estimator='MI',contrast=c,n=len(z),imputations=0,**v))
            states.append(dict(component=label,status='not_needed',reason='No missing covariates'));continue
        if any(original[k].notna().sum()<30 for k in missing):
            states.append(dict(component=label,status='not_estimable',reason='Fewer than 30 observed values for a missing covariate'));continue
        draws=[];failed=False
        for chain in range(chains):
            def compute_chain():
                mice=mi.new_mice(original,missing,list(original),seed+job*10000+chain)
                mice.update_all(5 if diagnostic else 20);rows=[]
                for draw in range(per):
                    mice.update_all(2 if diagnostic else 5);completed=mice.data.copy();mi.assert_observed_preserved(completed,original)
                    parameters=sensitivity_parameters(completed,original,missing)
                    rng=np.random.default_rng(seed+job*100000+chain*1000+draw)
                    uniforms={k:rng.random(len(original)) for k in missing}
                    for s in specs:
                        work,nclip=shifted(completed,original,parameters,s['shifts'],uniforms)
                        for contrast,v in mi.fit(work,y,model,adjusted).items():
                            rows.append(dict(chain=chain+1,draw=draw+1,scenario=s['name'],contrast=contrast,n_clipped=nclip,**v))
                    if (draw+1)%5==0:print(f'  chain {chain+1}/4 draw {draw+1}/{per}',flush=True)
                return rows
            try:draws.extend(checkpoint(out,label+'_chain'+str(chain+1),compute_chain))
            except (NotEstimable,np.linalg.LinAlgError) as e:
                states.append(dict(component=label,status='not_estimable',reason=str(e)));failed=True;break
        if failed:continue
        draws=pd.DataFrame(draws);draws['component']=label;all_draws.append(draws)
        for (scenario,contrast),v in draws.groupby(['scenario','contrast']):
            result=mi.pool(v.beta.to_numpy(),v.variance.to_numpy(),float(v.df_complete.min()))
            matrix=v.sort_values(['chain','draw']).beta.to_numpy().reshape(chains,per)
            result.update(mi.chain_diagnostics(matrix,result['se']))
            pooled.append(dict(component=label,outcome=y,model=model,adjustment='Centiloid' if adjusted else 'original_covariates',scenario=scenario,estimator='MI',contrast=contrast,n=len(z),imputations=imputations,mean_clipped=v.n_clipped.mean(),**result))
        states.append(dict(component=label,status='completed',reason='Assumption-based covariate sensitivity, not outcome/dropout/mortality imputation'))
        csv(out/'analysis_status.csv',states);csv(out/'scenario_definitions.csv',scenario_specs)
    result=pd.DataFrame(pooled)
    if len(result):
        # Variable-specific scenarios may be inactive in some targets. MAR fills only those inactive scenario/target combinations.
        expanded=[]
        for scenario in sorted(result.scenario.unique()):
            for component,group in result.groupby('component'):
                chosen=group.loc[group.scenario.eq(scenario)]
                if chosen.empty:chosen=group.loc[group.scenario.eq('MAR')].assign(scenario=scenario,scenario_inactive=True)
                else:chosen=chosen.assign(scenario_inactive=False)
                expanded.append(chosen)
        result=pd.concat(expanded,ignore_index=True)
        result=pd.concat([mi.family_adjust(group) for _,group in result.groupby('scenario')],ignore_index=True)
        mar=result.loc[result.scenario.eq('MAR'),['component','contrast','beta']].rename(columns={'beta':'MAR_beta'})
        result=result.merge(mar,on=['component','contrast'],validate='many_to_one')
        result['direction_changed']=np.sign(result.beta)!=np.sign(result.MAR_beta)
        result['interval_contains_zero']=result.ci_low.le(0)&result.ci_high.ge(0)
    csv(out/'covariate_MNAR_results.csv',result);csv(out/'scenario_definitions.csv',scenario_specs);csv(out/'analysis_status.csv',states)
    csv(out/'draw_diagnostics.csv',pd.concat(all_draws,ignore_index=True) if all_draws else [])
    dump(out/'specification.json',dict(imputations=imputations,chains=4,diagnostic_only=diagnostic,
        interpretation='Post-imputation pattern-mixture stress scenarios; offsets do not feed back into the MAR FCS chain. They are not estimated from data or externally established as plausible.',
        continuous='Missing values shifted by +/-0.5 or +/-1 observed conditional residual SD, bounded by observed covariate support; clipping reported.',
        binary='Missing binary draws odds-tilted using an observed-data logistic conditional probability, OR 0.5 or 2, preserving observed values.',
        multiplicity='Separate four-test BH families per scenario and adjustment; no correction across the sensitivity-scenario search. Report the complete grid.',
        limitation='Outcome-observed targets only. No missing MRI/outcomes/death status imputed; this does not solve longitudinal MNAR.'))
