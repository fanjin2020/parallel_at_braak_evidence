"""Covariate-only MI for the original A/T contrasts and A-positive tau slopes."""
import argparse
import json
import platform
import inspect
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy.stats import t as t_dist, norm
import statsmodels
import statsmodels.api as sm
from statsmodels.imputation.mice import MICEData
from statsmodels.stats.multitest import multipletests

COV=['age_at_tau','sex_male','education_years','apoe4_carrier']
STAGES=['A-T-','A+T-','A+T+']
OUTCOMES=['hippocampus_icv','cdrsb']

def rid(x):
    return str(int(float(x))) if pd.notna(x) else None

def validate_master(d):
    if d.rid.isna().any() or d.rid.duplicated().any():raise ValueError('Duplicate/missing participant IDs')
    if not d.at_stage.isin(STAGES+['A-T+']).all():raise ValueError('Invalid A/T classification')
    tau_date=pd.to_datetime(d.tau_pet_date,errors='coerce')
    amy_date=pd.to_datetime(d.amy_pet_date,errors='coerce')
    if tau_date.isna().any() or amy_date.isna().any():raise ValueError('Missing/invalid selected PET dates')
    gap=(tau_date-amy_date).dt.days.abs()
    if not gap.eq(d.at_pair_days).all() or not gap.between(0,180).all():raise ValueError('PET date/alignment mismatch')
    if not np.isfinite(d.tau_temporal_suvr).all():raise ValueError('Missing/nonfinite selected tau SUVR')
    if not d.tau_temporal_suvr.ge(1.34).eq(d.at_stage.str.endswith('T+')).all():raise ValueError('Tau threshold and T classification disagree')
    for k in ['sex_male','apoe4_carrier']:
        if not d[k].dropna().isin([0,1]).all():raise ValueError('Invalid binary covariate: '+k)
    for k in COV+OUTCOMES:
        if not np.isfinite(d[k].dropna()).all():raise ValueError('Nonfinite observed values: '+k)
    if not d.cdrsb.dropna().between(0,18).all():raise ValueError('CDR-SB outside 0-18')

def add_centiloid(d,path):
    raw=pd.read_csv(path,usecols=['RID','SCANDATE','qc_flag','TRACER','AMYLOID_STATUS_COMPOSITE_REF','CENTILOIDS'],low_memory=False)
    raw['rid']=raw.RID.map(rid);raw['date']=pd.to_datetime(raw.SCANDATE,errors='coerce').dt.normalize()
    for c in ['qc_flag','AMYLOID_STATUS_COMPOSITE_REF','CENTILOIDS']:raw[c]=pd.to_numeric(raw[c],errors='coerce')
    groups=dict(tuple(raw.loc[raw.qc_flag.eq(2)&raw.date.notna()].groupby(['rid','date'])))
    values=[];status=[]
    for _,p in d.iterrows():
        date=pd.to_datetime(p.amy_pet_date,errors='coerce')
        g=groups.get((p.rid,date.normalize())) if pd.notna(date) else None
        value=np.nan;s='no_exact_passed_scan'
        if g is not None:
            u=g[['CENTILOIDS','TRACER','AMYLOID_STATUS_COMPOSITE_REF']].drop_duplicates()
            if len(u)!=1:s='ambiguous_same_day_scan'
            elif u.iloc[0].AMYLOID_STATUS_COMPOSITE_REF!=int(p.at_stage.startswith('A+')):s='classification_mismatch'
            elif not np.isfinite(u.iloc[0].CENTILOIDS):s='centiloid_missing'
            else:value=float(u.iloc[0].CENTILOIDS);s='matched'
        values.append(value);status.append(s)
    d=d.copy();d['centiloid10']=np.asarray(values)/10
    return d,pd.Series(status).value_counts().rename_axis('status').reset_index(name='n')

def eligible(d,outcome,model,adjusted):
    keep=d.at_stage.isin(STAGES)&d.at_pair_days.between(0,180)&d[outcome+'_days_from_tau'].between(0,180)&d[outcome].notna()
    if outcome=='hippocampus_icv':keep&=~d.hippocampus_outlier_flag.astype(str).str.lower().isin(['true','1','1.0'])
    if model=='within_Aplus':keep&=d.at_stage.isin(['A+T-','A+T+'])
    keep&=d.tau_temporal_suvr.notna()
    if adjusted:keep&=d.centiloid10.notna()
    z=d.loc[keep].copy().reset_index(drop=True)
    z['tau_per_01']=z.tau_temporal_suvr/.1
    z['group_AplusTminus']=z.at_stage.eq('A+T-').astype(float)
    z['group_AplusTplus']=z.at_stage.eq('A+T+').astype(float)
    return z

def fit(z,y,model,adjusted):
    exposure=['group_AplusTminus','group_AplusTplus'] if model=='AT' else ['tau_per_01']
    columns=exposure+COV+(['centiloid10'] if adjusted else [])
    X=sm.add_constant(z[columns],has_constant='add').astype(float)
    if len(z)<=len(columns)+10 or np.linalg.matrix_rank(X)!=X.shape[1]:raise ValueError('Rank-deficient or under-supported design')
    f=sm.OLS(z[y],X,missing='raise').fit(cov_type='HC3')
    specs={'A+T-_vs_A-T-':{'group_AplusTminus':1},'A+T+_vs_A-T-':{'group_AplusTplus':1},'A+T+_vs_A+T-':{'group_AplusTplus':1,'group_AplusTminus':-1}} if model=='AT' else {'tau_per_0.1_in_Aplus':{'tau_per_01':1}}
    result={}
    for name,weights in specs.items():
        v=np.array([weights.get(k,0) for k in X.columns]);q=float(v@f.params);u=float(v@f.cov_params()@v)
        if not np.isfinite(q) or not np.isfinite(u) or u<=0:raise ValueError('Invalid HC3 contrast')
        se=np.sqrt(u);result[name]=dict(beta=q,variance=u,se=se,ci_low=q-1.959963984540054*se,ci_high=q+1.959963984540054*se,p=2*norm.sf(abs(q/se)),df_complete=float(f.df_resid))
    return result

def pool(estimates,variances,df_complete):
    m=len(estimates)
    if m<2:raise ValueError('At least two imputations required')
    q=float(np.mean(estimates));u=float(np.mean(variances));b=float(np.var(estimates,ddof=1));total=u+(1+1/m)*b
    if not np.isfinite(total) or total<=0:raise ValueError('Invalid pooled variance')
    lam=(1+1/m)*b/total
    old=(m-1)/lam**2 if lam>0 else np.inf
    observed=(df_complete+1)/(df_complete+3)*df_complete*(1-lam)
    df=1/(1/old+1/observed) if observed>0 else np.nan
    if not np.isfinite(df) or df<=0:raise ValueError('Invalid Barnard-Rubin degrees of freedom')
    se=np.sqrt(total);critical=t_dist.ppf(.975,df);mcse=np.sqrt(b/m)
    return dict(beta=q,se=se,ci_low=q-critical*se,ci_high=q+critical*se,p=2*t_dist.sf(abs(q/se),df),
        within_variance=u,between_variance=b,total_variance=total,df=df,missing_information_lambda=lam,
        beta_mcse=mcse,mcse_over_se=mcse/se,mcse_review=mcse/se>.1)

def family_adjust(rows):
    df=pd.DataFrame(rows)
    if df.empty:return df
    df['family']=np.where(df.contrast.isin(['A+T-_vs_A-T-','A+T+_vs_A-T-']),'four_AT_contrasts','two_direct_contrasts_and_two_Aplus_slopes')
    df['planned_family_size']=4;df['q']=np.nan
    for _,idx in df.groupby(['adjustment','estimator','family']).groups.items():
        p=df.loc[idx,'p'].to_numpy(float)
        if len(p)>4:raise ValueError('Duplicate planned family member')
        q=multipletests(np.r_[np.nan_to_num(p,nan=1.),np.ones(4-len(p))],method='fdr_bh')[1][:len(p)]
        df.loc[idx,'q']=np.where(np.isfinite(p),q,np.nan)
    return df

def chain_diagnostics(values, pooled_se):
    """Descriptive chain checks and batch-mean Monte Carlo error, not a convergence proof."""
    v=np.asarray(values,float)
    if v.ndim!=2 or not np.isfinite(v).all():raise ValueError('Invalid chain draws')
    c,n=v.shape;flat=v.ravel();naive=float(np.std(flat,ddof=1)/np.sqrt(flat.size))
    within=float(np.mean(np.var(v,axis=1,ddof=1))) if n>1 else np.nan
    between=float(n*np.var(v.mean(axis=1),ddof=1)) if c>1 else np.nan
    rhat=float(np.sqrt(((n-1)/n*within+between/n)/within)) if c>1 and n>1 and within>0 else np.nan
    batch_size=max(2,int(np.sqrt(n)));nb=n//batch_size;batch_se=np.nan
    if nb>=2:
        means=v[:,:nb*batch_size].reshape(c,nb,batch_size).mean(axis=2)
        batch_se=float(np.sqrt(np.sum(np.var(means,axis=1,ddof=1)/nb)/c**2))
    mcse=max(naive,batch_se) if np.isfinite(batch_se) else naive
    ratio=mcse/pooled_se
    review=c<2 or n<20 or not np.isfinite(rhat) or rhat>1.1 or ratio>.1
    return dict(chains=c,draws_per_chain=n,rhat_descriptive=rhat,beta_mcse_iid=naive,
        beta_mcse_batch=batch_se,beta_mcse=mcse,mcse_over_se=ratio,mcse_review=bool(review))

def new_mice(x,missing,variables,seed):
    # statsmodels 0.14 uses NumPy's global RNG; newer versions accept an explicit RNG.
    kwargs=dict(perturbation_method='boot',k_pmm=5)
    if 'rng' in inspect.signature(MICEData).parameters:kwargs['rng']=np.random.default_rng(seed)
    else:np.random.seed(seed)
    mice=MICEData(x,**kwargs)
    for k in missing:
        predictors=[v for v in variables if v!=k and x[v].nunique(dropna=True)>1]
        mice.set_imputer(k,formula=' + '.join(predictors),model_class=sm.OLS,perturbation_method='boot',k_pmm=5)
    return mice

def assert_observed_preserved(completed,original):
    for k in original:
        mask=original[k].notna()
        np.testing.assert_allclose(completed.loc[mask,k],original.loc[mask,k],rtol=0,atol=0)

def current_run_reference(root):
    """Collect independent fits from this run for a numerical MI-implementation audit."""
    root=Path(root)
    primary_path='results/adni/adni_results.csv'
    extension_path='results/extensions/cross_sectional/extension_results.csv'
    clinical_path='revision/clinical/cross_sectional_results.csv'
    primary=pd.read_csv(root/primary_path)
    extension=pd.read_csv(root/extension_path)
    clinical=pd.read_csv(root/clinical_path)
    rows=[]
    def one(frame,where):
        if len(frame)!=1:raise ValueError('Expected one independent reference row: '+where)
        return frame.iloc[0]
    def append(row,y,model,adjustment,contrast,source,original=False):
        rows.append(dict(outcome=y,model=model,adjustment=adjustment,contrast=contrast,
            n=row['n_complete_case'] if original else row['n'],beta=row['beta'],
            ci_low=row['ci95_lower'] if 'ci95_lower' in row else row['ci_low'],
            ci_high=row['ci95_upper'] if 'ci95_upper' in row else row['ci_high'],
            ci_critical_value=1.96 if original else 1.959963984540054,
            reference_source=source))
    for y in OUTCOMES:
        for contrast in ['A+T-_vs_A-T-','A+T+_vs_A-T-']:
            row=one(primary.loc[primary.analysis.eq('confirmatory_180day')&
                primary.outcome.eq(y)&primary.contrast.eq(contrast.replace('_vs_',' vs '))],y+' '+contrast)
            append(row,y,'AT','original_covariates',contrast,primary_path,original=True)
        for analysis,model,contrast in [
            ('direct_A+T+_vs_A+T-','AT','A+T+_vs_A+T-'),
            ('continuous_tau_per_0.1_in_Apositive','within_Aplus','tau_per_0.1_in_Aplus')]:
            row=one(extension.loc[extension.outcome.eq(y)&extension.population.eq('all_eligible')&
                extension.analysis.eq(analysis)],y+' '+analysis)
            append(row,y,model,'original_covariates',contrast,extension_path)
        for analysis,model,contrast in [
            ('A+T-_vs_A-T-_CL1','AT','A+T-_vs_A-T-'),
            ('A+T+_vs_A-T-_CL1','AT','A+T+_vs_A-T-'),
            ('A+T+_vs_A+T-_CL1','AT','A+T+_vs_A+T-'),
            ('Aplus_tau_CL1','within_Aplus','tau_per_0.1_in_Aplus')]:
            row=one(clinical.loc[clinical.outcome.eq(y)&clinical.analysis.eq(analysis)],y+' '+analysis)
            append(row,y,model,'Centiloid',contrast,clinical_path)
    result=pd.DataFrame(rows)
    if len(result)!=16 or result.duplicated(['outcome','model','adjustment','contrast']).any():
        raise ValueError('Independent complete-case reference must contain 16 unique contrasts')
    return result

def check_reference(reference,expected,y,model,adjusted,n):
    rows=[]
    subset=expected.loc[expected.outcome.eq(y)&expected.model.eq(model)&expected.adjustment.eq('Centiloid' if adjusted else 'original_covariates')]
    if subset.contrast.duplicated().any() or set(subset.contrast)!=set(reference):raise ValueError('Reference contrast coverage mismatch or duplicate contrast')
    if 'ci_critical_value' not in subset:raise ValueError('Reference must declare ci_critical_value; do not infer historical precision from a failed fit')
    for r in subset.itertuples():
        v=reference[r.contrast]
        critical=float(r.ci_critical_value)
        if critical not in (1.96,1.959963984540054):raise ValueError('Unsupported declared reference CI critical value')
        if critical==1.96 and not (model=='AT' and not adjusted and r.contrast in ['A+T-_vs_A-T-','A+T+_vs_A-T-']):
            raise ValueError('Rounded 1.96 critical value is documented only for the original four A/T contrasts')
        if not np.isfinite([r.n,r.beta,r.ci_low,r.ci_high]).all() or r.n!=int(r.n) or r.ci_high<=r.ci_low:
            raise ValueError('Invalid reference count or interval')
        reference_se=(r.ci_high-r.ci_low)/(2*critical)
        comparable_low=v['beta']-critical*v['se'];comparable_high=v['beta']+critical*v['se']
        checks=dict(n_match=int(r.n)==n,
            beta_match=np.isclose(v['beta'],r.beta,rtol=1e-6,atol=1e-10),
            hc3_se_match=np.isclose(v['se'],reference_se,rtol=1e-6,atol=1e-12),
            reference_center_match=np.isclose((r.ci_low+r.ci_high)/2,r.beta,rtol=1e-6,atol=1e-10),
            ci_low_match=np.isclose(comparable_low,r.ci_low,rtol=1e-6,atol=1e-10),
            ci_high_match=np.isclose(comparable_high,r.ci_high,rtol=1e-6,atol=1e-10))
        rows.append(dict(outcome=y,model=model,adjustment=r.adjustment,contrast=r.contrast,n_observed=n,n_reference=r.n,
            beta_observed=v['beta'],beta_reference=r.beta,hc3_se_observed=v['se'],hc3_se_reference=reference_se,
            reference_source=getattr(r,'reference_source','supplied audit reference'),
            reference_se_source='reference CI width / (2 * declared critical value)',ci_critical_value=critical,
            ci_low_reference=r.ci_low,ci_high_reference=r.ci_high,ci_low_exact_normal=v['ci_low'],ci_high_exact_normal=v['ci_high'],
            ci_low_reference_scale=comparable_low,ci_high_reference_scale=comparable_high,
            ci_low_difference=comparable_low-r.ci_low,ci_high_difference=comparable_high-r.ci_high,
            **{k:bool(value) for k,value in checks.items()},passed=bool(all(checks.values())),
            reason='; '.join(k for k,value in checks.items() if not value)))
    return rows

def analyze(d,out,imputations=100,burnin=20,thin=5,seed=20261024,chains=4,expected=None):
    if chains<1 or imputations%chains or imputations//chains<2:raise ValueError('Imputations must divide evenly into chains with at least two retained draws per chain')
    out.mkdir(parents=True,exist_ok=True);rows=[];diagnostics=[];draws=[];states=[];traces=[];cohorts=[];warning_rows=[];checks=[];label=None
    def save():
        family_adjust(rows).to_csv(out/'MI_and_complete_case_results.csv',index=False)
        pd.DataFrame(states).to_csv(out/'analysis_status.csv',index=False)
        pd.DataFrame(diagnostics).to_csv(out/'covariate_diagnostics.csv',index=False)
        pd.DataFrame(draws).to_csv(out/'imputation_estimates.csv',index=False)
        pd.DataFrame(traces).to_csv(out/'chain_trace.csv',index=False)
        pd.DataFrame(cohorts).to_csv(out/'cohort_counts.csv',index=False)
        pd.DataFrame(warning_rows,columns=['component','chain','message']).to_csv(out/'fitting_warnings.csv',index=False)
        pd.DataFrame(checks).to_csv(out/'complete_case_reference_check.csv',index=False)
        (out/'progress.json').write_text(json.dumps(dict(completed_targets=len(states),planned_targets=8,active_target=label)),encoding='utf-8')
    job=0
    for adjusted in [False,True]:
      for model in ['AT','within_Aplus']:
       for y in OUTCOMES:
        job+=1;label=f'{model}_{y}_'+('Centiloid' if adjusted else 'original_covariates')
        print(f'[{job}/8] {label}',flush=True)
        try:
            z=eligible(d,y,model,adjusted)
            fixed=[y,'tau_per_01']+(['group_AplusTminus','group_AplusTplus'] if model=='AT' else [])+(['centiloid10'] if adjusted else [])
            if model=='AT' and set(z.at_stage)!=set(STAGES):raise ValueError('An A/T comparison group is absent')
            variables=COV+fixed;x=z[variables].astype(float)
            if not all(np.isfinite(x[k].dropna()).all() for k in x):raise ValueError('Nonfinite input')
            missing=[k for k in COV if x[k].isna().any()]
            if any(x[k].notna().sum()<30 for k in missing):raise ValueError('Fewer than 30 observed donors for an imputed covariate')
            cc=z.dropna(subset=variables);reference=fit(cc,y,model,adjusted)
            if expected is not None:
                checked=check_reference(reference,expected,y,model,adjusted,len(cc));checks.extend(checked)
                if not all(k['passed'] for k in checked):raise ValueError('Complete-case fit differs from its independent reference; inspect the numerical audit before interpreting MI')
            context=dict(outcome=y,model=model,adjustment='Centiloid' if adjusted else 'original_covariates',n_target=len(z),n_complete=len(cc),n_imputed_people=int(x[COV].isna().any(axis=1).sum()),target='Outcome-observed PET-aligned '+('three-state' if model=='AT' else 'A-positive')+' population')
            context['missing_covariate_fraction']=context['n_imputed_people']/len(z)
            for stage in STAGES:
                cohorts.append(dict(component=label,at_stage=stage,n_target=int(z.at_stage.eq(stage).sum()),n_complete=int(cc.at_stage.eq(stage).sum())))
            for name,v in reference.items():rows.append(dict(**context,contrast=name,estimator='complete_case',imputations=0,**v))
            for k in variables:diagnostics.append(dict(component=label,variable=k,n=len(x),missing=int(x[k].isna().sum()),observed_min=x[k].min(),observed_max=x[k].max()))
            if not missing:
                for name,v in reference.items():rows.append(dict(**context,contrast=name,estimator='MI',imputations=0,**v))
                states.append(dict(component=label,status='not_needed',reason='No missing covariates; MI row repeats complete-case estimate'));save();continue
            effects={k:[] for k in reference};variances={k:[] for k in reference};dfc=None
            per_chain=imputations//chains
            for chain in range(chains):
              mice=new_mice(x,missing,variables,seed+job*10000+chain)
              with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always')
                for cycle in range(1,burnin+thin*per_chain+1):
                    mice.update_all(1);completed=mice.data
                    assert_observed_preserved(completed,x)
                    for k in ['sex_male','apoe4_carrier']:
                        if not completed[k].isin([0,1]).all():raise ValueError('PMM did not preserve binary support')
                    for k in missing:
                        imp=completed.loc[x[k].isna(),k]
                        traces.append(dict(component=label,chain=chain+1,cycle=cycle,phase='burnin' if cycle<=burnin else 'sampling',variable=k,imputed_mean=imp.mean(),imputed_sd=imp.std(),imputed_min=imp.min(),imputed_max=imp.max()))
                    if cycle<=burnin or (cycle-burnin)%thin:continue
                    draw=(cycle-burnin)//thin;j=chain*per_chain+draw
                    estimates=fit(completed,y,model,adjusted)
                    for name,v in estimates.items():
                        effects[name].append(v['beta']);variances[name].append(v['variance']);dfc=v['df_complete']
                        draws.append(dict(component=label,imputation=j,chain=chain+1,draw=draw,contrast=name,beta=v['beta'],variance=v['variance']))
                    for k in missing:diagnostics.append(dict(component=label,variable=k,imputation=j,chain=chain+1,completed_mean=completed[k].mean(),completed_sd=completed[k].std(),imputed_mean=completed.loc[x[k].isna(),k].mean()))
                    if j%5==0:print(f'  {j}/{imputations} imputations; chain {chain+1}/{chains}',flush=True);save()
                warning_rows.extend(dict(component=label,chain=chain+1,message=str(w.message)) for w in caught)
            for name in reference:
                pooled=pool(effects[name],variances[name],dfc)
                pooled.update(chain_diagnostics(np.array(effects[name]).reshape(chains,per_chain),pooled['se']))
                rows.append(dict(**context,contrast=name,estimator='MI',imputations=imputations,**pooled))
            states.append(dict(component=label,status='completed',reason='All requested imputations fitted; covariates only'))
        except Exception as e:
            states.append(dict(component=label,status='failed',reason=f'{type(e).__name__}: {e}'))
        save()
    (out/'methods.txt').write_text('Outcome-observed, PET-aligned targets; original and Centiloid-complete populations are separate.\nOnly age, sex, education and APOE covariates are imputed using bootstrap predictive mean matching. Outcome, tau, A/T classification, Centiloid and absent MRI are not imputed.\nA/T group indicators and observed outcomes are included in the three-state imputation models. Within-A+ slopes use A-positive imputation targets. Original exploratory direct contrasts are obtained from three-state fits.\nHC3 within-imputation covariance; Rubin pooling with Barnard-Rubin finite-sample degrees of freedom. BH families contain the planned four A/T tests or two direct contrasts plus two slopes, separately by adjustment and estimator. These are revision sensitivities, not new confirmatory tests.\nMAR assumptions and imputation-model adequacy remain necessary. MCSE is reported; this does not resolve outcome missingness, MRI missingness, mortality or MNAR.\n',encoding='utf-8')
    with (out/'methods.txt').open('a',encoding='utf-8') as f:
        f.write(f'{chains} separately seeded chains; {imputations} total retained imputations; {burnin} burn-in cycles per chain; {thin} cycles between retained draws.\n')
        f.write('All observed values are checked unchanged after every cycle. Traces include burn-in and retained-sampling periods. Classical between-chain R-hat and batch-mean MCSE are descriptive checks, not proofs of convergence; short chains or MCSE/SE > 0.1 are flagged.\n')
        f.write('The four-test families used here are stated explicitly and need not reproduce the q values of differently grouped historical sensitivity tables. Current complete-case estimates use exact-normal intervals; MI pooling uses finite-sample t intervals.\n')
        if expected is not None:
            f.write('The complete-case audit compares n, beta, HC3 SE and interval endpoints against independent reference estimates, using their declared critical value. Reference SE is reconstructed from the interval width. Reference estimates are not used to fit or pool models.\n')
        else:
            f.write('No historical reference estimates were supplied. Complete-case models were refitted from current-run data; numerical agreement with an archived manuscript was not tested by this step.\n')
    return states

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--project-root',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    audit=p.add_mutually_exclusive_group()
    audit.add_argument('--reference-estimates',type=Path,help='Optional audit reference; never used to fit models')
    audit.add_argument('--audit-current-run',action='store_true',help='Audit against independently fitted A/T and tau models from this run')
    p.add_argument('--imputations',type=int,default=100);p.add_argument('--chains',type=int,default=4);p.add_argument('--burnin',type=int,default=20);p.add_argument('--thin',type=int,default=5);p.add_argument('--seed',type=int,default=20261024)
    a=p.parse_args()
    if a.chains<1 or a.imputations%a.chains or a.imputations//a.chains<2 or a.burnin<1 or a.thin<1:p.error('Invalid MI iteration counts; use at least two retained draws per chain')
    if a.output_dir.exists() and any(a.output_dir.iterdir()):p.error('Choose a new output directory')
    d=pd.read_csv(a.project_root/'data/private/adni/adni_pet_aligned_master.csv',dtype={'rid':str})
    d.rid=d.rid.map(rid)
    for k in COV+['tau_temporal_suvr','at_pair_days',*OUTCOMES,*[y+'_days_from_tau' for y in OUTCOMES]]:d[k]=pd.to_numeric(d[k],errors='raise')
    validate_master(d)
    d,matching=add_centiloid(d,a.project_root/'data/private/adni/raw/UCBERKELEY_AMY_6MM.csv')
    a.output_dir.mkdir(parents=True,exist_ok=True);matching.to_csv(a.output_dir/'centiloid_matching_counts.csv',index=False)
    expected=current_run_reference(a.project_root) if a.audit_current_run else (pd.read_csv(a.reference_estimates) if a.reference_estimates else None)
    (a.output_dir/'reference_audit_status.json').write_text(json.dumps(dict(
        requested=expected is not None,
        role='Independent numerical check only; reference values are not fitted-model inputs',
        source='current-run fits' if a.audit_current_run else ('supplied reference' if expected is not None else None),
        status='requested' if expected is not None else 'not_requested'),indent=2),encoding='utf-8')
    states=analyze(d,a.output_dir,a.imputations,a.burnin,a.thin,a.seed,a.chains,expected)
    if expected is not None:
        audit_rows=pd.read_csv(a.output_dir/'complete_case_reference_check.csv')
        (a.output_dir/'reference_audit_status.json').write_text(json.dumps(dict(requested=True,
            source='current-run fits' if a.audit_current_run else 'supplied reference',
            checks=len(audit_rows),passed=int(audit_rows.passed.astype(str).str.lower().eq('true').sum()),
            status='completed' if len(audit_rows)==len(expected) and audit_rows.passed.astype(str).str.lower().eq('true').all() else 'failed'),indent=2),encoding='utf-8')
    (a.output_dir/'environment.json').write_text(json.dumps(dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,statsmodels=statsmodels.__version__,imputations=a.imputations,chains=a.chains,burnin=a.burnin,thin=a.thin,seed=a.seed),indent=2),encoding='utf-8')
    if any(s['status']=='failed' for s in states):raise SystemExit('One or more MI targets failed; inspect analysis_status.csv')

if __name__=='__main__':main()
