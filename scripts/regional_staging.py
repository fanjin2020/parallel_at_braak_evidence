"""Outcome-independent exploratory PET staging in the current processing pipeline."""
import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import norm
import json


def mixture_cutoff(values):
    x=np.asarray(values,float);x=x[np.isfinite(x)]
    if len(x)<30 or np.std(x)<1e-5:raise ValueError('Too few varying values for two-component mixture')
    fits=[];floor=max(np.var(x)*1e-3,1e-6)
    for quantiles in ((.2,.8),(.35,.9),(.1,.6)):
        mu=np.quantile(x,quantiles);var=np.repeat(np.var(x),2);pi=np.array([.6,.4]);old=-np.inf
        for iteration in range(1000):
            z=np.column_stack([np.log(pi[j])+norm.logpdf(x,mu[j],np.sqrt(var[j])) for j in (0,1)])
            den=logsumexp(z,axis=1);ll=den.sum();w=np.exp(z-den[:,None]);counts=w.sum(axis=0)
            if (counts<5).any():break
            pi=counts/len(x);mu=(w*x[:,None]).sum(axis=0)/counts
            var=np.maximum((w*(x[:,None]-mu)**2).sum(axis=0)/counts,floor)
            if abs(ll-old)<1e-7:break
            old=ll
        else:continue
        if (counts<5).any():continue
        j=np.argsort(mu);mu,var,pi=mu[j],var[j],pi[j]
        fits.append(dict(cutoff=float(mu[0]+2*np.sqrt(var[0])),n=len(x),mean_lower=mu[0],sd_lower=np.sqrt(var[0]),mean_upper=mu[1],sd_upper=np.sqrt(var[1]),fraction_lower=pi[0],log_likelihood=ll,
                         bic_two=5*np.log(len(x))-2*ll,bic_one=2*np.log(len(x))-2*norm.logpdf(x,x.mean(),x.std()).sum()))
    if not fits:raise ValueError('Mixture did not converge with adequately populated components')
    return max(fits,key=lambda r:r['log_likelihood'])


def thresholds(data):
    m=mixture_cutoff(data.tau_entorhinal);n=mixture_cutoff(data.tau_neocortex)
    elevated=data.loc[data.tau_entorhinal.ge(m['cutoff'])&data.tau_neocortex.ge(n['cutoff'])]
    high=mixture_cutoff(elevated.tau_neocortex)
    if high['cutoff']<=n['cutoff']:raise ValueError('High-neocortical threshold not above moderate threshold')
    return dict(MTL=m,neocortex_moderate=n,neocortex_high=high)


def derive(data,out,bootstrap=200,seed=20261017):
    d=data.dropna(subset=['tau_entorhinal','tau_neocortex'])
    result=thresholds(d);rng=np.random.default_rng(seed);draws=[]
    for b in range(bootstrap):
        try:
            fit=thresholds(d.iloc[rng.integers(0,len(d),len(d))])
            draws.extend(dict(replicate=b+1,region=k,cutoff=v['cutoff'],status='ok') for k,v in fit.items())
        except ValueError as e:draws.append(dict(replicate=b+1,status='failed',reason=str(e)))
    draws=pd.DataFrame(draws).reindex(columns=['replicate','region','cutoff','status','reason']);draws.to_csv(out/'stage_threshold_bootstrap.csv',index=False)
    rows=[]
    for k,v in result.items():
        boot=draws.loc[draws.status.eq('ok') & draws.region.eq(k),'cutoff']
        rows.append(dict(region=k,**v,bootstrap_valid=len(boot),bootstrap_requested=bootstrap,
                         bootstrap_ci_low=boot.quantile(.025) if len(boot)>=.9*bootstrap else np.nan,
                         bootstrap_ci_high=boot.quantile(.975) if len(boot)>=.9*bootstrap else np.nan,
                         mixture_supported_by_BIC=v['bic_two']<v['bic_one']))
    pd.DataFrame(rows).to_csv(out/'stage_thresholds.csv',index=False)
    specification=dict(mtl_column='tau_entorhinal',neocortex_column='tau_neocortex',mtl_cutoff=result['MTL']['cutoff'],
                       neocortex_moderate_cutoff=result['neocortex_moderate']['cutoff'],neocortex_high_cutoff=result['neocortex_high']['cutoff'],
                       source='AA 2024 framework https://doi.org/10.1002/alz.13859; outcome-independent mixture approach adapted from Yushkevich 2024 supplement https://doi.org/10.1002/alz.14161',
                       pipeline_compatibility='Cutoffs newly estimated in the current Berkeley FTP non-PVC table; external numerical thresholds were not transferred',
                       classification='Exploratory AA-inspired operationalization, NOT externally validated AA stages',
                       derivation='Two-component 1D Gaussian mixtures, lower-component mean + 2 SD; high threshold fitted to MTL+ and neocortex+ subset; no clinical outcomes used',
                       ROI='MTL: bilateral entorhinal; neocortex: volume-weighted bilateral CTX parcels excluding entorhinal and parahippocampal; inferior cerebellar reference from existing table',
                       caveats='Different atlas and reference from Yushkevich: their BrainColor/cerebellar-gray thresholds 1.239/1.191/1.457 are incompatible. Conditional outcome CIs omit estimated-threshold uncertainty; unstable bootstrap or weak mixture separation must be reported.')
    (out/'exploratory_stage_config.json').write_text(json.dumps(specification,indent=2),encoding='utf-8')
    return specification
