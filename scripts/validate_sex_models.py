"""Independent arithmetic and complete-family checks for the added sex tests."""
from pathlib import Path
import numpy as np
import pandas as pd
from analysis_utils import dump

MODELS={'baseline','scan_date','joint_CA_joint_quantile','joint_NG_joint_quantile'}

def bh(p):
    p=np.asarray(p,dtype=float)
    if not np.isfinite(p).all() or ((p<0)|(p>1)).any():raise ValueError('Invalid interaction P values')
    order=np.argsort(p,kind='stable');sorted_q=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
    q=np.empty(len(p));q[order]=np.clip(sorted_q,0,1);return q

def validate(work):
    work=Path(work);folder=work/'sex_interaction';checks=[]
    def check(name,passed,value=None):
        checks.append(dict(check=name,passed=bool(passed),value=value))
        if not passed:raise ValueError('Sex-interaction validation failed: '+name)
    status=pd.read_csv(folder/'analysis_status.csv')
    check('Eight complete fits',len(status)==8 and status.status.eq('completed').all(),len(status))
    diag=pd.read_csv(folder/'model_diagnostics.csv')
    check('QC212 and 71 donors in all fits',diag.n_samples.eq(212).all() and diag.n_donors.eq(71).all())
    feature=None
    for name,n in [('feature_sex_interactions',22011),('GO_sex_interactions',4826)]:
        d=pd.read_csv(folder/(name+'.csv'),dtype={'id':str})
        check(name+' complete model set',set(d.model)==MODELS)
        for model,z in d.groupby('model'):
            check(name+' '+model+' full family',len(z)==n and not z.id.duplicated().any(),len(z))
            for prefix in ['female_late','male_late','interaction_late','omnibus']:
                delta=float(np.max(np.abs(z[prefix+'_q']-bh(z[prefix+'_p']))))
                check(name+' '+model+' '+prefix+' BH',delta<1e-10,delta)
            delta=float(np.max(np.abs(z.interaction_late_effect-(z.male_late_effect-z.female_late_effect))))
            check(name+' '+model+' contrast identity',delta<1e-10,delta)
            for prefix in ['female_late','male_late','interaction_late']:
                delta=float(np.max(np.abs(z[prefix+'_effect']-(z[prefix+'_ci_low']+z[prefix+'_ci_high'])/2)))
                check(name+' '+model+' '+prefix+' CI center',delta<1e-10 and (z[prefix+'_ci_high']>=z[prefix+'_ci_low']).all(),delta)
            check(name+' '+model+' omnibus df',z.omnibus_df.eq(3).all())
        if name=='feature_sex_interactions':feature=d
    targets=pd.read_csv(folder/'target_sex_interactions.csv',dtype={'id':str})
    check('16 non-Y targets in each model',len(targets)==64 and targets.groupby('model').size().eq(16).all() and targets.target_symbol.nunique()==16 and not targets.target_symbol.eq('UTY').any())
    merged=targets.merge(feature,on=['id','model'],suffixes=('_target','_full'),validate='one_to_one')
    check('All selected targets found in full family',len(merged)==64)
    for col in ['interaction_late_effect','interaction_late_ci_low','interaction_late_ci_high','interaction_late_p','interaction_late_q','omnibus_F','omnibus_p','omnibus_q']:
        check('Target/full identity: '+col,np.array_equal(merged[col+'_target'],merged[col+'_full']))
    coverage=pd.read_csv(folder/'cohort_coverage.csv')
    check('Eight sex-by-stage cells',len(coverage)==8 and coverage.donors.ge(3).all() and coverage.donors.sum()==71 and coverage.arrays.sum()==212)
    report=dict(passed=True,checks=checks,validation_scope='Independent full-family BH arithmetic, contrasts, confidence bounds, cohort coverage and target extraction; not independent raw-CEL preprocessing.')
    dump(work/'sex_validation.json',report);return report

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--work-dir',type=Path,required=True);a=p.parse_args()
    r=validate(a.work_dir);print('SEX_VALIDATION_OK',len(r['checks']))
