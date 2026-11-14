"""Cross-sectional amyloid, regional tau, MRI and observation sensitivities."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import longitudinal as core
import analyze_tau_incremental as inc
import cross_sectional as cross
from pet_thresholds import within_band

COV = 'age_at_tau + sex_male + education_years + apoe4_carrier'
TEMPORAL_NEO = ['FUSIFORM', 'INFERIORTEMPORAL', 'MIDDLETEMPORAL']
MTL = ['ENTORHINAL', 'PARAHIPPOCAMPAL']


def output_dir(path):
    p = Path(path).resolve()
    if p.exists() and any(p.iterdir()):
        raise ValueError('Choose a new empty output directory')
    p.mkdir(parents=True, exist_ok=True)
    return p


def write(out, name, rows):
    pd.DataFrame(rows).to_csv(out / (name + '.csv'), index=False, encoding='utf-8-sig')


def population_description(df):
    """Aggregate scan metadata and fixed-threshold support; no model fitting."""
    descriptions = []
    for group, g in [('all', df)] + list(df.groupby('at_stage')):
        for variable in ['race', 'ethnic', 'site', 'phase', 'amyloid_tracer', 'mri_field', 'mri_qc', 'mri_version']:
            for level, n in g[variable].fillna('Unavailable').value_counts().items():
                descriptions.append(dict(group=group, variable=variable, level=level, n=n, denominator=len(g)))
    tau_values = pd.to_numeric(df.tau_temporal_suvr, errors='coerce')
    if not np.isfinite(tau_values.dropna()).all() or not tau_values.notna().any():
        raise ValueError('No finite tau values for threshold support')
    bins = np.arange(np.floor(tau_values.min()/.05)*.05, np.ceil(tau_values.max()/.05)*.05+.050001, .05)
    hist = []
    for stage, g in df.groupby('at_stage'):
        values = pd.to_numeric(g.tau_temporal_suvr, errors='coerce')
        counts, edges = np.histogram(values.dropna(), bins=bins)
        if counts.sum() != values.notna().sum():
            raise ValueError('Histogram omits observed tau values')
        hist.extend(dict(at_stage=stage, lower=edges[i], upper=edges[i+1], n=int(n), threshold=1.34)
                    for i, n in enumerate(counts))
        for width in [.02, .05, .10]:
            descriptions.append(dict(group=stage, variable='tau_distance_from_1.34', level=f'within_{width}',
                                     n=int(within_band(values, 1.34, width).sum()), denominator=len(g)))
    return pd.DataFrame(descriptions), pd.DataFrame(hist)


def weighted_roi(row, regions):
    v = np.array([pd.to_numeric(row.get('CTX_'+r+'_VOLUME'), errors='coerce') for r in regions], float)
    s = np.array([pd.to_numeric(row.get('CTX_'+r+'_SUVR'), errors='coerce') for r in regions], float)
    if not np.isfinite(v).all() or not np.isfinite(s).all() or (v <= 0).any():
        return np.nan
    return np.sum(v*s)/np.sum(v)


def unique_row(g, signature):
    return g.iloc[0] if len(g) and len(g[signature].drop_duplicates()) == 1 else None


def load_data(root):
    raw = root/'data/private/adni/raw'
    master_path = root/'data/private/adni/adni_pet_aligned_master.csv'
    people, _ = core.load_master(master_path)
    master = pd.read_csv(master_path, dtype={'rid':str})
    master.rid = master.rid.map(core.rid)
    cl, matching = inc.match_centiloids(people, raw/'UCBERKELEY_AMY_6MM.csv')
    for k in ['centiloids', 'amyloid_tracer', 'centiloid_match_status']:
        master[k] = master.rid.map({r:v[k] for r,v in cl.items()})
    mri = pd.read_csv(raw/'UCSFFSX7.csv', low_memory=False)
    mri['rid'] = mri.RID.map(core.rid)
    mri['date'] = inc.dates(mri.EXAMDATE)
    for k in ['ST29SV', 'ST88SV', 'ST10CV']:
        mri[k] = pd.to_numeric(mri[k], errors='coerce')
    mri['ratio'] = (mri.ST29SV+mri.ST88SV)/mri.ST10CV
    mri_groups = dict(tuple(mri.groupby('rid')))
    tau = pd.read_csv(raw/'UCBERKELEY_TAU_6MM.csv', low_memory=False)
    tau['rid'], tau['date'] = tau.RID.map(core.rid), inc.dates(tau.SCANDATE)
    tau = tau.loc[pd.to_numeric(tau.qc_flag, errors='coerce').eq(2) & tau.TRACER.eq('FTP')]
    tau_groups = dict(tuple(tau.groupby(['rid','date'])))
    pvc = pd.read_csv(raw/'UCBERKELEY_TAUPVC_6MM.csv', low_memory=False)
    pvc['rid'], pvc['date'] = pvc.RID.map(core.rid), inc.dates(pvc.SCANDATE)
    pvc = pvc.loc[pvc.TRACER.eq('FTP')]
    pvc_groups = dict(tuple(pvc.groupby(['rid','date'])))
    neo_regions = [k[4:-5] for k in tau if k.startswith('CTX_') and k.endswith('_SUVR') and not k.startswith(('CTX_LH_','CTX_RH_')) and k[4:-5] not in MTL]
    supplemental = []
    for p in people:
        row = dict(rid=p['rid'], hip_mm3=np.nan, icv_mm3=np.nan, mri_date=pd.NaT,
                   mri_field='Unavailable', mri_qc='Unavailable', mri_version='Unavailable',
                   mri_3t_available=False, mri_non3t_only=False, mri_match='no_valid_aligned_mri',
                   tau_pvc=np.nan, tau_entorhinal=np.nan, tau_mtl=np.nan,
                   tau_temporal_neo=np.nan, tau_neocortex=np.nan, site='Unavailable', phase='Unavailable')
        g = mri_groups.get(p['rid'], pd.DataFrame())
        if len(g):
            g = g.loc[g.date.notna() & np.isfinite(g.ratio) & g.ratio.gt(0)].copy()
            g['gap'] = (g.date-pd.Timestamp(p['anchor'])).abs().dt.days
            g = g.loc[g.gap.le(180)].sort_values(['gap','date'])
            row['mri_3t_available'] = bool(g.FIELD_STRENGTH.eq('3T').any())
            row['mri_non3t_only'] = bool(len(g) and not row['mri_3t_available'])
            if len(g):
                c = unique_row(g.loc[g.date.eq(g.iloc[0].date)], ['ST29SV','ST88SV','ST10CV','FIELD_STRENGTH'])
                row['mri_match'] = 'ambiguous_nearest'
                if c is not None and np.isfinite(c[['ST29SV','ST88SV','ST10CV']].astype(float)).all():
                    row.update(mri_match='master_mismatch', mri_field=str(c.FIELD_STRENGTH),
                               mri_qc='Not_recorded' if pd.isna(c.OVERALLQC) else str(c.OVERALLQC), mri_version=str(c.FSVER))
                    if p['hippocampus_icv'] is not None and np.isclose(c.ratio,p['hippocampus_icv'],rtol=1e-9,atol=1e-12) and c.gap==p['hippocampus_icv_days_from_tau'] and (c[['ST29SV','ST88SV','ST10CV']].astype(float)>0).all():
                        row.update(hip_mm3=c.ST29SV+c.ST88SV, icv_mm3=c.ST10CV, mri_date=c.date, mri_match='matched')
        key = (p['rid'], pd.Timestamp(p['anchor']))
        g = tau_groups.get(key, pd.DataFrame())
        if len(g):
            sig = ['LONIUID','META_TEMPORAL_SUVR','SITEID','ORIGPROT']+[k for k in g if k.endswith(('_SUVR','_VOLUME'))]
            c = unique_row(g, list(dict.fromkeys(sig)))
            if c is not None and np.isclose(c.META_TEMPORAL_SUVR,p['tau_temporal_suvr'],rtol=1e-9):
                row.update(tau_entorhinal=c.CTX_ENTORHINAL_SUVR,tau_mtl=weighted_roi(c,MTL),
                           tau_temporal_neo=weighted_roi(c,TEMPORAL_NEO),tau_neocortex=weighted_roi(c,neo_regions),
                           site=str(c.SITEID),phase=str(c.ORIGPROT))
                pg = pvc_groups.get(key,pd.DataFrame())
                if len(pg):
                    pg = pg.loc[pg.LONIUID.astype(str).eq(str(c.LONIUID))]
                    pc = unique_row(pg,['META_TEMPORAL_SUVR'])
                    if pc is not None:
                        row['tau_pvc'] = pd.to_numeric(pc.META_TEMPORAL_SUVR,errors='coerce')
        supplemental.append(row)
    master = master.merge(pd.DataFrame(supplemental),on='rid',validate='one_to_one')
    adsl = pd.read_csv(raw/'ADSL.csv',low_memory=False)
    adsl['rid'] = adsl.SUBJID.map(core.rid)
    if adsl.rid.duplicated().any():
        raise ValueError('ADSL participant keys are not unique')
    for col in ['RACE','ETHNIC']:
        master[col.lower()] = master.rid.map(adsl.set_index('rid')[col]).fillna('Unavailable')
    master['centiloid10'] = master.centiloids/10
    master['tau_per_01'] = master.tau_temporal_suvr/.1
    master['icv_litre'] = master.icv_mm3/1e6
    master['hip_ml'] = master.hip_mm3/1000
    master['amyloid_positive'] = master.at_stage.str.startswith('A+').astype(int)
    return master,people,matching


def fit_effect(df, outcome, rhs, term, weights=None, covariance='HC3'):
    x = patsy.dmatrix(rhs,df,return_type='dataframe',NA_action='raise')
    if len(df)<=x.shape[1]+10 or np.linalg.matrix_rank(x)!=x.shape[1]:
        raise ValueError('Insufficient sample or rank-deficient design')
    if weights is None:
        model = sm.OLS(df[outcome],x)
    else:
        model = sm.WLS(df[outcome],x,weights=weights)
    f = model.fit(cov_type=covariance)
    if term not in x:
        raise ValueError('Requested term absent')
    ci = f.conf_int().loc[term]
    sd = df[outcome].std(ddof=1)
    return dict(n=len(df),beta=float(f.params[term]),ci_low=float(ci.iloc[0]),ci_high=float(ci.iloc[1]),
                p=float(f.pvalues[term]),beta_outcome_sd=float(f.params[term]/sd),
                ci_low_outcome_sd=float(ci.iloc[0]/sd),ci_high_outcome_sd=float(ci.iloc[1]/sd),
                formula=outcome+' ~ '+rhs,condition_number=float(np.linalg.cond(x)))


def smd_table(df, selected, label):
    rows=[]
    fields=cross.COVARIATES+['tau_temporal_suvr','centiloids']
    for category in ['at_stage','race','ethnic','phase']:
        for v in sorted(df[category].dropna().unique()):
            df=df.copy(); key=category+':'+str(v); df[key]=df[category].eq(v).astype(float); fields.append(key)
    for v in fields:
        a=pd.to_numeric(df.loc[selected,v],errors='coerce').dropna()
        b=pd.to_numeric(df.loc[~selected,v],errors='coerce').dropna()
        scale=np.sqrt((a.var(ddof=1)+b.var(ddof=1))/2)
        rows.append(dict(comparison=label,variable=v,selected_n=len(a),other_n=len(b),
                         selected_mean=a.mean(),other_mean=b.mean(),smd=(a.mean()-b.mean())/scale if scale>0 else np.nan))
    return rows


def observation_weights(df, selected, predictors):
    x=patsy.dmatrix(predictors,df,return_type='dataframe',NA_action='raise')
    if np.linalg.matrix_rank(x)!=x.shape[1] or selected.nunique()!=2:
        raise ValueError('Observation model unidentified')
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        f=sm.GLM(selected.astype(int),x,family=sm.families.Binomial()).fit(maxiter=200)
    p=np.asarray(f.predict(x))
    if not f.converged or any('Separation' in type(z.message).__name__ for z in w) or (p<1e-5).any() or (p>1-1e-5).any():
        raise ValueError('Observation model convergence/positivity failure')
    weights=selected.mean()/p[selected]
    cut=np.quantile(weights,[.01,.99]); trimmed=np.clip(weights,*cut)
    return trimmed,dict(n_target=len(df),n_observed=int(selected.sum()),p_min=p.min(),p_max=p.max(),
                        raw_weight_max=weights.max(),trim_low=cut[0],trim_high=cut[1],
                        effective_n=trimmed.sum()**2/(trimmed@trimmed),
                        limitation='Measured-predictor MAR sensitivity in covariate-complete target; not a correction for missing predictors or MNAR; HC3 conditions on estimated weights')


def biological_stage(df, specification):
    req=['mtl_column','neocortex_column','mtl_cutoff','neocortex_moderate_cutoff','neocortex_high_cutoff','source','pipeline_compatibility']
    if any(k not in specification or specification[k] in (None,'') for k in req):
        raise ValueError('Staging needs sourced, pipeline-compatible cutoffs; no default thresholds')
    s=specification; m=df[s['mtl_column']]; n=df[s['neocortex_column']]
    if s['neocortex_high_cutoff']<=s['neocortex_moderate_cutoff']:
        raise ValueError('High cutoff must exceed moderate cutoff')
    result=pd.Series('Unavailable',index=df.index)
    available=m.notna()&n.notna()&df.amyloid_positive.eq(1)
    mt=m.ge(s['mtl_cutoff']); mod=n.ge(s['neocortex_moderate_cutoff']); high=n.ge(s['neocortex_high_cutoff'])
    result.loc[available&~mt&~mod]='A'
    result.loc[available&mt&~mod]='B'
    result.loc[available&mt&mod&~high]='C'
    result.loc[available&mt&high]='D'
    result.loc[available&~mt&mod]='discordant_neocortical_without_MTL'
    result.loc[df.amyloid_positive.eq(0)]='amyloid_negative'
    return result


def run(root,out,stage_config=None,data_derived_staging=False,stage_bootstrap=200):
    df,people,matching=load_data(root)
    rows,failures,missing,balance,weights_out,flow=[],[],[],[],[],[]
    def add(data,y,rhs,term,analysis,family,weights=None):
        try:
            result=fit_effect(data,y,rhs,term,weights)
            if y in ('hippocampus_icv','hip_ml'):
                ref=df.loc[df.at_stage.eq('A-T-')].dropna(subset=cross.COVARIATES+['centiloids'])
                if y=='hippocampus_icv':ref=cross.eligible(ref,y)
                else:ref=cross.eligible(ref,'hippocampus_icv').loc[lambda x:x.mri_match.eq('matched')]
                reference=ref[y].mean()
                result.update(reference_mean=reference,percent_reference=100*result['beta']/reference,
                              percent_reference_ci_low=100*result['ci_low']/reference,
                              percent_reference_ci_high=100*result['ci_high']/reference,
                              percentage_note='Percent of A-T- outcome-eligible CL-complete reference mean, not within-person percentage change')
            rows.append(dict(outcome=y,analysis=analysis,family=family,**result))
        except (ValueError,np.linalg.LinAlgError) as e:
            failures.append(dict(component=analysis,outcome=y,status='not_estimable',reason=str(e)))
    for y in cross.OUTCOMES:
        original=cross.eligible(df,y)
        target=df.loc[df.at_stage.isin(cross.STAGES)]
        balance+=smd_table(target,target.rid.isin(original.rid),y+'_complete_vs_other_three_state')
        balance+=smd_table(df,df.rid.isin(original.rid),y+'_complete_vs_other_all_PET')
        criteria={'outcome_absent':target[y].isna(),'outside_or_missing_180d':~target[y+'_days_from_tau'].between(0,180),
                  'missing_covariate':target[cross.COVARIATES].isna().any(axis=1)}
        if y=='hippocampus_icv': criteria['outlier_flag']=target.hippocampus_outlier_flag.astype(str).str.lower().isin(['true','1','1.0'])
        for name,mask in criteria.items(): missing.append(dict(outcome=y,reason=name,n=int(mask.sum()),denominator=len(target),note='Overlapping reasons; no inferred clinical cause'))
        common=original.dropna(subset=['centiloids'])
        flow.extend(dict(outcome=y,step=k,n=n) for k,n in [('all_PET_aligned',len(df)),('three_confirmatory_AT_states',len(target)),('Aminus_Tplus_excluded_by_design',len(df)-len(target)),('original_complete_case',len(original)),('common_CL_complete_case',len(common)),('missing_CL_in_original_complete_case',len(original)-len(common))])
        for adjusted in (False,True):
            cov=COV+(' + centiloid10' if adjusted else '')
            for ref,stage in [('A-T-','A+T-'),('A-T-','A+T+'),('A+T-','A+T+')]:
                term=f"C(at_stage, Treatment(reference='{ref}'))[T.{stage}]"
                rhs=f"C(at_stage, Treatment(reference='{ref}')) + "+cov
                add(common,y,rhs,term,f'{stage}_vs_{ref}_CL{int(adjusted)}','AT_CL'+str(int(adjusted)))
            a=common.loc[common.amyloid_positive.eq(1)]
            add(a,y,'tau_per_01 + '+cov,'tau_per_01','Aplus_tau_CL'+str(int(adjusted)),'Aplus_CL'+str(int(adjusted)))
        full=df.loc[df[y].notna()&df[y+'_days_from_tau'].between(0,180)].dropna(subset=cross.COVARIATES+['centiloids'])
        if y=='hippocampus_icv': full=full.loc[~full.hippocampus_outlier_flag.astype(str).str.lower().isin(['true','1','1.0'])]
        add(full,y,'tau_per_01 + centiloid10 + '+COV,'tau_per_01','full_cohort_continuous_tau','full_cohort_tau')
        add(full,y,'tau_per_01 * amyloid_positive + centiloid10 + '+COV,'tau_per_01:amyloid_positive','full_cohort_tau_by_A_interaction','full_cohort_interaction')
        regional=common.loc[common.amyloid_positive.eq(1)].dropna(subset=['tau_pvc','tau_entorhinal','tau_mtl','tau_temporal_neo','tau_neocortex'])
        for var in ['tau_temporal_suvr','tau_pvc','tau_entorhinal','tau_mtl','tau_temporal_neo','tau_neocortex']:
            d=regional.copy(); d['exposure']=d[var]/.1
            add(d,y,'exposure + centiloid10 + '+COV,'exposure',var+'_per01','regional_tau')
        site=common.loc[common.amyloid_positive.eq(1)&common.site.ne('Unavailable')]
        counts=site.site.value_counts(); site=site.loc[site.site.isin(counts[counts>=5].index)]
        add(site,y,'tau_per_01 + centiloid10 + '+COV+' + C(site)','tau_per_01','site_fixed_effect_sites_n_ge5','site_sensitivity')
        target=target.dropna(subset=cross.COVARIATES+['centiloids','tau_temporal_suvr']).copy()
        selected=target.rid.isin(common.rid)
        rhs='C(at_stage) + tau_per_01 + centiloid10 + '+COV
        try:
            w,diagnostic=observation_weights(target,selected,rhs)
            observed=target.loc[selected]
            add(observed,y,"C(at_stage, Treatment(reference='A-T-')) + centiloid10 + "+COV,
                "C(at_stage, Treatment(reference='A-T-'))[T.A+T+]",'observation_weighted_A+T+_vs_A-T-','observation_weighting',w)
            weights_out.append(dict(outcome=y,**diagnostic))
        except (ValueError,np.linalg.LinAlgError) as e: failures.append(dict(component='observation_weighting',outcome=y,status='not_estimable',reason=str(e)))
    h=cross.eligible(df,'hippocampus_icv').dropna(subset=['hip_ml','icv_litre','centiloids'])
    for scope,d in [('matched_MRI',h),('matched_3T',h.loc[h.mri_field.eq('3T')]),('QC_pass_only',h.loc[h.mri_qc.eq('Pass')])]:
        a=d.loc[d.amyloid_positive.eq(1)]
        add(a,'hippocampus_icv','tau_per_01 + centiloid10 + '+COV,'tau_per_01',scope+'_ratio','MRI_adjustment')
        add(a,'hip_ml','tau_per_01 + centiloid10 + icv_litre + '+COV,'tau_per_01',scope+'_ICV_regression','MRI_adjustment')
    version_sample=h.loc[h.amyloid_positive.eq(1)&~h.mri_version.isin(['Unavailable','nan',''])]
    for y,extra in [('hippocampus_icv',''),('hip_ml',' + icv_litre')]:
        for adjusted in (False,True):
            rhs='tau_per_01 + centiloid10 + '+COV+extra+(' + C(mri_version)' if adjusted else '')
            add(version_sample,y,rhs,'tau_per_01','FreeSurfer_version_'+('adjusted' if adjusted else 'same_sample_reference'),'MRI_processing_version')
    balance+=smd_table(df,df.mri_3t_available,'aligned_3T_available_vs_other')
    qc=[]
    for (stage,quality,version),g in df.groupby(['at_stage','mri_qc','mri_version']):
        qc.append(dict(at_stage=stage,quality=quality,freesurfer_version=version,n=len(g),
                       n_original_hip_model=len(cross.eligible(g,'hippocampus_icv')),n_CL_complete_matched=int(g.rid.isin(h.rid).sum())))
    write(out,'MRI_QC_coverage',qc)
    descriptive=[]
    from scipy.stats import t as student_t
    for stage,g in df.groupby('at_stage'):
        for y in cross.OUTCOMES:
            d=g.loc[g[y+'_days_from_tau'].between(0,180)]
            if y=='hippocampus_icv':d=d.loc[~d.hippocampus_outlier_flag.astype(str).str.lower().isin(['true','1','1.0'])]
            x=d[y].dropna();n=len(x);mean=x.mean();se=x.std(ddof=1)/np.sqrt(n) if n else np.nan
            margin=student_t.ppf(.975,n-1)*se if n>1 else np.nan
            descriptive.append(dict(at_stage=stage,outcome=y,n=n,mean=mean,sd=x.std(ddof=1),ci_low=mean-margin,ci_high=mean+margin,
                                    interpretation='Unadjusted descriptive mean with nominal t interval, including A-T+; no small-group confirmatory inference'))
    write(out,'four_state_descriptive_intervals',descriptive)
    counts=[]
    for stage,g in df.groupby('at_stage'):
        counts.append(dict(at_stage=stage,n=len(g),aligned_3t=int(g.mri_3t_available.sum()),aligned_non3t_only=int(g.mri_non3t_only.sum()),
                           original_selected_3t=int(g.mri_field.eq('3T').sum()),original_selected_15t=int(g.mri_field.eq('1.5T').sum()),
                           original_hip_model_n=len(cross.eligible(g,'hippocampus_icv'))))
    descriptions,hist=population_description(df)
    config=None
    if stage_config:config=json.loads(Path(stage_config).read_text(encoding='utf-8-sig'))
    elif data_derived_staging:
        try:
            from regional_staging import derive
            config=derive(df,out,stage_bootstrap)
        except ValueError as e:failures.append(dict(component='exploratory_staging',status='not_estimable',reason=str(e)))
    if config:
        df['biological_stage']=biological_stage(df,config)
        write(out,'biological_stage_counts',df.groupby(['at_stage','biological_stage']).size().rename('n').reset_index())
        (out/'biological_stage_specification.json').write_text(json.dumps(config,indent=2),encoding='utf-8')
        for y in cross.OUTCOMES:
            d=cross.eligible(df,y)
            d=d.loc[d.biological_stage.isin(['A','B','C','D'])].dropna(subset=['centiloids'])
            for stage in ['B','C','D']:
                add(d,y,"C(biological_stage, Treatment(reference='A')) + centiloid10 + "+COV,
                    f"C(biological_stage, Treatment(reference='A'))[T.{stage}]",'biological_'+stage+'_vs_A','biological_staging')
    else:
        failures.append(dict(component='AA_biological_staging',status='requires_validated_cutoffs',reason='Regional continuous tau analyses provided; no universal pipeline-independent staging cutoff assumed'))
    result=pd.DataFrame(rows)
    if len(result):
        result['q']=np.nan
        for family,g in result.groupby('family'):
            result.loc[g.index,'q']=multipletests(g.p,method='fdr_bh')[1]
        result['inference']='Revision exploratory; BH within named analysis family; not across all analysis choices'
    for name,x in [('cross_sectional_results',result),('complete_case_comparisons',balance),('missingness_reasons',missing),
                   ('observation_weight_diagnostics',weights_out),('MRI_field_strength',counts),('population_description',descriptions),
                   ('tau_threshold_histogram',hist),('input_matching',matching),('cohort_sample_flow',flow),('analysis_status',failures)]: write(out,name,x)
    spec=dict(roi_temporal_neocortical=TEMPORAL_NEO,roi_MTL=MTL,roi_neocortex='All bilateral CTX parcels excluding ENTORHINAL and PARAHIPPOCAMPAL',
              roi_weight='Volume-weighted; every component required',pvc='Same RID/date/LONIUID as non-PVC QC=2 FTP scan; PVC table lacks own qc_flag',
              original_MRI='Original nearest scan preserved; 3T and pass-only subsets reported separately',
              sources=['https://adni.loni.usc.edu/quick-start-guide-asset/pet.html','https://doi.org/10.1002/alz.13859'])
    (out/'analysis_specification.json').write_text(json.dumps(spec,indent=2),encoding='utf-8')
    print('Clinical sensitivities:',len(result),'models;',len(failures),'items needing inspection',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--stage-config',type=Path)
    p.add_argument('--data-derived-staging',action='store_true',help='Exploratory in-pipeline staging, not externally validated AA stages')
    p.add_argument('--stage-bootstrap',type=int,default=200)
    a=p.parse_args()
    if a.stage_config and a.data_derived_staging:p.error('Choose sourced thresholds or data-derived exploratory staging')
    run(a.project_root.resolve(),output_dir(a.output_dir),a.stage_config,a.data_derived_staging,a.stage_bootstrap)
