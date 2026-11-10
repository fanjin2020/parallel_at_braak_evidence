"""Clinical-stage-adjusted Gaussian and fractional-logit comparisons."""
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels
import longitudinal as core
import longitudinal_sensitivity as diag
VERSION = '2026-09-15.stage-bounded.1'
RHS = diag.RHS + ' + C(diagnosis, Treatment(reference="CN")) * time_years'

def prepare(master, cdr, dx):
    eligible, long, _ = core.make_cohort(master, cdr)
    diagnoses = {}
    for p in eligible:
        d, _ = core.nearest(dx.get(p['rid'], {}), p['anchor'], 180)
        diagnoses[p['rid']] = {1: 'CN', 2: 'MCI', 3: 'Dementia'}[d[1]] if d else 'Unavailable'
    df = diag.design_frame(long)
    df['diagnosis'] = df.rid.map(diagnoses)
    missing = df.loc[df.diagnosis.eq('Unavailable'), 'rid'].nunique()
    df = df.loc[df.diagnosis.ne('Unavailable')].copy().reset_index(drop=True)
    return (df, dict(n_baseline_eligible=len(eligible), n_with_followup=len({p['rid'] for p in long}), n_missing_diagnosis_excluded=int(missing), n_analysis_people=int(df.rid.nunique()), n_analysis_observations=len(df)))

def fit_comparisons(df):
    models = []
    contrasts = []
    coefficients = []
    calibration = []
    for label, rhs, fractional in [('gaussian_original', diag.RHS, False), ('gaussian_stage', RHS, False), ('fractional_logit_original', diag.RHS, True), ('fractional_logit_stage', RHS, True)]:
        fit, x, diagnostics = diag.fit_gee(df, rhs, fractional)
        if not np.isfinite(fit.cov_params()).all().all():
            raise ValueError('Nonfinite robust covariance')
        target = 'time_years:tau01'
        ci = fit.conf_int().loc[target]
        models.append(dict(model=label, **diagnostics, time_tau_coefficient=float(fit.params[target]), coefficient_ci_low=float(ci.iloc[0]), coefficient_ci_high=float(ci.iloc[1]), coefficient_p=float(fit.pvalues[target]), coefficient_scale='logit mean fraction' if fractional else 'CDR-SB points/year per 0.1 SUVR', fractional_mean_not_binomial_trials=fractional, formula='cdr_fraction ~ ' + rhs if fractional else 'cdr ~ ' + rhs))
        contrasts.extend(diag.standardized_change(fit, x.design_info, df, fractional, label))
        ci_all = fit.conf_int()
        coefficients.extend((dict(model=label, term=t, beta=float(fit.params[t]), se=float(fit.bse[t]), ci_low=float(ci_all.loc[t].iloc[0]), ci_high=float(ci_all.loc[t].iloc[1])) for t in fit.params.index))
        work = df.copy()
        work['fitted'] = np.asarray(fit.fittedvalues) * (18 if fractional else 1)
        work['period'] = np.select([work.time_years.eq(0), work.time_years.le(1), work.time_years.le(2)], ['baseline', '(0,1]y', '(1,2]y'], default='>2y')
        for (diagnosis, period), g in work.groupby(['diagnosis', 'period']):
            calibration.append(dict(model=label, diagnosis=diagnosis, period=period, n_people=int(g.rid.nunique()), n_observations=len(g), observed_mean=float(g.cdr.mean()), fitted_mean=float(g.fitted.mean()), mean_residual=float((g.cdr - g.fitted).mean()), fitted_below_zero=int((g.fitted < 0).sum()), fitted_above_18=int((g.fitted > 18).sum())))
    return (models, contrasts, coefficients, calibration)

def run(args):
    root = args.project_root.resolve()
    out = args.output_dir.resolve()
    if out == root or out == root / 'data' or root / 'data' in out.parents:
        raise ValueError('Do not put output into input data directory')
    if out.exists() and any(out.iterdir()):
        raise ValueError('Choose a NEW empty output directory')
    raw = root / 'data/private/adni/raw'
    files = [root / 'data/private/adni/adni_pet_aligned_master.csv', raw / 'CDR.csv', raw / 'DXSUM.csv']
    master, checks = core.load_master(files[0])
    ids = {r['rid'] for r in master}
    cdr, q1 = core.load_visits(files[1], ids, 'cdr')
    dx, q2 = core.load_visits(files[2], ids, 'dx')
    df, flow = prepare(master, cdr, dx)
    models, contrasts, coef, calibration = fit_comparisons(df)
    people = df.drop_duplicates('rid')
    low, high = np.quantile(people.tau_temporal_suvr, [0.25, 0.75])
    support = []
    for name, g in people.groupby('diagnosis'):
        a, b = np.quantile(g.tau_temporal_suvr, [0.05, 0.95])
        support.append(dict(diagnosis=name, n_people=len(g), tau_min=float(g.tau_temporal_suvr.min()), tau_max=float(g.tau_temporal_suvr.max()), tau_p05=float(a), tau_p95=float(b), overall_tau_q25=float(low), overall_tau_q75=float(high), q25_outside_stage_p05_p95=bool(not a <= low <= b), q75_outside_stage_p05_p95=bool(not a <= high <= b)))
    comparison = []
    for year in (1, 2):
        selected = {r['model']: r for r in contrasts if r['year'] == year}
        g = selected['gaussian_stage']['estimate']
        f = selected['fractional_logit_stage']['estimate']
        comparison.append(dict(year=year, gaussian_stage_change=g, fractional_stage_change=f, difference=f - g, same_direction=bool(g * f > 0), note='Descriptive comparison of dependent model estimates; not a test of their difference'))
    out.mkdir(parents=True, exist_ok=True)
    for name, rows in [('stage_bounded_models', models), ('standardized_stage_change', contrasts), ('stage_model_comparison', comparison), ('stage_model_coefficients', coef), ('stage_time_calibration', calibration), ('tau_support_by_diagnosis', support), ('sample_flow', [flow])]:
        core.write_table(out / (name + '.csv'), rows)
    print(pd.DataFrame(models)[['model', 'n_participants', 'fitted_outside_0_18_fraction']].to_string(index=False))
    print(pd.DataFrame(contrasts)[['model', 'year', 'estimate', 'ci_low', 'ci_high']].to_string(index=False))
    print('Aggregate outputs:', out)
if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root', type=Path, required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    run(p.parse_args())
