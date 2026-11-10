"""Calendar opportunity, participant resampling and outcome-model sensitivity."""

from __future__ import annotations
import argparse
import json
import math
import warnings
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy.special import expit
from scipy.stats import chi2
import statsmodels
import statsmodels.api as sm
import patsy
import longitudinal as base
VERSION = '2026-09-15.1'
RHS = 'time_years * (tau01 + age10 + sex_male + education_c + apoe4_carrier)'
TARGET = 'time_years:tau01'

def design_frame(long):
    df = pd.DataFrame(long).copy()
    df['tau01'] = (df.tau_temporal_suvr - 1.34) / 0.1
    df['age10'] = (df.age_at_tau - 75) / 10
    df['education_c'] = df.education_years - 16
    return df

def cutoff_date(raw, supplied=None):
    if supplied:
        value = base.iso_date(supplied)
        if value is None:
            raise ValueError('Invalid --as-of; use YYYY-MM-DD')
        return (value, 'explicit_cli_override')
    dates = {base.iso_date(r['data_downloaded_date']) for r in base.read_csv(raw / 'DATA_DOWNLOADED_DATE.csv', ['data_downloaded_date'])}
    if None in dates or len(dates) != 1:
        raise ValueError('Download date missing/ambiguous; supply a verified --as-of YYYY-MM-DD')
    return (dates.pop(), 'DATA_DOWNLOADED_DATE.csv')

def followup_opportunity(eligible, cdr, raw, cutoff):
    """Counts are overlapping evidence, never inferred reasons for loss to follow-up."""
    adsl = {}
    duplicates = set()
    for r in base.read_csv(raw / 'ADSL.csv', ['SUBJID', 'DTHFL', 'DTHDT', 'DTHDTF', 'EOSDT', 'EOSDTF']):
        key = base.rid(r['SUBJID'])
        if key in adsl:
            duplicates.add(key)
        adsl[key] = r
    for key in duplicates:
        adsl.pop(key)
    counters = {'with_followup': Counter(), 'without_followup': Counter()}
    horizon_rows, participant_summaries = ([], [])
    for p in eligible:
        group = 'with_followup' if p['included'] else 'without_followup'
        c = counters[group]
        c['eligible_n'] += 1
        available_days = (cutoff - p['anchor']).days
        c['cutoff_before_tau_pet'] += available_days < 0
        c['calendar_opportunity_less_than_12_months'] += available_days < 365.25
        c['calendar_opportunity_at_least_12_months'] += available_days >= 365.25
        c['calendar_opportunity_at_least_24_months'] += available_days >= 2 * 365.25
        c['calendar_opportunity_at_least_36_months'] += available_days >= 3 * 365.25
        visits = cdr.get(p['rid'], {})
        future_dates = [d for d, v in visits.items() if len(v) == 1 and d > p['anchor']]
        c['any_valid_cdr_beyond_3y'] += any(((d - p['anchor']).days > 3 * 365.25 for d in future_dates))
        c['any_valid_cdr_after_download_date'] += any((d > cutoff for d in visits))
        a = adsl.get(p['rid'])
        c['adsl_matched_unique'] += a is not None
        c['adsl_unmatched_or_duplicate'] += a is None
        if a:
            death_flag = a['DTHFL'].strip().lower() in ('yes', 'y', '1', '1.0')
            c['death_flag_present_timing_not_assumed'] += death_flag
            c['death_flag_blank_unknown_not_alive'] += a['DTHFL'].strip() == ''
            for field, flag in [('DTHDT', 'DTHDTF'), ('EOSDT', 'EOSDTF')]:
                d = base.iso_date(a[field])
                marked = bool(a[flag].strip())
                c[field + '_date_recorded'] += d is not None
                c[field + '_date_flag_nonblank'] += marked
                c[field + '_invalid_nonblank_date'] += bool(a[field].strip()) and d is None
                if d is not None and (not marked):
                    c[field + '_unflagged_on_or_before_pet'] += d <= p['anchor']
                    c[field + '_unflagged_after_pet_within_3y'] += 0 < (d - p['anchor']).days <= 3 * 365.25
                    c[field + '_unflagged_before_12_months'] += 0 < (d - p['anchor']).days < 365.25
                c[field + '_date_after_download'] += d is not None and d > cutoff
            eos = base.iso_date(a['EOSDT'])
            if eos is not None and (not a['EOSDTF'].strip()):
                c['valid_cdr_after_unflagged_EOSDT'] += any((d > eos for d in visits))
        participant_summaries.append(dict(group=group, available_days=available_days, included=p['included']))
    for months in (12, 24, 36):
        people = [p for p in participant_summaries if p['available_days'] >= months / 12 * 365.25]
        horizon_rows.append(dict(calendar_opportunity_months=months, n=len(people), with_primary_qualifying_followup=sum((p['included'] for p in people)), without_primary_qualifying_followup=sum((not p['included'] for p in people))))
    records = [dict(group=g, check=k, n=v) for g, c in counters.items() for k, v in sorted(c.items())]
    return (records, horizon_rows)

def cluster_sufficient_stats(x, y, groups):
    unique, codes = np.unique(groups, return_inverse=True)
    xx = np.stack([x[codes == i].T @ x[codes == i] for i in range(len(unique))])
    xy = np.stack([x[codes == i].T @ y[codes == i] for i in range(len(unique))])
    return (xx, xy)

def checked_solve(a, b):
    if np.linalg.matrix_rank(a) != a.shape[0]:
        return None
    value = np.linalg.solve(a, b)
    return value if np.isfinite(value).all() else None

def cluster_resampling(x, y, groups, target_index, reference_beta, bootstrap, seed):
    """Exact OLS point estimates = Gaussian independence GEE point estimates.

    Resample participants, retaining their entire cluster including its size.
    Repeated sampled participants contribute multiplicities to X'X and X'y.
    """
    xx, xy = cluster_sufficient_stats(x, y, groups)
    total_xx, total_xy = (xx.sum(axis=0), xy.sum(axis=0))
    full = checked_solve(total_xx, total_xy)
    if full is None or not np.isclose(full[target_index], reference_beta, atol=1e-10, rtol=1e-09):
        raise ValueError('OLS sufficient statistics failed to reproduce primary GEE point estimate')
    loo = []
    for i in range(len(xx)):
        beta = checked_solve(total_xx - xx[i], total_xy - xy[i])
        if beta is not None:
            loo.append(float(beta[target_index]))
    if not loo:
        raise ValueError('All leave-one-participant-out fits failed')
    loo = np.asarray(loo)
    change = np.abs(loo - reference_beta)
    summary = [dict(method='leave_one_participant_out', n_requested=len(xx), n_valid=len(loo), reference_beta=reference_beta, min_beta=float(loo.min()), max_beta=float(loo.max()), median_beta=float(np.median(loo)), max_absolute_change=float(change.max()), max_relative_change=float(change.max() / abs(reference_beta)) if reference_beta else '', sign_reversal_count=int((loo * reference_beta < 0).sum()), note='diagnostic only; no participants removed from primary model')]
    rng = np.random.default_rng(seed)
    estimates = []
    for b in range(bootstrap):
        weights = np.bincount(rng.integers(0, len(xx), size=len(xx)), minlength=len(xx))
        beta = checked_solve(np.einsum('i,ijk->jk', weights, xx), weights @ xy)
        estimates.append(dict(replicate=b + 1, status='ok' if beta is not None else 'rank_or_numeric_failure', beta=float(beta[target_index]) if beta is not None else ''))
        if (b + 1) % 500 == 0:
            print(f'Participant bootstrap: {b + 1}/{bootstrap}', flush=True)
    good = np.array([e['beta'] for e in estimates if e['status'] == 'ok'])
    if not len(good):
        raise ValueError('All participant bootstrap estimates failed')
    ci = np.quantile(good, [0.025, 0.975])
    summary.append(dict(method='participant_pairs_bootstrap', n_requested=bootstrap, n_valid=len(good), reference_beta=reference_beta, bootstrap_se=float(good.std(ddof=1)) if len(good) > 1 else '', percentile_ci_low=float(ci[0]), percentile_ci_high=float(ci[1]), fraction_positive=float((good > 0).mean()), note='percentile interval; fraction_positive is NOT a p-value; unstratified clusters'))
    return (summary, estimates)

def fit_gee(df, rhs, fractional=False):
    x = patsy.dmatrix(rhs, df, return_type='dataframe')
    if np.linalg.matrix_rank(x.to_numpy()) != x.shape[1]:
        raise ValueError('Rank-deficient sensitivity design')
    y = df.cdr.to_numpy() / (18 if fractional else 1)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        model = sm.GEE(y, x, groups=df.rid.to_numpy(), family=sm.families.Binomial() if fractional else sm.families.Gaussian(), cov_struct=sm.cov_struct.Independence())
        fit = model.fit(maxiter=200, cov_type='robust')
    if not fit.converged or not np.isfinite(fit.params).all() or (not np.isfinite(fit.bse).all()):
        raise ValueError('Sensitivity fit did not converge to finite estimates')
    diag = dict(n_participants=int(df.rid.nunique()), n_observations=len(df), converged=True, condition_number=float(np.linalg.cond(x)), warning_types=';'.join(sorted({type(w.message).__name__ for w in caught})))
    fitted = fit.fittedvalues.to_numpy() if hasattr(fit.fittedvalues, 'to_numpy') else np.asarray(fit.fittedvalues)
    fitted = fitted * (18 if fractional else 1)
    diag.update(fitted_min=float(fitted.min()), fitted_max=float(fitted.max()), fitted_outside_0_18_fraction=float(((fitted < 0) | (fitted > 18)).mean()), rmse_in_sample=float(np.sqrt(np.mean((df.cdr.to_numpy() - fitted) ** 2))))
    return (fit, x, diag)

def robust_joint_test(fit, terms, label):
    idx = [list(fit.params.index).index(t) for t in terms]
    b = np.asarray(fit.params)[idx]
    v = np.asarray(fit.cov_params())[np.ix_(idx, idx)]
    if np.linalg.matrix_rank(v) < len(idx):
        return dict(test=label, status='singular_covariance', df=len(idx), terms=';'.join(terms))
    w = float(b @ np.linalg.solve(v, b))
    return dict(test=label, status='ok', df=len(idx), wald_chi2=w, p_value=float(chi2.sf(w, len(idx))), terms=';'.join(terms))

def standardized_change(fit, info, df, fractional, label):
    """IQR tau contrast in modeled change since baseline, normalized per 0.1 SUVR.

    Fractional result is not a constant annual rate. Delta-method intervals.
    Standardize over one record per participant, not repeated-visit frequencies.
    """
    people = df.drop_duplicates('rid').copy()
    low, high = np.quantile(people.tau_temporal_suvr, [0.25, 0.75])
    units = (high - low) / 0.1
    if units <= 0:
        raise ValueError('No tau IQR variation')
    result = []
    b, v = (np.asarray(fit.params), np.asarray(fit.cov_params()))
    for year in (1, 2):
        estimate = 0.0
        gradient = np.zeros_like(b)
        for t, tau, sign in [(year, high, 1), (0, high, -1), (year, low, -1), (0, low, 1)]:
            prediction = people.copy()
            prediction['time_years'] = t
            prediction['tau01'] = (tau - 1.34) / 0.1
            x = np.asarray(patsy.build_design_matrices([info], prediction)[0])
            if fractional:
                mu = expit(x @ b)
                estimate += sign * 18 * mu.mean() / units
                gradient += sign * 18 * ((mu * (1 - mu))[:, None] * x).mean(axis=0) / units
            else:
                estimate += sign * (x @ b).mean() / units
                gradient += sign * x.mean(axis=0) / units
        se = math.sqrt(max(0.0, float(gradient @ v @ gradient)))
        result.append(dict(model=label, year=year, estimate=float(estimate), ci_low=float(estimate - 1.96 * se), ci_high=float(estimate + 1.96 * se), tau_q25=float(low), tau_q75=float(high), n_standardized_people=len(people), n_people_with_observed_span_at_least_year=int((df.groupby('rid').time_years.max() >= year).sum()), units='CDR-SB modeled change contrast since baseline; IQR tau contrast normalized per 0.1 SUVR', note='not a time-constant annual slope; in-sample standardized association'))
    return result

def calibration_rows(df, fit, label, fractional=False):
    work = df.copy()
    work['fitted'] = np.asarray(fit.fittedvalues) * (18 if fractional else 1)
    work['period'] = np.select([work.time_years.eq(0), work.time_years.le(1), work.time_years.le(2)], ['baseline', '(0,1]y', '(1,2]y'], default='>2y')
    output = []
    for (stage, period), g in work.groupby(['at_stage', 'period'], observed=True):
        output.append(dict(model=label, stage=stage, period=period, n_observations=len(g), n_participants=int(g.rid.nunique()), observed_mean=float(g.cdr.mean()), fitted_mean=float(g.fitted.mean()), observed_zero_fraction=float(g.cdr.eq(0).mean()), mean_residual=float((g.cdr - g.fitted).mean())))
    return output

def run(args):
    root, out = (args.project_root.resolve(), args.output_dir.resolve())
    if out.exists() and any(out.iterdir()):
        raise ValueError('Choose a NEW empty output directory; existing results are not overwritten')
    if out == root or out == root / 'data' or root / 'data' in out.parents:
        raise ValueError('Do not put diagnostics into original input data directory')
    if args.bootstrap < 100:
        raise ValueError('--bootstrap must be >=100; use 2000 for the manuscript run')
    raw = root / 'data/private/adni/raw'
    master, checks = base.load_master(root / 'data/private/adni/adni_pet_aligned_master.csv')
    cdr, quality = base.load_visits(raw / 'CDR.csv', {r['rid'] for r in master}, 'cdr')
    eligible, long, flow = base.make_cohort(master, cdr)
    df = design_frame(long)
    cutoff, cutoff_source = cutoff_date(raw, args.as_of)
    opportunity, horizons = followup_opportunity(eligible, cdr, raw, cutoff)
    primary, x, diag = fit_gee(df, RHS)
    beta = float(primary.params[TARGET])
    resampling, bootstrap = cluster_resampling(x.to_numpy(), df.cdr.to_numpy(), df.rid.to_numpy(), list(x.columns).index(TARGET), beta, args.bootstrap, args.seed)
    model_diags = [dict(model='primary_gaussian_linear_time', status='ok', **diag)]
    tests = []
    contrasts = standardized_change(primary, x.design_info, df, False, 'primary_gaussian_linear_time')
    calibration = calibration_rows(df, primary, 'primary_gaussian_linear_time')
    coefficients = []
    failures = []
    for label, rhs, fractional in [('gaussian_quadratic_time', RHS + ' + I(time_years ** 2) + I(time_years ** 2):tau01', False), ('fractional_logit_linear_time', RHS, True)]:
        try:
            fit, xx, d = fit_gee(df, rhs, fractional)
            model_diags.append(dict(model=label, status='ok', **d))
            ci = fit.conf_int()
            coefficients.extend((dict(model=label, term=t, beta=float(fit.params[t]), se=float(fit.bse[t]), ci_low=float(ci.loc[t, 0]), ci_high=float(ci.loc[t, 1]), p_value=float(fit.pvalues[t]), scale='logit_mean_CDR_fraction' if fractional else 'CDR-SB') for t in fit.params.index))
            if not fractional:
                tests.append(robust_joint_test(fit, ['I(time_years ** 2)', 'I(time_years ** 2):tau01'], 'added_quadratic_time_terms'))
                tests.append(robust_joint_test(fit, [TARGET, 'I(time_years ** 2):tau01'], 'tau_time_association_in_quadratic_model'))
            contrasts.extend(standardized_change(fit, xx.design_info, df, fractional, label))
            calibration.extend(calibration_rows(df, fit, label, fractional))
        except (ValueError, np.linalg.LinAlgError, FloatingPointError) as e:
            failures.append(dict(component=label, error_type=type(e).__name__, status='failed'))
            model_diags.append(dict(model=label, status='failed', error_type=type(e).__name__))
    opportunity_effects = []
    for months in (12, 24, 36):
        chosen = [p for p in eligible if (cutoff - p['anchor']).days >= months / 12 * 365.25]
        ids = {p['rid'] for p in chosen if p['included']}
        subset = df.loc[df.rid.isin(ids)].copy()
        row = dict(min_calendar_opportunity_months=months, eligible_with_baseline_n=len(chosen), analyzed_n=len(ids), without_qualifying_followup_n=sum((not p['included'] for p in chosen)))
        if len(ids) < 30:
            row.update(status='insufficient_participants')
        else:
            f, _, _ = fit_gee(subset, RHS)
            ci = f.conf_int().loc[TARGET]
            row.update(status='ok', beta=float(f.params[TARGET]), ci_low=float(ci.iloc[0]), ci_high=float(ci.iloc[1]), p_value=float(f.pvalues[TARGET]))
        opportunity_effects.append(row)
    out.mkdir(parents=True, exist_ok=True)
    tables = {'followup_opportunity': opportunity, 'calendar_opportunity_summary': horizons, 'calendar_opportunity_sensitivity': opportunity_effects, 'cluster_resampling_summary': resampling, 'model_diagnostics': model_diags, 'quadratic_time_tests': tests, 'standardized_change_contrasts': contrasts, 'sensitivity_model_coefficients': coefficients, 'time_bin_calibration': calibration}
    for name, table in tables.items():
        base.write_table(out / (name + '.csv'), table)
    paths = [root / 'data/private/adni/adni_pet_aligned_master.csv', raw / 'CDR.csv', raw / 'ADSL.csv', raw / 'DATA_DOWNLOADED_DATE.csv', Path(base.__file__), Path(__file__)]
    print('Completed:', out, flush=True)
    print(json.dumps(resampling, indent=2), flush=True)
    if failures or resampling[1]['n_valid'] / args.bootstrap < 0.95:
        raise SystemExit(2)

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root', type=Path, required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--bootstrap', type=int, default=2000)
    p.add_argument('--seed', type=int, default=20260915)
    p.add_argument('--as-of', default=None, help='Verified administrative cutoff override, YYYY-MM-DD; otherwise download-date file')
    run(p.parse_args())
if __name__ == '__main__':
    main()
