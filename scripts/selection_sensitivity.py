"""Follow-up availability weighting with full-process participant bootstrap."""
from pathlib import Path
from datetime import datetime, timezone
import argparse, json, warnings
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
import longitudinal as core
VERSION = '2026-09-15.selection.1'
RHS = 'time_years * (tau01 + age10 + sex_male + education_c + apoe4_carrier)'
SEL = ['age_at_tau', 'sex_male', 'education_years', 'apoe4_carrier', 'tau_temporal_suvr', 'baseline_cdrsb', 'baseline_offset_days', 'opportunity_years']

def selection_weights(people):
    x = people[SEL].astype(float)
    sd = x.std(ddof=0)
    if (sd <= 0).any() or not np.isfinite(x).all().all():
        raise ValueError('Invalid or constant selection predictors')
    x = sm.add_constant((x - x.mean()) / sd)
    if np.linalg.matrix_rank(x) != x.shape[1]:
        raise ValueError('Selection model rank deficient')
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter('always')
        fit = sm.GLM(people.included.astype(int), x, family=sm.families.Binomial()).fit(maxiter=200)
    if not fit.converged or any(('Separation' in type(w.message).__name__ for w in ws)):
        raise ValueError('Selection model did not converge or has separation')
    p = np.asarray(fit.predict(x))
    if not np.isfinite(p).all() or p.min() <= 1e-06 or p.max() >= 1 - 1e-06:
        raise ValueError('Near-deterministic selection probabilities; do not silently clip')
    selected = people.included.to_numpy(bool)
    w = selected.mean() / p[selected]
    lo, hi = np.quantile(w, [0.01, 0.99])
    return (p, w, np.clip(w, lo, hi), fit)

def outcome(long, weights=None, stage=False):
    df = long.copy().reset_index(drop=True)
    df['tau01'] = (df.tau_temporal_suvr - 1.34) / 0.1
    df['age10'] = (df.age_at_tau - 75) / 10
    df['education_c'] = df.education_years - 16
    formula = 'cdr ~ ' + RHS
    if stage:
        df = df.loc[df.diagnosis.isin(['CN', 'MCI', 'Dementia'])].copy()
        formula += ' + C(diagnosis, Treatment(reference="CN")) * time_years'
    if df.rid.nunique() < 30:
        raise ValueError('Insufficient outcome participants')
    w = None if weights is None else df.rid.map(weights).to_numpy(float)
    model = smf.gee(formula, groups='rid', data=df, weights=w, family=sm.families.Gaussian(), cov_struct=sm.cov_struct.Independence())
    if np.linalg.matrix_rank(model.exog) != model.exog.shape[1]:
        raise ValueError('Outcome model rank deficient')
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter('always')
        fit = model.fit(maxiter=200, cov_type='robust')
    if not fit.converged or not np.isfinite(fit.params).all() or (not np.isfinite(fit.bse).all()):
        raise ValueError('Invalid outcome model')
    term = 'time_years:tau01'
    ci = fit.conf_int().loc[term]
    return dict(n_participants=int(df.rid.nunique()), n_observations=len(df), beta=float(fit.params[term]), ci_low=float(ci.iloc[0]), ci_high=float(ci.iloc[1]), p_value=float(fit.pvalues[term]), condition_number=float(np.linalg.cond(model.exog)), fitted_outside_0_18_fraction=float(((fit.fittedvalues < 0) | (fit.fittedvalues > 18)).mean()), warnings=';'.join(sorted({type(w.message).__name__ for w in ws})), formula=formula)

def ess(w):
    w = np.asarray(w)
    return float(w.sum() ** 2 / (w @ w))

def probability_summary(probabilities, included):
    p = np.asarray(probabilities, dtype=float)
    selected = np.asarray(included, dtype=bool)
    if p.ndim != 1 or p.shape != selected.shape or not np.isfinite(p).all() or ((p < 0) | (p > 1)).any():
        raise ValueError('Invalid selection probabilities or inclusion mask')
    rows = []
    for cohort, mask in [('all_eligible', np.ones(len(p), dtype=bool)),
                         ('with_followup', selected), ('without_followup', ~selected)]:
        v = p[mask]
        rows.append(dict(cohort=cohort, n=len(v),
                         minimum=float(v.min()) if len(v) else '',
                         median=float(np.median(v)) if len(v) else '',
                         maximum=float(v.max()) if len(v) else ''))
    return rows

def run(args):
    root = args.project_root.resolve()
    out = args.output_dir.resolve()
    if out == root or out == root / 'data' or root / 'data' in out.parents:
        raise ValueError('Output must not be in input data directory')
    if out.exists() and any(out.iterdir()):
        raise ValueError('Choose a NEW empty output directory')
    raw = root / 'data/private/adni/raw'
    files = [root / 'data/private/adni/adni_pet_aligned_master.csv', raw / 'CDR.csv', raw / 'DXSUM.csv', raw / 'DATA_DOWNLOADED_DATE.csv']
    dates = {core.iso_date(r['data_downloaded_date']) for r in core.read_csv(files[3], ['data_downloaded_date'])}
    if None in dates or len(dates) != 1:
        raise ValueError('Ambiguous download cutoff')
    cutoff = dates.pop()
    master, checks = core.load_master(files[0])
    ids = {r['rid'] for r in master}
    cdr, q1 = core.load_visits(files[1], ids, 'cdr')
    dx, q2 = core.load_visits(files[2], ids, 'dx')
    eligible, long, _ = core.make_cohort(master, cdr)
    for p in eligible:
        d, _ = core.nearest(dx.get(p['rid'], {}), p['anchor'], 180)
        p['diagnosis'] = {1: 'CN', 2: 'MCI', 3: 'Dementia'}[d[1]] if d else 'Unavailable'
        p['opportunity_years'] = min(3.0, (cutoff - p['anchor']).days / 365.25)
        if p['opportunity_years'] <= 0:
            raise ValueError('No calendar follow-up opportunity')
    people = pd.DataFrame(eligible)
    df = pd.DataFrame(long)
    df['diagnosis'] = df.rid.map(people.set_index('rid').diagnosis)
    results = [dict(analysis='unweighted', **outcome(df))]
    prob, w, wt, fit = selection_weights(people)
    selected = people.loc[people.included].copy()
    weight_summaries = []
    balances = []
    for label, weight in [('ipw_untruncated', w), ('ipw_truncated_1_99', wt)]:
        mapping = dict(zip(selected.rid, weight))
        results.append(dict(analysis=label, **outcome(df, mapping), inference_note='GEE CI treats estimated weights as fixed; use full-process bootstrap too'))
        weight_summaries.append(dict(analysis=label, n=len(weight), ess=ess(weight), minimum=float(weight.min()), median=float(np.median(weight)), maximum=float(weight.max()), p01=float(np.quantile(weight, 0.01)), p99=float(np.quantile(weight, 0.99)), n_changed_by_truncation=int(np.sum(w != weight))))
    for name in SEL:
        target = people[name].mean()
        scale = people[name].std(ddof=1)
        balances.append(dict(variable=name, target_mean=target, unweighted_mean=selected[name].mean(), weighted_mean=np.average(selected[name], weights=w), truncated_weighted_mean=np.average(selected[name], weights=wt), unweighted_standardized_difference=(selected[name].mean() - target) / scale, weighted_standardized_difference=(np.average(selected[name], weights=w) - target) / scale, truncated_standardized_difference=(np.average(selected[name], weights=wt) - target) / scale))
    results.append(dict(analysis='clinical_stage_adjusted', **outcome(df, stage=True), inference_note='Conditional clinical-stage association; not the same estimand as unadjusted model'))
    out.mkdir(parents=True, exist_ok=True)
    core.write_table(out / 'selection_probability_summary.csv', probability_summary(prob, people.included))
    boot = []
    rng = np.random.default_rng(args.seed)
    groups = {rid: g for rid, g in df.groupby('rid')}
    for b in range(args.bootstrap):
        sampled = people.iloc[rng.integers(0, len(people), len(people))].copy().reset_index(drop=True)
        blocks = []
        for j, p in sampled.iterrows():
            old = p['rid']
            sampled.loc[j, 'rid'] = str(j)
            if p['included']:
                block = groups[old].copy()
                block['rid'] = str(j)
                blocks.append(block)
        try:
            _, bw, bwt, _ = selection_weights(sampled)
            bl = pd.concat(blocks, ignore_index=True)
            chosen = sampled.loc[sampled.included, 'rid']
            for label, ww in [('ipw_untruncated', bw), ('ipw_truncated_1_99', bwt)]:
                try:
                    boot.append(dict(replicate=b + 1, analysis=label, status='ok', beta=outcome(bl, dict(zip(chosen, ww)))['beta']))
                except (ValueError, np.linalg.LinAlgError) as e:
                    boot.append(dict(replicate=b + 1, analysis=label, status=type(e).__name__))
        except (ValueError, np.linalg.LinAlgError) as e:
            for label in ['ipw_untruncated', 'ipw_truncated_1_99']:
                boot.append(dict(replicate=b + 1, analysis=label, status=type(e).__name__))
        if (b + 1) % 100 == 0:
            print('Full-process bootstrap', b + 1, '/', args.bootstrap, flush=True)
    bs = []
    for label in ['ipw_untruncated', 'ipw_truncated_1_99']:
        values = [r['beta'] for r in boot if r['analysis'] == label and r['status'] == 'ok']
        valid = len(values) >= max(200, 0.95 * args.bootstrap)
        row = dict(analysis=label, requested=args.bootstrap, valid=len(values), interval_reportable=bool(valid))
        if valid:
            row.update(ci_low=float(np.quantile(values, 0.025)), ci_high=float(np.quantile(values, 0.975)))
        bs.append(row)
    for name, rows in [('selection_sensitivity_results', results), ('weight_diagnostics', weight_summaries), ('balance_against_eligible', balances), ('bootstrap_summary', bs)]:
        core.write_table(out / (name + '.csv'), rows)
    core.write_table(out / 'diagnosis_counts.csv', [dict(diagnosis=d, n_eligible=len(g), n_followup=int(g.included.sum())) for d, g in people.groupby('diagnosis')])
    core.write_table(out / 'selection_model_coefficients.csv', [dict(term=t, beta=float(fit.params[t]), se=float(fit.bse[t])) for t in fit.params.index])
    print(pd.DataFrame(results)[['analysis', 'n_participants', 'beta', 'ci_low', 'ci_high']].to_string(index=False))
    print('Aggregate outputs:', out)
if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root', type=Path, required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--bootstrap', type=int, default=1000)
    p.add_argument('--seed', type=int, default=20260915)
    a = p.parse_args()
    if a.bootstrap < 0:
        p.error('--bootstrap must be nonnegative')
    run(a)
