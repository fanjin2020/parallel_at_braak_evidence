"""Exploratory A+ tau and subsequent CDR-SB associations."""

from __future__ import annotations
import argparse
import csv
import json
import math
import warnings
from collections import Counter, defaultdict
from datetime import date, datetime, timezone
from pathlib import Path
VERSION = '2026-09-15.1'
COVS = ['age_at_tau', 'sex_male', 'education_years', 'apoe4_carrier']
STAGES = ['A-T-', 'A-T+', 'A+T-', 'A+T+']
MISSING = {'', 'na', 'nan', 'none', 'null'}

def numeric(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (ValueError, TypeError):
        return None

def iso_date(x):
    try:
        return date.fromisoformat(x.strip()[:10])
    except (ValueError, TypeError, AttributeError):
        return None

def rid(x):
    v = numeric(x)
    if v is None or v < 0 or (not v.is_integer()):
        raise ValueError('Missing or invalid RID; values withheld from log')
    return str(int(v))

def read_csv(path, required):
    with path.open(encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames or []
        if len(header) != len(set(header)) or not set(required).issubset(header):
            raise ValueError(f'{path.name}: duplicate or missing required columns {sorted(set(required) - set(header))}')
        yield from reader

def load_master(path):
    required = ['rid', 'at_stage', 'tau_pet_date', 'amy_pet_date', 'at_pair_days', 'tau_temporal_suvr', 'cdrsb', 'cdrsb_days_from_tau', 'hippocampus_icv', 'hippocampus_icv_days_from_tau', 'hippocampus_outlier_flag'] + COVS
    records = list(read_csv(path, required))
    seen = set()
    for r in records:
        r['rid'] = rid(r['rid'])
        if r['rid'] in seen:
            raise ValueError('Duplicate master RID')
        seen.add(r['rid'])
        r['anchor'] = iso_date(r['tau_pet_date'])
        amy = iso_date(r['amy_pet_date'])
        if r['anchor'] is None or amy is None:
            raise ValueError('Invalid master PET date')
        for col in COVS + ['at_pair_days', 'tau_temporal_suvr', 'cdrsb', 'cdrsb_days_from_tau', 'hippocampus_icv', 'hippocampus_icv_days_from_tau']:
            original = r[col]
            r[col] = numeric(original)
            if r[col] is None and str(original).strip().lower() not in MISSING:
                raise ValueError(f'Invalid numeric master field: {col}')
        if r['at_pair_days'] != abs((amy - r['anchor']).days) or not 0 <= r['at_pair_days'] <= 180:
            raise ValueError('Invalid PET pairing window')
        if r['at_stage'] not in STAGES or r['tau_temporal_suvr'] is None:
            raise ValueError('Invalid A/T stage or tau value')
        if (r['tau_temporal_suvr'] >= 1.34) != r['at_stage'].endswith('T+'):
            raise ValueError('Tau threshold and stage disagree')
        for col in ['sex_male', 'apoe4_carrier']:
            if r[col] not in (None, 0, 1):
                raise ValueError(f'Invalid binary covariate {col}')
        r['hip_flag'] = str(r['hippocampus_outlier_flag']).lower() in ('true', '1', '1.0')
    return (records, [])

def load_visits(path, cohort, kind):
    dc, vc = ('VISDATE', 'CDRSB') if kind == 'cdr' else ('EXAMDATE', 'DIAGNOSIS')
    grouped = defaultdict(lambda: defaultdict(set))
    q = Counter()
    for row in read_csv(path, ['RID', dc, vc, 'HAS_QC_ERROR']):
        q['source_rows'] += 1
        who = rid(row['RID'])
        if who not in cohort:
            continue
        q['cohort_rows'] += 1
        qc = row['HAS_QC_ERROR'].strip().lower()
        if qc in ('1', '1.0', 'has qc error'):
            q['explicit_qc_error_excluded'] += 1
            continue
        if qc not in ('', '0', '0.0', 'does not have qc error or qc error has been approved'):
            raise ValueError(f'{path.name}: unrecognized QC label; inspect dictionary')
        if not qc:
            q['qc_blank_retained_not_assumed_pass'] += 1
        d = iso_date(row[dc])
        v = numeric(row[vc])
        if kind == 'dx' and v is None:
            v = {'CN': 1, 'MCI': 2, 'Dementia': 3}.get(row[vc].strip())
            if v is None and row[vc].strip().lower() not in MISSING:
                raise ValueError('Unrecognized DIAGNOSIS label; do not guess mapping')
        if d is None or v is None:
            q['missing_invalid_date_or_value_excluded'] += 1
            continue
        if kind == 'cdr' and (not 0 <= v <= 18) or (kind == 'dx' and v not in (1, 2, 3)):
            q['out_of_range_excluded'] += 1
            continue
        if v in grouped[who][d]:
            q['identical_value_same_day_deduplicated'] += 1
        grouped[who][d].add(v)
    q['conflicting_participant_dates'] = sum((len(v) > 1 for ds in grouped.values() for v in ds.values()))
    return (grouped, [dict(source=path.name, check=k, count=v) for k, v in sorted(q.items())])

def nearest(visits, anchor, window, pre_only=False):
    candidates = [(abs((d - anchor).days), d, values) for d, values in visits.items() if abs((d - anchor).days) <= window and (not pre_only or d <= anchor)]
    if not candidates:
        return (None, 'no_visit')
    _, d, values = min(candidates, key=lambda z: (z[0], z[1]))
    if len(values) != 1:
        return (None, 'conflicting_nearest_date')
    return ((d, next(iter(values))), 'available')

def make_cohort(master, cdr, window=180, pre_only=False, min_followups=1):
    flow = Counter()
    eligible, long = ([], [])
    for r in master:
        if not r['at_stage'].startswith('A+'):
            continue
        flow['01_A_positive'] += 1
        if r['at_pair_days'] > window:
            continue
        flow['02_pair_window'] += 1
        if any((r[c] is None for c in COVS)):
            continue
        flow['03_complete_covariates'] += 1
        visits = cdr.get(r['rid'], {})
        b, status = nearest(visits, r['anchor'], window, pre_only)
        if b is None:
            flow['excluded_baseline_' + status] += 1
            continue
        flow['04_unambiguous_baseline'] += 1
        bd, bv = b
        later = sorted(((d, next(iter(v))) for d, v in visits.items() if len(v) == 1 and d > max(bd, r['anchor']) and ((d - r['anchor']).days <= 3 * 365.25)))
        included = len(later) >= min_followups
        p = {**r, 'baseline_cdrsb': bv, 'baseline_offset_days': (bd - r['anchor']).days, 'n_followups': len(later), 'included': included}
        eligible.append(p)
        if not included:
            continue
        flow['05_with_required_followup'] += 1
        for d, v in [(bd, bv)] + later:
            long.append({**p, 'cdr': v, 'time_years': (d - bd).days / 365.25})
    return (eligible, long, flow)

def fit_model(long, label):
    import numpy as np
    import pandas as pd
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    df = pd.DataFrame(long)
    if df.empty or df.rid.nunique() < 30:
        return (dict(analysis=label, status='insufficient_participants'), [])
    df['tau01'] = (df.tau_temporal_suvr - 1.34) / 0.1
    df['age10'] = (df.age_at_tau - 75) / 10
    df['education_c'] = df.education_years - 16
    formula = 'cdr ~ time_years * (tau01 + age10 + sex_male + education_c + apoe4_carrier)'
    model = smf.gee(formula, groups='rid', data=df, family=sm.families.Gaussian(), cov_struct=sm.cov_struct.Independence())
    if np.linalg.matrix_rank(model.exog) != model.exog.shape[1]:
        return (dict(analysis=label, status='rank_deficient'), [])
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter('always')
        result = model.fit(maxiter=200, cov_type='robust')
    valid = bool(result.converged and np.isfinite(result.params).all() and np.isfinite(result.bse).all())
    target = 'time_years:tau01'
    ci = result.conf_int()
    summary = dict(analysis=label, status='ok' if valid else 'invalid_fit', n_participants=int(df.rid.nunique()), n_observations=len(df), n_AplusTminus=int(df.loc[df.at_stage.eq('A+T-'), 'rid'].nunique()), n_AplusTplus=int(df.loc[df.at_stage.eq('A+T+'), 'rid'].nunique()), formula=formula, working_correlation='independence', covariance='participant_cluster_robust', warning_types=';'.join(sorted({type(w.message).__name__ for w in ws})), condition_number=float(np.linalg.cond(model.exog)), observed_zero_fraction=float((df.cdr == 0).mean()), fitted_outside_0_18_fraction=float(((result.fittedvalues < 0) | (result.fittedvalues > 18)).mean()), units='CDR-SB points/year per 0.1 tau SUVR')
    if valid:
        summary.update(beta=float(result.params[target]), ci_low=float(ci.loc[target, 0]), ci_high=float(ci.loc[target, 1]), p_value=float(result.pvalues[target]))
    coefficients = [dict(analysis=label, term=t, beta=float(result.params[t]), se=float(result.bse[t]), ci_low=float(ci.loc[t, 0]), ci_high=float(ci.loc[t, 1]), p_value=float(result.pvalues[t])) for t in result.params.index] if valid else []
    return (summary, coefficients)

def attrition_summary(eligible):
    import numpy as np
    result = []
    for variable in COVS + ['tau_temporal_suvr', 'baseline_cdrsb', 'baseline_offset_days']:
        a = np.array([p[variable] for p in eligible if p['included']], dtype=float)
        b = np.array([p[variable] for p in eligible if not p['included']], dtype=float)
        scale = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2) if len(a) > 1 and len(b) > 1 else float('nan')
        result.append(dict(variable=variable, with_followup_n=len(a), without_followup_n=len(b), with_followup_mean=float(a.mean()) if len(a) else '', without_followup_mean=float(b.mean()) if len(b) else '', standardized_mean_difference=float((a.mean() - b.mean()) / scale) if scale > 0 else ''))
    return result

def write_table(path, records):
    keys = list(dict.fromkeys((k for row in records for k in row)))
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(records)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    root, out = (args.project_root.resolve(), args.output_dir.resolve())
    if out.exists() and any(out.iterdir()):
        parser.error('Output directory is not empty; choose a NEW output directory')
    if out == root or root / 'data' == out or root / 'data' in out.parents:
        parser.error('Do not put outputs in the input data directory')
    import numpy, pandas, scipy, statsmodels
    files = [root / 'data/private/adni/adni_pet_aligned_master.csv', root / 'data/private/adni/raw/CDR.csv', root / 'data/private/adni/raw/DXSUM.csv']
    master, checks = load_master(files[0])
    ids = {r['rid'] for r in master}
    cdr, q1 = load_visits(files[1], ids, 'cdr')
    dx, q2 = load_visits(files[2], ids, 'dx')
    diagnoses = Counter()
    for r in master:
        d, status = nearest(dx.get(r['rid'], {}), r['anchor'], 180)
        category = {1: 'CN', 2: 'MCI', 3: 'Dementia'}[d[1]] if d else status
        diagnoses[r['at_stage'], category] += 1
    results, coefficients, flows, coverages = ([], [], [], [])
    designs = [('primary_Aplus_180d', 180, False, 1), ('sensitivity_prePET_baseline', 180, True, 1), ('sensitivity_90d', 90, False, 1), ('sensitivity_2_followups', 180, False, 2)]
    attrition = []
    for label, window, pre, min_visits in designs:
        eligible, long, flow = make_cohort(master, cdr, window, pre, min_visits)
        flows.extend((dict(analysis=label, step=k, n=v) for k, v in sorted(flow.items())))
        result, coef = fit_model(long, label)
        results.append(result)
        coefficients.extend(coef)
        for stage in ['A+T-', 'A+T+']:
            people = [p for p in eligible if p['at_stage'] == stage and p['included']]
            coverages.append(dict(analysis=label, stage=stage, n=len(people), baseline_after_pet=sum((p['baseline_offset_days'] > 0 for p in people)), one_followup=sum((p['n_followups'] == 1 for p in people)), two_or_more_followups=sum((p['n_followups'] >= 2 for p in people))))
        if label.startswith('primary'):
            attrition = attrition_summary(eligible)
    out.mkdir(parents=True, exist_ok=True)
    write_table(out / 'diagnosis_at_tau.csv', [dict(stage=s, diagnosis=d, n=n) for (s, d), n in sorted(diagnoses.items())])
    write_table(out / 'sample_flow.csv', flows)
    write_table(out / 'followup_coverage.csv', coverages)
    write_table(out / 'attrition_comparison.csv', attrition)
    write_table(out / 'longitudinal_results.csv', results)
    write_table(out / 'model_coefficients.csv', coefficients)
    print('Completed aggregate outputs:', out)
    for r in results:
        print(r['analysis'], r['status'], 'N=', r.get('n_participants'), 'beta=', r.get('beta'), '95% CI=', r.get('ci_low'), r.get('ci_high'), 'p=', r.get('p_value'))
    if any((r['status'] != 'ok' for r in results)):
        raise SystemExit(2)
if __name__ == '__main__':
    main()
