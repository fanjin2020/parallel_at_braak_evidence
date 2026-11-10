"""Clinical baseline and follow-up comparison tables from source data."""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
VERSION = '2026-09-15.clinical-tables.1'
VARS = ['age_at_tau', 'sex_male', 'education_years', 'apoe4_carrier', 'tau_temporal_suvr', 'cdrsb']
LABELS = {'age_at_tau': 'Age, years', 'sex_male': 'Male sex', 'education_years': 'Education, years', 'apoe4_carrier': 'APOE ε4 carrier', 'tau_temporal_suvr': 'Tau SUVR', 'cdrsb': 'CDR-SB', 'baseline_cdrsb': 'Baseline CDR-SB', 'baseline_offset_days': 'Baseline CDR-SB offset, days'}
STAGES = ['A-T-', 'A-T+', 'A+T-', 'A+T+']

def read(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def write(path, rows):
    fields = list(dict.fromkeys((k for row in rows for k in row)))
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

def summarize(records, var):
    a = np.array([r[var] for r in records if r.get(var) is not None], dtype=float)
    row = {'n_total': len(records), 'n_nonmissing': len(a), 'n_missing': len(records) - len(a)}
    if len(a):
        row.update(mean=float(a.mean()), sd=float(a.std(ddof=1)) if len(a) > 1 else None, median=float(np.median(a)), q25=float(np.quantile(a, 0.25)), q75=float(np.quantile(a, 0.75)))
        if var in ['sex_male', 'apoe4_carrier']:
            row.update(n_positive=int((a == 1).sum()), proportion=float(a.mean()))
    return row

def number(v):
    return None if v is None or v == '' else float(v)

def display(s, var):
    n = int(s['n_nonmissing'])
    if not n:
        return 'Not available (n=0)'
    if var in ['sex_male', 'apoe4_carrier']:
        return f"{int(float(s['n_positive']))}/{n} ({float(s['proportion']) * 100:.1f}%)"
    if var in ['cdrsb', 'baseline_cdrsb', 'tau_temporal_suvr']:
        digits = 3 if var == 'tau_temporal_suvr' else 1
        return f"{float(s['median']):.{digits}f} ({float(s['q25']):.{digits}f}–{float(s['q75']):.{digits}f}); n={n}"
    return f"{float(s['mean']):.1f} ({float(s['sd']):.1f}); n={n}"

def dx_summary(records):
    counts = {k: sum((r['diagnosis_at_tau'] == k for r in records)) for k in ['CN', 'MCI', 'Dementia', 'no_visit', 'conflicting_nearest_date']}
    counts['known'] = sum((counts[k] for k in ['CN', 'MCI', 'Dementia']))
    counts['missing'] = len(records) - counts['known']
    assert sum((counts[k] for k in ['CN', 'MCI', 'Dementia', 'no_visit', 'conflicting_nearest_date'])) == len(records)
    return counts

def smd(a, b):
    a, b = (np.asarray(a, dtype=float), np.asarray(b, dtype=float))
    if len(a) < 2 or len(b) < 2:
        return None
    scale = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    if scale == 0:
        return 0.0 if a.mean() == b.mean() else None
    return float((a.mean() - b.mean()) / scale)

def build(root, out):
    import longitudinal as core
    root, out = (Path(root).resolve(), Path(out).resolve())
    if out == root or root / 'data' == out or root / 'data' in out.parents:
        raise ValueError('Output must not be an input directory')
    if out.exists() and any(out.iterdir()):
        raise ValueError('Choose a new, empty output directory')
    files = [root / 'data/private/adni/adni_pet_aligned_master.csv', root / 'data/private/adni/raw/CDR.csv', root / 'data/private/adni/raw/DXSUM.csv']
    master, checks = core.load_master(files[0])
    ids = {r['rid'] for r in master}
    cdr, qc_cdr = core.load_visits(files[1], ids, 'cdr')
    dx, qc_dx = core.load_visits(files[2], ids, 'dx')
    for r in master:
        visit, status = core.nearest(dx.get(r['rid'], {}), r['anchor'], 180)
        r['diagnosis_at_tau'] = {1: 'CN', 2: 'MCI', 3: 'Dementia'}[visit[1]] if visit else status
    eligible, long, flow = core.make_cohort(master, cdr)
    followed = [r for r in eligible if r['included']]
    unfollowed = [r for r in eligible if not r['included']]
    baseline = []
    for stage in STAGES:
        group = [r for r in master if r['at_stage'] == stage]
        for var in VARS:
            baseline.append({'group': stage, 'variable': var, **summarize(group, var)})
    old_attr = core.attrition_summary(eligible)
    groups = {stage: [r for r in master if r['at_stage'] == stage] for stage in STAGES}
    long_groups = {'Baseline eligible': eligible, 'With follow-up': followed, 'Without follow-up': unfollowed}
    summaries = []
    diagnoses = []
    missing = []
    for name, records in {**groups, **long_groups}.items():
        variables = VARS if name in STAGES else core.COVS + ['tau_temporal_suvr', 'baseline_cdrsb', 'baseline_offset_days']
        for v in variables:
            s = summarize(records, v)
            summaries.append({'group': name, 'variable': v, **s})
            missing.append({'group': name, 'variable': v, 'n_total': len(records), 'n_available': s['n_nonmissing'], 'n_missing': s['n_missing']})
        dd = dx_summary(records)
        for k in ['CN', 'MCI', 'Dementia', 'no_visit', 'conflicting_nearest_date']:
            diagnoses.append({'group': name, 'diagnosis': k, 'n': dd[k], 'n_total': len(records), 'n_known_diagnosis': dd['known']})
        missing.append({'group': name, 'variable': 'diagnosis_at_tau', 'n_total': len(records), 'n_available': dd['known'], 'n_missing': dd['missing']})
    index = {(r['group'], r['variable']): r for r in summaries}
    for r in baseline:
        index[r['group'], r['variable']] = r

    def panel(gps, variables, include_models=False):
        rows = [{'Characteristic': 'Participants, n', **{k: str(len(v)) for k, v in gps.items()}}]
        for var in variables:
            rows.append({'Characteristic': LABELS[var], **{k: display(index[k, var], var) for k in gps}})
        for label in ['CN', 'MCI', 'Dementia']:
            cells = {}
            for k, recs in gps.items():
                d = dx_summary(recs)
                cells[k] = f"{d[label]}/{d['known']} ({100 * d[label] / d['known']:.1f}%)" if d['known'] else 'Not available'
            rows.append({'Characteristic': 'Diagnosis: ' + label, **cells})
        rows.append({'Characteristic': 'Diagnosis unavailable, n', **{k: str(dx_summary(v)['missing']) for k, v in gps.items()}})
        if include_models:
            for outcome, label in [('hippocampus_icv', 'Hippocampus/ICV model, n'), ('cdrsb', 'CDR-SB model, n')]:
                counts = {}
                for k, recs in gps.items():
                    if k == 'A-T+':
                        counts[k] = 'Descriptive only'
                        continue
                    n = sum((all((r[c] is not None for c in core.COVS + [outcome])) and r[outcome + '_days_from_tau'] is not None and (0 <= r[outcome + '_days_from_tau'] <= 180) and (not (outcome == 'hippocampus_icv' and r['hip_flag'])) for r in recs))
                    counts[k] = str(n)
                rows.append({'Characteristic': label, **counts})
        return rows
    a = panel(groups, VARS, True)
    b = panel(long_groups, core.COVS + ['tau_temporal_suvr', 'baseline_cdrsb'])
    for stage in ['A+T-', 'A+T+']:
        b.insert(1 if stage == 'A+T-' else 2, {'Characteristic': stage + ', n', **{k: str(sum((r['at_stage'] == stage for r in v))) for k, v in long_groups.items()}})
    attr_full = []
    attr_display = []
    for old in old_attr:
        v = old['variable']
        sa = index['With follow-up', v]
        sb = index['Without follow-up', v]
        attr_full.append({**old, 'with_followup_sd': sa['sd'], 'without_followup_sd': sb['sd'], 'with_followup_missing': sa['n_missing'], 'without_followup_missing': sb['n_missing'], 'source': 'baseline_summary'})

        def attr_cell(s):
            if v in ['sex_male', 'apoe4_carrier']:
                return display(s, v)
            digits = 3 if v == 'tau_temporal_suvr' else 2
            return f"{s['mean']:.{digits}f} ({s['sd']:.{digits}f}); n={s['n_nonmissing']}"
        attr_display.append({'Characteristic': LABELS[v], 'With follow-up': attr_cell(sa), 'Without follow-up': attr_cell(sb), 'SMD': f"{float(old['standardized_mean_difference']):.3f}", 'Missing, with/without': f"{sa['n_missing']}/{sb['n_missing']}"})
    for label in ['CN', 'MCI', 'Dementia', 'Unavailable']:
        if label == 'Unavailable':
            va = [int(r['diagnosis_at_tau'] not in ['CN', 'MCI', 'Dementia']) for r in followed]
            vb = [int(r['diagnosis_at_tau'] not in ['CN', 'MCI', 'Dementia']) for r in unfollowed]
        else:
            va = [int(r['diagnosis_at_tau'] == label) for r in followed if r['diagnosis_at_tau'] in ['CN', 'MCI', 'Dementia']]
            vb = [int(r['diagnosis_at_tau'] == label) for r in unfollowed if r['diagnosis_at_tau'] in ['CN', 'MCI', 'Dementia']]
        value = smd(va, vb)
        attr_full.append({'variable': 'diagnosis_' + label, 'with_followup_n': len(va), 'without_followup_n': len(vb), 'with_followup_mean': float(np.mean(va)) if va else None, 'without_followup_mean': float(np.mean(vb)) if vb else None, 'standardized_mean_difference': value, 'source': 'diagnosis_summary'})
        attr_display.append({'Characteristic': 'Diagnosis: ' + label, 'With follow-up': f'{sum(va)}/{len(va)} ({100 * np.mean(va):.1f}%)' if va else 'Not available', 'Without follow-up': f'{sum(vb)}/{len(vb)} ({100 * np.mean(vb):.1f}%)' if vb else 'Not available', 'SMD': f'{value:.3f}' if value is not None else 'Not estimable', 'Missing, with/without': f"{dx_summary(followed)['missing']}/{dx_summary(unfollowed)['missing']}" if label != 'Unavailable' else '0/0'})
    notes = {'table1': 'Age and education are mean (SD); tau SUVR and CDR-SB are median (IQR). Binary variables and diagnoses are n/available N (%). Continuous summaries show available n. Diagnosis is the nearest eligible assessment within ±180 days of tau PET; equal-distance ties favor the earlier date. Conflicting nearest-date diagnoses are treated as unavailable, not replaced by a more distant visit. Descriptive denominators are variable-specific and differ from model complete-case sets. A−T+ is descriptive only.', 'table1b': 'The longitudinal baseline-eligible cohort requires complete model covariates and an eligible CDR-SB baseline. Eligible follow-up occurs after both baseline and tau PET and within three years of tau PET. Baseline CDR-SB is selected with the existing longitudinal analysis rules. This table does not introduce a new model.', 's9': 'Continuous variables are mean (SD), with available n; binary variables and diagnoses are n/available N (%). Signed SMD is (mean with follow-up − mean without follow-up) divided by the square root of the average sample variances, retaining the original attrition convention for continuous and binary variables. Diagnosis categories use one-versus-rest indicators among known diagnoses; diagnosis unavailability uses all participants. SMD is descriptive, not a significance test. Offsets are signed days from tau PET. '}
    out.mkdir(parents=True, exist_ok=True)
    outputs = {'table1a_clinical_by_at.csv': a, 'table1b_longitudinal_cohorts.csv': b, 'table_s9_followup_comparison.csv': attr_display, 'clinical_baseline_numeric.csv': summaries, 'diagnosis_by_analysis_cohort.csv': diagnoses, 'missingness_by_cohort.csv': missing, 'attrition_comparison_extended_numeric.csv': attr_full}
    for name, rows in outputs.items():
        write(out / name, rows)
    payload = {'version': VERSION, 'notes': notes, 'table1a': a, 'table1b': b, 'table_s9': attr_display, 'diagnosis': diagnoses, 'summary': summaries}
    (out / 'clinical_tables.json').write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding='utf-8')
    print('Clinical tables written to', out)

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root', type=Path, required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    a = p.parse_args()
    build(a.project_root, a.output_dir)
if __name__ == '__main__':
    main()
