"""Exploratory A/T contrasts and continuous-tau associations."""
from __future__ import annotations
import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
import statsmodels
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
COVARIATES = ['age_at_tau', 'sex_male', 'education_years', 'apoe4_carrier']
STAGES = ['A-T-', 'A+T-', 'A+T+']
OUTCOMES = ['hippocampus_icv', 'cdrsb']
ALL_STAGES = ['A-T-', 'A-T+', 'A+T-', 'A+T+']

def eligible(master, outcome, window=180, min_age=None):
    needed = ['rid', 'at_stage', 'at_pair_days', outcome, outcome + '_days_from_tau', *COVARIATES]
    keep = master.at_stage.isin(STAGES) & master.at_pair_days.between(0, window)
    keep &= master[outcome + '_days_from_tau'].between(0, window)
    keep &= master[needed].notna().all(axis=1)
    if outcome == 'hippocampus_icv':
        flagged = master.hippocampus_outlier_flag.astype(str).str.lower().isin(['true', '1'])
        keep &= ~flagged
    if min_age is not None:
        keep &= master.age_at_tau.ge(min_age)
    return master.loc[keep].copy()

def fit(formula, table):
    model = smf.ols(formula, data=table, missing='raise')
    x = model.exog
    if not np.isfinite(x).all() or np.linalg.matrix_rank(x) != x.shape[1]:
        raise ValueError('Nonfinite or rank-deficient design')
    if len(table) <= x.shape[1] + 2:
        raise ValueError('Insufficient residual degrees of freedom')
    return model.fit(cov_type='HC3')

def effect(result, vector):
    test = result.t_test(vector)
    interval = np.asarray(test.conf_int()).ravel()
    return dict(beta=float(np.asarray(test.effect).item()), standard_error=float(np.asarray(test.sd).item()), ci95_lower=float(interval[0]), ci95_upper=float(interval[1]), p_value=float(np.asarray(test.pvalue).item()))

def analyze(master, min_age=None):
    rows = []
    for outcome in OUTCOMES:
        table = eligible(master, outcome, min_age=min_age)
        result = fit(f"{outcome} ~ C(at_stage, Treatment(reference='A-T-')) + " + ' + '.join(COVARIATES), table)
        names = list(result.params.index)
        vector = np.zeros(len(names))
        for stage, weight in [('A+T+', 1), ('A+T-', -1)]:
            vector[next((i for i, name in enumerate(names) if f'[T.{stage}]' in name))] = weight
        stats = effect(result, vector)
        check = fit(f"{outcome} ~ C(at_stage, Treatment(reference='A+T-')) + " + ' + '.join(COVARIATES), table)
        term = next((n for n in check.params.index if '[T.A+T+]' in n))
        np.testing.assert_allclose([stats['beta'], stats['standard_error']], [check.params[term], check.bse[term]], rtol=1e-08, atol=1e-10)
        rows.append(dict(analysis='direct_A+T+_vs_A+T-', outcome=outcome, n=len(table), n_AminusTminus=int(table.at_stage.eq('A-T-').sum()), n_AplusTminus=int(table.at_stage.eq('A+T-').sum()), n_AplusTplus=int(table.at_stage.eq('A+T+').sum()), **stats))
        positive = table.loc[table.at_stage.isin(['A+T-', 'A+T+'])].dropna(subset=['tau_temporal_suvr']).copy()
        positive['tau_per_01'] = positive.tau_temporal_suvr / 0.1
        model = fit(f'{outcome} ~ tau_per_01 + ' + ' + '.join(COVARIATES), positive)
        v = np.zeros(len(model.params))
        v[list(model.params.index).index('tau_per_01')] = 1
        stats = effect(model, v)
        check = fit(f'{outcome} ~ tau_temporal_suvr + ' + ' + '.join(COVARIATES), positive)
        np.testing.assert_allclose([stats['beta'], stats['standard_error']], [check.params['tau_temporal_suvr'] * 0.1, check.bse['tau_temporal_suvr'] * 0.1], rtol=1e-08, atol=1e-10)
        rows.append(dict(analysis='continuous_tau_per_0.1_in_Apositive', outcome=outcome, n=len(positive), n_AminusTminus=0, n_AplusTminus=int(positive.at_stage.eq('A+T-').sum()), n_AplusTplus=int(positive.at_stage.eq('A+T+').sum()), tau_min=positive.tau_temporal_suvr.min(), tau_max=positive.tau_temporal_suvr.max(), **stats))
    output = pd.DataFrame(rows)
    output['fdr_bh'] = multipletests(output.p_value, method='fdr_bh')[1]
    output['population'] = 'all_eligible' if min_age is None else f'age_ge_{min_age}'
    output['fdr_family'] = 'four exploratory extension tests: ' + output.population
    return output

def describe(master):
    rows = []
    variables = ['age_at_tau', 'sex_male', 'education_years', 'apoe4_carrier', 'tau_temporal_suvr', 'hippocampus_icv', 'cdrsb', 'fdg_metaroi']
    for stage in ['Total', *ALL_STAGES]:
        table = master if stage == 'Total' else master.loc[master.at_stage.eq(stage)]
        for var in variables:
            values = table[var].dropna()
            row = dict(at_stage=stage, variable=var, n_total=len(table), n_nonmissing=len(values), n_missing=len(table) - len(values))
            if var in ['sex_male', 'apoe4_carrier']:
                row.update(n_positive=int(values.eq(1).sum()), proportion=values.mean())
            else:
                row.update(mean=values.mean(), sd=values.std(), median=values.median(), q25=values.quantile(0.25), q75=values.quantile(0.75), minimum=values.min(), maximum=values.max())
            rows.append(row)
    return pd.DataFrame(rows)

def validate(master):
    if master.rid.isna().any() or master.rid.duplicated().any():
        raise ValueError('Missing or duplicate participant IDs')
    for c in COVARIATES + ['at_pair_days', 'tau_temporal_suvr', *OUTCOMES]:
        if not np.isfinite(master[c].dropna()).all():
            raise ValueError(f'Nonfinite {c}')
    for c in ['sex_male', 'apoe4_carrier']:
        if not master[c].dropna().isin([0, 1]).all():
            raise ValueError(f'Invalid binary {c}')
    days = (pd.to_datetime(master.amy_pet_date) - pd.to_datetime(master.tau_pet_date)).dt.days.abs()
    np.testing.assert_array_equal(days, master.at_pair_days)
    np.testing.assert_array_equal(master.tau_temporal_suvr.ge(1.34), master.at_stage.str.endswith('T+'))
    return None

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--master', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    master = pd.read_csv(args.master)
    validate(master)
    primary = analyze(master)
    sensitivity = analyze(master, min_age=65)
    baseline = describe(master)
    age = []
    for name, table in [('full_cohort', master), ('three_group_cohort', master.loc[master.at_stage.isin(STAGES)])]:
        age.append(dict(population=name, n=len(table), age_missing=int(table.age_at_tau.isna().sum()), age_under65=int(table.age_at_tau.lt(65).sum()), age_ge65=int(table.age_at_tau.ge(65).sum()), age_min=table.age_at_tau.min(), age_max=table.age_at_tau.max()))
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    for name, table in [('extension_results', pd.concat([primary, sensitivity], ignore_index=True)), ('clinical_baseline', baseline), ('age_coverage', pd.DataFrame(age))]:
        table.to_csv(out / (name + '.csv'), index=False)
    print(pd.concat([primary, sensitivity]).to_string(index=False))
    print(pd.DataFrame(age).to_string(index=False))
    print('Saved aggregate outputs:', out)
if __name__ == '__main__':
    main()
