"""Natural-spline and influence sensitivity analyses for continuous tau."""

from __future__ import annotations
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
import patsy
import statsmodels
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from cross_sectional import COVARIATES, OUTCOMES, eligible, fit, validate

def model_fit(y, x):
    if not np.isfinite(x).all() or np.linalg.matrix_rank(x) != x.shape[1]:
        raise ValueError('Invalid or rank-deficient spline design')
    if len(y) <= x.shape[1] + 2:
        raise ValueError('Insufficient residual degrees of freedom')
    return sm.OLS(y, x).fit(cov_type='HC3')

def spline_design(tau, covariates):
    low, high = (float(tau.min()), float(tau.max()))
    knots = np.quantile(tau, [1 / 3, 2 / 3])
    if not low < knots[0] < knots[1] < high:
        raise ValueError('Insufficient distinct tau values for fixed spline')
    specification = 'cr(tau, knots=knots, lower_bound=low, upper_bound=high) - 1'
    basis = patsy.dmatrix(specification, dict(tau=tau, knots=knots, low=low, high=high))
    linear = np.column_stack([np.ones(len(tau)), tau / 0.1])
    projection = np.linalg.lstsq(linear, np.asarray(basis), rcond=None)[0]
    residual = np.asarray(basis) - linear @ projection
    _, values, vt = np.linalg.svd(residual, full_matrices=False)
    rank = int(np.sum(values > values[0] * 1e-10))
    if rank != 2:
        raise ValueError(f'Expected two nonlinear terms, obtained {rank}')
    rotation = vt[:rank].T
    nonlinear = residual @ rotation
    x = np.column_stack([linear, covariates, nonlinear])
    return (x, basis.design_info, projection, rotation, knots, low, high)

def summarize_slope(result):
    name = 'tau_per_01'
    ci = result.conf_int().loc[name]
    return dict(beta=float(result.params[name]), se=float(result.bse[name]), ci95_lower=float(ci.iloc[0]), ci95_upper=float(ci.iloc[1]), p_value=float(result.pvalues[name]))

def diagnose(master, outcome, min_age):
    table = eligible(master, outcome, min_age=min_age)
    table = table.loc[table.at_stage.isin(['A+T-', 'A+T+'])].dropna(subset=['tau_temporal_suvr']).copy()
    table['tau_per_01'] = table.tau_temporal_suvr / 0.1
    population = 'all_eligible' if min_age is None else 'age_ge_65'
    formula = f'{outcome} ~ tau_per_01 + ' + ' + '.join(COVARIATES)
    linear = fit(formula, table)
    original = summarize_slope(linear)
    tau = table.tau_temporal_suvr.to_numpy(float)
    y = table[outcome].to_numpy(float)
    cov = table[COVARIATES].to_numpy(float)
    x, info, projection, rotation, knots, low, high = spline_design(tau, cov)
    spline = model_fit(y, x)
    restriction = np.zeros((2, x.shape[1]))
    restriction[:, -2:] = np.eye(2)
    test = spline.wald_test(restriction, use_f=False, scalar=True)
    b = np.asarray(patsy.build_design_matrices([info], {'tau': tau})[0])
    alternate = model_fit(y, np.column_stack([b, cov]))
    np.testing.assert_allclose(spline.fittedvalues, alternate.fittedvalues, rtol=1e-07, atol=1e-09)
    if np.sum(spline.resid ** 2) > np.sum(linear.resid ** 2) + 1e-08:
        raise AssertionError('Spline does not contain the linear model')
    nl = dict(population=population, outcome=outcome, n=len(table), tau_boundary_low=low, tau_knot1=float(knots[0]), tau_knot2=float(knots[1]), tau_boundary_high=high, nonlinear_df=2, nonlinear_wald_chi2=float(test.statistic), nonlinear_p=float(test.pvalue), linear_r_squared=float(linear.rsquared), spline_r_squared=float(spline.rsquared))
    influence = linear.get_influence()
    cooks = np.asarray(influence.cooks_distance[0])
    leverage = np.asarray(influence.hat_matrix_diag)
    student = np.asarray(influence.resid_studentized_internal)
    if not np.isfinite(np.column_stack([cooks, leverage, student])).all():
        raise ValueError('Nonfinite influence statistics')
    n, p = linear.model.exog.shape
    cook_flag, lev_flag, resid_flag = (cooks > 4 / n, leverage > 2 * p / n, np.abs(student) > 3)
    flag = cook_flag | lev_flag | resid_flag
    reduced = table.loc[~flag].copy()
    reduced_result = fit(formula, reduced)
    sens = summarize_slope(reduced_result)
    diagnostic = dict(population=population, outcome=outcome, n=n, n_parameters=p, n_cook_flag=int(cook_flag.sum()), n_leverage_flag=int(lev_flag.sum()), n_studentized_flag=int(resid_flag.sum()), n_union_flag=int(flag.sum()), n_retained=len(reduced), max_cook_d=float(cooks.max()), max_leverage=float(leverage.max()), max_abs_studentized=float(np.abs(student).max()), cook_threshold=4 / n, leverage_threshold=2 * p / n, studentized_threshold=3, tau_min_retained=float(reduced.tau_temporal_suvr.min()), tau_max_retained=float(reduced.tau_temporal_suvr.max()), direction_same=bool(np.sign(original['beta']) == np.sign(sens['beta'])), beta_change_percent=100 * (sens['beta'] - original['beta']) / abs(original['beta']) if original['beta'] != 0 else np.nan)
    diagnostic.update({'full_' + k: v for k, v in original.items()})
    diagnostic.update({'sensitivity_' + k: v for k, v in sens.items()})
    grid = np.linspace(*np.quantile(tau, [0.05, 0.95]), 101)
    g_linear = np.column_stack([np.ones(len(grid)), grid / 0.1])
    gb = np.asarray(patsy.build_design_matrices([info], {'tau': grid})[0])
    gcov = np.repeat(cov.mean(axis=0)[None, :], len(grid), axis=0)
    gx = np.column_stack([g_linear, gcov, (gb - g_linear @ projection) @ rotation])
    prediction = spline.get_prediction(gx).summary_frame(alpha=0.05)
    new = pd.DataFrame(gcov, columns=COVARIATES)
    new['tau_per_01'] = grid / 0.1
    lp = linear.get_prediction(new).summary_frame(alpha=0.05)
    curve = pd.DataFrame(dict(population=population, outcome=outcome, tau_suvr=grid, spline_mean=prediction['mean'].to_numpy(), spline_ci95_lower=prediction['mean_ci_lower'].to_numpy(), spline_ci95_upper=prediction['mean_ci_upper'].to_numpy(), linear_mean=lp['mean'].to_numpy(), linear_ci95_lower=lp['mean_ci_lower'].to_numpy(), linear_ci95_upper=lp['mean_ci_upper'].to_numpy()))
    curve['interpretation'] = 'adjusted mean at covariate means; pointwise CI; 5th-95th tau percentiles'
    return (nl, diagnostic, curve)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--master', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    master = pd.read_csv(args.master)
    validate(master)
    nl, inf, curves = ([], [], [])
    for age in [None, 65]:
        for outcome in OUTCOMES:
            a, b, c = diagnose(master, outcome, age)
            nl.append(a)
            inf.append(b)
            curves.append(c)
    nl = pd.DataFrame(nl)
    inf = pd.DataFrame(inf)
    nl['nonlinear_fdr_bh'] = multipletests(nl.nonlinear_p, method='fdr_bh')[1]
    nl['fdr_family'] = 'four nonlinear tests: two outcomes by two age populations'
    inf['sensitivity_fdr_bh'] = multipletests(inf.sensitivity_p_value, method='fdr_bh')[1]
    inf['sensitivity_fdr_family'] = 'four influence-refit slopes; sensitivity only'
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    nl.to_csv(out / 'nonlinearity_results.csv', index=False)
    inf.to_csv(out / 'influence_sensitivity.csv', index=False)
    pd.concat(curves, ignore_index=True).to_csv(out / 'adjusted_tau_curves.csv', index=False)
    print(nl.to_string(index=False))
    print(inf[['population', 'outcome', 'n', 'n_union_flag', 'n_retained', 'full_beta', 'sensitivity_beta', 'beta_change_percent', 'direction_same']].to_string(index=False))
    print('Diagnostics written to:', out)
if __name__ == '__main__':
    main()
