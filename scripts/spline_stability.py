"""Spline shape comparisons using fixed full-sample knots and a common tau range."""


from __future__ import annotations
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
import patsy
import statsmodels
from statsmodels.stats.multitest import multipletests
from cross_sectional import COVARIATES, OUTCOMES, eligible, fit, validate
from tau_sensitivity import spline_design, model_fit

def contrast_predictions(model, design):
    value = design @ np.asarray(model.params)
    variance = np.einsum('ij,jk,ik->i', design, np.asarray(model.cov_params()), design)
    if np.any(variance < -1e-10):
        raise ValueError('Negative contrast variance')
    se = np.sqrt(np.maximum(variance, 0))
    return (value, value - 1.959963984540054 * se, value + 1.959963984540054 * se)

def nonlinear(model):
    r = np.zeros((2, len(model.params)))
    r[:, -2:] = np.eye(2)
    test = model.wald_test(r, use_f=False, scalar=True)
    return (float(test.statistic), float(test.pvalue))

def analyze(master, outcome, age):
    table = eligible(master, outcome, min_age=age)
    table = table.loc[table.at_stage.isin(['A+T-', 'A+T+'])].dropna(subset=['tau_temporal_suvr']).copy()
    table['tau_per_01'] = table.tau_temporal_suvr / 0.1
    linear = fit(f'{outcome} ~ tau_per_01 + ' + ' + '.join(COVARIATES), table)
    influence = linear.get_influence()
    cooks = np.asarray(influence.cooks_distance[0])
    leverage = np.asarray(influence.hat_matrix_diag)
    resid = np.asarray(influence.resid_studentized_internal)
    if not np.isfinite(np.column_stack([cooks, leverage, resid])).all():
        raise ValueError('Nonfinite influence diagnostics')
    n, p = linear.model.exog.shape
    flagged = (cooks > 4 / n) | (leverage > 2 * p / n) | (abs(resid) > 3)
    keep = ~flagged
    tau = table.tau_temporal_suvr.to_numpy(float)
    y = table[outcome].to_numpy(float)
    cov = table[COVARIATES].to_numpy(float)
    x, info, projection, rotation, knots, low, high = spline_design(tau, cov)
    full = model_fit(y, x)
    retained = model_fit(y[keep], x[keep])
    b = np.asarray(patsy.build_design_matrices([info], {'tau': tau[keep]})[0])
    alternate = model_fit(y[keep], np.column_stack([b, cov[keep]]))
    np.testing.assert_allclose(retained.fittedvalues, alternate.fittedvalues, rtol=1e-07, atol=1e-09)
    lower = max(np.quantile(tau, 0.05), np.quantile(tau[keep], 0.05))
    upper = min(np.quantile(tau, 0.95), np.quantile(tau[keep], 0.95))
    if not lower < upper:
        raise ValueError('Empty common central tau range')
    grid = np.linspace(lower, upper, 101)
    reference = float(np.clip(np.median(tau), lower, upper))

    def design(values):
        values = np.asarray(values)
        lin = np.column_stack([np.ones(len(values)), values / 0.1])
        bas = np.asarray(patsy.build_design_matrices([info], {'tau': values})[0])
        return np.column_stack([lin, np.repeat(cov.mean(axis=0)[None, :], len(values), axis=0), (bas - lin @ projection) @ rotation])
    gx = design(grid)
    refx = design([reference])
    centered = gx - refx
    pop = 'all_eligible' if age is None else 'age_ge_65'
    curves = pd.DataFrame(dict(population=pop, outcome=outcome, tau_suvr=grid, reference_tau=reference))
    for name, model in [('full', full), ('retained', retained)]:
        for label, mat in [('mean', gx), ('contrast_to_reference', centered)]:
            val, lo, hi = contrast_predictions(model, mat)
            curves[f'{name}_{label}'] = val
            curves[f'{name}_{label}_ci95_lower'] = lo
            curves[f'{name}_{label}_ci95_upper'] = hi
    curves['retained_minus_full_mean'] = curves.retained_mean - curves.full_mean
    curves['retained_minus_full_shape'] = curves.retained_contrast_to_reference - curves.full_contrast_to_reference
    fstat, fp = nonlinear(full)
    rstat, rp = nonlinear(retained)
    summary = dict(population=pop, outcome=outcome, n_full=n, n_flagged=int(flagged.sum()), n_retained=int(keep.sum()), boundary_low=low, knot1=float(knots[0]), knot2=float(knots[1]), boundary_high=high, common_tau_low=float(lower), common_tau_high=float(upper), reference_tau=reference, full_nonlinear_chi2=fstat, full_nonlinear_p=fp, retained_nonlinear_chi2=rstat, retained_nonlinear_p=rp, max_abs_mean_difference=float(curves.retained_minus_full_mean.abs().max()), max_abs_shape_difference=float(curves.retained_minus_full_shape.abs().max()), rms_shape_difference=float(np.sqrt(np.mean(curves.retained_minus_full_shape ** 2))), full_outcome_sd=float(np.std(y, ddof=1)))
    summary['max_shape_difference_in_outcome_sd'] = summary['max_abs_shape_difference'] / summary['full_outcome_sd']
    for name, model in [('full', full), ('retained', retained)]:
        est, lo, hi = contrast_predictions(model, gx[-1:] - gx[:1])
        summary.update({f'{name}_high_minus_low': float(est[0]), f'{name}_high_minus_low_ci95_lower': float(lo[0]), f'{name}_high_minus_low_ci95_upper': float(hi[0])})
    summary['interpretation'] = 'descriptive shape comparison; separate pointwise CIs, no test of difference between overlapping samples'
    if curves.tau_suvr.min() < tau[keep].min() or curves.tau_suvr.max() > tau[keep].max():
        raise AssertionError('Grid exceeds retained support')
    return (summary, curves)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--master', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    master = pd.read_csv(args.master)
    validate(master)
    rows, curves = ([], [])
    for age in [None, 65]:
        for outcome in OUTCOMES:
            row, curve = analyze(master, outcome, age)
            rows.append(row)
            curves.append(curve)
    summary = pd.DataFrame(rows)
    for kind in ['full', 'retained']:
        summary[f'{kind}_nonlinear_fdr_bh'] = multipletests(summary[f'{kind}_nonlinear_p'], method='fdr_bh')[1]
    summary['fdr_note'] = 'full: four full-sample tests; retained: separate four-test sensitivity family'
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out / 'spline_stability_summary.csv', index=False)
    pd.concat(curves, ignore_index=True).to_csv(out / 'spline_common_range_comparison.csv', index=False)
    scripts = [Path(__file__), Path(__file__).with_name('tau_sensitivity.py'), Path(__file__).with_name('cross_sectional.py')]
    print(summary[['population', 'outcome', 'n_full', 'n_flagged', 'n_retained', 'full_nonlinear_fdr_bh', 'retained_nonlinear_fdr_bh', 'max_shape_difference_in_outcome_sd']].to_string(index=False))
    print('Saved aggregate stability outputs:', out)
if __name__ == '__main__':
    main()
