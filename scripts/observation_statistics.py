"""Numerical guards and reporting counts; no statistical model fitting."""
import numpy as np
import pandas as pd


def dispersion(values):
    x = np.asarray(values, dtype=float)
    valid = len(x) > 1 and np.isfinite(x).all()
    sd = float(np.std(x, ddof=1)) if valid else np.nan
    span = float(np.ptp(x)) if valid else np.nan
    scale = max(np.finfo(float).tiny, float(np.max(np.abs(x)))) if valid else np.nan
    tolerance = 64 * np.finfo(float).eps * scale if valid else np.nan
    if not valid:
        status = 'not_applicable_nonfinite_or_insufficient'
    elif span == 0:
        status = 'not_applicable_constant'
    elif span <= tolerance or sd <= tolerance:
        status = 'not_applicable_near_constant'
    else:
        status = 'computed'
    return dict(target_sd=sd, target_range=span, numerical_tolerance=tolerance,
                smd_status=status, smd_applicable=status == 'computed')


def rule_value(value):
    if pd.isna(value):
        return None
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, float, np.integer, np.floating)) and value in (0, 1):
        return bool(value)
    if isinstance(value, str) and value.lower() in ('true', 'false'):
        return value.lower() == 'true'
    raise ValueError('Invalid probability-rule flag: ' + repr(value))


def bootstrap_counts(draws):
    if draws.duplicated(['replicate', 'model']).any():
        raise ValueError('Duplicate bootstrap replicate/model rows')
    rows = []
    for model, g in draws.groupby('model', sort=False):
        estimates = pd.to_numeric(g.get('estimate', pd.Series(np.nan, index=g.index)), errors='coerce')
        fit_valid = g.status.eq('estimated') & np.isfinite(estimates)
        applicable = model != 'unweighted'
        flags = [rule_value(v) for v in g.get('original_rule_passed', pd.Series(None, index=g.index))] if applicable else []
        passed = sum(v is True for v in flags)
        failed = sum(v is False for v in flags)
        rows.append(dict(model=model, bootstrap_attempted=int(g.replicate.nunique()),
            bootstrap_fit_valid=int(fit_valid.sum()), original_rule_applicable=applicable,
            bootstrap_original_rule_checked=passed + failed,
            bootstrap_original_rule_passed=passed, bootstrap_original_rule_failed=failed,
            bootstrap_original_rule_unknown=sum(v is None for v in flags),
            bootstrap_original_rule_not_applicable=0 if applicable else len(g),
            inference_policy='additional_exploratory_numerical_support' if applicable else 'unweighted_reference'))
    return rows


def annotate_effects(effects, draws):
    out = effects.copy()
    summary = pd.DataFrame(bootstrap_counts(draws))
    if out.model.duplicated().any() or set(out.model) != set(summary.model):
        raise ValueError('Point estimates and bootstrap model names differ')
    if 'source_status' not in out:
        out['source_status'] = out.status
    summary = summary.set_index('model')
    for row in out.itertuples():
        counts = summary.loc[row.model]
        if int(row.bootstrap_valid) != counts.bootstrap_fit_valid or int(row.bootstrap_attempted) != counts.bootstrap_attempted:
            raise ValueError('Bootstrap counts disagree with recorded inference: ' + row.model)
    for column in summary.columns:
        if column != 'bootstrap_attempted':
            out[column] = out.model.map(summary[column])
    out['original_rule_passed'] = pd.Series(
        [pd.NA if m == 'unweighted' else rule_value(v) for m,v in zip(out.model,out.original_rule_passed)],
        index=out.index,dtype='boolean')
    out['status'] = out.status.replace({'bootstrap_stable': 'estimated_with_bootstrap',
                                        'bootstrap_unstable': 'bootstrap_insufficient_valid'})
    return out, summary.reset_index()
