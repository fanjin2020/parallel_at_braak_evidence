"""Inclusive intervals around PET cutoffs expressed as decimal SUVR values."""
from decimal import Decimal, InvalidOperation

import numpy as np
import pandas as pd


def band_bounds(center, width):
    try:
        center, width = Decimal(str(center)), Decimal(str(width))
    except InvalidOperation as exc:
        raise ValueError('Cutoff and width must be finite numbers') from exc
    if not center.is_finite() or not width.is_finite() or width < 0:
        raise ValueError('Cutoff must be finite and width must be finite and nonnegative')
    lower, upper = float(center - width), float(center + width)
    if not np.isfinite([lower, upper]).all():
        raise ValueError('Cutoff interval cannot be represented as finite numbers')
    return lower, upper


def within_band(values, center, width):
    """Include both decimal endpoints without rounding the observations."""
    lower, upper = band_bounds(center, width)
    values = pd.to_numeric(pd.Series(values, copy=False), errors='coerce')
    return (np.isfinite(values) & values.ge(lower) & values.le(upper)).fillna(False)


def outside_band(values, center, width):
    """Exclude both endpoints and all missing/nonfinite observations."""
    values = pd.to_numeric(pd.Series(values, copy=False), errors='coerce')
    return (np.isfinite(values) & ~within_band(values, center, width)).fillna(False)
