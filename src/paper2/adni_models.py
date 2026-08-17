"""Pre-specified ADNI regression models used by ``scripts/analyze_adni.py``."""

from __future__ import annotations

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests


STAGES = ["A-T-", "A+T-", "A+T+"]
COVARIATES = ["age_at_tau", "sex_male", "education_years", "apoe4_carrier"]


def _complete_case(data: pd.DataFrame, outcome: str, day_column: str, window_days: int = 180) -> pd.DataFrame:
    required = ["at_stage", outcome, day_column, *COVARIATES]
    table = data.loc[data["at_stage"].isin(STAGES), required].dropna().copy()
    table = table.loc[table[day_column].between(0, window_days)].copy()
    if outcome == "hippocampus_icv" and "hippocampus_outlier_flag" in data.columns:
        outlier = (
            data.loc[table.index, "hippocampus_outlier_flag"]
            .astype("string").str.strip().str.casefold().isin(["true", "1"])
        )
        table = table.loc[~outlier].copy()
    table["at_stage"] = pd.Categorical(table["at_stage"], categories=STAGES)
    return table


def _stage_rows(result, outcome: str, family: str, n: int) -> list[dict]:
    rows: list[dict] = []
    for stage in STAGES[1:]:
        term = next(name for name in result.params.index if f"[T.{stage}]" in name)
        estimate = float(result.params[term])
        standard_error = float(result.bse[term])
        rows.append(
            {
                "outcome": outcome,
                "family": family,
                "contrast": f"{stage} vs A-T-",
                "n_complete_case": n,
                "beta": estimate,
                "hc3_standard_error": standard_error,
                "ci95_lower": estimate - 1.96 * standard_error,
                "ci95_upper": estimate + 1.96 * standard_error,
                "p_value": float(result.pvalues[term]),
            }
        )
    return rows


def fit_stage_model(data: pd.DataFrame, outcome: str, family: str, day_column: str, window_days: int = 180) -> pd.DataFrame:
    table = _complete_case(data, outcome, day_column, window_days=window_days)
    formula = f"{outcome} ~ C(at_stage, Treatment(reference='A-T-')) + " + " + ".join(COVARIATES)
    result = smf.ols(formula, data=table).fit(cov_type="HC3")
    output = pd.DataFrame(_stage_rows(result, outcome, family, len(table)))
    output["fdr_bh"] = multipletests(output["p_value"], method="fdr_bh")[1]
    return output


def fit_apoe_interaction(data: pd.DataFrame, outcome: str, day_column: str) -> pd.DataFrame:
    table = _complete_case(data, outcome, day_column)
    formula = (
        f"{outcome} ~ C(at_stage, Treatment(reference='A-T-')) * apoe4_carrier "
        "+ age_at_tau + sex_male + education_years"
    )
    result = smf.ols(formula, data=table).fit(cov_type="HC3")
    terms = [name for name in result.params.index if "at_stage" in name and "apoe4_carrier" in name]
    if len(terms) != 2:
        raise ValueError(f"Unexpected APOE interaction terms for {outcome}: {terms}")
    restriction = np.zeros((2, len(result.params)))
    for row, term in enumerate(terms):
        restriction[row, list(result.params.index).index(term)] = 1
    test = result.wald_test(restriction, scalar=False)
    return pd.DataFrame(
        [{
            "outcome": outcome,
            "n_complete_case": len(table),
            "joint_chi_square": float(np.asarray(test.statistic).squeeze()),
            "degrees_of_freedom": 2,
            "p_value": float(test.pvalue),
        }]
    )


def fit_apoe_stratified_effects(data: pd.DataFrame, outcome: str, day_column: str) -> pd.DataFrame:
    """Obtain simple A/T effects from the same HC3 interaction model.

    The estimates are unchanged from the former implementation. The output now
    separates the full interaction-model sample size from the actual APOE ε4
    stratum size, so the two quantities cannot be confused in a table or plot.
    """
    table = _complete_case(data, outcome, day_column)
    formula = (
        f"{outcome} ~ C(at_stage, Treatment(reference='A-T-')) * apoe4_carrier "
        "+ age_at_tau + sex_male + education_years"
    )
    result = smf.ols(formula, data=table).fit(cov_type="HC3")
    parameter_names = list(result.params.index)
    covariance = result.cov_params()
    total_n = len(table)
    rows: list[dict] = []

    for stage in STAGES[1:]:
        stage_term = next(name for name in parameter_names if f"[T.{stage}]" in name and ":apoe4_carrier" not in name)
        interaction_term = next(name for name in parameter_names if f"[T.{stage}]" in name and ":apoe4_carrier" in name)
        for status, apoe4_value, weights in (
            ("noncarrier", 0, {stage_term: 1.0}),
            ("carrier", 1, {stage_term: 1.0, interaction_term: 1.0}),
        ):
            vector = pd.Series(0.0, index=parameter_names)
            for term, value in weights.items():
                vector.loc[term] = value
            estimate = float(vector @ result.params)
            standard_error = float(np.sqrt(vector @ covariance @ vector))
            p_value = float(2.0 * norm.sf(abs(estimate / standard_error)))
            rows.append(
                {
                    "outcome": outcome,
                    "contrast": f"{stage} vs A-T-",
                    "apoe4_status": status,
                    "n_total_interaction_model": total_n,
                    "n_apoe4_stratum": int((table["apoe4_carrier"] == apoe4_value).sum()),
                    "beta": estimate,
                    "hc3_standard_error": standard_error,
                    "ci95_lower": estimate - 1.96 * standard_error,
                    "ci95_upper": estimate + 1.96 * standard_error,
                    "p_value": p_value,
                }
            )
    return pd.DataFrame(rows)
