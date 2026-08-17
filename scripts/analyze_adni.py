#!/usr/bin/env python3
"""Run every ADNI analysis reported in the manuscript from the prepared master table.

中文说明：
本脚本只使用 prepare_adni.py 生成的 adni_pet_aligned_master.csv，不再读取原始
ADNI 文件。它在一个文件中完成论文报告的全部 ADNI 统计分析：180 天主分析、
90 天窗口敏感性分析、FDG 探索性分析、APOE ε4 交互/分层、影响点删除分析，
以及 CDR-SB 的 Huber 回归与分组 bootstrap。

为避免不必要的审计文件和受试者级数据输出，本脚本只写出两个聚合结果文件：
  results/adni/adni_results.csv       所有模型的效应量和检验结果
  results/adni/adni_model_counts.csv  每个模型实际使用的人数及敏感性摘要

The script intentionally produces no RID-level diagnostics, no frozen/legacy files,
and no intermediate tables.  The published 90-day analysis constrains BOTH the
amyloid--tau pairing and the outcome--tau distance to 90 days.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import yaml
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from paper2.adni_models import (  # noqa: E402
    COVARIATES,
    STAGES,
    fit_apoe_interaction,
    fit_apoe_stratified_effects,
    fit_stage_model,
)


OUTCOME_DAYS = {
    "hippocampus_icv": "hippocampus_icv_days_from_tau",
    "cdrsb": "cdrsb_days_from_tau",
    "fdg_metaroi": "fdg_metaroi_days_from_tau",
}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="Path to config/project.yaml")
    parser.add_argument("--bootstrap-reps", type=int, default=2000,
                        help="Number of within-stage Huber bootstrap replicates (default: 2000)")
    parser.add_argument("--seed", type=int, default=20260715,
                        help="Random seed for the Huber bootstrap (default: 20260715)")
    return parser.parse_args()


def resolve(root: Path, configured_path: str) -> Path:
    """Allow either an absolute path or a project-relative YAML path."""
    path = Path(configured_path)
    return path if path.is_absolute() else root / path


def apply_bh(table: pd.DataFrame, family: str) -> pd.DataFrame:
    """Apply BH adjustment to one explicitly declared family of P values."""
    table = table.copy()
    table["fdr_bh"] = multipletests(table["p_value"], method="fdr_bh")[1]
    table["fdr_family"] = family
    return table


def analysis_table(master: pd.DataFrame, outcome: str, window_days: int) -> pd.DataFrame:
    """Return the exact complete-case records for one outcome/window model.

    The two time rules are deliberately both present: at_pair_days is the
    amyloid--tau pairing distance and OUTCOME_DAYS[outcome] is the distance from
    the structural, clinical, or FDG outcome to the tau PET scan.
    """
    outcome_days = OUTCOME_DAYS[outcome]
    required = ["rid", "at_stage", "at_pair_days", outcome, outcome_days, *COVARIATES]
    missing = [column for column in required if column not in master.columns]
    if missing:
        raise ValueError("ADNI master is missing required columns: " + ", ".join(missing))

    table = master.loc[master["at_stage"].isin(STAGES), required].copy()
    table = table.loc[table["at_pair_days"].between(0, window_days)]
    table = table.loc[table[outcome_days].between(0, window_days)].dropna().copy()
    if outcome == "hippocampus_icv" and "hippocampus_outlier_flag" in master.columns:
        flagged = master.loc[table.index, "hippocampus_outlier_flag"].astype("string")
        flagged = flagged.str.strip().str.casefold().isin({"true", "1"})
        table = table.loc[~flagged].copy()
    table["at_stage"] = pd.Categorical(table["at_stage"], categories=STAGES)
    return table


def stage_counts(table: pd.DataFrame, analysis: str, outcome: str,
                 window_days: int, n_removed: int = 0) -> list[dict]:
    """Make transparent aggregate sample counts; never export RID-level rows."""
    rows = [{
        "analysis": analysis,
        "outcome": outcome,
        "window_days": window_days,
        "at_stage": stage,
        "n": int((table["at_stage"] == stage).sum()),
        "n_removed": n_removed,
    } for stage in STAGES]
    rows.append({
        "analysis": analysis,
        "outcome": outcome,
        "window_days": window_days,
        "at_stage": "Total",
        "n": int(len(table)),
        "n_removed": n_removed,
    })
    return rows


def add_analysis(table: pd.DataFrame, analysis: str, family: str | None = None) -> pd.DataFrame:
    """Standardise columns from the compact modelling helpers."""
    table = table.copy()
    table.insert(0, "analysis", analysis)
    if family is not None:
        table["fdr_family"] = family
    return table


def formula(outcome: str) -> str:
    return f"{outcome} ~ C(at_stage, Treatment(reference='A-T-')) + " + " + ".join(COVARIATES)


def stage_term(parameter_names: list[str], stage: str) -> str:
    return next(name for name in parameter_names if f"[T.{stage}]" in name)


def influence_sensitivity(master: pd.DataFrame, outcome: str, window_days: int) -> tuple[pd.DataFrame, pd.DataFrame, list[dict]]:
    """Refit after removing observations flagged by the manuscript's fixed rule.

    A record is flagged when Cook's distance > 4/n, leverage > 2p/n, or the
    absolute internally studentised residual > 3.  Only the number removed is
    retained in the output; RID-level diagnostics are intentionally discarded.
    """
    table = analysis_table(master, outcome, window_days)
    fitted = smf.ols(formula(outcome), data=table).fit()
    diagnostic = fitted.get_influence().summary_frame()
    n, p = len(table), len(fitted.params)
    flagged = (
        (diagnostic["cooks_d"] > 4 / n)
        | (diagnostic["hat_diag"] > 2 * p / n)
        | (diagnostic["student_resid"].abs() > 3)
    )
    reduced = table.loc[~flagged].copy()
    # fit_stage_model applies the same outcome/covariate completeness checks again.
    effects = fit_stage_model(reduced, outcome, "influence", OUTCOME_DAYS[outcome], window_days)
    effects = add_analysis(effects, "influence_sensitivity")
    effects["n_removed"] = int(flagged.sum())
    summary = pd.DataFrame([{
        "analysis": "influence_sensitivity",
        "outcome": outcome,
        "window_days": window_days,
        "influence_rule": "CookD>4/n OR leverage>2p/n OR abs(studentized_residual)>3",
        "n_input": n,
        "n_removed": int(flagged.sum()),
        "n_retained": int(len(reduced)),
    }])
    return effects, summary, stage_counts(reduced, "influence_sensitivity", outcome, window_days, int(flagged.sum()))


def huber_design(table: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, list[str]]:
    """Construct the historical fixed-column Huber design matrix.

    The explicit RID ordering is essential: with a fixed random seed it gives
    exactly the same within-stage bootstrap draws as the frozen analysis.
    """
    ordered = table.sort_values("rid", key=lambda value: value.astype(str)).reset_index(drop=True).copy()
    stage = ordered["at_stage"].astype(str)
    columns = ["Intercept", "A+T-", "A+T+", "age_at_tau", "sex_male", "education_years", "apoe4_carrier"]
    matrix = np.column_stack([
        np.ones(len(ordered)),
        (stage == "A+T-").astype(float),
        (stage == "A+T+").astype(float),
        ordered["age_at_tau"].to_numpy(float),
        ordered["sex_male"].to_numpy(float),
        ordered["education_years"].to_numpy(float),
        ordered["apoe4_carrier"].to_numpy(float),
    ])
    return ordered, matrix, ordered["cdrsb"].to_numpy(float), columns


def robust_scale(residual: np.ndarray) -> float:
    """MAD scale rule used in the frozen Huber IRLS implementation."""
    median = np.median(residual)
    scale = np.median(np.abs(residual - median)) / 0.6744897501960817
    if not np.isfinite(scale) or scale < 1e-8:
        scale = np.sqrt(np.mean(residual ** 2))
    return float(max(scale, 1e-8))


def fit_huber_irls(matrix: np.ndarray, values: np.ndarray,
                   tuning_constant: float = 1.345) -> tuple[np.ndarray, bool]:
    """Exact historical Huber IRLS rule: 200 iterations and 1e-10 tolerance."""
    beta, _, rank, _ = np.linalg.lstsq(matrix, values, rcond=None)
    if rank < matrix.shape[1]:
        raise RuntimeError("Huber initial design matrix is rank deficient.")
    for _ in range(200):
        residual = values - matrix @ beta
        standardized = residual / robust_scale(residual)
        weights = np.minimum(1.0, tuning_constant / np.maximum(np.abs(standardized), 1e-12))
        root_weights = np.sqrt(weights)
        updated, _, rank, _ = np.linalg.lstsq(matrix * root_weights[:, None], values * root_weights, rcond=None)
        if rank < matrix.shape[1]:
            raise RuntimeError("Huber weighted design matrix is rank deficient.")
        if np.max(np.abs(updated - beta)) < 1e-10:
            return updated, True
        beta = updated
    return beta, False


def empirical_two_sided_sign_p(values: np.ndarray) -> float:
    """Historical empirical sign P value used for the Huber bootstrap summary."""
    lower_tail = (float(np.sum(values <= 0.0)) + 1.0) / (len(values) + 1.0)
    upper_tail = (float(np.sum(values >= 0.0)) + 1.0) / (len(values) + 1.0)
    return float(min(1.0, 2.0 * min(lower_tail, upper_tail)))


def huber_bootstrap(master: pd.DataFrame, repetitions: int, seed: int,
                    window_days: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Reproduce the frozen CDR-SB Huber/bootstrap procedure without exporting draws.

    The historical procedure used a custom Huber IRLS rather than a package
    default.  Keeping it here is necessary for exact numerical reproduction;
    only the final aggregate summary is retained.
    """
    ordered, matrix, values, columns = huber_design(analysis_table(master, "cdrsb", window_days))
    beta, converged = fit_huber_irls(matrix, values)
    if not converged:
        raise RuntimeError("Huber IRLS did not converge for the original complete-case data.")

    rng = np.random.default_rng(seed)
    indices = {stage: np.flatnonzero(ordered["at_stage"].astype(str).to_numpy() == stage) for stage in STAGES}
    successful_draws: list[np.ndarray] = []
    for _ in range(repetitions):
        draw = np.concatenate([rng.choice(indices[stage], size=len(indices[stage]), replace=True) for stage in STAGES])
        try:
            boot_beta, converged = fit_huber_irls(matrix[draw, :], values[draw])
            if not converged or not np.all(np.isfinite(boot_beta)):
                raise RuntimeError("Huber bootstrap fit did not converge.")
            successful_draws.append(np.asarray(boot_beta, dtype=float))
        except (np.linalg.LinAlgError, RuntimeError, ValueError):
            continue
    if len(successful_draws) < math.ceil(repetitions * 0.95):
        raise RuntimeError("Fewer than 95% of Huber bootstrap replicates succeeded.")

    bootstrap = np.asarray(successful_draws, dtype=float)
    rows = []
    for stage in STAGES[1:]:
        index = columns.index(stage)
        estimates = bootstrap[:, index]
        rows.append({
            "analysis": "huber_bootstrap_sensitivity",
            "outcome": "cdrsb",
            "contrast": f"{stage} vs A-T-",
            "n_complete_case": int(len(ordered)),
            "beta": float(beta[index]),
            "ci95_lower": float(np.quantile(estimates, 0.025)),
            "ci95_upper": float(np.quantile(estimates, 0.975)),
            "bootstrap_empirical_two_sided_p": empirical_two_sided_sign_p(estimates),
            "bootstrap_repetitions_planned": repetitions,
            "bootstrap_repetitions_successful": int(len(estimates)),
            "bootstrap_repetitions_failed": int(repetitions - len(estimates)),
            "bootstrap_seed": seed,
            "huber_tuning_constant": 1.345,
        })
    summary = pd.DataFrame([{
        "analysis": "huber_bootstrap_sensitivity",
        "outcome": "cdrsb",
        "window_days": window_days,
        "n_input": int(len(ordered)),
        "bootstrap_repetitions_planned": repetitions,
        "bootstrap_repetitions_successful": int(len(successful_draws)),
        "bootstrap_repetitions_failed": int(repetitions - len(successful_draws)),
        "bootstrap_seed": seed,
        "huber_tuning_constant": 1.345,
    }])
    return pd.DataFrame(rows), summary


def main() -> None:
    args = arguments()
    config_path = Path(args.config).resolve()
    root = config_path.parents[1]
    with config_path.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)

    master_path = resolve(root, config["paths"]["adni_master"])
    output_dir = resolve(root, config["paths"].get("output_dir", "results")) / "adni"
    output_dir.mkdir(parents=True, exist_ok=True)
    master = pd.read_csv(master_path)
    primary_window = int(config["analysis"].get("primary_window_days", 180))
    strict_window = int(config["analysis"].get("sensitivity_window_days", 90))

    # Main models.  The FDR family is explicitly four tests across two outcomes.
    primary_data = master.loc[master["at_pair_days"].between(0, primary_window)].copy()
    primary = pd.concat([
        fit_stage_model(primary_data, outcome, "confirmatory", OUTCOME_DAYS[outcome], primary_window)
        for outcome in ("hippocampus_icv", "cdrsb")
    ], ignore_index=True)
    primary = add_analysis(apply_bh(primary, "four confirmatory contrasts"), "confirmatory_180day")

    fdg = fit_stage_model(primary_data, "fdg_metaroi", "exploratory_FDG", OUTCOME_DAYS["fdg_metaroi"], primary_window)
    fdg = add_analysis(apply_bh(fdg, "two exploratory FDG contrasts"), "exploratory_FDG")

    interaction = pd.concat([
        fit_apoe_interaction(primary_data, outcome, OUTCOME_DAYS[outcome])
        for outcome in ("hippocampus_icv", "cdrsb")
    ], ignore_index=True)
    interaction = add_analysis(apply_bh(interaction, "two APOE interaction tests"), "APOE_joint_interaction")

    stratified = pd.concat([
        fit_apoe_stratified_effects(primary_data, outcome, OUTCOME_DAYS[outcome])
        for outcome in ("hippocampus_icv", "cdrsb")
    ], ignore_index=True)
    stratified = add_analysis(stratified, "APOE_stratified_display", "display only; no additional FDR family")

    # Strict 90-day sensitivity: require both time distances to be <= 90 days.
    strict_data = master.loc[master["at_pair_days"].between(0, strict_window)].copy()
    strict = pd.concat([
        fit_stage_model(strict_data, outcome, "sensitivity_90day", OUTCOME_DAYS[outcome], strict_window)
        for outcome in ("hippocampus_icv", "cdrsb")
    ], ignore_index=True)
    strict = add_analysis(apply_bh(strict, "four 90-day contrasts"), "sensitivity_90day")

    influence_parts, influence_notes, count_rows = [], [], []
    for outcome in ("hippocampus_icv", "cdrsb"):
        effect, note, counts = influence_sensitivity(master, outcome, primary_window)
        influence_parts.append(effect)
        influence_notes.append(note)
        count_rows.extend(counts)
    influence = pd.concat(influence_parts, ignore_index=True)
    influence["fdr_bh"] = multipletests(influence["p_value"], method="fdr_bh")[1]
    influence["fdr_family"] = "four influence-sensitivity contrasts"

    huber, huber_note = huber_bootstrap(master, args.bootstrap_reps, args.seed, primary_window)

    results = pd.concat([primary, fdg, interaction, stratified, strict, influence, huber],
                        ignore_index=True, sort=False)
    results.to_csv(output_dir / "adni_results.csv", index=False)

    counts: list[dict] = []
    for outcome in ("hippocampus_icv", "cdrsb"):
        counts.extend(stage_counts(analysis_table(master, outcome, primary_window), "confirmatory_180day", outcome, primary_window))
        counts.extend(stage_counts(analysis_table(master, outcome, strict_window), "sensitivity_90day", outcome, strict_window))
    counts.extend(stage_counts(analysis_table(master, "fdg_metaroi", primary_window), "exploratory_FDG", "fdg_metaroi", primary_window))
    counts.extend(count_rows)
    notes = pd.concat(influence_notes + [huber_note], ignore_index=True, sort=False)
    pd.concat([pd.DataFrame(counts), notes], ignore_index=True, sort=False).to_csv(
        output_dir / "adni_model_counts.csv", index=False
    )
    print("ADNI analyses completed:", output_dir)
    print("Aggregate outputs: adni_results.csv and adni_model_counts.csv")


if __name__ == "__main__":
    main()
