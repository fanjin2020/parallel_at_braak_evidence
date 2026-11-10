"""Covariate imputation and bounded outcome scenarios with explicit targets."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm
from statsmodels.imputation.mice import MICEData
from scipy.stats import t as student_t
from scipy.special import expit
import pet_staging_covariates as adni
import cross_sectional as cross
import analyze_tau_incremental as inc
from pet_day_clinical_models import match_landmark
import longitudinal as core


def impute_covariates(root, out, m, seed):
    df, _, _ = adni.load_data(root)
    rows, diagnostic, status = [], [], []
    np.random.seed(seed)
    for outcome in cross.OUTCOMES:
        d = df.loc[df.at_stage.isin(["A-T-", "A+T-", "A+T+"])]
        d = d.loc[d[outcome].notna() & d[outcome+"_days_from_tau"].between(0, 180)]
        if outcome == "hippocampus_icv":
            d = d.loc[~d.hippocampus_outlier_flag.astype(str).str.lower().isin(["true", "1", "1.0"])]
        # No MRI, PET, outcome, or entirely unavailable variable is imputed.
        fixed = ["tau_per_01", "centiloid10", outcome]
        d = d.dropna(subset=fixed).copy()
        variables = [*cross.COVARIATES, *fixed]
        x = d[variables].astype(float).copy()
        missing = [c for c in cross.COVARIATES if x[c].isna().any()]
        for c in variables:
            diagnostic.append(dict(outcome=outcome, variable=c, n=len(x), missing=int(x[c].isna().sum())))
        if not missing:
            status.append(dict(component="MI_"+outcome, status="not_needed", reason="No missing demographic/genotype covariates in fixed-outcome target"))
            continue
        if any(x[c].notna().sum() < 30 for c in missing):
            status.append(dict(component="MI_"+outcome, status="insufficient_observed_covariates"))
            continue
        try:
            mice = MICEData(x, perturbation_method="boot", k_pmm=5)
            for c in missing:
                predictors = [k for k in variables if k != c]
                if len(predictors) >= 2:
                    mice.set_imputer(c, formula=" + ".join(predictors), model_class=sm.OLS,
                                     perturbation_method="boot", k_pmm=5)
            mice.update_all(10)
            estimates, variances = [], []
            for j in range(m):
                mice.update_all(5)
                z = mice.data
                design = patsy.dmatrix("tau_per_01 + centiloid10 + "+adni.COV, z, return_type="dataframe")
                fit = sm.OLS(z[outcome], design).fit(cov_type="HC3")
                if np.linalg.matrix_rank(design) != design.shape[1]:
                    raise ValueError("Imputed design rank deficient")
                estimates.append(float(fit.params["tau_per_01"]))
                variances.append(float(fit.cov_params().loc["tau_per_01", "tau_per_01"]))
                diagnostic.append(dict(outcome=outcome, variable="tau_coefficient", imputation=j+1,
                                       value=estimates[-1], variance=variances[-1]))
            q, within, between = np.mean(estimates), np.mean(variances), np.var(estimates, ddof=1)
            total = within+(1+1/m)*between
            dof = (m-1)*(1+within/((1+1/m)*between))**2 if between > 0 else np.inf
            crit = student_t.ppf(.975, dof)
            complete = d.dropna(subset=variables)
            reference = adni.fit_effect(complete, outcome, "tau_per_01 + centiloid10 + "+adni.COV, "tau_per_01")
            rows.append(dict(outcome=outcome, n_target=len(d), n_complete=len(complete), imputations=m,
                beta=q, ci_low=q-crit*np.sqrt(total), ci_high=q+crit*np.sqrt(total),
                total_variance=total, within_variance=within, between_variance=between,
                rubin_df=dof, complete_case_beta=reference["beta"], complete_case_ci_low=reference["ci_low"],
                complete_case_ci_high=reference["ci_high"],
                estimand="Full three-state outcome-observed, PET/CL-complete continuous-tau slope; not the primary A+ slope"))
            status.append(dict(component="MI_"+outcome, status="completed"))
        except (ValueError, np.linalg.LinAlgError) as e:
            status.append(dict(component="MI_"+outcome, status="failed", reason=str(e)))
    adni.write(out, "MI_pooled_results", rows)
    adni.write(out, "MI_diagnostics", diagnostic)
    return status


def mnar(root, out, status_file, bootstrap, seed):
    if not status_file:
        return [dict(component="MNAR_survivors", status="requires_data",
                     reason="Need rid, alive_at_window_end, status_date; no death or survival status inferred from missing visits")]
    survival = pd.read_csv(status_file, dtype={"rid": str})
    required = {"rid", "alive_at_window_end", "status_date"}
    if not required.issubset(survival) or survival.rid.duplicated().any():
        raise ValueError("Invalid survival table; one row per rid required")
    survival.rid = survival.rid.map(core.rid)
    survival["status_date"] = pd.to_datetime(survival.status_date, errors="raise")
    people, cdr, _, _ = match_landmark(root, 180)
    endpoint, _ = inc.endpoint_cohort(people, cdr, 2, 180, 180)
    target = pd.DataFrame(people)
    target = target.merge(survival, on="rid", validate="one_to_one")
    target = target.loc[target.alive_at_window_end.eq(1) &
                        target.status_date.ge(pd.to_datetime(target.anchor)+pd.to_timedelta(2*365.25+180, unit="D"))]
    target = target.merge(endpoint[["rid", "outcome_cdrsb"]], on="rid", how="left", validate="one_to_one")
    target = inc.scaled(target).reset_index(drop=True)
    if len(target) < 40 or target.outcome_cdrsb.notna().sum() < 40:
        return [dict(component="MNAR_survivors", status="insufficient_participants", n=len(target))]
    rhs = inc.MODEL_RHS["clinical_amyloid_mri_tau"].replace(" + visit_offset", "")
    rng = np.random.default_rng(seed)
    rows, failures = [], []
    for b in range(bootstrap+1):
        d = target if b == 0 else target.iloc[rng.integers(0, len(target), len(target))].reset_index(drop=True)
        try:
            observed = d.outcome_cdrsb.notna()
            x = patsy.dmatrix(rhs, d, return_type="dataframe")
            fit = inc.fit_endpoint(x.loc[observed], d.loc[observed, "outcome_cdrsb"], True)
            mu = inc.predict(fit, x, True)
            # Delta shifts on the logit-mean scale respect the 0..18 support.
            for delta in [-1., -.5, 0., .5, 1.]:
                completed = d.outcome_cdrsb.to_numpy().copy()
                completed[~observed] = 18*expit(np.asarray(x)[~observed]@fit.params+delta)
                f = inc.fit_endpoint(x, completed, True)
                z = d.copy()
                z["visit_offset"] = 0.
                effect = inc.tau_contrast(f, x.design_info, z, True)
                rows.append(dict(replicate=b, delta_logit=delta, n=len(d), n_observed=int(observed.sum()),
                    n_missing=int((~observed).sum()),
                    mean_missing_shift=np.mean(completed[~observed]-mu[~observed]) if (~observed).any() else 0.,
                    estimate=effect["difference_per_01"]))
        except (ValueError, np.linalg.LinAlgError, patsy.PatsyError) as e:
            failures.append(dict(component="MNAR_survivors", replicate=b, status="failed", reason=str(e).splitlines()[0]))
    adni.write(out, "MNAR_draws", rows)
    summary = []
    if rows:
        for delta, d in pd.DataFrame(rows).groupby("delta_logit"):
            draws = d.loc[d.replicate.gt(0), "estimate"]
            original = d.loc[d.replicate.eq(0)]
            summary.append(dict(delta_logit=delta, requested=bootstrap, valid=len(draws),
                estimate=original.estimate.iloc[0] if len(original) else np.nan,
                ci_low=draws.quantile(.025) if len(draws)>=.9*bootstrap else np.nan,
                ci_high=draws.quantile(.975) if len(draws)>=.9*bootstrap else np.nan))
    adni.write(out, "MNAR_summary", summary)
    return failures+[dict(component="MNAR_survivors", status="completed",
        reason="Conditional-mean completion plus participant bootstrap; explicit survivor target, delta not estimated from data")]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project-root", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--survival-status", type=Path)
    p.add_argument("--imputations", type=int, default=20)
    p.add_argument("--bootstrap", type=int, default=500)
    p.add_argument("--seed", type=int, default=20261020)
    a = p.parse_args()
    out = adni.output_dir(a.output_dir)
    status = impute_covariates(a.project_root, out, a.imputations, a.seed)
    status += mnar(a.project_root, out, a.survival_status, a.bootstrap, a.seed)
    adni.write(out, "analysis_status", status)
    (out/"methods.txt").write_text(
        "MI: demographic/genotype covariates only, predictive mean matching, outcome included as auxiliary predictor, "
        "Rubin pooling of HC3 estimates. Requires MAR; structural/outcome/PET/MRI and entirely absent variables not imputed.\n"
        "MNAR: documented survivors through endpoint window, fixed logit-mean shifts, two-year conditional CDR contrast. "
        "Conditional-mean completion is a sensitivity calculation, not recovered individual outcomes; bootstrap refits "
        "the observed-mean and completed-target models. No mortality model or identification of MNAR is claimed.\n",
        encoding="utf-8")


if __name__ == "__main__":
    main()
