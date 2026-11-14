"""Centiloid-adjusted tau associations and internal validation near two years.

Only aggregate results are written. Participant records stay in memory.
"""
from __future__ import annotations

import argparse
import json
import platform
import warnings
from collections import Counter
from datetime import timedelta
from pathlib import Path

import numpy as np
import pandas as pd
import patsy
import scipy
from scipy.special import expit
import statsmodels
import statsmodels.api as sm
import longitudinal as core

VERSION = "2026-10-17.1"
DEMOGRAPHIC = "age10 + sex_male + education_c + apoe4_carrier"
CLINICAL = DEMOGRAPHIC + " + baseline_cdrsb + C(diagnosis) + visit_offset"
MODEL_RHS = {
    "clinical": CLINICAL,
    "clinical_amyloid": CLINICAL + " + centiloid10",
    "clinical_amyloid_mri": CLINICAL + " + centiloid10 + hip1000",
    "clinical_amyloid_mri_tau": CLINICAL + " + centiloid10 + hip1000 + tau01",
}
PREDICTOR_COLUMNS = [*core.COVS, "baseline_cdrsb", "centiloids", "hippocampus_icv", "tau_temporal_suvr"]


def dates(series):
    text = series.astype("string").str.strip()
    out = pd.Series(pd.NaT, index=series.index, dtype="datetime64[ns]")
    for fmt in ("%Y-%m-%d", "%m/%d/%Y", "%m/%d/%y", "%Y/%m/%d"):
        out = out.fillna(pd.to_datetime(text, format=fmt, errors="coerce"))
    return out


def resolve(root, supplied, default):
    p = supplied or Path(default)
    return p.resolve() if p.is_absolute() else (root / p).resolve()


def save(out, name, rows):
    pd.DataFrame(rows).to_csv(out / (name + ".csv"), index=False, encoding="utf-8-sig")


def match_centiloids(master, path):
    cols = ["RID", "SCANDATE", "qc_flag", "TRACER", "AMYLOID_STATUS_COMPOSITE_REF", "CENTILOIDS"]
    raw = pd.read_csv(path, usecols=cols, low_memory=False)
    raw["rid"] = raw.RID.map(core.rid)
    raw["scan_date"] = dates(raw.SCANDATE)
    for c in ["qc_flag", "AMYLOID_STATUS_COMPOSITE_REF", "CENTILOIDS"]:
        raw[c] = pd.to_numeric(raw[c], errors="coerce")
    raw = raw.loc[raw.qc_flag.eq(2) & raw.scan_date.notna()]
    groups = {(r, d.date()): g for (r, d), g in raw.groupby(["rid", "scan_date"])}
    output, counts = {}, Counter()
    for p in master:
        key = (p["rid"], core.iso_date(p["amy_pet_date"]))
        g = groups.get(key)
        row = {"centiloids": np.nan, "amyloid_tracer": "Unavailable"}
        status = "no_exact_passed_scan"
        if g is not None:
            unique = g[["CENTILOIDS", "TRACER", "AMYLOID_STATUS_COMPOSITE_REF"]].drop_duplicates()
            if len(unique) != 1:
                status = "ambiguous_same_day_scan"
            elif unique.iloc[0].AMYLOID_STATUS_COMPOSITE_REF != int(p["at_stage"].startswith("A+")):
                status = "amyloid_classification_mismatch"
            else:
                scan = unique.iloc[0]
                # Tracer is scan metadata, even when its Centiloid value is absent.
                if pd.notna(scan.TRACER) and str(scan.TRACER).strip():
                    row["amyloid_tracer"] = str(scan.TRACER).strip()
                if np.isfinite(scan.CENTILOIDS):
                    status = "matched"
                    row["centiloids"] = float(scan.CENTILOIDS)
                else:
                    status = "centiloid_missing"
        row["centiloid_match_status"] = status
        output[p["rid"]] = row
        counts[status] += 1
    return output, [dict(component="centiloid_matching", status=k, n=v) for k, v in counts.items()]


def match_mri_dates(master, path):
    cols = ["RID", "EXAMDATE", "ST29SV", "ST88SV", "ST10CV"]
    raw = pd.read_csv(path, usecols=cols, low_memory=False)
    raw["rid"] = raw.RID.map(core.rid)
    raw["date"] = dates(raw.EXAMDATE)
    for c in cols[2:]:
        raw[c] = pd.to_numeric(raw[c], errors="coerce")
    raw["ratio"] = (raw.ST29SV + raw.ST88SV) / raw.ST10CV
    raw = raw.loc[np.isfinite(raw.ratio) & raw.ratio.gt(0) & raw.date.notna()]
    groups = dict(tuple(raw.groupby("rid")))
    output, counts = {}, Counter()
    for p in master:
        row = {"mri_date": None}
        status = "master_missing_or_flagged"
        if p["hippocampus_icv"] is not None and not p["hip_flag"]:
            g = groups.get(p["rid"])
            status = "no_aligned_mri"
            if g is not None:
                g = g.copy()
                g["gap"] = (g.date - pd.Timestamp(p["anchor"])).abs().dt.days
                g = g.loc[g.gap.le(180)].sort_values(["gap", "date"])
                if len(g):
                    chosen = g.loc[g.date.eq(g.iloc[0].date)]
                    status = "ambiguous_mri"
                    if len(chosen[cols[2:]].drop_duplicates()) == 1:
                        c = chosen.iloc[0]
                        status = "mri_master_mismatch"
                        if np.isclose(c.ratio, p["hippocampus_icv"], rtol=1e-9, atol=1e-12) and c.gap == p["hippocampus_icv_days_from_tau"]:
                            if (c[cols[2:]] > 0).all():
                                status = "matched"
                                row["mri_date"] = c.date.date()
                            else:
                                status = "nonpositive_source_volume"
        row["mri_match_status"] = status
        output[p["rid"]] = row
        counts[status] += 1
    return output, [dict(component="mri_matching", status=k, n=v) for k, v in counts.items()]


def enrich(eligible, amy, mri, dx):
    result = []
    for p in eligible:
        d, status = core.nearest(dx.get(p["rid"], {}), p["anchor"], 180)
        q = {**p, **amy[p["rid"]], **mri[p["rid"]]}
        q["diagnosis"] = {1: "CN", 2: "MCI", 3: "Dementia"}[d[1]] if d else "Unavailable"
        q["diagnosis_date"] = d[0] if d else None
        q["diagnosis_status"] = status
        q["baseline_date"] = q["anchor"] + timedelta(days=q["baseline_offset_days"])
        result.append(q)
    return result


def scaled(df):
    df = df.copy()
    df["tau01"] = (df.tau_temporal_suvr - 1.34) / .1
    df["centiloid10"] = df.centiloids / 10
    df["age10"] = (df.age_at_tau - 75) / 10
    df["education_c"] = df.education_years - 16
    df["hip1000"] = df.hippocampus_icv * 1000
    return df


def endpoint_cohort(people, cdr, year, window_days, landmark_days=180):
    """One observed visit per person after the specified predictor landmark."""
    rows, reasons = [], Counter()
    for p in people:
        reasons["01_baseline_eligible"] += 1
        if any(not np.isfinite(p[k]) if p[k] is not None else True for k in PREDICTOR_COLUMNS):
            reasons["excluded_missing_predictor"] += 1
            continue
        if p["diagnosis"] == "Unavailable" or p["mri_match_status"] != "matched":
            reasons["excluded_diagnosis_or_mri"] += 1
            continue
        landmark = p["anchor"] + timedelta(days=landmark_days)
        availability = [core.iso_date(p["amy_pet_date"]), p["baseline_date"], p["mri_date"], p["diagnosis_date"], p["anchor"]]
        if any(d is None or d > landmark for d in availability):
            reasons["excluded_predictor_after_landmark"] += 1
            continue
        reasons["02_complete_landmark_predictors"] += 1
        candidates = [(abs((d - p["anchor"]).days - year * 365.25), d, values)
                      for d, values in cdr.get(p["rid"], {}).items()
                      if d > landmark and abs((d - p["anchor"]).days - year * 365.25) <= window_days]
        if not candidates:
            reasons["excluded_no_visit_in_window"] += 1
            continue
        _, d, values = min(candidates, key=lambda x: (x[0], x[1]))
        if len(values) != 1:
            reasons["excluded_conflicting_nearest_outcome"] += 1
            continue
        actual = (d - p["anchor"]).days / 365.25
        rows.append({**p, "outcome_cdrsb": next(iter(values)), "actual_year": actual,
                     "visit_offset": actual - year, "year": year})
        reasons["03_observed_endpoint"] += 1
    df = pd.DataFrame(rows)
    if len(df):
        if df.rid.duplicated().any():
            raise ValueError("Endpoint cohort contains repeated participants")
        df = scaled(df)
    return df, [dict(year=year, step=k, n=v) for k, v in reasons.items()]


def fit_endpoint(x, y, fractional):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("Nonfinite model input")
    if np.linalg.matrix_rank(x) < x.shape[1]:
        raise ValueError("Rank deficient model")
    if len(y) <= x.shape[1] + 5:
        raise ValueError("Insufficient residual degrees of freedom")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        fit = (sm.GLM(y / 18, x, family=sm.families.Binomial()).fit(maxiter=200, cov_type="HC3")
               if fractional else sm.OLS(y, x).fit(cov_type="HC3"))
    if not getattr(fit, "converged", True) or not np.isfinite(fit.params).all() or not np.isfinite(fit.cov_params()).all():
        raise ValueError("Model did not converge to finite estimates")
    if any("PerfectSeparation" in type(w.message).__name__ for w in caught):
        raise ValueError("Perfect separation")
    return fit


def predict(fit, x, fractional):
    return np.asarray(fit.predict(np.asarray(x))) * (18 if fractional else 1)


def tau_contrast(fit, info, df, fractional):
    low, high = df.tau_temporal_suvr.quantile([.25, .75])
    if high <= low:
        raise ValueError("No tau IQR variation")
    xs = []
    for t in (low, high):
        z = df.copy()
        z["tau01"] = (t - 1.34) / .1
        z["visit_offset"] = 0.0
        xs.append(np.asarray(patsy.build_design_matrices([info], z)[0]))
    b = np.asarray(fit.params)
    if fractional:
        lo, hi = (expit(x @ b) for x in xs)
        estimate = 18 * (hi - lo).mean()
        gradient = 18 * ((hi * (1-hi))[:, None] * xs[1] - (lo * (1-lo))[:, None] * xs[0]).mean(axis=0)
    else:
        gradient = (xs[1] - xs[0]).mean(axis=0)
        estimate = gradient @ b
    se = np.sqrt(max(0., gradient @ np.asarray(fit.cov_params()) @ gradient))
    units = (high-low) / .1
    return dict(tau_q25=low, tau_q75=high, iqr_difference=estimate,
                iqr_ci_low=estimate-1.96*se, iqr_ci_high=estimate+1.96*se,
                difference_per_01=estimate/units, ci_low_per_01=(estimate-1.96*se)/units,
                ci_high_per_01=(estimate+1.96*se)/units,
                interpretation="Adjusted CDR-SB level contrast at target year, holding baseline severity fixed; not an annual slope or causal effect")


def endpoint_models(df, year):
    summaries, coefficients, contrasts, failures = [], [], [], []
    if len(df) < 40:
        return [], [], [], [dict(year=year, component="endpoint_models", status="insufficient_participants", n=len(df))]
    y = df.outcome_cdrsb.to_numpy()
    for family in ("gaussian", "fractional_logit"):
        fractional = family == "fractional_logit"
        for label, rhs in MODEL_RHS.items():
            try:
                x = patsy.dmatrix(rhs, df, return_type="dataframe")
                f = fit_endpoint(x, y, fractional)
                pred = predict(f, x, fractional)
                ci = np.asarray(f.conf_int())
                for j, term in enumerate(x.columns):
                    coefficients.append(dict(year=year, family=family, model=label, term=term,
                        beta=f.params[j], ci_low=ci[j, 0], ci_high=ci[j, 1], nominal_p=f.pvalues[j],
                        scale="logit mean CDR-SB fraction" if fractional else "CDR-SB points"))
                summaries.append(dict(year=year, family=family, model=label, n=len(df), n_parameters=x.shape[1],
                    formula="outcome_cdrsb ~ " + rhs, condition_number=np.linalg.cond(x),
                    in_sample_rmse=np.sqrt(np.mean((y-pred)**2)), fitted_min=pred.min(), fitted_max=pred.max(),
                    fitted_outside_0_18_fraction=np.mean((pred < 0) | (pred > 18)), covariance="HC3"))
                if label.endswith("_tau"):
                    contrasts.append(dict(year=year, family=family, n=len(df), **tau_contrast(f, x.design_info, df, fractional)))
            except (ValueError, np.linalg.LinAlgError) as e:
                failures.append(dict(year=year, family=family, component=label, status="failed", reason=str(e)))
    return summaries, coefficients, contrasts, failures


def metrics(y, pred):
    y, pred = np.asarray(y), np.asarray(pred)
    err = y-pred
    sst = np.sum((y-y.mean())**2)
    calibration = np.linalg.lstsq(np.column_stack([np.ones(len(y)), pred]), y, rcond=None)[0]
    return dict(rmse=np.sqrt(np.mean(err**2)), mae=np.mean(abs(err)), r2=1-np.sum(err**2)/sst if sst > 0 else np.nan,
                mean_error=np.mean(pred-y), calibration_intercept=calibration[0], calibration_slope=calibration[1],
                prediction_outside_0_18_fraction=np.mean((pred < 0) | (pred > 18)))


def paired_cv(df, folds, repeats, bootstrap, seed):
    """No tuning or imputation. One row per participant; all models share splits."""
    if df.rid.duplicated().any():
        raise ValueError("CV requires one endpoint per participant")
    y = df.outcome_cdrsb.to_numpy()
    labels = ["clinical_amyloid_mri", "clinical_amyloid_mri_tau"]
    rng = np.random.default_rng(seed)
    by_repeat, warnings_out = [], []
    predictions = {family: [] for family in ("gaussian", "fractional_logit")}
    for repeat in range(repeats):
        fold_id = np.empty(len(df), int)
        fold_id[rng.permutation(len(df))] = np.arange(len(df)) % folds
        for family in predictions:
            paired = np.full((len(df), 2), np.nan)
            valid = True
            for fold in range(folds):
                train, test = df.iloc[np.flatnonzero(fold_id != fold)], df.iloc[np.flatnonzero(fold_id == fold)]
                if set(train.rid) & set(test.rid):
                    raise AssertionError("Participant leakage")
                for j, label in enumerate(labels):
                    try:
                        x = patsy.dmatrix(MODEL_RHS[label], train, return_type="dataframe")
                        xt = patsy.build_design_matrices([x.design_info], test)[0]
                        f = fit_endpoint(x, train.outcome_cdrsb.to_numpy(), family == "fractional_logit")
                        paired[fold_id == fold, j] = predict(f, xt, family == "fractional_logit")
                    except (ValueError, patsy.PatsyError, np.linalg.LinAlgError) as e:
                        warnings_out.append(dict(component="cross_validation", family=family, repeat=repeat+1,
                                                 fold=fold+1, model=label, status="failed", reason=str(e).splitlines()[0]))
                        valid = False
            if not valid or not np.isfinite(paired).all():
                continue
            predictions[family].append(paired)
            vals = [metrics(y, paired[:, j]) for j in range(2)]
            by_repeat.append(dict(family=family, repeat=repeat+1, n=len(df),
                **{"base_"+k:v for k,v in vals[0].items()}, **{"tau_"+k:v for k,v in vals[1].items()},
                delta_rmse=vals[1]["rmse"]-vals[0]["rmse"], delta_mae=vals[1]["mae"]-vals[0]["mae"],
                delta_r2=vals[1]["r2"]-vals[0]["r2"]))
        if (repeat+1) % 5 == 0:
            print(f"Internal validation: {repeat+1}/{repeats} repeats", flush=True)
    summaries = []
    for family, values in predictions.items():
        n_valid = len(values)
        if n_valid != repeats:
            summaries.append(dict(family=family, status="incomplete_cv_do_not_interpret", requested_repeats=repeats,
                                  valid_repeats=n_valid, n=len(df)))
            continue
        a = np.asarray(values)  # repeat, participant, model
        squared = (a-y[None,:,None])**2
        absolute = abs(a-y[None,:,None])
        rng_b = np.random.default_rng(seed+701)
        draws = []
        for _ in range(bootstrap):
            idx = rng_b.integers(0, len(df), len(df))
            rmse = np.sqrt(squared[:,idx,:].mean(axis=1)).mean(axis=0)
            mae = absolute[:,idx,:].mean(axis=(0,1))
            draws.append([rmse[1]-rmse[0], mae[1]-mae[0]])
        interval = np.quantile(draws, [.025, .975], axis=0)
        selected = [r for r in by_repeat if r["family"] == family]
        summaries.append(dict(family=family, status="ok", n=len(df), folds=folds, valid_repeats=n_valid,
            **{k:np.mean([r[k] for r in selected]) for k in selected[0] if k not in ("family","repeat","n")},
            delta_rmse_ci_low=interval[0,0], delta_rmse_ci_high=interval[1,0],
            delta_mae_ci_low=interval[0,1], delta_mae_ci_high=interval[1,1],
            interval_method="Participant bootstrap of paired fixed out-of-fold losses; models NOT refitted in bootstrap",
            limitation="Conditional intervals omit model-refitting uncertainty; CV repeats are not independent samples; no external validation",
            interpretation="Negative delta RMSE/MAE favors added tau; positive delta R2 favors added tau"))
    return summaries, by_repeat, warnings_out


def amyloid_gee(long, amy):
    df = pd.DataFrame(long)
    df["centiloids"] = df.rid.map({k:v["centiloids"] for k,v in amy.items()})
    df = df.dropna(subset=["centiloids"]).copy()
    df = scaled(df)
    rows, failures = [], []
    for adjusted in (False, True):
        label = "same_sample_plus_centiloid" if adjusted else "same_sample_original"
        rhs = "time_years * (tau01 + " + DEMOGRAPHIC + (" + centiloid10" if adjusted else "") + ")"
        try:
            x = patsy.dmatrix(rhs, df, return_type="dataframe")
            if np.linalg.matrix_rank(x) != x.shape[1]:
                raise ValueError("Rank deficient GEE")
            f = sm.GEE(df.cdr, x, groups=df.rid, family=sm.families.Gaussian(), cov_struct=sm.cov_struct.Independence()).fit(maxiter=200)
            if not f.converged or not np.isfinite(f.cov_params()).all().all():
                raise ValueError("Invalid GEE fit")
            t = "time_years:tau01"
            ci = f.conf_int().loc[t]
            rows.append(dict(model=label, n_participants=df.rid.nunique(), n_observations=len(df),
                beta=f.params[t], ci_low=ci.iloc[0], ci_high=ci.iloc[1], nominal_p=f.pvalues[t],
                scale="CDR-SB points/year per 0.1 tau SUVR", formula="cdr ~ "+rhs,
                purpose="Replay of original association on identical Centiloid-complete sample; not prediction"))
        except (ValueError, np.linalg.LinAlgError) as e:
            failures.append(dict(component=label, status="failed", reason=str(e)))
    return rows, failures


def describe_followup(people, cdr, window):
    rows = []
    q25, q75 = np.quantile([p["tau_temporal_suvr"] for p in people], [.25,.75])
    for grouping in ("all", "diagnosis", "tau"):
        groups = {}
        for p in people:
            level = "all" if grouping == "all" else (p["diagnosis"] if grouping == "diagnosis" else
                    ("lower_quartile" if p["tau_temporal_suvr"] <= q25 else "upper_quartile" if p["tau_temporal_suvr"] >= q75 else "middle_half"))
            groups.setdefault(level, []).append(p)
        for level, ps in groups.items():
            spans, near = [], {1:0,2:0,3:0}
            for p in ps:
                times = [(d-p["anchor"]).days/365.25 for d,v in cdr.get(p["rid"],{}).items()
                         if len(v)==1 and d > max(p["anchor"],p["baseline_date"]) and (d-p["anchor"]).days <= 3*365.25]
                if times:
                    spans.append(max(times))
                for year in near:
                    near[year] += any(abs(t-year)*365.25 <= window for t in times)
            rows.append(dict(grouping=grouping, group=level, n_baseline=len(ps), n_with_followup=len(spans),
                last_visit_median_years=np.median(spans) if spans else np.nan,
                last_visit_q25=np.quantile(spans,.25) if spans else np.nan,
                last_visit_q75=np.quantile(spans,.75) if spans else np.nan,
                last_visit_max=max(spans) if spans else np.nan,
                n_near_1year=near[1], n_near_2years=near[2], n_near_3years=near[3], window_days=window,
                note="Time from tau PET; follow-up capped at 3 years; near-3-year window is right-truncated"))
    return rows


def endpoint_selection(people, df, year):
    selected = set(df.rid) if len(df) else set()
    table = pd.DataFrame(people)
    table["selected"] = table.rid.isin(selected)
    for diagnosis in ("CN", "MCI", "Dementia"):
        table["diagnosis_"+diagnosis] = table.diagnosis.eq(diagnosis).astype(float)
    output = []
    for field in [*PREDICTOR_COLUMNS, "diagnosis_CN", "diagnosis_MCI", "diagnosis_Dementia"]:
        a = pd.to_numeric(table.loc[table.selected, field],errors="coerce").dropna()
        b = pd.to_numeric(table.loc[~table.selected, field],errors="coerce").dropna()
        binary = field in ("sex_male", "apoe4_carrier") or field.startswith("diagnosis_")
        variance = (a.mean()*(1-a.mean()) + b.mean()*(1-b.mean()))/2 if binary else (a.var(ddof=1)+b.var(ddof=1))/2
        scale = np.sqrt(variance) if variance > 0 else np.nan
        output.append(dict(year=year,variable=field,endpoint_n=len(a),other_eligible_n=len(b),
            endpoint_mean=a.mean(),other_eligible_mean=b.mean(),standardized_difference=(a.mean()-b.mean())/scale))
    return output


def tau_support(df, year):
    if df.empty:
        return []
    low, high = df.tau_temporal_suvr.quantile([.25,.75])
    rows=[]
    for diagnosis,g in df.groupby("diagnosis"):
        p05,p95 = g.tau_temporal_suvr.quantile([.05,.95])
        rows.append(dict(year=year,diagnosis=diagnosis,n=len(g),tau_p05=p05,tau_p95=p95,
            overall_q25=low,overall_q75=high,q25_outside_stage_p05_p95=not p05<=low<=p95,
            q75_outside_stage_p05_p95=not p05<=high<=p95))
    return rows


def run(a):
    root = a.project_root.resolve()
    raw = resolve(root, a.raw_dir, "data/private/adni/raw")
    master_path = resolve(root, a.master, "data/private/adni/adni_pet_aligned_master.csv")
    out = resolve(root, a.output_dir, "results/revision_step1")
    if out == root or root/"data" == out or root/"data" in out.parents or out == raw or raw in out.parents:
        raise ValueError("Output cannot be an input directory")
    if out.exists() and any(out.iterdir()):
        raise ValueError("Output directory is not empty; choose a new directory")
    if not 2 <= a.folds <= 10 or a.repeats < 1 or a.bootstrap < 100 or not 1 <= a.window_days <= 180:
        raise ValueError("Invalid folds/repeats/bootstrap/window parameters")
    master, _ = core.load_master(master_path)
    ids = {p["rid"] for p in master}
    cdr, _ = core.load_visits(raw/"CDR.csv", ids, "cdr")
    dx, _ = core.load_visits(raw/"DXSUM.csv", ids, "dx")
    eligible, long, _ = core.make_cohort(master, cdr)
    amy, matching = match_centiloids(master, raw/"UCBERKELEY_AMY_6MM.csv")
    mri, mri_matching = match_mri_dates(master, raw/"UCSFFSX7.csv")
    people = enrich(eligible, amy, mri, dx)
    out.mkdir(parents=True)
    save(out, "input_matching", matching+mri_matching)
    save(out, "cohort_counts", [dict(master_n=len(master), baseline_eligible_n=len(eligible),
        original_followup_n=len({r['rid'] for r in long}), original_observations_n=len(long))])
    save(out, "followup_support", describe_followup(people,cdr,a.window_days))
    gee, failures = amyloid_gee(long, amy)
    save(out, "centiloid_adjusted_gee", gee)
    samples, models, coefs, contrasts, descriptions, selection, support = [], [], [], [], [], [], []
    main_endpoint = None
    for year in (1,2):
        df, flow = endpoint_cohort(people,cdr,year,a.window_days)
        samples.extend(flow)
        selection.extend(endpoint_selection(people,df,year))
        support.extend(tau_support(df,year))
        if len(df):
            for dx_name,g in df.groupby("diagnosis"):
                descriptions.append(dict(year=year, diagnosis=dx_name, n=len(g),
                    outcome_mean=g.outcome_cdrsb.mean(), baseline_mean=g.baseline_cdrsb.mean(),
                    tau_mean=g.tau_temporal_suvr.mean(), centiloid_mean=g.centiloids.mean(),
                    observed_year_min=g.actual_year.min(), observed_year_max=g.actual_year.max()))
        if year == 2:
            main_endpoint = df
        m,c,t,f = endpoint_models(df,year)
        models.extend(m); coefs.extend(c); contrasts.extend(t); failures.extend(f)
    for name, rows in [("endpoint_sample_flow",samples),("endpoint_models",models),("endpoint_coefficients",coefs),
                       ("endpoint_tau_contrasts",contrasts),("endpoint_description",descriptions),
                       ("endpoint_selection_comparison",selection),("tau_support_by_diagnosis",support)]:
        save(out,name,rows)
    cv_summary, cv_repeats = [], []
    if len(main_endpoint) >= a.min_cv_n:
        cv_summary,cv_repeats,f = paired_cv(main_endpoint,a.folds,a.repeats,a.bootstrap,a.seed)
        failures.extend(f)
    else:
        failures.append(dict(component="cross_validation",status="insufficient_participants", n=len(main_endpoint), minimum=a.min_cv_n))
    save(out,"prediction_comparison",cv_summary)
    save(out,"prediction_by_repeat",cv_repeats)
    save(out,"analysis_status",failures or [dict(component="all_requested_analyses",status="ok")])
    metadata = dict(script_version=VERSION, parameters={k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()},
        versions=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,
                      statsmodels=statsmodels.__version__,patsy=patsy.__version__),
        primary_endpoint="Observed CDR-SB nearest 2 years after tau PET, within requested window; one visit per participant",
        landmark="Predictor assessment window closes 180 days after tau PET; not a prediction available on the PET date",
        primary_prediction_comparison="clinical_amyloid_mri versus clinical_amyloid_mri_tau; fractional-logit primary, Gaussian sensitivity",
        secondary_endpoint="Observed visit near 1 year; descriptive/model sensitivity only; not cross-validated",
        sample="Complete-case A+ participants with observed endpoint; no correction for endpoint observation selection",
        hippocampus="Original aligned hippocampus/ICV ratio, scaled by 1000; not regression-based ICV adjustment",
        inference="All analyses added during revision, exploratory; P values nominal; intervals omit full model-selection uncertainty",
        cv_interval="Conditional paired participant bootstrap of fixed out-of-fold errors, not a full refit bootstrap",
        output_privacy="Aggregate outputs only; no participant-level predictions or identifiers written")
    (out/"analysis_specification.json").write_text(json.dumps(metadata,indent=2,ensure_ascii=False),encoding="utf-8")
    print("Completed:",out,flush=True)
    print("Two-year observed endpoint N:",len(main_endpoint),flush=True)
    if cv_summary:
        print(pd.DataFrame(cv_summary).to_string(index=False),flush=True)
    print("Inspect analysis_status.csv before interpreting results.",flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project-root",type=Path,required=True)
    p.add_argument("--master",type=Path)
    p.add_argument("--raw-dir",type=Path)
    p.add_argument("--output-dir",type=Path,default=Path("results/revision_step1"))
    p.add_argument("--window-days",type=int,default=180)
    p.add_argument("--folds",type=int,default=5)
    p.add_argument("--repeats",type=int,default=20)
    p.add_argument("--bootstrap",type=int,default=2000)
    p.add_argument("--min-cv-n",type=int,default=60)
    p.add_argument("--seed",type=int,default=20261017)
    run(p.parse_args())
