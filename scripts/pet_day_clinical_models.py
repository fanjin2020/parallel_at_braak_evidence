"""Clinical timing, common-sample support and classification sensitivities."""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import timedelta
import json
from pathlib import Path
import numpy as np
import pandas as pd
import patsy
from statsmodels.stats.multitest import multipletests
import longitudinal as core
import analyze_tau_incremental as inc
import longitudinal_model_comparison as longitudinal
import pet_staging_covariates as adni
import regional_staging as staging
import cross_sectional as cross


def choose_before(records, anchor, cutoff, signature, lookback=180):
    """Nearest record to PET within a fixed lookback and information cutoff."""
    if records is None or not len(records):
        return None, "no_record"
    d = records.loc[records.date.ge(pd.Timestamp(anchor-timedelta(days=lookback))) &
                    records.date.le(pd.Timestamp(cutoff))].copy()
    if not len(d):
        return None, "no_record_before_cutoff"
    d["gap"] = (d.date-pd.Timestamp(anchor)).abs().dt.days
    d = d.sort_values(["gap", "date"], kind="stable")
    same = d.loc[d.date.eq(d.iloc[0].date)]
    if len(same[signature].drop_duplicates()) != 1:
        return None, "conflicting_nearest_record"
    return same.iloc[0], "available"


def match_landmark(root, landmark_days):
    raw = root/"data/private/adni/raw"
    master, _ = core.load_master(root/"data/private/adni/adni_pet_aligned_master.csv")
    ids = {p["rid"] for p in master}
    cdr, _ = core.load_visits(raw/"CDR.csv", ids, "cdr")
    dx, _ = core.load_visits(raw/"DXSUM.csv", ids, "dx")
    amy = pd.read_csv(raw/"UCBERKELEY_AMY_6MM.csv", low_memory=False)
    amy["rid"], amy["date"] = amy.RID.map(core.rid), inc.dates(amy.SCANDATE)
    for c in ["qc_flag", "CENTILOIDS", "AMYLOID_STATUS_COMPOSITE_REF"]:
        amy[c] = pd.to_numeric(amy[c], errors="coerce")
    amy = amy.loc[amy.qc_flag.eq(2) & amy.CENTILOIDS.notna() &
                  amy.AMYLOID_STATUS_COMPOSITE_REF.isin([0, 1])]
    mri = pd.read_csv(raw/"UCSFFSX7.csv", low_memory=False)
    mri["rid"], mri["date"] = mri.RID.map(core.rid), inc.dates(mri.EXAMDATE)
    vols = ["ST29SV", "ST88SV", "ST10CV"]
    mri[vols] = mri[vols].apply(pd.to_numeric, errors="coerce")
    mri["ratio"] = (mri.ST29SV+mri.ST88SV)/mri.ST10CV
    mri = mri.loc[mri[vols].gt(0).all(axis=1) & mri.ratio.between(0, 1, inclusive="neither")]
    if "OVERALLQC" in mri:
        mri = mri.loc[~mri.OVERALLQC.astype(str).str.strip().str.lower().isin(["fail", "failed", "reject", "rejected"])]
    mri_signature = vols+[k for k in ["FSVER", "FIELD_STRENGTH", "OVERALLQC"] if k in mri]
    ag, mg = dict(tuple(amy.groupby("rid"))), dict(tuple(mri.groupby("rid")))
    people, flow, lags = [], Counter(), []
    for p in master:
        flow["01_original_PET_aligned_cohort"] += 1
        if any(p[k] is None for k in core.COVS):
            flow["excluded_missing_demographic_genotype"] += 1
            continue
        cutoff = p["anchor"]+timedelta(days=landmark_days)
        aa, reason = choose_before(ag.get(p["rid"]), p["anchor"], cutoff,
                                   ["CENTILOIDS", "AMYLOID_STATUS_COMPOSITE_REF", "TRACER"])
        if aa is None:
            flow["excluded_amyloid_"+reason] += 1
            continue
        if aa.AMYLOID_STATUS_COMPOSITE_REF != 1:
            flow["excluded_A_negative_at_cutoff"] += 1
            continue
        mm, reason = choose_before(mg.get(p["rid"]), p["anchor"], cutoff, mri_signature)
        if mm is None:
            flow["excluded_MRI_"+reason] += 1
            continue
        chosen = {}
        for label, source in [("baseline", cdr), ("diagnosis", dx)]:
            records = pd.DataFrame([dict(date=pd.Timestamp(d), value=v)
                                    for d, values in source.get(p["rid"], {}).items()
                                    for v in values], columns=["date", "value"])
            chosen[label], reason = choose_before(records, p["anchor"], cutoff, ["value"])
            if chosen[label] is None:
                flow["excluded_"+label+"_"+reason] += 1
                break
        if any(chosen.get(k) is None for k in ["baseline", "diagnosis"]):
            continue
        b, d = chosen["baseline"], chosen["diagnosis"]
        q = {**p, "amy_pet_date": aa.date.date().isoformat(), "centiloids": float(aa.CENTILOIDS),
             "hippocampus_icv": float(mm.ratio), "mri_date": mm.date.date(),
             "mri_match_status": "matched", "baseline_cdrsb": float(b.value),
             "baseline_date": b.date.date(), "diagnosis_date": d.date.date(),
             "diagnosis": {1: "CN", 2: "MCI", 3: "Dementia"}[int(d.value)],
             "baseline_offset_days": (b.date.date()-p["anchor"]).days,
             "at_stage": "A+T+" if p["tau_temporal_suvr"] >= 1.34 else "A+T-"}
        people.append(q)
        flow["02_complete_rematched_predictors"] += 1
        lags.append({k: (day-p["anchor"]).days for k, day in [
            ("amyloid", aa.date.date()), ("MRI", mm.date.date()),
            ("CDR", b.date.date()), ("diagnosis", d.date.date())]})
    lag_summary = []
    for k in ["amyloid", "MRI", "CDR", "diagnosis"]:
        v = [d[k] for d in lags]
        if v:
            lag_summary.append(dict(variable=k, n=len(v), minimum=min(v), median=np.median(v),
                                    maximum=max(v), landmark_days=landmark_days))
    return people, cdr, [dict(step=k, n=v, landmark_days=landmark_days) for k, v in flow.items()], lag_summary


def run_prediction(root, out, bootstrap, repeats, seed):
    # Actual outcome-visit timing is not available at the prediction cutoff.
    old = dict(inc.MODEL_RHS)
    inc.MODEL_RHS = {k: v.replace(" + visit_offset", "") for k, v in old.items()}
    statuses, flows, lags = [], [], []
    try:
        for landmark in [0, 180]:
            people, cdr, flow, lag = match_landmark(root, landmark)
            flows += flow
            lags += lag
            d, flow = inc.endpoint_cohort(people, cdr, 2, 180, landmark)
            flows += [dict(**r, landmark_days=landmark) for r in flow]
            if len(d) < 40:
                statuses.append(dict(component=f"prediction_{landmark}", status="insufficient_participants",
                                     n=len(d), reason="Implementation minimum 40; not evidence of no association"))
                continue
            prefix = f"landmark_{landmark}"
            model, coef, contrasts, failure = inc.endpoint_models(d, 2)
            for name, rows in [("models", model), ("coefficients", coef), ("contrasts", contrasts)]:
                adni.write(out, prefix+"_"+name, rows)
            summary, by_repeat, failure_cv = inc.paired_cv(d, 5, repeats, bootstrap, seed)
            refit, draws = longitudinal.refit_cv_uncertainty(d, bootstrap, 5, seed+17)
            for name, rows in [("CV", summary), ("CV_repeats", by_repeat),
                               ("refitted_CV", refit), ("refitted_CV_draws", draws)]:
                adni.write(out, prefix+"_"+name, rows)
            statuses += failure+failure_cv
            statuses.append(dict(component=prefix, status="completed", n=len(d)))
    finally:
        inc.MODEL_RHS = old
    adni.write(out, "landmark_flow", flows)
    adni.write(out, "predictor_date_lags", lags)
    return statuses


def common_support(root, out):
    people, long, cdr = longitudinal.prepare(root)
    common = long.dropna(subset=["centiloids", "hippocampus_icv", "baseline_cdrsb"])
    common = common.loc[common.diagnosis.ne("Unavailable") & common.mri_match_status.eq("matched")]
    ids = set(common.rid)
    chosen = [p for p in people if p["rid"] in ids]
    for name, table in longitudinal.followup_tables(chosen, cdr).items():
        adni.write(out, "common_"+name, table)
    support = []
    common_q25, common_q75 = common.drop_duplicates("rid").tau_temporal_suvr.quantile([.25, .75])
    for diagnosis, d in common.drop_duplicates("rid").groupby("diagnosis"):
        support.append(dict(diagnosis=diagnosis, n=len(d), tau_min=d.tau_temporal_suvr.min(),
                            tau_q25=d.tau_temporal_suvr.quantile(.25),
                            tau_median=d.tau_temporal_suvr.median(),
                            tau_q75=d.tau_temporal_suvr.quantile(.75), tau_max=d.tau_temporal_suvr.max(),
                            common_q25=common_q25, common_q75=common_q75,
                            common_q25_within_stage_range=d.tau_temporal_suvr.min()<=common_q25<=d.tau_temporal_suvr.max(),
                            common_q75_within_stage_range=d.tau_temporal_suvr.min()<=common_q75<=d.tau_temporal_suvr.max(),
                            fraction_within_common_IQR=d.tau_temporal_suvr.between(common_q25,common_q75).mean()))
    adni.write(out, "common_tau_support_by_stage", support)
    adni.write(out, "common_sample_size", [dict(n_people=len(ids), n_observations=len(common),
               expected_reference_people=185, expected_reference_observations=466,
               matches_reference=len(ids)==185 and len(common)==466)])
    adni.write(out, "common_model_observation_totals", [
        dict(observation_type="baseline", n_people=common.loc[common.time_years.eq(0), "rid"].nunique(),
             n_observations=int(common.time_years.eq(0).sum())),
        dict(observation_type="followup", n_people=common.loc[common.time_years.gt(0), "rid"].nunique(),
             n_observations=int(common.time_years.gt(0).sum()))])
    return dict(component="common_support", status="completed", n=len(ids))


def mri_qc_coverage(root, out):
    raw = pd.read_csv(root/"data/private/adni/raw/UCSFFSX7.csv", low_memory=False)
    fields = [k for k in raw if any(word in k.upper() for word in ["QC", "FSVER", "FIELD_STRENGTH"])]
    coverage = []
    for k in fields:
        counts = raw[k].fillna("Not_recorded").astype(str).value_counts()
        coverage += [dict(field=k, value=v, rows=int(n)) for v, n in counts.items()]
    adni.write(out, "MRI_available_QC_fields", coverage)
    names = ["FSVER", "FIELD_STRENGTH", "OVERALLQC"]
    if set(names).issubset(raw):
        d = raw[names].fillna("Not_recorded").astype(str)
        adni.write(out, "MRI_QC_by_version_field", d.groupby(names, dropna=False).size().rename("source_rows").reset_index())
    return dict(component="MRI_QC_coverage", status="completed",
                reason="Raw labels reported by software version; missing or unknown codes are not treated as pass or fail. Visual QC requires images.")


def cutoff_sensitivity(df, out):
    rows, states, errors = [], [], []
    for threshold in [1.30, 1.32, 1.34, 1.36, 1.38]:
        for band in [0., .02]:
            x = df.copy()
            x["at_stage"] = np.char.add(np.where(x.amyloid_positive.eq(1), "A+", "A-"), np.where(x.tau_temporal_suvr.ge(threshold), "T+", "T-"))
            keep = (x.tau_temporal_suvr-threshold).abs().gt(band) if band else pd.Series(True, index=x.index)
            for (old, new), g in pd.DataFrame({"old": df.at_stage, "new": x.at_stage, "keep": keep}).groupby(["old", "new"]):
                states.append(dict(threshold=threshold, exclusion_band=band, original=old, new=new,
                                   n=len(g), n_retained=int(g.keep.sum())))
            for y in cross.OUTCOMES:
                d = cross.eligible(x.loc[keep], y).dropna(subset=["centiloids"])
                rhs = "C(at_stage, Treatment(reference='A-T-')) + centiloid10 + "+adni.COV
                for group in ["A+T-", "A+T+"]:
                    term = f"C(at_stage, Treatment(reference='A-T-'))[T.{group}]"
                    try:
                        v = adni.fit_effect(d, y, rhs, term)
                        rows.append(dict(outcome=y, threshold=threshold, exclusion_band=band, contrast=group, **v))
                    except (ValueError, np.linalg.LinAlgError, patsy.PatsyError) as e:
                        errors.append(dict(component="cutoff", outcome=y, threshold=threshold,
                                           exclusion_band=band, contrast=group, status="not_estimable", reason=str(e)))
    if rows:
        t = pd.DataFrame(rows)
        planned = 5*2*len(cross.OUTCOMES)*2
        t["q_all_cutoff_sensitivity_tests"] = multipletests(
            np.r_[t.p, np.ones(planned-len(t))], method="fdr_bh")[1][:len(t)]
        adni.write(out, "cutoff_effects", t)
    adni.write(out, "cutoff_classification_changes", states)
    return errors


def stage_config(t):
    return dict(mtl_column="tau_entorhinal", neocortex_column="tau_neocortex",
                mtl_cutoff=t["MTL"]["cutoff"], neocortex_moderate_cutoff=t["neocortex_moderate"]["cutoff"],
                neocortex_high_cutoff=t["neocortex_high"]["cutoff"],
                source="Current-pipeline mixture", pipeline_compatibility="Internal exploratory estimates")


def stage_uncertainty(df, out, bootstrap, seed):
    base = df.dropna(subset=["tau_entorhinal", "tau_neocortex"]).copy()
    threshold = staging.thresholds(base)
    original = adni.biological_stage(base, stage_config(threshold))
    support = []
    for label, d in base.assign(regional_stage=original).groupby("regional_stage"):
        for variable in ["tau_temporal_suvr", "tau_entorhinal", "tau_neocortex"]:
            support.append(dict(regional_stage=label, variable=variable, n=len(d),
                minimum=d[variable].min(), q25=d[variable].quantile(.25),
                median=d[variable].median(), q75=d[variable].quantile(.75), maximum=d[variable].max()))
    adni.write(out, "regional_stage_tau_support", support)
    rng = np.random.default_rng(seed)
    effects, migrations, thresholds, errors = [], [], [], []
    for b in range(bootstrap+1):
        d = base.copy() if b == 0 else base.iloc[rng.integers(0, len(base), len(base))].copy().reset_index(drop=True)
        try:
            t = staging.thresholds(d)
            d["regional_stage"] = adni.biological_stage(d, stage_config(t))
            for k, v in t.items():
                thresholds.append(dict(replicate=b, region=k, cutoff=v["cutoff"]))
            moved = adni.biological_stage(base, stage_config(t))
            for (a, z), n in pd.DataFrame({"original": original, "new": moved}).value_counts().items():
                migrations.append(dict(replicate=b, original=a, new=z, n=int(n)))
            for y in cross.OUTCOMES:
                x = cross.eligible(d, y).dropna(subset=["centiloids"])
                x = x.loc[x.regional_stage.isin(["A", "B", "C", "D"])]
                rhs = "C(regional_stage, Treatment(reference='A')) + centiloid10 + "+adni.COV
                for group in ["B", "C", "D"]:
                    term = f"C(regional_stage, Treatment(reference='A'))[T.{group}]"
                    try:
                        f = adni.fit_effect(x, y, rhs, term)
                        effects.append(dict(replicate=b, outcome=y, contrast=group+"_vs_A", **f))
                    except (ValueError, np.linalg.LinAlgError, patsy.PatsyError) as e:
                        errors.append(dict(component="stage_outcome", replicate=b, outcome=y,
                                           contrast=group+"_vs_A", status="not_estimable", reason=str(e)))
        except (ValueError, np.linalg.LinAlgError, patsy.PatsyError) as e:
            errors.append(dict(component="stage_threshold", replicate=b, status="not_estimable", reason=str(e)))
        if b and b % 50 == 0:
            print("Threshold and outcome bootstrap", b, "/", bootstrap, flush=True)
    adni.write(out, "stage_effect_draws", effects)
    adni.write(out, "stage_threshold_draws", thresholds)
    adni.write(out, "stage_migrations", migrations)
    summary = []
    if effects:
        for (y, contrast), d in pd.DataFrame(effects).groupby(["outcome", "contrast"]):
            v = d.loc[d.replicate.gt(0), "beta"]
            estimate = d.loc[d.replicate.eq(0), "beta"]
            summary.append(dict(outcome=y, contrast=contrast, requested=bootstrap, valid=len(v),
                estimate=estimate.iloc[0] if len(estimate) else np.nan,
                ci_low=v.quantile(.025) if len(v)>=.9*bootstrap else np.nan,
                ci_high=v.quantile(.975) if len(v)>=.9*bootstrap else np.nan,
                status="ok" if len(v)>=.9*bootstrap else "too_many_failed_bootstraps"))
    adni.write(out, "stage_propagated_intervals", summary)
    return errors


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project-root", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--bootstrap", type=int, default=500)
    p.add_argument("--repeats", type=int, default=20)
    p.add_argument("--seed", type=int, default=20261020)
    a = p.parse_args()
    out = adni.output_dir(a.output_dir)
    status = []
    for name, function in [
        ("common_support", lambda: [common_support(a.project_root, out)]),
        ("MRI_QC_coverage", lambda: [mri_qc_coverage(a.project_root, out)]),
        ("prediction", lambda: run_prediction(a.project_root, out, a.bootstrap, a.repeats, a.seed)),
        ("cutoffs", lambda: cutoff_sensitivity(adni.load_data(a.project_root)[0], out)),
        ("regional_staging", lambda: stage_uncertainty(adni.load_data(a.project_root)[0], out, a.bootstrap, a.seed)),
    ]:
        try:
            status += function()
            status.append(dict(component=name, status="completed"))
        except (ValueError, KeyError, FileNotFoundError, np.linalg.LinAlgError, patsy.PatsyError) as e:
            status.append(dict(component=name, status="failed", reason=str(e)))
    adni.write(out, "analysis_status", status)
    (out/"methods.json").write_text(json.dumps(dict(
        predictor_window="PET-180 days through landmark; nearest PET-date valid record, earlier date breaks ties",
        landmarks_days=[0, 180], prediction_target="PET+2 years, nearest observed CDR within +/-180 days after landmark",
        prediction_covariates="Actual future visit offset excluded from all prospective CV models",
        static_covariates="Existing master demographics/genotype; collection-time availability not newly verified",
        MRI="Positive volumes, ratio in (0,1), explicit Fail/Failed/Reject/Rejected excluded. Missing/unknown QC not assumed pass; original post-PET outlier flag not transferred",
        cutoff_grid="Fixed 1.30/1.32/1.34/1.36/1.38 with 0/.02 exclusion bands; sensitivity, not validated replacements",
        staging="Each participant bootstrap re-estimates thresholds, reclassifies and fits outcomes; nominal exploratory intervals",
        scope="Original PET-aligned cohort reanalysis, not a new prospective recruitment cohort"), indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
