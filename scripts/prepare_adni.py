#!/usr/bin/env python3
# =============================================================================
# 中文说明
# -----------------------------------------------------------------------------
# 脚本名称：prepare_adni.py
# 作用：从 ADNI 的原始导出 CSV 文件构建以 tau PET 为时间锚点的分析主表。
#
# 输入文件（位于 config/project.yaml 的 paths.adni_raw_dir）：
#   1. UCBERKELEY_AMY_6MM.csv      ：amyloid PET 分类与扫描日期；
#   2. UCBERKELEY_TAU_6MM.csv      ：flortaucipir tau PET、时间 MetaROI SUVR；
#   3. UCSFFSX7.csv                ：左右海马体积和颅内容积；
#   4. CDR.csv                     ：Clinical Dementia Rating Sum of Boxes；
#   5. APOERES.csv                 ：APOE 基因型；
#   6. PTDEMOG.csv                 ：出生日期、性别、教育年限；
#   7. UCBERKELEYFDG_8mm.csv       ：FDG PET MetaROI（探索性结局）。
#
# 主要处理规则：
#   - 仅保留通过质量控制的 amyloid PET 与 FTP tau PET；
#   - 以 amyloid PET Core 分类确定 A+/A−；
#   - 以 temporal MetaROI SUVR >= 1.34 确定 T+/T−；
#   - 在预设窗口内选择 amyloid PET 与 tau PET 最接近的一对扫描；
#   - 以 tau PET 日期为时间锚点，选择窗口内最接近的 MRI、CDR-SB 和 FDG；
#   - 相同距离优先选择较早日期；同日记录若数值冲突，则该结局保留缺失；
#   - hippocampus/ICV = (左海马 + 右海马) / 颅内容积；
#   - 从 APOE 基因型生成 ε4 携带状态；
#   - 使用 tau PET 日期计算年龄，并选取最接近 tau PET 的教育年限记录；
#   - 按三组确认性队列的 hippocampus/ICV 分布标记预定义 MRI 异常值。
#
# 输出：
#   仅生成 paths.adni_master 指定的 adni_pet_aligned_master.csv。
#   输出表每位参与者仅保留一行，供后续 analyze_adni.py 直接建模。
#
# 本脚本不会：
#   - 输出 RID 级个体诊断、影响点或审计文件；
#   - 执行回归、FDR 校正、Huber 回归或 bootstrap；
#   - 将 Braak 转录组结果与 ADNI 个体数据进行匹配或整合。
# =============================================================================
#
# =============================================================================
# English documentation
# -----------------------------------------------------------------------------
# Script: prepare_adni.py
# Purpose: Build one tau-PET-anchored ADNI master table from the raw ADNI
#          export files. This table is the direct input for analyze_adni.py.
#
# Key rules:
#   - Amyloid status is taken from the PET Core composite classification.
#   - Tau positivity is defined as temporal MetaROI SUVR >= 1.34.
#   - Eligible amyloid and tau scans are paired within the pre-specified window.
#   - The tau PET scan is the temporal anchor for MRI, CDR-SB, and FDG MetaROI.
#   - The nearest measurement is selected; ties use the earlier date.
#   - Conflicting measurements on an otherwise identical selected date are set
#     to missing rather than resolved arbitrarily.
#   - Hippocampus/ICV is calculated from bilateral hippocampal volume divided
#     by intracranial volume.
#
# Output:
#   One participant-level master table at paths.adni_master. No RID-level
#   diagnostics, audit files, statistical models, or public release files are
#   created by this preparation script.
# =============================================================================

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="Project YAML file")
    return parser.parse_args()


def read_csv(folder: Path, name: str) -> pd.DataFrame:
    path = folder / name
    if not path.exists():
        raise FileNotFoundError(f"Missing ADNI export: {path}")
    return pd.read_csv(path, low_memory=False)


def need(data: pd.DataFrame, columns: set[str], label: str) -> None:
    missing = sorted(columns.difference(data.columns))
    if missing:
        raise ValueError(f"{label} is missing columns: {', '.join(missing)}")


def date(values: pd.Series) -> pd.Series:
    """Accept only complete calendar dates; do not guess incomplete dates."""
    text = values.astype("string").str.strip()
    out = pd.Series(pd.NaT, index=values.index, dtype="datetime64[ns]")
    for fmt in ("%Y-%m-%d", "%m/%d/%Y", "%m/%d/%y", "%Y/%m/%d"):
        out = out.fillna(pd.to_datetime(text, format=fmt, errors="coerce"))
    return out


def birth_date(values: pd.Series) -> pd.Series:
    """Convert ADNI DOB encodings to a fixed day solely for age calculation."""
    text = values.astype("string").str.strip().str.replace(r"^(\d{4})\.0+$", r"\1", regex=True)
    normal = pd.Series(pd.NA, index=values.index, dtype="string")
    normal.loc[text.str.fullmatch(r"\d{1,2}[/\-]\d{4}", na=False)] = text.str.replace(
        r"^(\d{1,2})[/\-](\d{4})$", r"\2-\1-15", regex=True
    )
    normal.loc[text.str.fullmatch(r"\d{4}[/\-]\d{1,2}", na=False)] = (
        text.str.replace("/", "-", regex=False) + "-15"
    )
    normal.loc[text.str.fullmatch(r"\d{4}", na=False)] = text + "-07-01"
    normal.loc[text.str.fullmatch(r"\d{1,2}/\d{1,2}/\d{4}", na=False)] = text
    out = pd.Series(pd.NaT, index=values.index, dtype="datetime64[ns]")
    for fmt in ("%Y-%m-%d", "%m/%d/%Y"):
        out = out.fillna(pd.to_datetime(normal, format=fmt, errors="coerce"))
    return out


def closest(pairs: pd.DataFrame, measures: pd.DataFrame, value: str, days: int, signature: list[str]) -> pd.DataFrame:
    """Select the nearest valid measurement; ties use the earlier date.

    If equal-date duplicate rows disagree in source values, leave the outcome
    missing rather than choosing arbitrarily.
    """
    records = []
    for rid, tau_date in pairs[["rid", "tau_pet_date"]].itertuples(index=False):
        x = measures.loc[measures.rid.eq(rid)].copy()
        if x.empty or pd.isna(tau_date):
            continue
        x["gap"] = (x.date - tau_date).abs().dt.days
        x = x.loc[x.gap.between(0, days)]
        if x.empty:
            continue
        x = x.loc[x.gap.eq(x.gap.min())]
        x = x.loc[x.date.eq(x.date.min())]
        if x[signature].drop_duplicates().shape[0] != 1:
            continue
        records.append({"rid": rid, value: x.iloc[0][value], f"{value}_days_from_tau": int(x.iloc[0].gap)})
    selected = pd.DataFrame(records, columns=["rid", value, f"{value}_days_from_tau"])
    return pairs.merge(selected, on="rid", how="left")


def apoe_carrier(raw: pd.DataFrame) -> pd.DataFrame:
    need(raw, {"RID", "GENOTYPE"}, "APOERES.csv")
    x = raw[["RID", "GENOTYPE"]].rename(columns={"RID": "rid", "GENOTYPE": "genotype"}).dropna()
    alleles = x.genotype.astype("string").str.strip().str.split("/")
    valid = alleles.map(lambda a: isinstance(a, list) and len(a) == 2 and set(a).issubset({"2", "3", "4"}))
    x = x.loc[valid].copy()
    x["genotype"] = alleles.loc[valid].map(lambda a: "/".join(sorted(a, key=int)))
    x = x.loc[x.groupby("rid").genotype.transform("nunique").eq(1)].drop_duplicates("rid")
    x["apoe4_carrier"] = x.genotype.str.contains("4", regex=False).astype(int)
    return x[["rid", "apoe4_carrier"]]


def demographics(raw: pd.DataFrame, pairs: pd.DataFrame) -> pd.DataFrame:
    need(raw, {"RID", "PTGENDER", "PTEDUCAT", "PTDOB", "VISDATE"}, "PTDEMOG.csv")
    x = raw[["RID", "PTGENDER", "PTEDUCAT", "PTDOB", "VISDATE"]].rename(columns={"RID": "rid", "VISDATE": "date"})
    x["date"], x["dob"] = date(x.date), birth_date(x.PTDOB)
    sex = x.PTGENDER.astype("string").str.upper().str.strip()
    x["sex_male"] = np.select([sex.isin(["M", "MALE", "1"]), sex.isin(["F", "FEMALE", "2"])], [1.0, 0.0], default=np.nan)
    x["education_years"] = pd.to_numeric(x.PTEDUCAT, errors="coerce")

    # Sex and DOB must be internally consistent for a participant.
    dob = x.dropna(subset=["dob"]).loc[lambda d: d.groupby("rid").dob.transform("nunique").eq(1)].drop_duplicates("rid").set_index("rid").dob
    male = x.dropna(subset=["sex_male"]).loc[lambda d: d.groupby("rid").sex_male.transform("nunique").eq(1)].drop_duplicates("rid").set_index("rid").sex_male

    # Education is visit-specific: select the closest date to tau PET.
    out = pairs[["rid", "tau_pet_date"]].merge(x[["rid", "date", "education_years"]], on="rid", how="left")
    out["gap"] = (out.date - out.tau_pet_date).abs().dt.days
    out = out.sort_values(["rid", "gap", "date"], na_position="last").drop_duplicates("rid")
    out["age_at_tau"] = ((out.tau_pet_date - out.rid.map(dob)).dt.days / 365.25).round(4)
    out["sex_male"] = out.rid.map(male)
    return out[["rid", "age_at_tau", "sex_male", "education_years"]]


def main():
    args = arguments()
    config_path = Path(args.config).resolve()
    root = config_path.parents[1]
    with config_path.open(encoding="utf-8") as f:
        config = yaml.safe_load(f)
    raw_dir = root / config["paths"]["adni_raw_dir"]
    output = root / config["paths"]["adni_master"]
    window = int(config["analysis"].get("primary_window_days", 180))
    tau_cutoff = float(config["analysis"].get("tau_positive_threshold", 1.34))

    amy, tau = read_csv(raw_dir, "UCBERKELEY_AMY_6MM.csv"), read_csv(raw_dir, "UCBERKELEY_TAU_6MM.csv")
    need(amy, {"RID", "SCANDATE", "qc_flag", "TRACER", "LONIUID", "AMYLOID_STATUS_COMPOSITE_REF", "CENTILOIDS"}, "amyloid PET")
    need(tau, {"RID", "SCANDATE", "qc_flag", "TRACER", "LONIUID", "META_TEMPORAL_SUVR"}, "tau PET")
    amy = amy.rename(columns={"RID": "rid", "SCANDATE": "amy_pet_date", "LONIUID": "amy_loniuid", "AMYLOID_STATUS_COMPOSITE_REF": "amyloid_positive"})
    amy["amy_pet_date"], amy["amyloid_positive"], amy["centiloids"] = date(amy.amy_pet_date), pd.to_numeric(amy.amyloid_positive, errors="coerce"), pd.to_numeric(amy.CENTILOIDS, errors="coerce")
    amy = amy.loc[amy.qc_flag.eq(2) & amy.amyloid_positive.isin([0, 1]) & amy.amy_pet_date.notna(), ["rid", "amy_pet_date", "TRACER", "amy_loniuid", "amyloid_positive", "centiloids"]]
    key = ["rid", "amy_pet_date", "TRACER", "amy_loniuid"]
    amy["signature"] = amy[["amyloid_positive", "centiloids"]].astype("string").fillna("<NA>").agg("|".join, axis=1)
    amy = amy.loc[amy.groupby(key, dropna=False).signature.transform("nunique").eq(1)].drop_duplicates(key).drop(columns="signature")

    tau = tau.rename(columns={"RID": "rid", "SCANDATE": "tau_pet_date", "LONIUID": "tau_loniuid", "META_TEMPORAL_SUVR": "tau_temporal_suvr"})
    tau["tau_pet_date"], tau["tau_temporal_suvr"] = date(tau.tau_pet_date), pd.to_numeric(tau.tau_temporal_suvr, errors="coerce")
    tau = tau.loc[tau.qc_flag.eq(2) & tau.TRACER.astype("string").eq("FTP") & tau.tau_pet_date.notna() & tau.tau_temporal_suvr.notna(), ["rid", "tau_pet_date", "TRACER", "tau_loniuid", "tau_temporal_suvr"]]
    tau = tau.loc[~tau.duplicated(["rid", "tau_pet_date", "TRACER", "tau_loniuid"], keep=False)]

    pairs = amy.merge(tau, on="rid")
    pairs["at_pair_days"] = (pairs.amy_pet_date - pairs.tau_pet_date).abs().dt.days
    pairs = pairs.loc[pairs.at_pair_days.between(0, window)].sort_values(["rid", "at_pair_days", "tau_pet_date", "amy_pet_date", "amy_loniuid", "tau_loniuid"]).drop_duplicates("rid")
    pairs["tau_positive"] = (pairs.tau_temporal_suvr >= tau_cutoff).astype(int)
    pairs["at_stage"] = np.select([(pairs.amyloid_positive == 0) & (pairs.tau_positive == 0), (pairs.amyloid_positive == 0) & (pairs.tau_positive == 1), (pairs.amyloid_positive == 1) & (pairs.tau_positive == 0), (pairs.amyloid_positive == 1) & (pairs.tau_positive == 1)], ["A-T-", "A-T+", "A+T-", "A+T+"], default="UNCLASSIFIED")
    if pairs.at_stage.eq("UNCLASSIFIED").any():
        raise ValueError("A PET pair could not be assigned to an A/T group.")

    mri = read_csv(raw_dir, "UCSFFSX7.csv")
    need(mri, {"RID", "EXAMDATE", "ST29SV", "ST88SV", "ST10CV"}, "UCSFFSX7.csv")
    mri = mri.rename(columns={"RID": "rid", "EXAMDATE": "date"})[["rid", "date", "ST29SV", "ST88SV", "ST10CV"]]
    mri["date"] = date(mri.date)
    for c in ["ST29SV", "ST88SV", "ST10CV"]: mri[c] = pd.to_numeric(mri[c], errors="coerce")
    mri["hippocampus_icv"] = (mri.ST29SV + mri.ST88SV) / mri.ST10CV
    mri = mri.loc[np.isfinite(mri.hippocampus_icv) & mri.hippocampus_icv.gt(0)]
    pairs = closest(pairs, mri, "hippocampus_icv", window, ["ST29SV", "ST88SV", "ST10CV"])

    cdr = read_csv(raw_dir, "CDR.csv")
    need(cdr, {"RID", "VISDATE", "CDRSB"}, "CDR.csv")
    cdr = cdr.rename(columns={"RID": "rid", "VISDATE": "date", "CDRSB": "cdrsb"})
    cdr["date"], cdr["cdrsb"] = date(cdr.date), pd.to_numeric(cdr.cdrsb, errors="coerce")
    pairs = closest(pairs, cdr.loc[cdr.cdrsb.between(0, 18)], "cdrsb", window, ["cdrsb"])

    fdg = read_csv(raw_dir, "UCBERKELEYFDG_8mm.csv")
    need(fdg, {"RID", "EXAMDATE", "ROINAME", "MEAN"}, "UCBERKELEYFDG_8mm.csv")
    fdg = fdg.rename(columns={"RID": "rid", "EXAMDATE": "date", "MEAN": "fdg_metaroi"})
    fdg = fdg.loc[fdg.ROINAME.astype("string").str.strip().eq("MetaROI")]
    fdg["date"], fdg["fdg_metaroi"] = date(fdg.date), pd.to_numeric(fdg.fdg_metaroi, errors="coerce")
    pairs = closest(pairs, fdg, "fdg_metaroi", window, ["fdg_metaroi"])

    pairs = pairs.merge(apoe_carrier(read_csv(raw_dir, "APOERES.csv")), on="rid", how="left")
    pairs = pairs.merge(demographics(read_csv(raw_dir, "PTDEMOG.csv"), pairs), on="rid", how="left")
    primary = pairs.loc[pairs.at_stage.isin(["A-T-", "A+T-", "A+T+"]), "hippocampus_icv"].dropna()
    mad = (primary - primary.median()).abs().median()
    pairs["hippocampus_outlier_flag"] = False if not np.isfinite(mad) or mad == 0 else ((0.6745 * (pairs.hippocampus_icv - primary.median()) / mad).abs() > 3.5)

    columns = ["rid", "at_stage", "amy_pet_date", "tau_pet_date", "at_pair_days", "tau_temporal_suvr", "age_at_tau", "sex_male", "education_years", "apoe4_carrier", "hippocampus_icv", "hippocampus_icv_days_from_tau", "hippocampus_outlier_flag", "cdrsb", "cdrsb_days_from_tau", "fdg_metaroi", "fdg_metaroi_days_from_tau"]
    output.parent.mkdir(parents=True, exist_ok=True)
    pairs.loc[:, columns].sort_values("rid").to_csv(output, index=False)
    print(f"Prepared ADNI master table: {output}")


if __name__ == "__main__":
    main()























