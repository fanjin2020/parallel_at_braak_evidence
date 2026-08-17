#!/usr/bin/env python3
"""Create compact manuscript source tables from the completed analysis outputs.

The script does not perform statistical modelling. It checks completed result
files, selects the rows used by the manuscript, and writes Table 1–5 source
tables to the configured manuscript-table directory.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml


def fail(message: str) -> None:
    raise RuntimeError(message)


def read_csv(path: Path, label: str) -> pd.DataFrame:
    if not path.is_file():
        fail(f"Missing {label}: {path}")
    return pd.read_csv(path)


def require_columns(table: pd.DataFrame, columns: list[str], label: str) -> None:
    missing = [column for column in columns if column not in table.columns]
    if missing:
        fail(f"{label} is missing required columns: {', '.join(missing)}")


def stable_pathways(pathways: pd.DataFrame) -> pd.DataFrame:
    # The R pathway script writes the conventional CAMERA column name NGenes.
    # Use one internal spelling below while accepting the source-file spelling.
    if "NGenes" in pathways.columns and "n_genes" not in pathways.columns:
        pathways = pathways.rename(columns={"NGenes": "n_genes"})
    require_columns(pathways, ["go_id", "full_fdr", "full_direction", "n_genes", "directionally_stable"], "ranked pathway output")
    flag = pathways["directionally_stable"].astype(str).str.strip().str.lower().isin({"true", "1", "yes"})
    selected = pathways.loc[flag].copy().sort_values(["full_fdr", "go_id"], kind="stable")
    if selected.shape[0] != 36:
        fail(f"Expected 36 directionally stable GO terms, found {selected.shape[0]}.")
    if selected["go_id"].duplicated().any():
        fail("Directionally stable GO terms are not unique by GO identifier.")
    return selected


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/project.yaml")
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    with config_path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}
    root = config_path.parent.parent
    paths = config.get("paths", {})

    def configured(key: str, default: str) -> Path:
        return (root / paths.get(key, default)).resolve()

    adni_results = read_csv(configured("adni_results", "results/adni/adni_results.csv"), "ADNI results")
    adni_counts = read_csv(configured("adni_model_counts", "results/adni/adni_model_counts.csv"), "ADNI model counts")
    robust_genes = read_csv(configured("gse_robust_genes", "results/gse131617/braak_robust_genes.csv"), "robust gene results")
    pathways = read_csv(configured("stable_pathways", "results/pathways/stable_pathways.csv"), "stable pathway results")
    gene_symbols = read_csv(configured("gse_gene_symbols", "data/contracts/gene_symbols.csv"), "gene-symbol contract")
    term_dictionary = read_csv(configured("go_bp_term_dictionary", "data/contracts/go_bp_term_dictionary.csv"), "GO term dictionary")
    output_dir = configured("manuscript_tables_dir", "results/manuscript_tables")
    output_dir.mkdir(parents=True, exist_ok=True)

    require_columns(adni_results, ["analysis", "outcome", "n_complete_case"], "ADNI results")
    require_columns(adni_counts, ["analysis", "outcome", "at_stage", "n"], "ADNI model counts")
    require_columns(term_dictionary, ["go_id", "term_name"], "GO term dictionary")
    require_columns(gene_symbols, ["feature_id", "entrez_id", "gene_symbol"], "gene-symbol contract")
    if term_dictionary["go_id"].duplicated().any():
        fail("GO term dictionary contains duplicate GO identifiers.")

    table1 = adni_counts.sort_values(["analysis", "outcome"], kind="stable")

    confirmatory = adni_results.loc[adni_results["analysis"].eq("confirmatory_180day")].copy()
    require_columns(confirmatory, ["contrast", "beta", "ci95_lower", "ci95_upper", "p_value", "fdr_bh"], "confirmatory ADNI results")
    if confirmatory.shape[0] != 4:
        fail(f"Expected four confirmatory rows, found {confirmatory.shape[0]}.")
    table2 = confirmatory.sort_values(["outcome", "contrast"], kind="stable")

    fdg = adni_results.loc[adni_results["analysis"].eq("exploratory_FDG")].copy()
    interaction = adni_results.loc[adni_results["analysis"].eq("APOE_joint_interaction")].copy()
    if fdg.shape[0] != 2 or interaction.shape[0] != 2:
        fail("Expected two FDG rows and two APOE interaction rows.")
    table3 = pd.concat([fdg, interaction], ignore_index=True, sort=False)

    require_columns(robust_genes, ["feature_id", "entrez_id", "log2_fold_change", "fdr_bh"], "robust gene results")
    if gene_symbols["feature_id"].duplicated().any():
        fail("Gene-symbol contract contains duplicate feature identifiers.")
    table4 = robust_genes.merge(
        gene_symbols.loc[:, ["feature_id", "gene_symbol"]],
        on="feature_id", how="left", validate="one_to_one"
    )
    if table4["gene_symbol"].isna().any():
        missing = ", ".join(table4.loc[table4["gene_symbol"].isna(), "feature_id"].astype(str).tolist())
        fail(f"Gene-symbol contract has no symbol for feature ID(s): {missing}")
    table4 = table4.sort_values(["log2_fold_change", "gene_symbol"], kind="stable")

    selected_pathways = stable_pathways(pathways)
    table5 = selected_pathways.merge(term_dictionary, on="go_id", how="left", validate="one_to_one")
    if table5["term_name"].isna().any():
        missing = ", ".join(table5.loc[table5["term_name"].isna(), "go_id"].tolist())
        fail(f"GO dictionary has no term name for: {missing}")
    table5 = table5.loc[:, ["go_id", "term_name", "full_direction", "n_genes", "full_fdr"]]

    outputs = {
        "table1_adni_analysis_sets.csv": table1,
        "table2_adni_confirmatory_associations.csv": table2,
        "table3_adni_exploratory_results.csv": table3,
        "table4_braak_robust_genes.csv": table4,
        "table5_stable_go_pathways.csv": table5,
    }
    for name, table in outputs.items():
        table.to_csv(output_dir / name, index=False)

    summary = pd.DataFrame(
        [
            {"item": "confirmatory ADNI rows", "value": len(table2)},
            {"item": "exploratory ADNI rows", "value": len(table3)},
            {"item": "robust genes", "value": len(table4)},
            {"item": "stable GO terms", "value": len(table5)},
        ]
    )
    summary.to_csv(output_dir / "reproduction_summary.csv", index=False)
    print(f"Manuscript source tables written to: {output_dir}")


if __name__ == "__main__":
    main()
