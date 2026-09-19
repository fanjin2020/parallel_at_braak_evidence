"""Select region-consistent GO terms, excluding obsolete GO:0090309."""
import pandas as pd

EXCLUDED_GO_IDS = {"GO:0090309"}
POLICY_VERSION = "2026-09-15"
POLICY_SOURCE = "https://amigo.geneontology.org/amigo/term/GO%3A0090309"

def submission_pathways(frame):
    required = {"go_id", "full_fdr", "full_direction", "directionally_stable"}
    if not required.issubset(frame.columns):
        raise ValueError(f"Missing GO columns: {sorted(required-set(frame.columns))}")
    if frame.go_id.duplicated().any():
        raise ValueError("Duplicate GO identifiers")
    flag = frame.directionally_stable.astype(str).str.lower().isin(["true", "1", "yes"])
    selected = frame.loc[flag & ~frame.go_id.isin(EXCLUDED_GO_IDS)].copy()
    if not (pd.to_numeric(selected.full_fdr) < .05).all():
        raise ValueError("Retained term fails original full-model FDR threshold")
    for col in ["direction_leave_out_EC", "direction_leave_out_FC", "direction_leave_out_TC"]:
        if col in selected and not selected[col].eq(selected.full_direction).all():
            raise ValueError(f"Direction mismatch: {col}")
    return selected.sort_values(["full_fdr", "go_id"], kind="stable")
