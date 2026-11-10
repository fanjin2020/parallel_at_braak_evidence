"""Plot model sensitivity and observed support; every panel exports its source table."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 9, "svg.fonttype": "none", "pdf.fonttype": 42,
    "axes.spines.top": False, "axes.spines.right": False})
COLORS = ["#444444", "#a86b38", "#447c9a", "#85599a"]


def save(fig, out, name, data, caption, captions):
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(out/(name+"."+ext), dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    data.to_csv(out/(name+"_source.csv"), index=False)
    captions[name] = caption


def run(root,save_hook=None):
    out = root/"figures"
    if save_hook is None:out.mkdir(exist_ok=True)
    captions = {}
    def emit(fig,out,name,data,caption,captions):
        if save_hook is None:save(fig,out,name,data,caption,captions)
        else:save_hook(fig,name,{'source':data},name,caption)
    def model_label(model):
        if model=='original':return 'Provided-matrix baseline'
        if model=='scan_date':return 'Scan date'
        kind,reference,mode=model.split('_',2)
        prefix='Scan date + ' if kind=='joint' else ''
        return prefix+reference+' ('+{'joint_quantile':'quantile','native_logscale':'native scale'}[mode]+')'
    file = root/"tissue/target_effects_full_background.csv"
    if file.exists():
        d = pd.read_csv(file)
        d = d.loc[d.gene_symbol.ne("UTY")]
        models = ["original", "scan_date", "composition_CA_joint_quantile", "joint_CA_joint_quantile"]
        labels = ["Provided-matrix baseline", "Scan date", "CA (quantile)", "Scan date + CA (quantile)"]
        genes = sorted(d.gene_symbol.unique())
        fig, ax = plt.subplots(figsize=(7.2, 8.4))
        for j, (model, label) in enumerate(zip(models, labels)):
            s = d.loc[d.model.eq(model)].set_index("gene_symbol").reindex(genes)
            y = np.arange(len(genes))+(j-1.5)*.17
            good = s[["log2FC", "ci_low", "ci_high"]].notna().all(axis=1)
            ax.errorbar(s.loc[good, "log2FC"], y[good], xerr=np.vstack([
                s.loc[good, "log2FC"]-s.loc[good, "ci_low"],
                s.loc[good, "ci_high"]-s.loc[good, "log2FC"]]), fmt="o", markersize=3,
                linewidth=.8, color=COLORS[j], label=label)
        ax.set_yticks(range(len(genes)), genes, fontstyle="italic")
        ax.invert_yaxis()
        ax.axvline(0, color="#777777", lw=.7)
        ax.set_xlabel("Late Braak versus Braak 0 log2 fold change")
        ax.legend(loc="lower center", bbox_to_anchor=(.5, 1.01), ncol=2, frameon=False)
        emit(fig, out, "Enhanced_01_tissue_models", d.loc[d.model.isin(models)],
             "Selected non-Y genes extracted after full-feature donor-blocked limma fitting. "
             "Points and bars: estimates and nominal moderated 95% CIs, conditional on original selection "
             "and estimated composition. All-feature q values are in the source table; CI overlap is not "
             "a formal test of model differences. Scan date is not RNA extraction batch.", captions)
    file = root/"tissue/GO_all_models.csv"
    if file.exists():
        d = pd.read_csv(file)
        summary = []
        for model, g in d.groupby("model", sort=False):
            for direction in ["Down", "Up"]:
                summary.append(dict(model=model, direction=direction, n=int(((g.q<.05)&g.direction.eq(direction)).sum()),
                                    tested=len(g)))
        s = pd.DataFrame(summary)
        planned = ["original", "scan_date"]
        for reference in ["CA", "NG"]:
            for mode in ["joint_quantile", "native_logscale"]:
                planned += ["composition_"+reference+"_"+mode, "joint_"+reference+"_"+mode]
        models = [m for m in planned if m in set(s.model)]
        fig, ax = plt.subplots(figsize=(8, max(4, .45*len(models))))
        for j, direction in enumerate(["Down", "Up"]):
            g = s.loc[s.direction.eq(direction)].set_index("model").reindex(models)
            ax.barh(np.arange(len(models))+(j-.5)*.32, g.n, height=.3,
                    color=["#447c9a", "#b78145"][j], label=direction+"-ranked")
            for k, value in enumerate(g.n):
                if value == 0:
                    ax.text(0, k+(j-.5)*.32, "0", va="center", ha="left",
                            color=["#447c9a", "#b78145"][j], fontsize=7)
        ax.set_yticks(range(len(models)), [model_label(m) for m in models])
        ax.invert_yaxis()
        ax.set_xlabel("GO sets with within-model BH q < 0.05")
        ax.set_xlim(0, max(1, float(s.n.max())*1.12))
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.legend(frameon=False)
        fig.tight_layout()
        emit(fig, out, "Enhanced_02_GO_models", s,
             "All eligible GO biological-process sets are retested in every model, with separate full-family BH adjustment. "
             "Counts are descriptive and are not independent evidence units. Down/Up refers to ranked association, not measured pathway activation.", captions)
    file = root/"clinical/common_observed_visit_periods.csv"
    if file.exists():
        d = pd.read_csv(file)
        s = d.loc[d.grouping.eq("all") & d.group.eq("all")]
        fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6), sharey=True)
        periods = ["0_to_1", "1_to_2", "over_2"]
        for ax, clock in zip(axes, ["baseline_CDR", "tau_PET"]):
            x = s.loc[s.clock.eq(clock)].set_index("period").reindex(periods)
            ax.bar(range(3), x.visits, color=COLORS[2], width=.6)
            for j, row in enumerate(x.itertuples()):
                if pd.notna(row.visits):
                    ax.text(j, row.visits, f"{int(row.people)} people", ha="center", va="bottom", fontsize=8)
            ax.set_xticks(range(3), ["0–1", "1–2", ">2"])
            ax.set_title("From "+clock.replace("_", " "))
            ax.set_xlabel("Years")
            ax.margins(y=.2)
        axes[0].set_ylabel("Observed follow-up assessments")
        fig.tight_layout()
        emit(fig, out, "Enhanced_04_common_followup_support", s,
             "Observed follow-up in the common-covariate trajectory sample. Bars count assessments; labels "
             "count unique participants within each interval. Participants can appear in several intervals. "
             "The two panels use different time origins. This is interval support, not an anniversary-window count.", captions)
    file = root/"clinical/stage_propagated_intervals.csv"
    if file.exists() and file.stat().st_size > 10:
        d = pd.read_csv(file)
        if len(d) and {"ci_low", "ci_high"}.issubset(d):
            outcomes = list(d.outcome.unique())
            fig, axes = plt.subplots(1, len(outcomes), figsize=(7.2, 3.8), squeeze=False)
            for i,(ax, outcome) in enumerate(zip(axes.flat, outcomes)):
                s = d.loc[d.outcome.eq(outcome)].dropna(subset=["estimate", "ci_low", "ci_high"])
                y = np.arange(len(s))
                scale=1e4 if outcome=='hippocampus_icv' else 1
                ax.hlines(y, s.ci_low*scale, s.ci_high*scale, color=COLORS[2])
                ax.plot(s.estimate*scale, y, "o", color=COLORS[2])
                ax.set_yticks(y, s.contrast.str.replace("_vs_", " vs ", regex=False))
                ax.invert_yaxis()
                ax.set_title(chr(97+i)+'  '+{"cdrsb": "CDR-SB", "hippocampus_icv": "Hippocampus/ICV"}.get(outcome, outcome),loc='left')
                ax.axvline(0, color="#777777", lw=.7)
                ax.set_xlabel("Adjusted stage difference"+(r" ($\times 10^{-4}$)" if scale!=1 else ""))
                ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
            fig.tight_layout()
            emit(fig, out, "Enhanced_03_stage_uncertainty", d,
                 "Exploratory regional-stage differences. Each donor/participant bootstrap draw re-estimates "
                 "the PET cutoffs, assigns stages and fits the outcome model. Bars show percentile 95% intervals; "
                 "intervals are withheld when fewer than 90% of requested refits succeed. No clinical-stage validation is claimed.", captions)
    if save_hook is None:(out/"captions.json").write_text(json.dumps(captions, indent=2), encoding="utf-8")


