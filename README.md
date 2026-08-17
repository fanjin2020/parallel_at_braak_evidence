# Parallel A/T–Braak Evidence

This repository contains the analysis code used for a study of PET-defined amyloid/tau (A/T) states in ADNI and late-Braak cortical transcriptomic patterns in GSE131617.

The two datasets answer different questions and are analysed independently. ADNI is used for cross-sectional clinical and imaging associations. GSE131617 is used to describe bulk-cortical expression and GO biological-process context associated with late Braak pathology. The repository does not match participants between datasets or use one dataset to validate the other.

## What this release reproduces

- PET-aligned ADNI master table and A/T groups.
- Confirmatory ADNI models for hippocampus/ICV and CDR-SB.
- Exploratory FDG and APOE interaction/stratified analyses.
- GSE131617 feature-level late-Braak differential expression.
- Cross-region direction checks and the robust gene set.
- Ranked GO biological-process analysis and stable pathway set.
- Compact manuscript tables derived from the above results.

The source data and analysis tables are reproduced by this code. Publication artwork is maintained separately from this repository; the scripts here do not regenerate the final formatted figures.

## Data access

ADNI participant-level data are controlled by the ADNI data-use agreement and are not included. Obtain the required ADNI files through the ADNI data-access process, then place them in `data/private/adni/raw/`.

GSE131617 is public. Place the series matrix and donor metadata in `data/public/gse131617/raw/` as specified in `config/project.example.yaml`. Derived GEO files are ignored by Git so that a local run does not add large files to a commit.

## Software

Python 3.10 or later is recommended. Install the Python dependencies in `environment/requirements.txt`.

The R pipeline requires R plus `yaml`, `readxl`, `DBI`, `RSQLite`, `limma`, `AnnotationDbi`, `org.Hs.eg.db`, `GO.db`, and `huex10sttranscriptcluster.db`. From the repository root, install missing R packages with:

```bash
Rscript environment/install_r_packages.R
```

For a strict reproduction of the pathway result, use the R/Bioconductor versions recorded by:

```bash
Rscript tools/check_r_environment.R --report validation/r_environment_check.csv
```

## Configuration

Copy the example configuration before running any script:

```bash
cp config/project.example.yaml config/project.yaml
```

On Windows PowerShell, use:

```powershell
Copy-Item config/project.example.yaml config/project.yaml
```

The local `config/project.yaml` is ignored by Git. Check or change paths and the `runtime.rscript` setting there.

## Run order

Run commands from the repository root. Each step stops when a required input is absent.

```bash
python scripts/prepare_adni.py --config config/project.yaml
python scripts/analyze_adni.py --config config/project.yaml

Rscript scripts/prepare_gse131617.R --config config/project.yaml
Rscript scripts/analyze_gse131617.R --config config/project.yaml
Rscript scripts/prepare_pathway_expression.R --config config/project.yaml
Rscript scripts/analyze_pathways.R --config config/project.yaml
Rscript scripts/summarize_uty_sensitivity.R --config config/project.yaml

python scripts/compile_results_tables.py --config config/project.yaml
```

Before a complete R run, the following lightweight checks are useful:

```bash
python tools/check_environment.py --config config/project.yaml
Rscript tools/check_r_environment.R --report validation/r_environment_check.csv
```

## Main outputs

| Location | Contents |
|---|---|
| `results/adni/adni_model_counts.csv` | ADNI analysis-set counts, including model-eligible counts. |
| `results/adni/adni_results.csv` | Confirmatory and exploratory ADNI estimates. |
| `results/gse131617/` | Feature results, representative genes, and region checks. |
| `results/pathways/` | Full ranked-pathway results and directionally stable terms. |
| `results/gse_sensitivity/` | Donor sex and UTY sensitivity summaries. |
| `results/manuscript_tables/` | Table 1 through Table 5 source tables. |

The APOE-stratified output reports both the total sample used by the interaction model and the actual sample size in each APOE ε4 stratum. These values have different meanings and should not be interchanged.

## Repository boundaries

Do not commit ADNI downloads, participant-level records, raw exports, local configuration, logs, or generated results. The supplied `.gitignore` excludes them. The contracts in `data/contracts/` are small, non-identifiable files required to document analysis inputs and labels.

## Citation and data acknowledgement

Please cite the final study when using this code. Any publication using ADNI data must also follow the current ADNI acknowledgement and publication requirements.

