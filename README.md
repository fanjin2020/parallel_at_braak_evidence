# PET A/T states and late-Braak transcriptomic context

Analysis code for *PET-defined amyloid and tau states and late-Braak transcriptomic context in Alzheimer’s disease*.

The ADNI analyses examine PET-defined amyloid/tau (A/T) states and continuous tau burden in relation to hippocampal volume and CDR-SB. Additional exploratory analyses examine subsequent CDR-SB change among amyloid-positive participants. The GSE131617 analyses examine cortical differential expression and GO biological processes associated with late Braak pathology. The two datasets are analyzed separately, without participant-level matching.

The workflow produces statistical results and source tables for the manuscript. It does not generate the submission Word documents or the final publication-layout figures.

## Code versions

The submitted manuscript cites the original code archive: [Zenodo, version 1.0.0](https://doi.org/10.5281/zenodo.21971630).

A subsequent code release is archived as [Zenodo, version 3.0.0](https://doi.org/10.5281/zenodo.22846279).

[This GitHub repository](https://github.com/fanjin2020/parallel_at_braak_evidence) contains the maintained code and may include changes made after the archived releases. To reproduce an archived version, use the files associated with its version-specific DOI. Record the Git commit or release tag when running code from GitHub.

## Project layout

```text
scripts/                     Python and R analysis scripts
config/project.example.yaml  Configuration template
data/contracts/              Field definitions and gene/GO dictionaries
environment/requirements.txt Python dependencies
run_analysis.sh              Full-workflow launcher for Linux
LICENSE                      Software license
```

All Python and R analysis scripts are stored directly in `scripts/`. If applying the V3 code update package to an existing project, retain the original `data/contracts/` directory and prepare the input data below.

## Input data

### ADNI

Obtain data through the [ADNI data-access process](https://adni.loni.usc.edu/). Participant-level data are not distributed in this repository. Place the following files in `data/private/adni/raw/`:

```text
UCBERKELEY_AMY_6MM.csv
UCBERKELEY_TAU_6MM.csv
UCSFFSX7.csv
CDR.csv
APOERES.csv
PTDEMOG.csv
UCBERKELEYFDG_8mm.csv
DXSUM.csv
ADSL.csv
DATA_DOWNLOADED_DATE.csv
```

`DATA_DOWNLOADED_DATE.csv` supplies the date used to assess calendar follow-up opportunity. It must contain a `data_downloaded_date` column in `YYYY-MM-DD` format, corresponding to the actual data download batch.

### GSE131617

Obtain the expression matrix and donor information from [GEO accession GSE131617](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131617). Place them in `data/public/gse131617/raw/`. The configuration template uses these filenames:

```text
GSE131617-GPL5175_series_matrix.txt.gz
GSE131617_subject_info.xlsx
```

Update the configuration if the downloaded filenames differ. Feature mapping also requires the HuEx annotation package listed below.

Keep the dictionaries in `data/contracts/`. Changing annotation versions can change feature mapping, pathway membership, and results.

## Software

The analysis records specify Python 3.12.13 and statsmodels 0.14.6. `environment/requirements.txt` uses version ranges rather than a complete environment lock; save the installed package versions with each run.

```bash
python -m pip install -r environment/requirements.txt
```

The R workflow uses the following core versions:

| Component | Version |
|---|---|
| R | 4.4.2 |
| Bioconductor | 3.20 |
| limma | 3.62.2 |
| AnnotationDbi | 1.68.0 |
| org.Hs.eg.db | 3.20.0 |
| GO.db | 3.20.0 |
| huex10sttranscriptcluster.db | 8.8.0 |

Other dependencies include yaml, readxl, DBI, RSQLite, and dependencies used by limma. From the repository root, use the `Rscript` executable for R 4.4.2:

```bash
Rscript scripts/install_r_packages.R
Rscript scripts/check_r_environment.R
```

Run the installer when preparing the environment. The version check stops if a checked core version differs. The full-workflow launcher checks package availability but does not replace this exact-version check.

## Configuration

Run the commands below from the repository root. Create a local configuration:

```bash
cp config/project.example.yaml config/project.yaml
```

Check the data paths and `runtime.rscript` in `config/project.yaml`. Relative input paths are resolved from the project root; absolute paths can be used for data stored elsewhere. Do not commit the local configuration.

## Run the full workflow

Activate the Python environment containing the dependencies, then run:

```bash
export RSCRIPT=/path/to/R-4.4.2/bin/Rscript
bash run_analysis.sh
```

Replace the example path with the actual R 4.4.2 executable. The launcher uses `python3` by default; set the `PYTHON` environment variable to use another executable.

The launcher checks the inputs and runs the analyses with 2,000 participant bootstrap replicates, 1,000 follow-up-weighting bootstrap replicates, and random seed `20260915`. Results are written to a new timestamped directory:

```text
results/paper_run_YYYYMMDD_HHMMSS/
```

To set the output directory and parameters explicitly:

```bash
python scripts/run_paper.py \
  --project-root . \
  --output-dir results/paper_run_v3_01 \
  --rscript /path/to/R-4.4.2/bin/Rscript \
  --bootstrap 2000 \
  --selection-bootstrap 1000 \
  --seed 20260915
```

Use a new or empty output directory for each run. Add `--dry-run` to check inputs and print the commands without fitting models. This mode does not complete the R package checks.

## Scripts

| Script | Purpose |
|---|---|
| `prepare_adni.py` | Prepare ADNI input tables and the PET-aligned analysis dataset |
| `analyze_adni.py`, `adni_models.py` | Fit the original A/T cross-sectional, FDG, APOE, and sensitivity models |
| `cross_sectional.py` | Estimate direct A+T+ versus A+T− contrasts and continuous-tau associations within A+ participants |
| `tau_sensitivity.py`, `spline_stability.py` | Examine nonlinearity, influential observations, and spline stability |
| `longitudinal.py` | Analyze tau burden and subsequent CDR-SB change in A+ participants |
| `clinical_tables.py` | Summarize baseline characteristics, longitudinal cohorts, and follow-up availability |
| `longitudinal_sensitivity.py` | Run participant resampling, calendar-opportunity checks, and outcome-model sensitivity analyses |
| `selection_sensitivity.py` | Estimate follow-up-availability weights, assess balance, and bootstrap the weighting workflow |
| `stage_sensitivity.py` | Adjust for clinical stage and compare Gaussian and fractional-logit models |
| `prepare_gse131617.R` | Prepare expression data, feature annotation, and donor information |
| `analyze_gse131617.R` | Fit donor-blocked differential-expression and leave-one-region-out models |
| `prepare_pathway_expression.R` | Construct the Entrez gene-level matrix for pathway analysis |
| `analyze_pathways.R` | Run ranked GO analyses and assess directional agreement across region deletions |
| `summarize_uty_sensitivity.R` | Summarize donor sex composition and the gene list after excluding UTY |
| `go_submission_policy.py` | Apply the GO reporting rules |
| `compile_results_tables.py` | Compile original-analysis source tables and provide the clinical-only entry point |
| `run_additional.py` | Run the additional clinical analyses in sequence |
| `run_paper.py` | Run the complete workflow |
| `collect_paper_results.py` | Collect manuscript source tables and create the results archive |
| `install_r_packages.R`, `check_r_environment.R` | Install R dependencies and check core versions |

The UTY summary is a gene-list check, not a male-only model refit. Excluding obsolete term `GO:0090309` from the reported GO list does not recalculate FDR for the original test family.

## Generated outputs

The paths below describe files generated by a completed run, not a required set of files in the GitHub download. The analysis can be run without precomputed results once the input data and dependencies are available.

Within a completed run directory:

| Relative path | Contents |
|---|---|
| `results/adni/` | Original ADNI model estimates and sample sizes |
| `results/gse131617/` | Feature-level differential expression and region-sensitivity results |
| `results/pathways/` | Pathway tests and reported terms |
| `results/gse_sensitivity/` | Donor sex and UTY summaries |
| `results/manuscript_tables/` | Compact source tables for the original analyses |
| `results/extensions/` | Additional cross-sectional, longitudinal, weighting, and clinical-stage analyses |
| `paper_results/` | Collected manuscript source tables and their mapping |
| `paper_results.zip` | Archive of the collected results |

The collector creates `paper_results/paper_table_sources.csv` to map manuscript items to source files. It collects the feature-level S2 results, donor summaries, clinical tables, model estimates, and sensitivity results from the completed run. Use the generated collection when checking table sources rather than assuming that any precomputed results included in a repository download are complete.

Some source filenames retain earlier table numbers. Use the mapping file and analysis names to identify the corresponding final manuscript tables. Weighted GEE confidence intervals and full-workflow bootstrap intervals are saved separately and should not be interchanged. Fractional-logit standardized changes are time-specific contrasts, not constant annual slopes.

To collect results from a completed run without refitting models:

```bash
python scripts/collect_paper_results.py \
  --run-dir results/paper_run_v3_01 \
  --output-dir results/paper_run_v3_01/paper_results_new
```

`--run-dir` must point to the full run directory containing `data/` and `results/`, not an extracted results-sharing archive. Choose an output path that does not conflict with an existing collection or ZIP file.

## Data sharing and license

Do not upload ADNI downloads, participant-level derived tables, local configuration files, virtual environments, or the full run directory. The full run directory contains restricted inputs and derivatives.

The collector limits its output to listed files and checks for explicit participant-identifier columns. Review the files against the applicable data-use agreement before sharing them; this check is not a substitute for that review. GSE131617 remains available through GEO. ADNI data use and acknowledgements must follow the ADNI requirements.

The code is distributed under the MIT license; see `LICENSE`. This software license does not grant rights to redistribute third-party data.
