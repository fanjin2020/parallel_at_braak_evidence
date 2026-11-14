# PET-defined amyloid/tau states and late-Braak transcriptomic context

V5 contains the source-data analysis and publication-export code for this study. It produces one numbered figure set and its supporting data tables. A fresh run does not require an earlier code release, an author's server directory or previously fitted models. The code does not generate or edit Word manuscripts.

## Reproduce from source data

Use Linux with Python 3.12. Full analysis also requires R 4.4.2 and the packages listed under [Software environment](#software-environment). Download or clone this repository, then run the commands below from its root directory:

~~~bash
git clone https://github.com/fanjin2020/parallel_at_braak_evidence.git
cd parallel_at_braak_evidence
~~~

### 1. Download the public inputs

~~~bash
bash download.sh \
  --cache-dir "$PWD/cache/public" \
  --public-cache "$PWD/cache/public/GSE174367"
~~~

`download.sh` calls `scripts/download_public.py`. Its default is `--only all`, which retrieves both groups of public inputs used by the workflow:

| Input | Files retrieved |
| --- | --- |
| GSE131617 | Series matrix, official file list and the CEL file for each sample in the discovery matrix |
| GSE44770 | Series matrix and GPL4372 platform annotation table |
| GSE106241 | Expression matrix and GEO SOFT metadata |
| GSE174367 | Single-nucleus H5 count matrix and cell metadata |
| Published supporting data | Donor supplement and the Sutton et al. cell-reference workbook |
| Gene and pathway annotations | Human gene information, human GO associations (GAF) and the matching GO ontology (OBO) |

These are the public research inputs required by this workflow, not every file deposited under each GEO accession. The downloader does not retrieve ADNI data, install software packages or run analyses. The CEL and single-nucleus files are large; allow sufficient storage for both downloaded inputs and analysis working files.

Existing complete files are reused and partial downloads can resume. Files are checked against recorded hashes when available; an existing snapshot with a conflicting identity is not overwritten. Each cache receives a `public_sources_manifest.json` recording its sources and file identities. Keep the manifests with the data.

The main and single-nucleus caches must be distinct directories; the layout above uses a subdirectory for the latter. To download only one group, use `--only main` or `--only single-nucleus`.

### 2. Obtain authorized ADNI inputs

Apply through the [ADNI data-access page](https://adni.loni.usc.edu/data-samples/adni-data/) and obtain the required tables under the ADNI data-use agreement. ADNI participant-level data and restricted derivatives are not included in this repository or downloaded by its scripts.

For the example run below, place your authorized tables in these local directories:

| Directory | Required filenames |
| --- | --- |
| `private_inputs/adni/primary_raw/` | `ADSL.csv`, `APOERES.csv`, `CDR.csv`, `DATA_DOWNLOADED_DATE.csv`, `DXSUM.csv`, `PTDEMOG.csv`, `UCBERKELEY_AMY_6MM.csv`, `UCBERKELEY_TAU_6MM.csv`, `UCBERKELEY_TAUPVC_6MM.csv`, `UCBERKELEYFDG_8mm.csv`, `UCSFFSX7.csv` |
| `private_inputs/adni/confounder_raw/` | `GDSCALE.csv`, `BACKMEDS.csv`, `MEDHIST.csv`, `MODHACH.csv`, `RECCMEDS.csv`, `INITHEALTH.csv` |
| `private_inputs/adni/event_tables/` | `STUDYSUM.csv`, `REGISTRY.csv`, `ADVERSE.csv`, `RECADV.csv`, `TREATDIS.csv`, `DATADIC.csv`, `DATA_DOWNLOADED_DATE.csv`; also `OTDISC.csv` when available |

Use the original source fields and consistent data snapshots. `DATA_DOWNLOADED_DATE.csv` records the source snapshot in its `data_downloaded_date` field; do not substitute the analysis execution date. Field names and transformations are documented in [the variable dictionary](docs/Variable_Dictionary.csv). Missing required tables stop the run; empty placeholder files are not substitutes.

Alternatively, replace `--confounder-raw-dir` in the command below with `--adnimerge2-dir /path/to/authorized/ADNIMERGE2`. This must be an already downloaded, unpacked package containing `DESCRIPTION` and the six corresponding `data/*.rda` tables. No package download or ADNI login is performed by the analysis.

`private_inputs/` is a local, ignored directory. Do not commit or upload it.

### 3. Run the analyses and final export

After installing the [software environment](#software-environment) and preparing both public and authorized ADNI inputs:

~~~bash
bash run.sh \
  --adni-raw-dir "$PWD/private_inputs/adni/primary_raw" \
  --confounder-raw-dir "$PWD/private_inputs/adni/confounder_raw" \
  --adni-tables "$PWD/private_inputs/adni/event_tables" \
  --cache-dir "$PWD/cache/public" \
  --public-cache "$PWD/cache/public/GSE174367"
~~~

The launcher prints the resolved input locations, fits the models and exports the publication files. `run.sh` is offline: it does not download data or install packages. Missing inputs or dependencies stop the run with an error. The environment name is not fixed; activate whichever Linux environment contains the required versions before running the command.

### Data versions and exact numerical reproduction

A fresh download can rebuild the analyses, but matching the manuscript's exact numbers requires the same ADNI release, public-data snapshots, annotations and software environment. Unversioned upstream gene and GO annotations can change. A new download of their latest release is not equivalent to the manuscript's original snapshot.

To use specific, already available annotation snapshots during public-data preparation, supply `--gene-info`, `--go-gaf` and `--go-ontology` to `download.sh`. The OBO release must match the release declared in the GAF. An unavailable historical snapshot is not silently replaced by another release. Run logs record input identities and installed software versions.

## Software environment

The analysis entry point requires Linux, Python 3.12, R 4.4.2, Bioconductor 3.20, `statsmodels` 0.14.6 and `oligo` 1.70.1. Bioconductor 3.20 belongs to the R 4.4 release series; see the [official version table](https://www.bioconductor.org/about/release-announcements/). The conda environment used for the manuscript was named `paper2-r442`, but that name and the authors' filesystem layout are not required.

In your Python 3.12 environment, install the Python packages before starting the offline analysis:

~~~bash
python3 -m pip install -r requirements.txt
~~~

In R 4.4.2, install the required CRAN and Bioconductor packages and their dependencies:

~~~r
stopifnot(getRversion() == "4.4.2")
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 3600)
install.packages(c("BiocManager", "yaml", "readxl", "jsonlite", "statmod", "dtangle"))
BiocManager::install(version = "3.20", ask = FALSE, update = FALSE)
BiocManager::install(
  c("BiocVersion", "limma", "edgeR", "Biobase", "AnnotationDbi", "org.Hs.eg.db",
    "huex10sttranscriptcluster.db", "oligo", "pd.huex.1.0.st.v2", "preprocessCore"),
  version = "3.20", ask = FALSE, update = FALSE
)
stopifnot(as.character(packageVersion("oligo")) == "1.70.1")
~~~

These installation commands require network access and may require system compilers and development libraries on Linux. They are separate from data downloading and analysis. The workflow checks required versions and packages before fitting and records the versions actually used; `requirements.txt` contains ranges for several dependencies rather than a complete environment lockfile.

## Reuse data already downloaded

Skip the download step and run the same analysis command with your existing input directories. Keep the public-cache manifests and matching source files together. No previous model results are needed.

For an existing source-data project with a supported directory layout or recorded input locations, the optional convenience argument is:

~~~bash
bash run.sh --data-root /path/to/existing_source_data_project
~~~

With `--data-root`, the launcher searches only that project and input locations recorded by its runs. It does not search neighboring projects or reuse their fitted estimates. Different candidate snapshots require an explicit choice. The individual path arguments in step 3 work independently of any previous project layout.

## Outputs

The visible deliverable is `runs/v5_<timestamp>/submission/`.

| Output | Content |
| --- | --- |
| `Figures/` | Figures 1–4 and S1–S24, one set in PDF, PNG and SVG (28 figures, 84 files) |
| `Supplementary_Data_1.zip` | Grouped aggregate table and figure sources |
| `Supplementary_Data_2.csv` | All QC-retained baseline microarray feature results, including nonsignificant results |
| `Supplementary_Data_3.zip` | Twelve survivor, single-nucleus and observation-process result tables |
| `Variable_Dictionary.csv` | Used source fields, model variables and exported fields |
| `SOURCE_INDEX.csv` | Source identifiers, fields, row counts and file hashes |
| `TABLE_SOURCE_MAP.csv`, `FIGURE_SOURCE_MAP.csv` | Published table/figure-to-source links |
| `PUBLICATION_SUMMARY.json` | Actual resampling counts, dictionary coverage and output scope |
| `SUBMISSION_SHA256SUMS.csv`, `READ_FIRST.txt` | Output-file identity and reading instructions |

The three image formats contain the same figures. Source records retain full precision. `_record_id` and `_source_row` identify each source inside the grouped CSV files. Old manuscript values are not used as model inputs or export gates.

`.work/` stores private fit checkpoints; `audit/` stores execution logs, input identity and completion records. These are local working files, not manuscript attachments or GitHub contents. They are kept for resuming or re-exporting. Use `--discard-work` only if no later re-export will be needed.

## Resume or regenerate outputs

Resume with the same inputs, code and output directory:

~~~bash
bash run.sh --data-root /path/to/existing_data_project \
  --output-dir /path/to/interrupted_run --resume
~~~

Export from an intact completed fitted workspace:

~~~bash
bash run.sh --export-only --from-run /path/to/completed_run
~~~

Export-only does not fit models or download data. It requires the original `.work`, `audit` and `COMPLETED.json`, including any linked workspaces. A submission ZIP alone cannot supply fitted model inputs. Export-only cannot repair outdated model estimates; run the full workflow after a statistical-code change.

This writes a new `runs/v5_export_YYYYMMDD_HHMMSS/submission` folder and prints its location. It redraws the same 28 numbered figures in PDF, SVG and PNG, regenerates the submission data and dictionary, and updates the output checksum manifest. Earlier fitted results and submission folders are not overwritten. A/T minus signs are formatted only in plot labels; CSV category keys and numeric results retain their original values.

## Repository contents

- `run.py`, `run.sh`: offline analysis and final export.
- `download.sh`, `scripts/download_public.py`: separate public-data retrieval.
- `scripts/`, `workflows/`: source preparation, statistical models, figures and export.
- `config/`: model settings, variable definitions and publication-source maps.
- `docs/Variable_Dictionary.csv`: documentation dictionary; its provenance is recorded alongside it.
- `requirements.txt`, `LICENSE`, `.gitignore`, `.gitattributes`.

The dictionary documents 2,265 field instances in 180 scopes: 449 source fields, 47 model variables and 1,769 exported fields. It contains 1,261 documented entries and 1,004 entries whose source semantics remain explicitly unverified. It contains no participant observations. Each run generates its own dictionary from the fields actually used.

Independent regression comparisons, synthetic model tests and manuscript comparisons are maintained outside this repository in the local `checks` folder. Neither model fitting nor export imports or runs that folder. The publication renderer builds only the retained panels, without creating and discarding older figure sets. Input-format guards, model estimability checks, scientific QC and restrictions on participant-level export remain part of the analysis.

Upload only this repository directory. Do not upload caches, local input bundles, runs, `.work`, `audit` or the separate `checks` folder. Ignore rules do not remove files already tracked on GitHub.

The MIT license covers the software, not third-party data.
