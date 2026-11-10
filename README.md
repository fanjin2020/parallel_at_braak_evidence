# PET A/T and late-Braak transcriptomic analyses — V4

Analysis code for *PET-defined amyloid and tau states and late-Braak transcriptomic context in Alzheimer’s disease*. The clinical and tissue analyses use separate populations.

The pipeline runs on Linux using local authorized ADNI inputs and existing public-data caches. It does not download data or install packages. This repository contains code, configuration, field metadata, and regression-test fingerprints; it contains no participant observations or packaged datasets.

## Repository contents

| Path | Purpose |
| --- | --- |
| `run.sh`, `run.py` | Offline entry point and workflow orchestration |
| `scripts/` | Input preparation, clinical and tissue models, validation, figures, and export |
| `workflows/tissue_survival/` | Event-record linkage, survivor outcome sensitivity, pseudobulk, and cell-type interaction analyses |
| `config/` | Analysis settings, source maps, field definitions, and reference fingerprints |
| `docs/Variable_Dictionary.csv` | Dictionary supplied with the checked manuscript export |
| `tests/` | Synthetic workflow and publication checks; no clinical records |
| `requirements.txt` | Python dependency ranges |
| `CODE_SHA256SUMS.csv` | File-integrity manifest checked before execution |
| `LICENSE` | Software license |

Upload the contents of this repository directory, including the dotfiles, rather than a parent directory containing local inputs. Keep the tests and integrity manifest: they check rounding, source coverage, bootstrap-source consistency, and export-only behavior.

## Variable dictionary

[Variable_Dictionary.csv](docs/Variable_Dictionary.csv) describes the input fields, model variables, and exported result fields used in the manuscript analyses. It is identical to the dictionary supplied with the checked export `v4_export_20261110_103922`.

The snapshot contains 2,265 field instances across 180 scopes: 449 input fields, 47 model variables, and 1,769 exported-result fields. A scope identifies an input source, model-variable group, or logical output; it is not necessarily an original data table. The row key is `(scope, field)`.

| Column | Meaning |
| --- | --- |
| `scope` | Source, model-variable group, or logical result |
| `field` | Exact field or variable name |
| `role` | `source_input`, `model_variable`, or `exported_result` |
| `definition` | Available definition or documented analytical use |
| `unit` | Recorded unit, where available |
| `coding_or_missing` | Available coding and missing-value information |
| `definition_status` | Documented definition or unverified source-field semantics |
| `code` | Script reference associated with the field |

There are 1,230 `documented` entries and 1,035 `source_definition_unverified` entries. A documented entry has a recorded definition; this label does not certify every original provider definition. Blank unit or coding fields do not establish that a variable is unitless or has no missing values.

[Variable_Dictionary.provenance.json](docs/Variable_Dictionary.provenance.json) records the source export, date, counts, and SHA256. This documentation snapshot is not a model input. Each completed export generates `submission/Variable_Dictionary.csv`; update the documentation snapshot and provenance together after checking a newer export.

The dictionary contains field metadata only, with no participant observations or observed-value examples. It does not grant access to ADNI data. The separate metadata lookup in `config/adni_source_dictionary_snapshot.csv` is used by the dictionary generator and must also be retained.

中文说明：字典仅列本文使用的输入字段、模型变量及导出结果字段，不列整个R包中的全部表；未核实的源语义保留标记。投稿附件与本地上传包使用同一份字典，不含ADNI个体记录。

## Environment and inputs

Use the analysis environment on Linux: Python 3.12, R 4.4.2, Bioconductor 3.20, and statsmodels 0.14.6. The analysis uses limma, edgeR, statmod, oligo, dtangle, and the additional packages checked by the workflow. Python ranges are listed in `requirements.txt`; actual versions are recorded during execution. The launcher limits numerical-library threads to one.

Provide local inputs explicitly when using a different server layout:

| Argument | Expected local content |
| --- | --- |
| `--adni-raw-dir` | ADSL, APOERES, CDR, DATA_DOWNLOADED_DATE, DXSUM, PTDEMOG, UCBERKELEY_AMY_6MM, UCBERKELEY_TAU_6MM, UCBERKELEY_TAUPVC_6MM, UCBERKELEYFDG_8mm, and UCSFFSX7 CSV files |
| `--adnimerge2-dir` | Already downloaded, unpacked ADNIMERGE2 package with DESCRIPTION and local R data files |
| `--confounder-raw-dir` | Alternative to the local R package: GDSCALE, BACKMEDS, MEDHIST, MODHACH, RECCMEDS, and INITHEALTH raw CSV files |
| `--cache-dir` | Main public cache with `public_sources_manifest.json`, source matrices, CEL files, annotations, GO snapshots, cell reference, and published donor supplement |
| `--public-cache` | GSE174367 H5 counts and cell metadata cache; distinct from the main public cache |
| `--adni-tables` | STUDYSUM, REGISTRY, ADVERSE, RECADV, TREATDIS, DATADIC, and DATA_DOWNLOADED_DATE CSV files; OTDISC is used when provided |
| `--private-bundle` | Alternative to `--adni-tables`: the existing local ADNI event-table archive |

Caches are verified against their recorded sizes and checksums. Missing or changed inputs stop the run. Public source identifiers are recorded in `scripts/public_sources.py` and the input cache manifest. Restricted ADNI data must be obtained under the user's own authorization and must stay outside the public repository.

## Run from local source data

Replace the example paths with directories on the server. Check the environment and inputs first:

```bash
conda activate paper2-r442
cd /path/to/parallel_at_braak_evidence_V4

bash run.sh --check \
  --adni-raw-dir /data/adni/primary_raw \
  --adnimerge2-dir /data/adni/ADNIMERGE2_local \
  --cache-dir /data/public/main_cache \
  --public-cache /data/public/GSE174367 \
  --adni-tables /data/adni/event_tables
```

Then run the same command without `--check`:

```bash
bash run.sh \
  --adni-raw-dir /data/adni/primary_raw \
  --adnimerge2-dir /data/adni/ADNIMERGE2_local \
  --cache-dir /data/public/main_cache \
  --public-cache /data/public/GSE174367 \
  --adni-tables /data/adni/event_tables
```

The default output is a new directory under `runs/`. Use `--output-dir` to choose another location. Explicit input paths avoid relying on the legacy server-directory defaults; existing server layouts remain supported.

For compatibility, `--prepared-inputs DIR` can reuse a prepared primary workspace, including its primary estimates, and refit the subsequent models. It is not a complete refit from raw source data, is disabled by default, and cannot be combined with primary raw-input options.

## Resume and export-only

To resume an interrupted run, use its original input arguments, the same `--output-dir`, and `--resume`. Inputs, code, environment, and completed-step files must pass integrity checks. A code change requires a new run directory; do not use `--resume` to combine different code versions. Individual bootstrap iterations are not all checkpointed.

To regenerate summaries, figures, and submission files from a completed V4 fit workspace:

```bash
bash run.sh --export-only --from-run /path/to/completed/v4_run
```

Without `--from-run`, export-only mode selects the latest completed fit run under this project's `runs/` and prints its path. It requires the saved `.work`, `audit`, and `COMPLETED.json`; the exported submission files alone are insufficient. It does not call R, refit models, download inputs, or change the old run. Outputs go to a new `runs/v4_export_*/` directory.

## Publication outputs

A completed export contains one figure set and the required source data:

| Output | Content |
| --- | --- |
| `Figures/` | Main Figures 1–4 and Supplementary Figures S1–S24, each as PDF, PNG, and SVG: 28 figures and 84 files |
| `Supplementary_Data_1.zip` | Grouped source tables for manuscript results and figures, including all-feature processing counts and Table 4 display values |
| `Supplementary_Data_2.csv` | Complete QC-retained baseline feature-level testing results |
| `Supplementary_Data_3.zip` | Twelve source tables covering survivor outcome sensitivity, single-nucleus analyses, OPC–INH interaction, and observation-process analyses |
| `Variable_Dictionary.csv` | Used input, model, and exported-result fields |
| `SOURCE_INDEX.csv` | Stable source keys, row counts, fields, and checksums |
| `TABLE_SOURCE_MAP.csv`, `FIGURE_SOURCE_MAP.csv` | Manuscript table/figure-to-source mapping |
| `PUBLICATION_SUMMARY.json` | Actual resampling counts, dictionary coverage, and validation status |
| `SUBMISSION_SHA256SUMS.csv`, `READ_FIRST.txt` | Export checksums and reading instructions |

Stable record keys are not display figure numbers. Data 1 retains `_record_id` and `_source_row`; the source index records each table's columns. The three image formats represent the same analysis, not separate figure sets.

All-feature counts use the complete fitted background rather than selected targets. Display rounding leaves full-precision estimates unchanged. Figure S8 is checked against the same threshold-bootstrap source as Table S5D.

Reference mismatches are reported in `results_to_review/`; the program does not alter estimates to match the manuscript. `COMPLETED.json` records workflow completion. A matching regression fingerprint does not certify Word contents: the program neither creates nor updates the manuscript documents.

Private model checkpoints are kept in `.work/`. They are removed after a successful matching export unless `--keep-work` is set; failed or mismatching runs retain them for inspection. Local `audit/` records include logs, input fingerprints, environment details, and step states. Do not upload entire run directories.

## Validation and access boundaries

```bash
python -m unittest discover -s tests -v
```

The tests use synthetic inputs. Passing them does not replace a full Linux/R analysis run. `config/result_reference.csv` contains regression-test fingerprints, not fitted-data inputs; its provenance is recorded separately. Code checksums must be updated when repository files change.

Git ignore rules exclude local inputs, run outputs, caches, archives, logs, credentials, and editor files. Only the specific configuration CSVs, code manifest, and public variable dictionary are allowed through the CSV rule. Ignore rules do not remove files already tracked in an existing remote repository.

The MIT license applies to repository software. It does not grant redistribution rights to ADNI or other third-party datasets; their access and reuse remain subject to the respective providers' conditions.
