# PET-defined amyloid/tau states and late-Braak transcriptomic context

V5 contains the analysis and publication-export code for this study. It fits models from local source data and produces one numbered figure set with its source tables. It does not generate or edit Word manuscripts.

## Run with data already available

Use Linux, Python 3.12, R 4.4.2 and the existing `paper2-r442` environment. Python dependencies are listed in `requirements.txt`; R analyses use Bioconductor 3.20, including limma, edgeR, oligo 1.70.1 and their dependencies, plus statmod, dtangle, readxl, jsonlite and yaml. Actual versions are recorded with each run.

~~~bash
conda activate paper2-r442
cd /home/research/git/repository/parallel_at_braak_evidence_V5
bash run.sh --data-root /home/research/git/repository/parallel_at_braak_evidence_V4
~~~

The launcher reads local raw ADNI tables, public-data caches and saved input locations from the selected project. It prints the resolved locations before fitting. With `--data-root`, it searches only that project and the input locations recorded by its previous runs. Without `--data-root`, it uses this repository's local input locations; it does not search neighboring projects. Multiple candidates with different recorded snapshots require an explicit path. It does not reuse earlier model estimates, move source files, download data or install packages. Missing inputs stop the run with the corresponding argument name.

For another layout, supply the directories explicitly:

~~~bash
bash run.sh \
  --adni-raw-dir /data/adni/primary_raw \
  --adnimerge2-dir /data/adni/ADNIMERGE2_local \
  --adni-tables /data/adni/event_tables \
  --cache-dir /data/public/main_cache \
  --public-cache /data/public/GSE174367
~~~

| Argument | Local input |
| --- | --- |
| `--adni-raw-dir` | ADSL, APOERES, CDR, DATA_DOWNLOADED_DATE, DXSUM, PTDEMOG, UCBERKELEY_AMY_6MM, UCBERKELEY_TAU_6MM, UCBERKELEY_TAUPVC_6MM, UCBERKELEYFDG_8mm and UCSFFSX7 CSV files |
| `--adnimerge2-dir` | Downloaded, unpacked ADNIMERGE2 package, with DESCRIPTION and local R data files |
| `--confounder-raw-dir` | Alternative to the R package: GDSCALE, BACKMEDS, MEDHIST, MODHACH, RECCMEDS and INITHEALTH CSV files |
| `--adni-tables` | STUDYSUM, REGISTRY, ADVERSE, RECADV, TREATDIS, DATADIC and DATA_DOWNLOADED_DATE; OTDISC when available |
| `--private-bundle` | Alternative to event-table directory: the existing private ADNI event archive |
| `--cache-dir` | Public cache with public_sources_manifest.json, source matrices, CEL files, matching annotation snapshots and reference workbook |
| `--public-cache` | GSE174367 H5 raw counts and cell metadata |

ADNI access remains subject to its data-use agreement. There is no ADNI downloader. Source data and restricted derivatives must remain outside the public repository.

## Download public inputs separately

Skip this section when the data are already available.

~~~bash
bash download.sh --cache-dir /data/public/main_cache \
                 --public-cache /data/public/GSE174367
~~~

This separate command retrieves GSE131617 expression and CEL files, GSE44770, GSE106241, GSE174367, the published donor supplement, cell-reference workbook, human gene annotation and GO annotation. It reuses existing complete files and resumes partial downloads. It never overwrites an existing snapshot with different bytes.

The main cache and single-nucleus cache must be different directories. The default places the single-nucleus files in a subdirectory. A shared directory, including a symbolic-link alias to the same directory, is rejected before downloading; an existing manifest of the other cache type is also protected.

Use `--only main` or `--only single-nucleus` to fetch one group. Optional `--go-gaf`, `--go-ontology` and `--gene-info` select already downloaded public snapshots. The OBO release must match the release declared by the GAF. If an older ontology projection is unavailable upstream, retain the original local snapshot; a different release is not an equivalent replacement.

Unversioned upstream annotations can change. Exact manuscript reproduction requires the recorded cache snapshot, not a new download of the latest annotation.

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

Export from an intact completed V4 or V5 fitted workspace:

~~~bash
bash run.sh --export-only --from-run /path/to/completed_run
~~~

Export-only does not fit models or download data. It requires the original `.work`, `audit` and `COMPLETED.json`, including any linked workspaces. A submission ZIP alone cannot supply fitted model inputs. Export-only cannot repair outdated model estimates; run the full workflow after a statistical-code change.

After a figure-label or export-only code change, reuse the latest completed model run in this project:

~~~bash
conda activate paper2-r442
cd /home/research/git/repository/parallel_at_braak_evidence_V5
bash run.sh --export-only --from-run "$PWD"
~~~

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
