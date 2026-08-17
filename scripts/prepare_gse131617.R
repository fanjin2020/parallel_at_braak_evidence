#!/usr/bin/env Rscript

# =============================================================================
# GSE131617 原始数据准备 / Prepare canonical GSE131617 analysis inputs
# =============================================================================
#
# 中文：本脚本只做一件事：把 GEO 原始 series matrix 与作者补充表转换为后续
# 差异表达分析唯一使用的四个标准输入文件。它不做差异表达、富集、图形、审计。
#
# English: This script converts the raw GEO series matrix and the accompanying
# donor spreadsheet into the four canonical inputs required for the downstream
# differential-expression analysis. It does not run models, enrichment, plots,
# audits, or create any historical/legacy output.
#
# Annotation is obtained from huex10sttranscriptcluster.db. GPL5175.annot.gz is
# deliberately NOT required because that GEO download endpoint is unreliable.
# =============================================================================

options(stringsAsFactors = FALSE)

stop_with <- function(message) stop(message, call. = FALSE)

argument_value <- function(name) {
  arguments <- commandArgs(trailingOnly = TRUE)
  position <- match(name, arguments)
  if (is.na(position) || position == length(arguments)) stop_with(paste("Missing argument", name))
  arguments[[position + 1L]]
}

need_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) stop_with(paste("Missing R package:", package))
}

project_path <- function(root, value) {
  if (is.null(value) || length(value) != 1L || is.na(value) || !nzchar(value)) return(NULL)
  if (grepl("^/", value) || grepl("^[A-Za-z]:[/\\\\]", value)) return(value)
  file.path(root, value)
}

configured_path <- function(config, root, name, default) {
  value <- project_path(root, config$paths[[name]])
  if (is.null(value)) file.path(root, default) else value
}

write_csv_gz <- function(data, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  connection <- gzfile(file, open = "wt")
  on.exit(close(connection), add = TRUE)
  utils::write.csv(data, connection, row.names = FALSE, na = "")
}

unquote_geo <- function(value) gsub('^"|"$', "", value)

parse_titles <- function(matrix_file) {
  connection <- gzfile(matrix_file, open = "rt")
  on.exit(close(connection), add = TRUE)
  title <- readLines(connection, warn = FALSE)
  title <- title[grepl("^!Sample_title", title)]
  if (length(title) != 1L) stop_with("Unable to find a unique !Sample_title line in the GEO matrix.")
  values <- unquote_geo(strsplit(title, "\\t")[[1L]][-1L])
  matches <- regexec("Braak_NFT_stage_(0|I-II|III-IV|V-VI)_subject-([0-9]+)_(TC|FC|EC)", values)
  pieces <- regmatches(values, matches)
  if (any(lengths(pieces) != 4L)) stop_with("Unexpected GEO sample-title format.")
  data.frame(
    braak_stage = vapply(pieces, `[[`, character(1), 2L),
    donor_number = vapply(pieces, `[[`, character(1), 3L),
    brain_region = vapply(pieces, `[[`, character(1), 4L),
    stringsAsFactors = FALSE
  )
}

donor_key <- function(stage, number) paste0("BN_", stage, "_", sprintf("%02d", as.integer(number)))

build_annotation <- function(feature_ids, database) {
  mapping <- suppressMessages(AnnotationDbi::select(
    database, keys = feature_ids, columns = "ENTREZID", keytype = "PROBEID"
  ))
  if (!all(c("PROBEID", "ENTREZID") %in% names(mapping))) stop_with("HuEx annotation package did not provide PROBEID and ENTREZID.")
  mapping$PROBEID <- trimws(as.character(mapping$PROBEID))
  mapping$ENTREZID <- trimws(as.character(mapping$ENTREZID))
  mapping <- mapping[mapping$PROBEID %in% feature_ids & grepl("^[0-9]+$", mapping$ENTREZID), c("PROBEID", "ENTREZID"), drop = FALSE]

  rows <- lapply(feature_ids, function(id) {
    entrez <- sort(unique(mapping$ENTREZID[mapping$PROBEID == id]))
    if (length(entrez) == 0L) {
      data.frame(feature_id = id, entrez_id = "", mapping_status = "unmapped_entrez")
    } else if (length(entrez) == 1L) {
      data.frame(feature_id = id, entrez_id = entrez, mapping_status = "unique_entrez_id")
    } else {
      data.frame(feature_id = id, entrez_id = paste(entrez, collapse = ";"), mapping_status = "ambiguous_multiple_entrez")
    }
  })
  annotation <- do.call(rbind, rows)
  if (nrow(annotation) != length(feature_ids) || anyDuplicated(annotation$feature_id)) stop_with("Annotation must have exactly one row per feature.")
  annotation
}

# ---- Read configuration and validate prerequisites --------------------------------
config_file <- normalizePath(argument_value("--config"), mustWork = TRUE)
root <- dirname(dirname(config_file))
need_package("yaml")
config <- yaml::read_yaml(config_file)

for (package in c("readxl", "AnnotationDbi", "huex10sttranscriptcluster.db")) need_package(package)

matrix_file <- configured_path(config, root, "gse_series_matrix", "data/public/gse131617/raw/GSE131617-GPL5175_series_matrix.txt.gz")
subject_file <- configured_path(config, root, "gse_subject_info", "data/public/gse131617/raw/GSE131617_Subject_info_Miyashita_TranslPsychiatry_2014_Supple_TS8.xlsx")
feature_expression_file <- configured_path(config, root, "gse_feature_expression", "data/public/gse131617/expression_feature.csv.gz")
annotation_file <- configured_path(config, root, "gse_feature_annotation", "data/public/gse131617/feature_annotation.csv")
entrez_expression_file <- configured_path(config, root, "gse_expression", "data/public/gse131617/expression_entrez.csv.gz")
manifest_file <- configured_path(config, root, "gse_manifest", "data/public/gse131617/sample_manifest.csv")

for (file in c(matrix_file, subject_file)) if (!file.exists(file)) stop_with(paste("Missing required GSE input:", file))

# ---- Read the already processed GEO expression values (no re-normalisation) -------
connection <- gzfile(matrix_file, open = "rt")
on.exit(close(connection), add = TRUE)
series <- utils::read.delim(connection, header = TRUE, sep = "\t", quote = "\"", comment.char = "!", check.names = FALSE)
if (ncol(series) != 214L || toupper(names(series)[[1L]]) != "ID_REF") stop_with("Expected ID_REF plus 213 GSM columns in the GEO matrix.")

titles <- parse_titles(matrix_file)
if (nrow(titles) != 213L) stop_with("Expected 213 GEO sample titles.")
titles$sample_id <- names(series)[-1L]

feature_ids <- trimws(as.character(series[[1L]]))
if (any(!nzchar(feature_ids)) || anyDuplicated(feature_ids)) stop_with("Feature IDs must be nonempty and unique.")
values <- as.matrix(series[, titles$sample_id, drop = FALSE])
storage.mode(values) <- "double"
if (anyNA(values) || any(!is.finite(values))) stop_with("Expression matrix contains missing or non-finite values.")
rownames(values) <- feature_ids

# Output 1: full feature-level matrix for the primary differential-expression model.
write_csv_gz(data.frame(feature_id = feature_ids, values, check.names = FALSE), feature_expression_file)

# Output 2: one mapping-status record for every feature.
chip_database <- get("huex10sttranscriptcluster.db", envir = asNamespace("huex10sttranscriptcluster.db"))
annotation <- build_annotation(feature_ids, chip_database)
dir.create(dirname(annotation_file), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(annotation, annotation_file, row.names = FALSE, na = "")

# Output 3: one representative feature per uniquely mapped Entrez identifier.
# Selection uses overall mean expression only; it never uses Braak effect size, P, or FDR.
unique_map <- annotation[annotation$mapping_status == "unique_entrez_id", c("feature_id", "entrez_id")]
unique_map$mean_expression <- rowMeans(values[unique_map$feature_id, , drop = FALSE])
unique_map <- unique_map[order(unique_map$entrez_id, -unique_map$mean_expression, unique_map$feature_id), ]
representative <- unique_map[!duplicated(unique_map$entrez_id), ]
if (anyDuplicated(representative$entrez_id)) stop_with("Entrez representative selection failed.")
write_csv_gz(data.frame(entrez_id = representative$entrez_id, values[representative$feature_id, , drop = FALSE], check.names = FALSE), entrez_expression_file)

# ---- Build donor-level manifest from Supplementary Table S8 ------------------------
header <- readxl::read_excel(subject_file, sheet = 1, range = "D4:K5", col_names = FALSE, .name_repair = "minimal")
required_header_terms <- c("Subject ID", "Gender", "NFT", "SP", "APOE", "AAD", "PMI")
if (any(!vapply(required_header_terms, grepl, logical(1), x = paste(unlist(header), collapse = " | "), fixed = TRUE))) stop_with("Unexpected S8 donor-table header in D4:K5.")
s8 <- readxl::read_excel(subject_file, sheet = 1, range = "D6:K100", col_names = FALSE, .name_repair = "minimal")
if (ncol(s8) != 8L) stop_with("Expected eight columns in S8 range D6:K100.")

donors <- data.frame(
  donor_key = trimws(as.character(s8[[1L]])),
  sex = trimws(as.character(s8[[3L]])),
  braak_s8 = trimws(as.character(s8[[4L]])),
  apoe = trimws(as.character(s8[[6L]])),
  age_at_death = suppressWarnings(as.numeric(s8[[7L]])),
  pmi_hours = suppressWarnings(as.numeric(s8[[8L]]))
)
donors <- donors[grepl("^BN_", donors$donor_key), ]
if (nrow(donors) != 71L || anyDuplicated(donors$donor_key)) stop_with("Expected exactly 71 unique BN_ donors in S8.")
donors$sex_male <- ifelse(donors$sex == "M", 1, ifelse(donors$sex == "F", 0, NA_real_))
donors$apoe4_carrier <- ifelse(grepl("4", donors$apoe, fixed = TRUE), 1, 0)

titles$donor_key <- donor_key(titles$braak_stage, titles$donor_number)
manifest <- merge(titles, donors[, c("donor_key", "braak_s8", "age_at_death", "sex_male", "pmi_hours", "apoe4_carrier")], by = "donor_key", all.x = TRUE, sort = FALSE)
manifest <- manifest[match(titles$sample_id, manifest$sample_id), c("sample_id", "donor_key", "braak_stage", "brain_region", "age_at_death", "sex_male", "pmi_hours", "apoe4_carrier")]
if (nrow(manifest) != 213L || anyNA(manifest) || anyDuplicated(manifest$sample_id)) stop_with("Prepared manifest is incomplete.")
if (!all(titles$braak_stage == donors$braak_s8[match(titles$donor_key, donors$donor_key)])) stop_with("GEO Braak labels and S8 NFT labels disagree.")
if (length(unique(manifest$donor_key)) != 71L || any(table(manifest$donor_key) != 3L)) stop_with("Expected three cortical samples per donor.")
dir.create(dirname(manifest_file), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(manifest, manifest_file, row.names = FALSE, na = "")

message("Prepared GSE131617 inputs: 22,011 features; 213 samples from 71 donors.")
