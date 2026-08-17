#!/usr/bin/env Rscript

# GO 通路分析的固定 Entrez 表达矩阵
# -----------------------------------------------------------------------------
# 中文：
# 本脚本只为 GO ranked-pathway 分析建立一个固定输入矩阵。
# 它从已准备的 HuEx feature 表达矩阵中，按 HuEx SQLite 注释规则保留：
#   - 仅能唯一映射到一个 Entrez ID 的 feature；
#   - SQLite 中不标记为 multiple 的 feature；
#   - 同一 Entrez ID 有多个候选时，保留全样本平均表达最高的 feature。
# 得到 14,894 个 Entrez ID，每行一个代表性 feature。这里不使用 Braak、P 值或 FDR，
# 因而该输入矩阵在通路检验前独立确定。
#
# English:
# This script builds the fixed Entrez-level expression input for pathway analysis.
# It selects one HuEx feature per uniquely mapped Entrez ID using all-sample mean
# expression only; Braak status, P values, and FDR are never used at this step.
# -----------------------------------------------------------------------------

options(stringsAsFactors = FALSE)

stop_with <- function(message) stop(message, call. = FALSE)

argument_value <- function(flag) {
  args <- commandArgs(trailingOnly = TRUE)
  position <- match(flag, args)
  if (is.na(position) || position == length(args)) stop_with(paste("Missing", flag))
  args[[position + 1L]]
}

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) stop_with(paste("Missing R package:", package))
}

configured_path <- function(config, root, name, default) {
  value <- config$paths[[name]]
  if (is.null(value) || length(value) != 1L || is.na(value) || !nzchar(value)) value <- default
  if (grepl("^/", value) || grepl("^[A-Za-z]:[/\\\\]", value)) return(value)
  file.path(root, value)
}

write_csv_gz <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  connection <- gzfile(path, open = "wt")
  on.exit(close(connection), add = TRUE)
  utils::write.csv(data, connection, row.names = FALSE, na = "", quote = TRUE)
}

require_package("yaml")
require_package("DBI")
require_package("RSQLite")
require_package("huex10sttranscriptcluster.db")

config_file <- normalizePath(argument_value("--config"), mustWork = TRUE)
root <- dirname(dirname(config_file))
config <- yaml::read_yaml(config_file)
feature_file <- configured_path(config, root, "gse_feature_expression", "data/public/gse131617/expression_feature.csv.gz")
output_file <- configured_path(config, root, "gse_pathway_expression", "data/public/gse131617/pathway_expression.csv.gz")
expected_n <- config$analysis$pathway_expected_entrez_n
if (is.null(expected_n) || length(expected_n) != 1L || is.na(expected_n)) expected_n <- 14894L
expected_n <- as.integer(expected_n)
if (!file.exists(feature_file)) stop_with(paste("Missing feature expression file:", feature_file))

# Read the prepared feature matrix.
feature_table <- utils::read.csv(gzfile(feature_file, open = "rt"), check.names = FALSE)
if (names(feature_table)[1L] != "feature_id" || anyDuplicated(feature_table$feature_id)) {
  stop_with("Feature-expression input must begin with unique feature_id values.")
}
feature_ids <- trimws(as.character(feature_table$feature_id))
expression <- as.matrix(feature_table[, -1L, drop = FALSE])
storage.mode(expression) <- "double"
rownames(expression) <- feature_ids
if (anyNA(expression) || any(!is.finite(expression))) stop_with("Feature-expression matrix contains missing or non-finite values.")

# Read the annotation database directly. The probes table preserves the manuscript rule
# that excludes mappings flagged as multiple, even where an Entrez identifier is present.
sqlite_file <- system.file("extdata", "huex10sttranscriptcluster.sqlite", package = "huex10sttranscriptcluster.db")
if (!nzchar(sqlite_file) || !file.exists(sqlite_file)) {
  candidates <- list.files(system.file(package = "huex10sttranscriptcluster.db"), pattern = "\\.sqlite$", recursive = TRUE, full.names = TRUE)
  if (length(candidates) != 1L) stop_with("Could not locate the HuEx SQLite annotation database.")
  sqlite_file <- candidates[[1L]]
}
connection <- DBI::dbConnect(RSQLite::SQLite(), sqlite_file)
on.exit(DBI::dbDisconnect(connection), add = TRUE)
fields <- DBI::dbListFields(connection, "probes")
find_field <- function(candidates) {
  index <- match(tolower(candidates), tolower(fields), nomatch = 0L)
  index <- index[index > 0L]
  if (length(index) == 0L) NA_character_ else fields[[index[[1L]]]]
}
probe_field <- find_field(c("probe_id", "probeset_id"))
entrez_field <- find_field(c("gene_id", "entrez_id"))
multiple_field <- find_field(c("is_multiple", "multiple"))
if (is.na(probe_field) || is.na(entrez_field)) stop_with("The HuEx database does not contain the required mapping fields.")
quote_field <- function(field) as.character(DBI::dbQuoteIdentifier(connection, field))
multiple_sql <- if (is.na(multiple_field)) "0 AS is_multiple" else paste0(quote_field(multiple_field), " AS is_multiple")
mapping <- DBI::dbGetQuery(connection, paste0(
  "SELECT ", quote_field(probe_field), " AS feature_id, ", quote_field(entrez_field),
  " AS entrez_id, ", multiple_sql, " FROM probes"
))
mapping$feature_id <- trimws(as.character(mapping$feature_id))
mapping$entrez_id <- trimws(as.character(mapping$entrez_id))
mapping$entrez_id[mapping$entrez_id %in% c("", "NA", "N/A", "NULL")] <- NA_character_
mapping <- mapping[mapping$feature_id %in% feature_ids, , drop = FALSE]

# A feature is eligible only if it has exactly one numeric Entrez ID and no multiple flag.
collapse_mapping <- function(feature_id) {
  rows <- mapping[mapping$feature_id == feature_id, , drop = FALSE]
  ids <- sort(unique(rows$entrez_id[!is.na(rows$entrez_id) & grepl("^[0-9]+$", rows$entrez_id)]))
  is_multiple <- nrow(rows) > 1L || any(suppressWarnings(as.integer(rows$is_multiple)) != 0L, na.rm = TRUE)
  if (length(ids) != 1L || is_multiple) return(NULL)
  data.frame(feature_id = feature_id, entrez_id = ids[[1L]], stringsAsFactors = FALSE)
}
unique_mapping <- do.call(rbind, lapply(feature_ids, collapse_mapping))
if (is.null(unique_mapping) || nrow(unique_mapping) == 0L) stop_with("No uniquely mapped Entrez features were found.")

# Deterministic representative rule. Tie-breaker: feature ID.
unique_mapping$mean_expression <- rowMeans(expression[unique_mapping$feature_id, , drop = FALSE])
unique_mapping <- unique_mapping[order(unique_mapping$entrez_id, -unique_mapping$mean_expression, unique_mapping$feature_id), , drop = FALSE]
representative <- unique_mapping[!duplicated(unique_mapping$entrez_id), , drop = FALSE]
if (nrow(representative) != expected_n) {
  stop_with(paste("Expected", expected_n, "Entrez IDs but obtained", nrow(representative), ". Check the locked HuEx annotation package."))
}

pathway_expression <- data.frame(
  entrez_id = representative$entrez_id,
  expression[representative$feature_id, , drop = FALSE],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write_csv_gz(pathway_expression, output_file)
message("Prepared pathway expression matrix: ", normalizePath(output_file))
message("Unique Entrez IDs: ", nrow(pathway_expression))
