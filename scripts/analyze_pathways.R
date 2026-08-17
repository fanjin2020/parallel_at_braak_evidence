#!/usr/bin/env Rscript

# GO biological-process ranked-pathway analysis
# -----------------------------------------------------------------------------
# 中文：
# 本脚本对固定的 14,894 Entrez 基因背景进行 GO biological-process CAMERA-PR 分析。
# 模型、协变量及留一脑区方向核验与主基因分析一致。一个通路必须同时满足：
# 全模型 FDR < 0.05，且留出 EC、FC、TC 后的上/下调方向均一致，才称为方向稳定。
# 输出完整通路结果和最终稳定通路表；不会改动 17 个稳健基因结果。
#
# English:
# This script performs CAMERA-PR on the fixed Entrez universe. A GO biological
# process is retained only when full-model FDR is below 0.05 and the direction
# agrees in all three leave-one-region-out analyses. It does not alter the
# feature-level robust-gene analysis.
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

# Fit the same donor-blocked model used for feature-level differential expression,
# then convert moderated t statistics to limma z scores for CAMERA-PR.
fit_z_scores <- function(expression, manifest) {
  design <- stats::model.matrix(
    ~ 0 + braak_stage + brain_region + age_at_death + sex_male + pmi_hours + apoe4_carrier,
    data = manifest
  )
  colnames(design) <- make.names(colnames(design), unique = TRUE)
  if (qr(design)$rank != ncol(design)) stop_with("Pathway design matrix is rank deficient.")
  donor_correlation <- limma::duplicateCorrelation(expression, design, block = manifest$donor_key)$consensus
  fit <- limma::lmFit(expression, design, block = manifest$donor_key, correlation = donor_correlation)
  contrast <- limma::makeContrasts(late_vs_zero = braak_stageV.VI - braak_stage0, levels = design)
  fit <- limma::eBayes(limma::contrasts.fit(fit, contrast), robust = TRUE)
  limma::zscoreT(fit$t[, "late_vs_zero"], df = fit$df.total)
}

build_go_sets <- function(entrez_ids, min_size, max_size) {
  organism_db <- get("org.Hs.eg.db", envir = asNamespace("org.Hs.eg.db"))
  annotation <- AnnotationDbi::select(
    organism_db, keys = entrez_ids, keytype = "ENTREZID", columns = c("GO", "ONTOLOGY")
  )
  annotation <- annotation[annotation$ONTOLOGY == "BP" & !is.na(annotation$GO), c("ENTREZID", "GO")]
  sets <- split(as.character(annotation$ENTREZID), annotation$GO)
  sets <- lapply(sets, function(ids) intersect(unique(ids), entrez_ids))
  sets[lengths(sets) >= min_size & lengths(sets) <= max_size]
}

run_camera <- function(z_scores, sets, inter_gene_correlation) {
  result <- limma::cameraPR(
    statistic = z_scores,
    index = sets,
    use.ranks = FALSE,
    inter.gene.cor = inter_gene_correlation,
    sort = TRUE
  )
  result$go_id <- rownames(result)
  result
}

require_package("yaml")
require_package("limma")
require_package("AnnotationDbi")
require_package("org.Hs.eg.db")
require_package("GO.db")

config_file <- normalizePath(argument_value("--config"), mustWork = TRUE)
root <- dirname(dirname(config_file))
config <- yaml::read_yaml(config_file)
expression_file <- configured_path(config, root, "gse_pathway_expression", "data/public/gse131617/pathway_expression.csv.gz")
manifest_file <- configured_path(config, root, "gse_manifest", "data/public/gse131617/sample_manifest.csv")
results_dir <- configured_path(config, root, "pathway_results_dir", "results/pathways")
alpha_value <- config$analysis$fdr_alpha
if (is.null(alpha_value) || length(alpha_value) != 1L || is.na(alpha_value)) alpha_value <- 0.05
alpha <- as.numeric(alpha_value)
min_size_value <- config$analysis$pathway_min_genes
max_size_value <- config$analysis$pathway_max_genes
correlation_value <- config$analysis$camera_inter_gene_correlation
if (is.null(min_size_value) || length(min_size_value) != 1L || is.na(min_size_value)) min_size_value <- 10L
if (is.null(max_size_value) || length(max_size_value) != 1L || is.na(max_size_value)) max_size_value <- 300L
if (is.null(correlation_value) || length(correlation_value) != 1L || is.na(correlation_value)) correlation_value <- 0.01
min_size <- as.integer(min_size_value)
max_size <- as.integer(max_size_value)
inter_gene_correlation <- as.numeric(correlation_value)
for (file in c(expression_file, manifest_file)) if (!file.exists(file)) stop_with(paste("Missing input:", file))
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

expression_table <- utils::read.csv(gzfile(expression_file), check.names = FALSE)
if (names(expression_table)[1L] != "entrez_id" || anyDuplicated(expression_table$entrez_id)) {
  stop_with("Pathway-expression input must begin with unique entrez_id values.")
}
if (nrow(expression_table) != 14894L) stop_with("The fixed pathway input must contain 14,894 Entrez IDs.")
expression <- as.matrix(expression_table[, -1L, drop = FALSE])
storage.mode(expression) <- "double"
rownames(expression) <- as.character(expression_table$entrez_id)
if (anyNA(expression) || any(!is.finite(expression))) stop_with("Pathway-expression matrix contains missing or non-finite values.")

manifest <- utils::read.csv(manifest_file, check.names = FALSE)
needed <- c("sample_id", "donor_key", "braak_stage", "brain_region", "age_at_death", "sex_male", "pmi_hours", "apoe4_carrier")
if (length(setdiff(needed, names(manifest))) > 0L) stop_with("Sample manifest lacks a required column.")
if (!setequal(colnames(expression), manifest$sample_id)) stop_with("Expression and manifest sample IDs do not match.")
manifest <- manifest[match(colnames(expression), manifest$sample_id), , drop = FALSE]
manifest$braak_stage <- factor(manifest$braak_stage, levels = c("0", "I-II", "III-IV", "V-VI"))
manifest$brain_region <- factor(manifest$brain_region, levels = c("EC", "FC", "TC"))
manifest$donor_key <- factor(manifest$donor_key)
for (column in c("age_at_death", "sex_male", "pmi_hours", "apoe4_carrier")) manifest[[column]] <- as.numeric(manifest[[column]])
if (nrow(manifest) != 213L || length(unique(manifest$donor_key)) != 71L ||
    anyNA(manifest[, c("age_at_death", "sex_male", "pmi_hours", "apoe4_carrier")])) {
  stop_with("Manifest does not satisfy the required 213-sample, 71-donor complete-covariate structure.")
}

sets <- build_go_sets(rownames(expression), min_size, max_size)
full <- run_camera(fit_z_scores(expression, manifest), sets, inter_gene_correlation)
names(full)[names(full) == "Direction"] <- "full_direction"
names(full)[names(full) == "PValue"] <- "full_p_value"
names(full)[names(full) == "FDR"] <- "full_fdr"

# Direction-only leave-one-region-out checks: they are not separate FDR screens.
for (region in levels(manifest$brain_region)) {
  keep <- manifest$brain_region != region
  leave_out <- run_camera(
    fit_z_scores(expression[, keep, drop = FALSE], droplevels(manifest[keep, , drop = FALSE])),
    sets, inter_gene_correlation
  )
  leave_out <- leave_out[, c("go_id", "Direction")]
  names(leave_out)[2L] <- paste0("direction_leave_out_", region)
  full <- merge(full, leave_out, by = "go_id", all.x = TRUE, sort = FALSE)
}
full$directionally_stable <- with(full,
  full_fdr < alpha & full_direction == direction_leave_out_EC &
    full_direction == direction_leave_out_FC & full_direction == direction_leave_out_TC
)
full <- full[order(full$full_fdr, full$go_id), , drop = FALSE]
stable <- full[full$directionally_stable, , drop = FALSE]

utils::write.csv(full, file.path(results_dir, "ranked_pathways.csv"), row.names = FALSE)
utils::write.csv(stable, file.path(results_dir, "stable_pathways.csv"), row.names = FALSE)
message("Ranked-pathway analysis completed: ", normalizePath(results_dir))
message("Directionally stable GO biological-process terms: ", nrow(stable))
