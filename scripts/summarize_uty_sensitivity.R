#!/usr/bin/env Rscript

# GSE131617 donor-sex and UTY sensitivity summary
# -----------------------------------------------------------------------------
# 中文：
# 本脚本不重新进行差异表达分析，也不产生性别分层差异表达结论。
# 它仅完成论文中必要的解释边界检查：
# 1. 在供体层面汇总各 Braak 分期的生物学性别；
# 2. 将最终 17 个稳健 feature 映射为基因符号；
# 3. 标记唯一的 Y 染色体基因 UTY；
# 4. 记录移除 UTY 后仍有 16 个非 UTY 稳健基因。
# 主 limma 模型的 sex_male 协变量已在 analyze_gse131617.R 中纳入。
#
# English:
# This script does not rerun differential expression or claim sex-stratified
# effects. It documents donor sex by Braak stage and identifies UTY as the
# single Y-chromosome gene among the 17 directionally robust genes.
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

require_package("yaml")

config_file <- normalizePath(argument_value("--config"), mustWork = TRUE)
root <- dirname(dirname(config_file))
config <- yaml::read_yaml(config_file)
manifest_file <- configured_path(config, root, "gse_manifest", "data/public/gse131617/sample_manifest.csv")
robust_gene_file <- configured_path(config, root, "gse_robust_genes", "results/gse131617/braak_robust_genes.csv")
symbol_file <- configured_path(config, root, "gse_gene_symbols", "data/contracts/gene_symbols.csv")
results_dir <- configured_path(config, root, "gse_sensitivity_dir", "results/gse_sensitivity")
for (file in c(manifest_file, robust_gene_file, symbol_file)) if (!file.exists(file)) stop_with(paste("Missing input:", file))
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Each donor has EC, FC, and TC samples; count each donor only once.
manifest <- utils::read.csv(manifest_file, check.names = FALSE)
needed_manifest <- c("donor_key", "braak_stage", "brain_region", "sex_male")
if (length(setdiff(needed_manifest, names(manifest))) > 0L) stop_with("Sample manifest lacks a required column.")
donors <- manifest[!duplicated(manifest$donor_key), c("donor_key", "braak_stage", "sex_male"), drop = FALSE]
if (nrow(donors) != 71L) stop_with("Expected 71 donors.")
donors$sex_label <- ifelse(donors$sex_male == 1, "male", ifelse(donors$sex_male == 0, "female", "missing_or_other"))
stage_order <- c("0", "I-II", "III-IV", "V-VI")
donors$braak_stage <- factor(donors$braak_stage, levels = stage_order)
sex_by_stage <- as.data.frame(with(donors, table(braak_stage, sex_label)), stringsAsFactors = FALSE)
names(sex_by_stage) <- c("braak_stage", "sex_label", "n_donors")
sex_by_stage$braak_stage <- as.character(sex_by_stage$braak_stage)

# Join the fixed gene-symbol dictionary to the reproduced 17-gene result.
robust <- utils::read.csv(robust_gene_file, check.names = FALSE, colClasses = "character")
symbols <- utils::read.csv(symbol_file, check.names = FALSE, colClasses = "character")
needed_robust <- c("feature_id", "entrez_id", "log2_fold_change", "fdr_bh", "primary_direction", "direction_robust")
needed_symbols <- c("feature_id", "entrez_id", "gene_symbol")
if (length(setdiff(needed_robust, names(robust))) > 0L) stop_with("Robust-gene result lacks a required column.")
if (length(setdiff(needed_symbols, names(symbols))) > 0L) stop_with("Gene-symbol dictionary lacks a required column.")

summary <- merge(robust, symbols[, needed_symbols, drop = FALSE], by = c("feature_id", "entrez_id"), all.x = TRUE, sort = FALSE)
if (nrow(summary) != 17L || anyNA(summary$gene_symbol) || sum(summary$gene_symbol == "UTY") != 1L) {
  stop_with("Expected 17 mapped robust genes containing exactly one UTY entry.")
}
summary$sex_chromosome_interpretation <- ifelse(
  summary$gene_symbol == "UTY",
  "Y-chromosome gene; interpret separately from non-sex-chromosome molecular context",
  "non-UTY robust gene"
)
summary$retained_after_excluding_uty <- summary$gene_symbol != "UTY"
summary <- summary[order(summary$gene_symbol), , drop = FALSE]

# Only two concise outputs are needed for the manuscript/supplement.
utils::write.csv(sex_by_stage, file.path(results_dir, "donor_sex_by_braak.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(results_dir, "uty_sensitivity.csv"), row.names = FALSE)
message("UTY sensitivity summary completed: ", normalizePath(results_dir))
message("Robust genes: 17; non-UTY genes: ", sum(summary$retained_after_excluding_uty))
