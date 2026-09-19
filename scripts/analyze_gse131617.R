# Fit donor-blocked late-Braak differential-expression models.
#!/usr/bin/env Rscript 

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

fit_braak_model <- function(expression, manifest) {
  design <- stats::model.matrix(
    ~ 0 + braak_stage + brain_region + age_at_death + sex_male + pmi_hours + apoe4_carrier,
    data = manifest
  )
  colnames(design) <- make.names(colnames(design), unique = TRUE)
  if (qr(design)$rank != ncol(design)) stop_with("Design matrix is rank deficient.")

  donor_correlation <- limma::duplicateCorrelation(expression, design, block = manifest$donor_key)$consensus
  fit <- limma::lmFit(expression, design, block = manifest$donor_key, correlation = donor_correlation)
  contrast <- limma::makeContrasts(late_vs_zero = braak_stageV.VI - braak_stage0, levels = design)
  fit <- limma::eBayes(limma::contrasts.fit(fit, contrast), robust = TRUE)
  result <- limma::topTable(fit, coef = "late_vs_zero", number = Inf, sort.by = "none")

  data.frame(
    feature_id = rownames(result),
    log2_fold_change = result$logFC,
    average_expression = result$AveExpr,
    moderated_t = result$t,
    p_value = result$P.Value,
    fdr_bh = stats::p.adjust(result$P.Value, method = "BH"),
    B = result$B,
    check.names = FALSE
  )
}

select_representatives <- function(primary, annotation, alpha) {
  merged <- merge(primary, annotation, by = "feature_id", all.x = TRUE, sort = FALSE)
  if (nrow(merged) != nrow(primary)) stop_with("Annotation changed the number of features.")

  eligible <- merged[
    merged$mapping_status == "unique_entrez_id" &
      grepl("^[0-9]+$", merged$entrez_id) &
      merged$fdr_bh < alpha,
    , drop = FALSE
  ]
  eligible$abs_effect <- abs(eligible$log2_fold_change)
  eligible <- eligible[order(eligible$entrez_id, eligible$fdr_bh, eligible$p_value,
                             -eligible$abs_effect, eligible$feature_id), , drop = FALSE]
  eligible <- eligible[!duplicated(eligible$entrez_id), , drop = FALSE]
  eligible$primary_direction <- ifelse(eligible$log2_fold_change > 0, "up", "down")
  eligible[, c("feature_id", "entrez_id", "log2_fold_change", "average_expression",
               "moderated_t", "p_value", "fdr_bh", "B", "primary_direction"), drop = FALSE]
}

require_package("limma")
require_package("yaml")

config_file <- normalizePath(argument_value("--config"), mustWork = TRUE)
root <- dirname(dirname(config_file))
config <- yaml::read_yaml(config_file)

expression_file <- configured_path(config, root, "gse_feature_expression", "data/public/gse131617/expression_feature.csv.gz")
annotation_file <- configured_path(config, root, "gse_feature_annotation", "data/public/gse131617/feature_annotation.csv")
manifest_file <- configured_path(config, root, "gse_manifest", "data/public/gse131617/sample_manifest.csv")
results_dir <- configured_path(config, root, "gse_results_dir", "results/gse131617")
alpha_value <- config$analysis$fdr_alpha
if (is.null(alpha_value) || length(alpha_value) != 1L || is.na(alpha_value)) alpha_value <- 0.05
alpha <- as.numeric(alpha_value)
for (file in c(expression_file, annotation_file, manifest_file)) {
  if (!file.exists(file)) stop_with(paste("Missing input:", file))
}
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

expression_table <- utils::read.csv(gzfile(expression_file), check.names = FALSE)
if (names(expression_table)[1L] != "feature_id" || anyDuplicated(expression_table$feature_id)) {
  stop_with("Expression table must begin with unique feature_id values.")
}
expression <- as.matrix(expression_table[, -1L, drop = FALSE])
storage.mode(expression) <- "double"
rownames(expression) <- expression_table$feature_id
if (anyNA(expression) || any(!is.finite(expression))) stop_with("Expression matrix contains missing or non-finite values.")

annotation <- utils::read.csv(annotation_file, colClasses = "character", check.names = FALSE)
manifest <- utils::read.csv(manifest_file, check.names = FALSE)
needed <- c("sample_id", "donor_key", "braak_stage", "brain_region", "age_at_death", "sex_male", "pmi_hours", "apoe4_carrier")
if (length(setdiff(needed, names(manifest))) > 0L) stop_with("Sample manifest lacks a required column.")
if (!setequal(colnames(expression), manifest$sample_id)) stop_with("Expression and manifest sample IDs do not match.")
manifest <- manifest[match(colnames(expression), manifest$sample_id), , drop = FALSE]

manifest$braak_stage <- factor(manifest$braak_stage, levels = c("0", "I-II", "III-IV", "V-VI"))
manifest$brain_region <- factor(manifest$brain_region, levels = c("EC", "FC", "TC"))
manifest$donor_key <- factor(manifest$donor_key)
for (column in c("age_at_death", "sex_male", "pmi_hours", "apoe4_carrier")) manifest[[column]] <- as.numeric(manifest[[column]])
if (anyNA(manifest[, c("age_at_death", "sex_male", "pmi_hours", "apoe4_carrier")])) stop_with("Manifest covariates contain missing values.")

primary <- fit_braak_model(expression, manifest)

regional <- list()
for (region in levels(manifest$brain_region)) {
  keep <- manifest$brain_region != region
  temporary <- fit_braak_model(expression[, keep, drop = FALSE], droplevels(manifest[keep, , drop = FALSE]))
  regional[[region]] <- temporary[, c("feature_id", "log2_fold_change")]
  names(regional[[region]])[2L] <- paste0("log2_fold_change_leave_out_", region)
}
sensitivity <- Reduce(function(left, right) merge(left, right, by = "feature_id", all = TRUE, sort = FALSE),
                      c(list(primary), regional))
sensitivity$direction_robust <- with(sensitivity,
  sign(log2_fold_change) == sign(log2_fold_change_leave_out_EC) &
    sign(log2_fold_change) == sign(log2_fold_change_leave_out_FC) &
    sign(log2_fold_change) == sign(log2_fold_change_leave_out_TC)
)

representatives <- select_representatives(primary, annotation, alpha)
representatives <- merge(representatives,
  sensitivity[, c("feature_id", "log2_fold_change_leave_out_EC", "log2_fold_change_leave_out_FC",
                  "log2_fold_change_leave_out_TC", "direction_robust")],
  by = "feature_id", all.x = TRUE, sort = FALSE
)
robust_genes <- representatives[representatives$direction_robust %in% TRUE, , drop = FALSE]

utils::write.csv(primary, file.path(results_dir, "braak_feature_results.csv"), row.names = FALSE)
utils::write.csv(sensitivity, file.path(results_dir, "braak_region_sensitivity.csv"), row.names = FALSE)
utils::write.csv(representatives, file.path(results_dir, "braak_representative_genes.csv"), row.names = FALSE)
utils::write.csv(robust_genes, file.path(results_dir, "braak_robust_genes.csv"), row.names = FALSE)

message("GSE131617 analysis completed: ", normalizePath(results_dir))
message("Directionally robust genes: ", nrow(robust_genes))

