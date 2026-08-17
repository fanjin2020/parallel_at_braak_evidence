#!/usr/bin/env Rscript

# Verify the R/Bioconductor snapshot required for strict pathway reproduction.
# This script uses only base R before it tests package availability.

args <- commandArgs(trailingOnly = TRUE)
report_path <- "validation/r_environment_lock_check.csv"
if (length(args) > 0L) {
  report_flag <- which(args == "--report")
  if (length(report_flag) == 1L && report_flag < length(args)) {
    report_path <- args[[report_flag + 1L]]
  }
}

locked <- c(
  R = "4.4.2",
  limma = "3.62.2",
  AnnotationDbi = "1.68.0",
  "org.Hs.eg.db" = "3.20.0",
  "GO.db" = "3.20.0",
  "huex10sttranscriptcluster.db" = "8.8.0"
)

operational <- c("BiocManager", "yaml", "readxl", "DBI", "RSQLite")

detected_r <- paste(R.version$major, R.version$minor, sep = ".")
records <- data.frame(
  component = "R",
  expected = locked[["R"]],
  detected = detected_r,
  required_for_strict_pathway_reproduction = TRUE,
  stringsAsFactors = FALSE
)

for (package_name in names(locked)[names(locked) != "R"]) {
  installed <- requireNamespace(package_name, quietly = TRUE)
  detected <- if (installed) as.character(utils::packageVersion(package_name)) else NA_character_
  records <- rbind(
    records,
    data.frame(
      component = package_name,
      expected = locked[[package_name]],
      detected = detected,
      required_for_strict_pathway_reproduction = TRUE,
      stringsAsFactors = FALSE
    )
  )
}

for (package_name in operational) {
  installed <- requireNamespace(package_name, quietly = TRUE)
  detected <- if (installed) as.character(utils::packageVersion(package_name)) else NA_character_
  records <- rbind(
    records,
    data.frame(
      component = package_name,
      expected = "installed (version recorded)",
      detected = detected,
      required_for_strict_pathway_reproduction = FALSE,
      stringsAsFactors = FALSE
    )
  )
}

locked_rows <- records$required_for_strict_pathway_reproduction
records$status <- ifelse(
  locked_rows,
  ifelse(!is.na(records$detected) & records$detected == records$expected, "OK", "MISMATCH"),
  ifelse(!is.na(records$detected), "RECORDED", "MISSING")
)

report_dir <- dirname(report_path)
if (!identical(report_dir, ".")) {
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
}
utils::write.csv(records, report_path, row.names = FALSE, na = "")

print(records, row.names = FALSE)
cat("\nR environment report written to:", normalizePath(report_path, mustWork = FALSE), "\n")

if (any(records$status[locked_rows] != "OK")) {
  stop("Locked R/Bioconductor environment check failed. Do not run strict pathway reproduction.", call. = FALSE)
}
