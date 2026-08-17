#!/usr/bin/env Rscript

# Install the R/Bioconductor environment locked for strict pathway reproduction.
# Run this script using R 4.4.2 in an isolated environment.

if (paste(R.version$major, R.version$minor, sep = ".") != "4.4.2") {
  stop("This installer requires R 4.4.2. Create/use the paper2-r442 environment first.", call. = FALSE)
}

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(BioC_mirror = "https://bioconductor.org")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(version = "3.20", ask = FALSE, update = FALSE)
options(repos = BiocManager::repositories(version = "3.20"))

BiocManager::install(
  c(
    "limma",
    "AnnotationDbi",
    "org.Hs.eg.db",
    "GO.db",
    "huex10sttranscriptcluster.db"
  ),
  version = "3.20",
  ask = FALSE,
  update = FALSE,
  force = TRUE
)

install.packages(
  c("yaml", "readxl", "DBI", "RSQLite"),
  repos = "https://cloud.r-project.org"
)

source("tools/check_r_environment.R", local = TRUE)
