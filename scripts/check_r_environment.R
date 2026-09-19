#!/usr/bin/env Rscript

reference <- c(
  R = "4.4.2", limma = "3.62.2", AnnotationDbi = "1.68.0",
  "org.Hs.eg.db" = "3.20.0", "GO.db" = "3.20.0",
  "huex10sttranscriptcluster.db" = "8.8.0"
)
packages <- c(names(reference)[-1L], "yaml", "readxl", "DBI", "RSQLite")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop(paste("Missing packages:", paste(missing, collapse = ", ")), call. = FALSE)
}
versions <- c(R = as.character(getRversion()),
              setNames(vapply(packages, function(p) as.character(utils::packageVersion(p)), ""), packages))
print(data.frame(component = names(versions), version = unname(versions)), row.names = FALSE)
different <- names(reference)[versions[names(reference)] != reference]
if (length(different)) {
  stop(paste("Package versions differ from the manuscript environment:",
             paste(different, collapse = ", "),
             ". Use the original environment for numerical reproduction."), call. = FALSE)
}
