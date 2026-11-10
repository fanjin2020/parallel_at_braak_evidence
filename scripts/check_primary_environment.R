#!/usr/bin/env Rscript
if(getRversion()!='4.4.2')stop('Use R 4.4.2 for this analysis environment')
cran <- c('yaml','readxl','jsonlite','statmod','dtangle')
bioc <- c('BiocVersion','limma','Biobase','AnnotationDbi','org.Hs.eg.db','huex10sttranscriptcluster.db','oligo','pd.huex.1.0.st.v2','preprocessCore')
missing <- c(cran,bioc)[!vapply(c(cran,bioc),requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)){
 stop('Unavailable R packages: ',paste(missing,collapse=', '),
      '\nR library paths: ',paste(.libPaths(),collapse='; '),call.=FALSE)
}
if(!startsWith(as.character(packageVersion('BiocVersion')),'3.20.'))stop('Use Bioconductor 3.20')
if(as.character(packageVersion('oligo'))!='1.70.1')stop('The CEL workflow requires oligo 1.70.1')
for(p in c(cran,bioc))cat(p,as.character(packageVersion(p)),'\n')
