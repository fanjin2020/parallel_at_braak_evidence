#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
args <- commandArgs(TRUE)
i <- match("--args-file",args)
if(!is.na(i)) args <- readLines(args[i+1L],encoding="UTF-8",warn=FALSE)
arg <- function(k,default=NULL) {
  j <- match(k,args)
  if(is.na(j)) { if(is.null(default)) stop("Missing ",k); return(default) }
  if(j==length(args)) stop("Missing value for ",k)
  args[j+1L]
}
script <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1L])
source(file.path(dirname(normalizePath(script,mustWork=TRUE)),"cel_qc_helpers.R"))
root <- normalizePath(arg("--project-root"),mustWork=TRUE)
arrays <- normalizePath(arg("--cel-dir"),mustWork=TRUE)
out <- arg("--output-dir")
if(dir.exists(out) && length(list.files(out,all.files=TRUE,no..=TRUE)))
  stop("CEL output directory is not empty; choose a new output directory")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
out <- normalizePath(out,mustWork=TRUE)
write_table <- function(x,n) write.csv(x,file.path(out,paste0(n,".csv")),row.names=FALSE,na="NA")
set_state <- function(state,reason) {
  write_table(data.frame(status=state,reason=reason),"processing_status")
  cat("CEL_QC:",state,reason,"\n"); flush.console()
}

main <- function() {
  for(p in c("oligo","pd.huex.1.0.st.v2","Biobase","preprocessCore","limma","statmod"))
    if(!requireNamespace(p,quietly=TRUE)) stop("Missing R package: ",p)
  if(as.character(packageVersion("oligo"))!="1.70.1") stop("Core-probe preflight is verified for oligo 1.70.1; review other versions before use")
  allow <- "--allow-qc-exclusion" %in% args
  maximum_fraction <- as.numeric(arg("--max-excluded-fraction","0.05"))
  base <- file.path(root,"data/public/gse131617")
  targets <- read.csv(file.path(root,"data/contracts/gene_symbols.csv"),colClasses="character")
  cel_validate_targets(targets)
  m <- read.csv(file.path(base,"sample_manifest.csv"))
  files <- read.csv(file.path(arrays,"downloaded_arrays.csv"))
  if(!all(c("sample_id","donor_key","brain_region","braak_stage") %in% names(m)) ||
     !all(c("sample_id","filename") %in% names(files))) stop("Missing manifest/download columns")
  stopifnot(!anyNA(m$sample_id),!anyDuplicated(m$sample_id),!anyDuplicated(files$sample_id))
  files <- files[match(m$sample_id,files$sample_id),,drop=FALSE]
  stopifnot(!anyNA(files$sample_id),all(file.exists(files$filename)))
  write_table(files,"input_arrays")
  capture.output(sessionInfo(),file=file.path(out,"R_session_info.txt"))
  writeLines(c("Post hoc raw-CEL computational QC sensitivity; not a prespecified confirmatory analysis.",
    "Core intensities are checked array-by-array before joint normalization.",
    "Technical failure: nonfinite/negative raw core input or nonfinite/nonpositive corrected core values.",
    "Background-correction execution errors stop the run without an array-exclusion decision.",
    "A high zero fraction is recorded but is not by itself an automatic exclusion rule.",
    paste0("Explicit technical exclusion permission: ",allow,"; maximum excluded fraction: ",maximum_fraction),
    "Excluded input files are preserved. No pseudocounts, imputation, feature filtering, or target-only fitting is used.",
    "RLE review flags: absolute median RLE >0.2 or RLE IQR > median +3 scaled MAD; used only in separately labelled sensitivity models.",
    "Provided-original and provided-QC-pass fits separate sample-exclusion effects from processing differences.",
    "All models share the same intersected full feature background; gene representatives are fixed using all provided original samples.",
    "The 17-target contract must be unique, fully matched to the background, and have finite effect estimates, intervals, P and q values.",
    "RLE/donor subset fits use the QC-pass joint RMA matrix; these subsets are not independently renormalized.",
    "GO membership is intersected with fixed gene representatives (10-300 genes); each model has its own complete tested-set BH family.",
    "Missing RIN, experimental batch or cause of death is not replaced by CEL QC."),file.path(out,"methods.txt"))
  set_state("reading_arrays",paste(nrow(m),"manifest-matched CEL files"))
  raw <- oligo::read.celfiles(files$filename,pkgname="pd.huex.1.0.st.v2")
  Biobase::sampleNames(raw) <- as.character(m$sample_id)
  info <- getFromNamespace("stArrayPmInfo","oligo")(raw,target="core")
  if(!all(c("fid","fsetid") %in% names(info)) || anyNA(info$fid) || anyNA(info$fsetid) ||
     min(info$fid)<1L || max(info$fid)>nrow(Biobase::exprs(raw))) stop("Invalid core-probe mapping")
  qc_rows <- vector("list",nrow(m))
  for(j in seq_len(nrow(m))) {
    cat("Pre-normalization QC:",j,"/",nrow(m),m$sample_id[j],"\n"); flush.console()
    v <- Biobase::exprs(raw)[info$fid,j]
    column_qc <- tryCatch(cel_check_column(v,preprocessCore::rma.background.correct),
      error=function(e) stop("Pre-normalization QC stopped at ",m$sample_id[j],": ",
        conditionMessage(e),call.=FALSE))
    qc_rows[[j]] <- cbind(m[j,c("sample_id","donor_key","brain_region","braak_stage"),drop=FALSE],
      column_qc)
    write_table(do.call(rbind,qc_rows[seq_len(j)]),"technical_array_qc")
  }
  qc <- do.call(rbind,qc_rows)
  write_table(qc[!qc$technical_pass,,drop=FALSE],"technical_exclusions")
  pass <- cel_exclusion_guard(qc,allow,maximum_fraction)
  retained <- m[pass,,drop=FALSE]
  write_table(retained,"retained_sample_manifest")
  write_table(as.data.frame(table(brain_region=m$brain_region,braak_stage=m$braak_stage,
    technical_pass=pass)),"technical_qc_group_counts")
  set_state("RMA",paste(sum(pass),"retained arrays;",sum(!pass),"documented technical exclusions"))
  raw <- raw[,pass]
  rm(v,info); gc()
  es <- oligo::rma(raw,target="core")
  y <- Biobase::exprs(es)
  saveRDS(y,file.path(out,"rma_core_expression_UNFILTERED.rds"),compress=FALSE)
  write_table(data.frame(sample_id=colnames(y),n_features=nrow(y),n_nonfinite=colSums(!is.finite(y))),"rma_finiteness")
  y <- cel_require_finite(y,as.character(retained$sample_id))
  rm(raw,es); gc()
  con <- gzfile(file.path(out,"rma_core_expression.csv.gz"),"wt")
  tryCatch(write.csv(data.frame(feature_id=rownames(y),y,check.names=FALSE),con,row.names=FALSE),finally=close(con))
  pc <- prcomp(t(y),center=TRUE,scale.=FALSE)
  rle <- sweep(y,1L,apply(y,1L,median),"-")
  quality <- data.frame(sample_id=retained$sample_id,rle_median=apply(rle,2L,median),
    rle_IQR=apply(rle,2L,IQR),PC1=pc$x[,1L],PC2=pc$x[,2L])
  mad_iqr <- mad(quality$rle_IQR,constant=1.4826)
  quality$flag_rle_spread <- if(mad_iqr>0) quality$rle_IQR>median(quality$rle_IQR)+3*mad_iqr else FALSE
  quality$flag_rle_center <- abs(quality$rle_median)>.2
  quality$flag_review <- quality$flag_rle_spread | quality$flag_rle_center
  write_table(quality,"array_quality")
  pdf(file.path(out,"array_quality_report.pdf"),width=11,height=7)
  tryCatch({
    par(mfrow=c(1,2))
    plot(quality$rle_median,ylab="Median gene-level RLE",xlab="Retained array index",pch=16,
      col=ifelse(quality$flag_review,"firebrick","grey40")); abline(h=c(-.2,.2),lty=2)
    plot(quality$rle_IQR,ylab="Gene-level RLE IQR",xlab="Retained array index",pch=16,
      col=ifelse(quality$flag_review,"firebrick","grey40")); abline(h=median(quality$rle_IQR)+3*mad_iqr,lty=2)
    par(mfrow=c(1,1))
    plot(quality$PC1,quality$PC2,col=as.integer(factor(retained$brain_region)),
      pch=as.integer(factor(retained$braak_stage))+14,xlab="PC1",ylab="PC2",main="Technical-QC-pass RMA core expression")
    legend("topright",legend=levels(factor(retained$brain_region)),col=seq_along(levels(factor(retained$brain_region))),pch=16)
  },finally=dev.off())
  rm(rle,pc); gc()
  old <- read.csv(gzfile(file.path(base,"expression_feature.csv.gz")),check.names=FALSE)
  old_y <- as.matrix(old[,-1L]); storage.mode(old_y) <- "double"; rownames(old_y) <- old[[1L]]
  old_y <- cel_require_finite(old_y,as.character(m$sample_id))
  common <- intersect(rownames(y),rownames(old_y))
  if(length(common)<1000L) stop("Insufficient common features for full-background comparison")
  old_y <- old_y[common,,drop=FALSE]; new_y <- y[common,,drop=FALSE]
  old_pair <- old_y[,retained$sample_id,drop=FALSE]
  write_table(data.frame(sample_id=retained$sample_id,n_common_features=length(common),
    pearson=vapply(seq_len(ncol(new_y)),function(j)cor(old_pair[,j],new_y[,j]),numeric(1)),
    median_difference=apply(new_y-old_pair,2L,median)),"processed_matrix_agreement")
  cel_validate_targets(targets,common)
  annotation <- read.csv(file.path(base,"feature_annotation.csv"),colClasses="character")
  annotation <- annotation[annotation$mapping_status=="unique_entrez_id" & annotation$feature_id %in% common,,drop=FALSE]
  annotation$average <- rowMeans(old_y[annotation$feature_id,,drop=FALSE])
  annotation <- annotation[order(annotation$entrez_id,-annotation$average,annotation$feature_id),,drop=FALSE]
  annotation <- annotation[!duplicated(annotation$entrez_id),,drop=FALSE]
  write_table(annotation,"fixed_gene_representatives")
  go <- arg("--go-dir","")
  sets <- NULL
  if(nzchar(go)) {
    members <- read.csv(file.path(go,"go_gene_membership.csv"),colClasses="character")
    sets <- lapply(split(members$entrez_id,members$go_id),function(ids)intersect(unique(ids),annotation$entrez_id))
    sizes <- lengths(sets)
    write_table(data.frame(go_id=names(sets),n_measurable_genes=sizes,tested=sizes>=10L & sizes<=300L),"GO_coverage")
    sets <- sets[sizes>=10L & sizes<=300L]
    if(!length(sets)) stop("No measurable GO sets under the specified size limits")
  }
  scopes <- list(original_provided=list(ids=as.character(m$sample_id),matrices="provided"),
    technical_qc_pass=list(ids=as.character(retained$sample_id),matrices=c("provided","reprocessed")),
    exclude_rle_flagged=list(ids=as.character(quality$sample_id[!quality$flag_review]),matrices=c("provided","reprocessed")))
  excluded_donors <- unique(m$donor_key[!pass])
  if(length(excluded_donors)) scopes$exclude_technical_failure_donors <- list(
    ids=as.character(m$sample_id[!m$donor_key %in% excluded_donors]),matrices=c("provided","reprocessed"))
  fits <- go_rows <- states <- cohorts <- list()
  save_states <- function() write_table(do.call(rbind,states),"analysis_status")
  for(scope in names(scopes)) for(matrix_name in scopes[[scope]]$matrices) {
    label <- paste(scope,matrix_name,sep="_")
    ids <- scopes[[scope]]$ids
    d <- m[match(ids,m$sample_id),,drop=FALSE]
    cohorts[[label]] <- data.frame(model=label,n_samples=length(ids),n_donors=length(unique(d$donor_key)))
    write_table(do.call(rbind,cohorts),"model_cohorts")
    write_table(d,paste0("samples_",label))
    cat("Full-background model:",label,"n =",nrow(d),"\n"); flush.console()
    result <- tryCatch({
      if(nrow(d)<3L) cel_not_estimable("Fewer than three samples in the QC subset")
      Y <- if(matrix_name=="provided") old_y[,ids,drop=FALSE] else new_y[,ids,drop=FALSE]
      f <- cel_fit_full_background(Y,d); tt <- f$table
      all <- data.frame(feature_id=rownames(tt),log2FC=tt$logFC,ci_low=tt$CI.L,ci_high=tt$CI.R,
        p=tt$P.Value,q=p.adjust(tt$P.Value,"BH"))
      write_table(all,paste0("features_",label))
      target_result <- cel_extract_targets(all,targets)
      target_result$model <- label; target_result$n_samples <- nrow(d)
      target_result$n_donors <- length(unique(d$donor_key)); target_result$n_tested <- nrow(Y)
      fits[[label]] <- target_result
      write_table(do.call(rbind,fits),"target_processing_comparison")
      if(!is.null(sets)) {
        G <- Y[annotation$feature_id,,drop=FALSE]; rownames(G) <- annotation$entrez_id
        tg <- cel_fit_full_background(G,d)$table
        z <- sign(tg$t)*qnorm(pmax(tg$P.Value,.Machine$double.xmin)/2,lower.tail=FALSE); names(z) <- rownames(tg)
        gs <- limma::cameraPR(z,lapply(sets,function(ids)match(ids,names(z))),inter.gene.cor=.01,use.ranks=FALSE,sort=FALSE)
        if(nrow(gs)!=length(sets) || any(!is.finite(gs$PValue))) stop("Incomplete/nonfinite GO inference")
        go_rows[[label]] <- data.frame(model=label,go_id=rownames(gs),n_genes=gs$NGenes,
          direction=gs$Direction,p=gs$PValue,q=p.adjust(gs$PValue,"BH"))
        write_table(do.call(rbind,go_rows),"GO_processing_comparison")
      }
      data.frame(component=label,status="completed",reason="")
    },cel_not_estimable=function(e) data.frame(component=label,status="not_estimable",reason=conditionMessage(e)),
      error=function(e) data.frame(component=label,status="failed",reason=conditionMessage(e)))
    states[[label]] <- result
    save_states()
  }
  statuses <- do.call(rbind,states)
  if(any(statuses$status=="failed")) stop("One or more models failed; inspect analysis_status.csv")
  if(!length(fits)) stop("No full-background QC model was estimable")
  set_state(if(all(statuses$status=="completed")) "completed" else "partial",
    paste("Technical exclusions:",sum(!pass),"; inspect model status and exclusion records before interpretation"))
  capture.output(sessionInfo(),file=file.path(out,"R_session_info.txt"))
}
tryCatch(main(),error=function(e) {
  set_state("failed",conditionMessage(e))
  stop(conditionMessage(e),call.=FALSE)
})
