cel_check_column <- function(v, background_fun) {
  finite <- is.finite(v)
  s <- data.frame(n_core_probes=length(v), n_raw_nonfinite=sum(!finite),
    n_raw_negative=sum(v < 0 & finite, na.rm=TRUE), n_raw_zero=sum(v == 0, na.rm=TRUE),
    zero_fraction=sum(v == 0, na.rm=TRUE)/length(v),
    raw_min=if(any(finite)) min(v[finite]) else NA_real_,
    raw_median=if(any(finite)) median(v[finite]) else NA_real_,
    raw_max=if(any(finite)) max(v[finite]) else NA_real_,
    n_background_nonfinite=NA_integer_, n_background_nonpositive=NA_integer_,
    technical_pass=FALSE, reason="", background_warning="", stringsAsFactors=FALSE)
  if(s$n_raw_nonfinite > 0 || s$n_raw_negative > 0) {
    s$reason <- "nonfinite_or_negative_core_input"
    return(s)
  }
  warnings <- character()
  b <- tryCatch(withCallingHandlers(background_fun(matrix(as.double(v), ncol=1L)),
    warning=function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error=function(e) e)
  s$background_warning <- paste(unique(warnings), collapse=" | ")
  if(inherits(b,"error")) {
    stop("Background correction failed; no exclusion decision made: ",
      conditionMessage(b),call.=FALSE)
  }
  if(!is.numeric(b) || !identical(dim(b), c(length(v),1L))) {
    stop("Background correction returned unexpected dimensions/type; no exclusion decision made")
  }
  s$n_background_nonfinite <- sum(!is.finite(b))
  s$n_background_nonpositive <- sum(is.finite(b) & b <= 0, na.rm=TRUE)
  s$technical_pass <- s$n_background_nonfinite == 0 && s$n_background_nonpositive == 0
  s$reason <- if(s$technical_pass) "pass" else "nonfinite_or_nonpositive_background"
  s
}

cel_exclusion_guard <- function(qc, allow=FALSE, maximum_fraction=.05) {
  if(!is.numeric(maximum_fraction) || length(maximum_fraction)!=1L ||
     !is.finite(maximum_fraction) || maximum_fraction<0 || maximum_fraction>=1)
    stop("maximum_fraction must be in [0, 1)")
  if(!nrow(qc) || anyNA(qc$technical_pass)) stop("Incomplete technical QC decisions")
  failed <- !qc$technical_pass
  if(any(failed) && !allow)
    stop("Technical QC failures recorded. Review technical_array_qc.csv; explicit --allow-qc-exclusion is required to continue")
  if(mean(failed)>maximum_fraction)
    stop("Technical exclusions exceed the allowed fraction; review inputs/environment before proceeding")
  if(sum(!failed)<3L) stop("Fewer than three technically usable arrays")
  !failed
}

cel_require_finite <- function(y, expected_samples) {
  if(!is.matrix(y) || !is.numeric(y) || !nrow(y) || any(!is.finite(y)))
    stop("Expression matrix contains nonfinite values or is invalid; saved unfiltered output must be reviewed")
  if(is.null(rownames(y)) || anyDuplicated(rownames(y)) || anyDuplicated(colnames(y)) ||
     !setequal(colnames(y),expected_samples)) stop("Expression identifiers do not match the analysis samples")
  y[,expected_samples,drop=FALSE]
}

cel_validate_targets <- function(targets, feature_ids=NULL, expected_n=17L) {
  required <- c("feature_id","gene_symbol")
  if(!is.data.frame(targets) || !all(required %in% names(targets)))
    stop("Target contract requires feature_id and gene_symbol")
  if(nrow(targets)!=expected_n)
    stop("Target contract must contain exactly ",expected_n," targets; found ",nrow(targets))
  for(k in required) {
    ids <- as.character(targets[[k]])
    if(anyNA(ids) || any(!nzchar(trimws(ids))) || any(ids!=trimws(ids)))
      stop("Missing, empty or whitespace-padded target ",k)
    if(anyDuplicated(ids)) stop("Duplicate target ",k)
  }
  if(!is.null(feature_ids)) {
    missing <- setdiff(targets$feature_id,feature_ids)
    if(length(missing)) stop("Targets absent from the full feature background: ",paste(missing,collapse=", "))
  }
  invisible(TRUE)
}

cel_extract_targets <- function(all, targets, expected_n=17L) {
  stats <- c("log2FC","ci_low","ci_high","p","q")
  if(!is.data.frame(all) || !all(c("feature_id",stats) %in% names(all)))
    stop("Full-feature results are missing required columns")
  if(anyNA(all$feature_id) || any(!nzchar(trimws(as.character(all$feature_id)))) ||
     anyDuplicated(all$feature_id)) stop("Invalid or duplicate full-feature identifiers")
  cel_validate_targets(targets,all$feature_id,expected_n)
  if(any(stats %in% names(targets))) stop("Target contract contains reserved statistical columns")
  result <- merge(all,targets,by="feature_id",all.y=TRUE,sort=FALSE)
  if(nrow(result)!=expected_n || anyDuplicated(result$feature_id) ||
     !setequal(result$feature_id,targets$feature_id)) stop("Incomplete or duplicated target results")
  if(!all(vapply(result[,stats,drop=FALSE],is.numeric,logical(1))) ||
     any(!is.finite(as.matrix(result[,stats,drop=FALSE]))))
    stop("Target results contain missing, nonnumeric or nonfinite statistics")
  if(any(result$p<0 | result$p>1 | result$q<0 | result$q>1))
    stop("Target P or q values are outside [0, 1]")
  if(any(result$ci_low>result$log2FC | result$ci_high<result$log2FC))
    stop("Target confidence intervals do not contain their estimates")
  result
}

cel_not_estimable <- function(reason) {
  stop(structure(list(message=reason,call=NULL),class=c("cel_not_estimable","error","condition")))
}

cel_fit_full_background <- function(Y,d) {
  if(!identical(colnames(Y),as.character(d$sample_id))) stop("Model matrix/sample metadata are not aligned")
  if(any(!is.finite(Y))) stop("Nonfinite expression reached model fitting")
  required <- c("sample_id","donor_key","braak_stage","brain_region","age_at_death","sex_male","pmi_hours","apoe4_carrier")
  if(!all(required %in% names(d)) || anyNA(d[,required])) stop("Missing model metadata; no implicit complete-case filtering permitted")
  d$braak_stage <- factor(d$braak_stage,levels=c("0","I-II","III-IV","V-VI"))
  d$brain_region <- factor(d$brain_region)
  if(anyNA(d$braak_stage)) stop("Unexpected Braak stage")
  x <- model.matrix(~0+braak_stage+brain_region+age_at_death+sex_male+pmi_hours+apoe4_carrier,d)
  if(qr(x)$rank<ncol(x) || nrow(x)<=ncol(x)+5L) cel_not_estimable("Rank-deficient model or insufficient residual observations")
  if(length(unique(d$donor_key))<3L || !anyDuplicated(d$donor_key)) cel_not_estimable("Insufficient repeated donor observations")
  rho <- limma::duplicateCorrelation(Y,x,block=d$donor_key)$consensus
  if(!is.finite(rho)) stop("Nonfinite estimated donor correlation")
  co <- matrix(0,ncol(x),1L)
  hi <- match("braak_stageV-VI",colnames(x)); lo <- match("braak_stage0",colnames(x))
  if(anyNA(c(hi,lo))) cel_not_estimable("Late/reference Braak contrast unavailable")
  co[hi,1L] <- 1; co[lo,1L] <- -1
  f <- limma::eBayes(limma::contrasts.fit(limma::lmFit(Y,x,block=d$donor_key,correlation=rho),co),robust=TRUE)
  tt <- limma::topTable(f,number=Inf,sort.by="none",confint=TRUE)
  if(nrow(tt)!=nrow(Y) || any(!is.finite(tt$P.Value))) stop("Incomplete/nonfinite full-background inference")
  list(table=tt,correlation=rho,n_samples=nrow(d),n_donors=length(unique(d$donor_key)))
}
