options(stringsAsFactors=FALSE)

sex_interaction_design <- function(d, extra=character()) {
  required <- c('sample_id','donor_key','braak_stage','brain_region','sex_male',
                'age_at_death','pmi_hours','apoe4_carrier',extra)
  if(!all(required %in% names(d)) || anyNA(d[,required])) stop('Missing interaction metadata')
  if(anyDuplicated(d$sample_id) || any(!d$sex_male %in% 0:1)) stop('Invalid sample/sex metadata')
  stages <- c('0','I-II','III-IV','V-VI')
  d$braak_stage <- factor(d$braak_stage,levels=stages)
  d$brain_region <- factor(d$brain_region,levels=c('EC','FC','TC'))
  if(anyNA(d$braak_stage) || anyNA(d$brain_region)) stop('Unknown Braak stage or brain region')
  for(k in c('braak_stage','sex_male','age_at_death','pmi_hours','apoe4_carrier'))
    if(any(vapply(split(d[[k]],d$donor_key),function(v)length(unique(v))!=1L,logical(1))))
      stop('Inconsistent donor metadata: ',k)
  donor <- d[!duplicated(d$donor_key),]
  counts <- table(factor(donor$sex_male,levels=0:1),factor(donor$braak_stage,levels=stages))
  if(any(counts<3L)) not_estimable('Fewer than three donors in a sex-by-Braak cell; the full interaction is not estimated')
  d$sex_braak <- factor(paste0(ifelse(d$sex_male==1,'male','female'),'__',d$braak_stage),
                       levels=as.vector(outer(c('female','male'),stages,paste,sep='__')))
  cov <- c('age_at_death','pmi_hours','apoe4_carrier',extra)
  dropped <- cov[!vapply(cov,function(k)length(unique(d[[k]]))>1L,logical(1))]
  cov <- setdiff(cov,dropped)
  X <- model.matrix(reformulate(c('sex_braak','brain_region',cov),intercept=FALSE),d)
  colnames(X) <- make.names(colnames(X))
  if(qr(X)$rank<ncol(X) || nrow(X)<=ncol(X)+5L)
    not_estimable('Rank-deficient or under-supported sex-interaction design; nonconstant covariates were not dropped')
  contrast <- function(sex,stage) {
    keys <- make.names(paste0('sex_braak',sex,'__',c(stage,'0')))
    if(!all(keys %in% colnames(X))) stop('Missing sex-by-stage coefficient')
    z <- setNames(rep(0,ncol(X)),colnames(X));z[keys] <- c(1,-1);z
  }
  female <- contrast('female','V-VI');male <- contrast('male','V-VI')
  H <- sapply(stages[-1],function(stage)contrast('male',stage)-contrast('female',stage))
  colnames(H) <- c('interaction_early','interaction_middle','interaction_late')
  C <- cbind(female_late=female,male_late=male,interaction_late=male-female)
  stopifnot(qr(H)$rank==3L)
  list(X=X,C=C,H=H,d=d,dropped=dropped,counts=counts)
}

sex_interaction_fit <- function(Y,d,extra=character()) {
  if(!identical(colnames(Y),as.character(d$sample_id)) || any(!is.finite(Y)) ||
     is.null(rownames(Y)) || anyDuplicated(rownames(Y))) stop('Expression/metadata mismatch')
  design <- sex_interaction_design(d,extra);X <- design$X
  if(!anyDuplicated(d$donor_key)) not_estimable('No repeated donor observations')
  rho <- limma::duplicateCorrelation(Y,X,block=d$donor_key)$consensus
  if(!is.finite(rho)) stop('Nonfinite donor correlation')
  fit <- limma::lmFit(Y,X,block=d$donor_key,correlation=rho)
  one <- limma::eBayes(limma::contrasts.fit(fit,design$C),robust=TRUE)
  result <- data.frame(id=rownames(Y))
  for(j in seq_len(ncol(design$C))) {
    z <- limma::topTable(one,coef=j,number=Inf,sort.by='none',confint=TRUE)
    stopifnot(identical(rownames(z),rownames(Y)),all(is.finite(z$P.Value)))
    v <- data.frame(effect=z$logFC,ci_low=z$CI.L,ci_high=z$CI.R,t=z$t,p=z$P.Value,
                    q=p.adjust(z$P.Value,'BH'))
    names(v) <- paste0(colnames(design$C)[j],'_',names(v));result <- cbind(result,v)
  }
  joint <- limma::eBayes(limma::contrasts.fit(fit,design$H),robust=TRUE)
  z <- limma::topTable(joint,coef=1:3,number=Inf,sort.by='none')
  stopifnot(identical(rownames(z),rownames(Y)),all(is.finite(z$P.Value)))
  result$omnibus_F <- z$F;result$omnibus_p <- z$P.Value;result$omnibus_q <- p.adjust(z$P.Value,'BH')
  result$omnibus_df <- 3L
  counts <- design$counts
  diag <- data.frame(n_samples=ncol(Y),n_donors=length(unique(d$donor_key)),n_tests=nrow(Y),
    n_parameters=ncol(X),residual_df=nrow(X)-ncol(X),donor_correlation=rho,
    condition_number=kappa(X),female_zero_donors=counts[1,1],male_zero_donors=counts[2,1],
    female_late_donors=counts[1,4],male_late_donors=counts[2,4],
    constant_covariates_omitted=paste(design$dropped,collapse=';'))
  list(effects=result,diagnostic=diag)
}

