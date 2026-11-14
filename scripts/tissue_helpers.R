options(stringsAsFactors=FALSE)

not_estimable <- function(message) {
  stop(structure(list(message=message,call=NULL),class=c('not_estimable','error','condition')))
}

review_design <- function(d, extra=character(), regional=FALSE, spline_age=FALSE) {
  required <- c('sample_id','donor_key','braak_stage','brain_region',
                'age_at_death','sex_male','pmi_hours','apoe4_carrier',extra)
  if(!all(required %in% names(d)) || anyNA(d[,required])) stop('Missing model metadata')
  d$braak_stage <- factor(d$braak_stage,levels=c('0','I-II','III-IV','V-VI'))
  d$brain_region <- factor(d$brain_region,levels=c('EC','FC','TC'))
  if(anyNA(d$braak_stage)||anyNA(d$brain_region)) stop('Unknown stage or region')
  if(anyDuplicated(d$sample_id)) stop('Duplicate sample identifiers')
  for(k in c('braak_stage','age_at_death','sex_male','pmi_hours','apoe4_carrier'))
    if(any(vapply(split(d[[k]],d$donor_key),function(z)length(unique(z))!=1L,logical(1))))
      stop('Inconsistent donor metadata: ',k)
  for(g in c('0','V-VI')) if(length(unique(d$donor_key[d$braak_stage==g]))<3L)
    not_estimable('Fewer than three donors in a contrast group')
  d <- droplevels(d)
  cov <- c('age_at_death','sex_male','pmi_hours','apoe4_carrier',extra)
  cov <- cov[vapply(cov,function(k)length(unique(d[[k]]))>1L,logical(1))]
  if(spline_age) cov[cov=='age_at_death'] <- 'splines::ns(age_at_death, df=3)'
  if(regional) {
    d$cell <- factor(paste(d$brain_region,d$braak_stage,sep='__'))
    terms <- c('cell',cov)
  } else terms <- c('braak_stage',if(nlevels(d$brain_region)>1L)'brain_region',cov)
  X <- model.matrix(reformulate(terms,intercept=FALSE),d)
  colnames(X) <- make.names(colnames(X))
  if(qr(X)$rank<ncol(X)||nrow(X)<=ncol(X)+5L)
    not_estimable('Rank-deficient or under-supported design; requested covariates not dropped')
  contrast <- function(positive,negative) {
    keys <- make.names(c(positive,negative))
    if(!all(keys %in% colnames(X))) not_estimable('A requested stage-region cell is absent')
    z <- rep(0,ncol(X));z[match(keys,colnames(X))] <- c(1,-1);z
  }
  C <- if(regional) sapply(c('EC','FC','TC'),function(r)
    contrast(paste0('cell',r,'__V-VI'),paste0('cell',r,'__0'))) else
    matrix(contrast('braak_stageV-VI','braak_stage0'),ncol=1,dimnames=list(NULL,'late_vs_zero'))
  list(X=X,C=C,d=d)
}

review_fit <- function(Y, d, extra=character(), regional=FALSE, spline_age=FALSE) {
  if(!identical(colnames(Y),as.character(d$sample_id))||any(!is.finite(Y))||
     is.null(rownames(Y))||anyDuplicated(rownames(Y))) stop('Expression/metadata mismatch')
  design <- review_design(d,extra,regional,spline_age);X <- design$X;d <- design$d
  if(!anyDuplicated(d$donor_key)) not_estimable('No repeated donor observations')
  rho <- limma::duplicateCorrelation(Y,X,block=d$donor_key)$consensus
  if(!is.finite(rho)) stop('Nonfinite donor correlation')
  fit <- limma::lmFit(Y,X,block=d$donor_key,correlation=rho)
  moderated <- limma::eBayes(limma::contrasts.fit(fit,design$C),robust=TRUE)
  effects <- lapply(seq_len(ncol(design$C)),function(j) {
    z <- limma::topTable(moderated,coef=j,number=Inf,sort.by='none',confint=TRUE)
    if(nrow(z)!=nrow(Y)||any(!is.finite(z$P.Value))) stop('Incomplete background tests')
    data.frame(id=rownames(z),region=colnames(design$C)[j],effect=z$logFC,
      ci_low=z$CI.L,ci_high=z$CI.R,p=z$P.Value,q=p.adjust(z$P.Value,'BH'))
  })
  interaction <- NULL
  if(regional) {
    C <- design$C
    H <- cbind(EC_minus_FC=C[,1]-C[,2],EC_minus_TC=C[,1]-C[,3])
    h <- limma::eBayes(limma::contrasts.fit(fit,H),robust=TRUE)
    z <- limma::topTable(h,coef=1:2,number=Inf,sort.by='none')
    if(nrow(z)!=nrow(Y)||any(!is.finite(z$P.Value))) stop('Incomplete interaction family')
    interaction <- data.frame(id=rownames(z),F=z$F,p=z$P.Value,q=p.adjust(z$P.Value,'BH'))
  }
  list(effects=do.call(rbind,effects),interaction=interaction,
       diagnostic=data.frame(n_samples=ncol(Y),n_donors=length(unique(d$donor_key)),
         n_zero=length(unique(d$donor_key[d$braak_stage=='0'])),
         n_late=length(unique(d$donor_key[d$braak_stage=='V-VI'])),
         n_tests=nrow(Y),n_parameters=ncol(X),residual_df=nrow(X)-ncol(X),
         donor_correlation=rho,condition_number=kappa(X)))
}

review_scores <- function(G,sets) {
  sd_gene <- apply(G,1,sd)
  if(any(!is.finite(sd_gene))||any(sd_gene<=0)) stop('Invalid pathway gene SD')
  Z <- (G-rowMeans(G))/sd_gene
  if(any(!unlist(sets,use.names=FALSE)%in%rownames(G))) stop('GO member outside fixed universe')
  scores <- t(vapply(sets,function(ids)colMeans(Z[ids,,drop=FALSE]),numeric(ncol(G))))
  rownames(scores)<-names(sets);colnames(scores)<-colnames(G);scores
}

review_targets <- function(effects,targets) {
  if(!all(targets$feature_id%in%effects$id)) stop('A selected target is absent')
  merge(targets,effects,by.x='feature_id',by.y='id',sort=FALSE)
}

