# Exploratory OPC-minus-INH difference in Braak effects, with paired donors.

fit_pair <- function(a, b, m, plus_pmi=FALSE) {
  for (p in c('edgeR','limma','statmod'))
    if (!requireNamespace(p, quietly=TRUE)) stop('Missing package: ', p, '; use the existing paper2-r442 environment')
  if (!identical(dim(a),dim(b)) || !identical(rownames(a),rownames(b)) ||
      !identical(colnames(a),colnames(b)) || !identical(colnames(a),as.character(m$donor_id)))
    stop('Paired count matrices or donor order differ')
  for (x in list(a,b)) if (anyNA(x) || any(!is.finite(x)) || any(x<0) || any(x!=round(x)))
    stop('Expected finite nonnegative integer pseudobulk counts')
  m$late <- as.integer(m$braak>=5); m$batch <- factor(m$batch)
  low <- sum(m$late==0); high <- sum(m$late==1)
  base <- list(n_low=low,n_high=high,n_donors=nrow(m),n_pseudobulks=2L*nrow(m),
               status='not_estimable',reason='',block_correlation=NA_real_,genes_tested=NA_integer_,
               independent_design_df=NA_integer_)
  if (min(low,high)<6) {base$reason <- 'Fewer than six paired independent donors in a stage group'; return(base)}
  formula <- if (plus_pmi) ~late+age+sex_male+batch+rin+pmi else ~late+age+sex_male+batch+rin
  X0 <- tryCatch(model.matrix(formula,m), error=function(e)e)
  if (inherits(X0,'error')) {base$reason <- conditionMessage(X0); return(base)}
  base$independent_design_df <- nrow(m)-ncol(X0)
  if (nrow(X0)!=nrow(m) || qr(X0)$rank!=ncol(X0) || base$independent_design_df<5) {
    base$reason <- 'Incomplete or rank-deficient donor design, or fewer than five donor-level residual df'
    return(base)
  }
  da <- edgeR::DGEList(counts=a); db <- edgeR::DGEList(counts=b)
  expressed <- edgeR::filterByExpr(da,design=X0) & edgeR::filterByExpr(db,design=X0)
  if (sum(expressed)<100) {base$reason <- 'Fewer than 100 features expressed in both cell types'; return(base)}
  da <- edgeR::calcNormFactors(da[expressed,,keep.lib.sizes=FALSE],method='TMM')
  db <- edgeR::calcNormFactors(db[expressed,,keep.lib.sizes=FALSE],method='TMM')
  joined <- cbind(da$counts,db$counts)
  colnames(joined) <- c(paste0(m$donor_id,'__OPC'),paste0(m$donor_id,'__INH'))
  y <- edgeR::DGEList(counts=joined,
                     lib.size=c(da$samples$lib.size,db$samples$lib.size),
                     norm.factors=c(da$samples$norm.factors,db$samples$norm.factors))
  colnames(y) <- c(paste0(m$donor_id,'__OPC'),paste0(m$donor_id,'__INH'))
  zero <- matrix(0,nrow(X0),ncol(X0))
  X <- rbind(cbind(X0,zero),cbind(zero,X0))
  colnames(X) <- c(paste0('OPC__',make.names(colnames(X0))),paste0('INH__',make.names(colnames(X0))))
  rownames(X) <- colnames(y)
  block <- factor(rep(m$donor_id,2))
  v <- limma::voom(y,X,plot=FALSE)
  rho <- limma::duplicateCorrelation(v,X,block=block)$consensus.correlation
  valid_rho <- function(x) length(x)==1 && is.finite(x) && x > -1+1e-6 && x < 1-1e-6
  if (!valid_rho(rho)) stop('Unusable within-donor consensus correlation; no independent-sample fallback')
  v <- limma::voom(y,X,block=block,correlation=rho,plot=FALSE)
  rho <- limma::duplicateCorrelation(v,X,block=block)$consensus.correlation
  if (!valid_rho(rho)) stop('Unusable updated within-donor correlation')
  v <- limma::voom(y,X,block=block,correlation=rho,plot=FALSE)
  fit <- limma::lmFit(v,X,block=block,correlation=rho)
  contrast <- rep(0,ncol(X)); names(contrast)<-colnames(X)
  contrast['OPC__late']<-1; contrast['INH__late']<--1
  if (length(contrast)!=ncol(X)) stop('Stage coefficients not found in interaction design')
  fit <- limma::eBayes(limma::contrasts.fit(fit,contrasts=matrix(contrast,ncol=1,dimnames=list(colnames(X),'OPC_minus_INH'))),robust=TRUE)
  tab <- limma::topTable(fit,coef=1,number=Inf,sort.by='none',confint=TRUE)
  base$status<-'estimated';base$block_correlation<-rho;base$genes_tested<-nrow(tab)
  base$table<-tab;base$t<-setNames(fit$t[,1],rownames(fit$t));base$formula<-paste(deparse(formula),collapse='')
  base
}


run_interaction <- function(src,go,out) {
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  write <- function(x,name) write.csv(x,file.path(out,paste0(name,'.csv')),row.names=FALSE,na='')
  donors<-read.csv(file.path(src,'donors.csv'),stringsAsFactors=FALSE)
  genes<-read.csv(file.path(src,'genes.csv'),stringsAsFactors=FALSE)
  cells<-read.csv(file.path(src,'donor_cell_counts.csv'),check.names=FALSE,stringsAsFactors=FALSE)
  if (anyDuplicated(donors$donor_id) || anyDuplicated(genes$gene_id) || anyDuplicated(cells$donor_id)) stop('Duplicate input keys')
  cells<-cells[match(donors$donor_id,cells$donor_id),,drop=FALSE]
  if (anyNA(cells$donor_id)) stop('Donor cell counts not aligned')
  read_counts <- function(kind) {
    x<-as.matrix(read.csv(gzfile(file.path(src,paste0('counts_',kind,'.csv.gz'))),row.names=1,check.names=FALSE))
    storage.mode(x)<-'double'
    if (!identical(colnames(x),donors$donor_id) || !identical(rownames(x),genes$gene_id)) stop('Count input order mismatch')
    x
  }
  a<-read_counts('OPC');b<-read_counts('INH')
  members<-read.csv(file.path(go,'go_gene_membership.csv'),colClasses='character')
  dictionary<-read.csv(file.path(go,'human_gene_dictionary.csv'),colClasses='character')
  terms<-sort(unique(members$go_id)); if(length(terms)!=4826L) stop('Expected the full retained 4826-set GO family')
  bad<-duplicated(dictionary$gene_symbol)|duplicated(dictionary$gene_symbol,fromLast=TRUE)
  dictionary<-dictionary[!bad,,drop=FALSE]
  unique_symbol<-!(duplicated(genes$gene_symbol)|duplicated(genes$gene_symbol,fromLast=TRUE))
  mapping<-dictionary$entrez_id[match(genes$gene_symbol,dictionary$gene_symbol)];mapping[!unique_symbol]<-NA_character_
  mapping[duplicated(mapping)|duplicated(mapping,fromLast=TRUE)]<-NA_character_;names(mapping)<-genes$gene_id
  sets<-split(members$entrez_id,members$go_id)
  targets<-c('MTHFD2','OXTR','HOPX','CXCL10','IRF1','CD14','TAP1','PTPN6','TNFAIP2','ICAM1','HAMP','ITPKC','JAK3','GMFG','HCK','CEBPB')
  statuses<-list();targets_out<-list();go_out<-list()
  for (model in c('age_sex_batch_RIN','plus_PMI')) {
    cat('Paired OPC–INH interaction:',model,'\n')
    keep<-cells$OPC>=20 & cells$INH>=20 & donors$braak %in% c(1,2,5,6)
    needed<-c('age','sex_male','batch','rin');if(model=='plus_PMI') needed<-c(needed,'pmi')
    keep<-keep & complete.cases(donors[,needed,drop=FALSE])
    z<-fit_pair(a[,keep,drop=FALSE],b[,keep,drop=FALSE],donors[keep,,drop=FALSE],model=='plus_PMI')
    statuses[[model]]<-data.frame(pair='OPC_minus_INH',model=model,
      n_low=z$n_low,n_high=z$n_high,n_donors=z$n_donors,n_pseudobulks=z$n_pseudobulks,
      status=z$status,reason=z$reason,genes_tested=z$genes_tested,
      independent_design_df=z$independent_design_df,block_correlation=z$block_correlation,
      contrast='(V-VI minus I-II in OPC) minus (V-VI minus I-II in INH)')
    r<-data.frame(pair='OPC_minus_INH',model=model,gene_symbol=targets,status='not_estimable',
      interaction_logFC=NA_real_,CI_L=NA_real_,CI_R=NA_real_,PValue=NA_real_,FDR_genomewide=NA_real_)
    p<-data.frame(pair='OPC_minus_INH',model=model,go_id=terms,n_genes=0L,status='model_not_estimable',
      direction=NA_character_,PValue=NA_real_)
    if (z$status=='estimated') {
      for (i in seq_along(targets)) {
        ids<-which(genes$gene_symbol==targets[i])
        if(length(ids)==0) {r$status[i]<-'absent_from_input';next}
        if(length(ids)>1) {r$status[i]<-'ambiguous_gene_symbol';next}
        j<-match(genes$gene_id[ids],rownames(z$table))
        if(is.na(j)) {r$status[i]<-'filtered_low_expression';next}
        r$status[i]<-'tested';t<-z$table[j,]
        r$interaction_logFC[i]<-t$logFC;r$CI_L[i]<-t$CI.L;r$CI_R[i]<-t$CI.R
        r$PValue[i]<-t$P.Value;r$FDR_genomewide[i]<-t$adj.P.Val
      }
      rank<-z$t; gid<-mapping[names(rank)];ok<-!is.na(gid);rank<-rank[ok];names(rank)<-gid[ok]
      index<-lapply(sets,function(x)which(names(rank)%in%x));sizes<-lengths(index);valid<-sizes>=10 & sizes<=300
      p$n_genes<-sizes[match(p$go_id,names(sizes))];p$status<-'insufficient_measured_members'
      if(any(valid)) {
        cam<-limma::cameraPR(rank,index[valid],inter.gene.cor=.01,use.ranks=FALSE,sort=FALSE)
        ix<-match(p$go_id,rownames(cam));good<-!is.na(ix)
        p$status[good]<-'tested';p$direction[good]<-cam$Direction[ix[good]];p$PValue[good]<-cam$PValue[ix[good]]
      }
    }
    r$q_fixed_targets<-p.adjust(ifelse(is.na(r$PValue),1,r$PValue),'BH',n=16)
    p$q_fixed_GO_family<-p.adjust(ifelse(is.na(p$PValue),1,p$PValue),'BH',n=4826)
    targets_out[[model]]<-r;go_out[[model]]<-p
    write(do.call(rbind,statuses),'GSE174367_interaction_status')
    cat(model,z$status,'paired donors=',z$n_donors,'\n')
  }
  r<-do.call(rbind,targets_out);p<-do.call(rbind,go_out)
  r$q_all_interaction_models<-p.adjust(ifelse(is.na(r$PValue),1,r$PValue),'BH',n=32)
  p$q_all_interaction_models<-p.adjust(ifelse(is.na(p$PValue),1,p$PValue),'BH',n=9652)
  write(r,'GSE174367_interaction_targets');write(p,'GSE174367_interaction_GO')
  writeLines(capture.output(sessionInfo()),file.path(out,'GSE174367_interaction_R_sessionInfo.txt'))
}

if (sys.nframe()==0L) {
  args<-commandArgs(TRUE)
  if(length(args)!=3) stop('Expected pseudobulk_directory go_annotation_directory output_directory')
  run_interaction(args[1],args[2],args[3])
}
