#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
a<-commandArgs(TRUE);i<-match('--args-file',a);if(!is.na(i))a<-readLines(a[i+1],warn=FALSE,encoding='UTF-8')
arg<-function(k,default=NULL){j<-match(k,a);if(is.na(j))default else a[j+1]}
need<-function(k){x<-arg(k);if(is.null(x))stop('Missing ',k);x}
for(p in c('limma','statmod'))if(!requireNamespace(p,quietly=TRUE))stop('Missing ',p)
input<-need('--input-dir');out<-need('--output-dir');dir.create(out,recursive=TRUE,showWarnings=FALSE)
m<-read.csv(file.path(input,'manifest.csv'),check.names=FALSE)
v<-read.csv(gzfile(file.path(input,'expression_gene.csv.gz')),check.names=FALSE)
y<-as.matrix(v[,-1]);storage.mode(y)<-'double';rownames(y)<-v[[1]]
stopifnot(setequal(colnames(y),m$sample_id),!anyDuplicated(m$donor_key),!anyDuplicated(rownames(y)),all(is.finite(y)))
m<-m[match(colnames(y),m$sample_id),]
targets<-read.csv(need('--targets'),colClasses='character')
stopifnot(all(c('entrez_id','gene_symbol','discovery_direction')%in%names(targets)),!anyDuplicated(targets$entrez_id),
 all(targets$discovery_direction%in%c('up','down')))
# UTY is a separate sensitivity, consistent with the revised non-Y target family.
targets<-targets[targets$gene_symbol!='UTY',]
models<-list(age_sex_APOE=c('age_at_death','sex_male','apoe4_count'),
 age_sex_APOE_RNAquality=c('age_at_death','sex_male','apoe4_count','rna_quality_category'))
results<-list();statuses<-list();specs<-list();pathways<-list()
go<-arg('--go-dir');membership<-if(!is.null(go))read.csv(file.path(go,'go_gene_membership.csv'),colClasses='character')else NULL
for(label in names(models)){
 cat('External model:',label,'\n');flush.console()
 tryCatch({
  cov<-models[[label]];ok<-complete.cases(m[,c('braak_group',cov)])
  d<-m[ok,];Y<-y[,ok,drop=FALSE]
  d$braak_group<-factor(d$braak_group,levels=c('0','I-II','III-IV','V-VI'))
  d$rna_quality_category<-factor(d$rna_quality_category)
  if(any(table(d$braak_group)[c('0','V-VI')]<5))stop('Fewer than five donors in contrast group')
  varying<-cov[vapply(cov,function(k)length(unique(d[[k]]))>1,logical(1))]
  X<-model.matrix(reformulate(c('braak_group',varying),intercept=FALSE),d)
  if(qr(X)$rank<ncol(X)||nrow(X)<=ncol(X)+5)stop('Nonestimable external design')
  fit<-limma::lmFit(Y,X);ct<-matrix(0,ncol(X),1)
  ct[match('braak_groupV-VI',colnames(X)),1]<-1;ct[match('braak_group0',colnames(X)),1]<- -1
  fit<-limma::eBayes(limma::contrasts.fit(fit,ct),robust=TRUE)
  z<-limma::topTable(fit,number=Inf,sort.by='none',confint=TRUE)
  all<-data.frame(entrez_id=rownames(z),log2FC=z$logFC,ci_low=z$CI.L,ci_high=z$CI.R,t=z$t,p=z$P.Value,q_genome=p.adjust(z$P.Value,'BH'))
  write.csv(all,file.path(out,paste0('genes_',label,'.csv')),row.names=FALSE)
  r<-merge(targets,all,by='entrez_id',all.x=TRUE);r$model<-label
  r$q_target_family<-p.adjust(r$p,'BH',n=nrow(targets))
  r$direction_agrees<-sign(r$log2FC)==ifelse(r$discovery_direction=='up',1,-1)
  r$available<-!is.na(r$p);results[[label]]<-r
  specs[[label]]<-data.frame(model=label,n_samples=nrow(d),n_zero=sum(d$braak_group=='0'),n_late=sum(d$braak_group=='V-VI'),
      n_genes=nrow(Y),excluded_missing_covariates=sum(!ok),covariates=paste(varying,collapse=';'),
      interpretation='Independent Braak contrast in inferior temporal cortex; not validation of all three discovery regions or ADNI')
  if(!is.null(membership)){
   sets<-lapply(split(membership$entrez_id,membership$go_id),function(g)intersect(g,all$entrez_id));sets<-sets[lengths(sets)>=10&lengths(sets)<=300]
   ranks<-sign(all$t)*qnorm(pmax(all$p,.Machine$double.xmin)/2,lower.tail=FALSE);names(ranks)<-all$entrez_id
   g<-limma::cameraPR(ranks,lapply(sets,function(ids)match(ids,names(ranks))),inter.gene.cor=.01,use.ranks=FALSE,sort=FALSE)
   pathways[[label]]<-data.frame(model=label,go_id=rownames(g),n_genes=g$NGenes,direction=g$Direction,p=g$PValue,q=p.adjust(g$PValue,'BH'))
  }
  statuses[[label]]<-data.frame(component=label,status='completed',reason='')
 },error=function(e){statuses[[label]]<<-data.frame(component=label,status='not_estimable',reason=conditionMessage(e))})
}
write.csv(do.call(rbind,statuses),file.path(out,'analysis_status.csv'),row.names=FALSE)
if(!length(results))stop('No external model was estimable')
r<-do.call(rbind,results)
# Report per-model target families and an additional family across both planned models.
r$q_targets_across_models<-p.adjust(r$p,'BH',n=nrow(targets)*length(models))
write.csv(r,file.path(out,'external_target_results.csv'),row.names=FALSE)
write.csv(do.call(rbind,specs),file.path(out,'model_specifications.csv'),row.names=FALSE)
if(length(pathways))write.csv(do.call(rbind,pathways),file.path(out,'external_GO_all_models.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'))
writeLines(c('External cohort and normalization were selected without viewing target effects.',
 'All four Braak groups retained; V-VI versus 0 contrast. Age, sex and APOE count included; published RNA-quality category is an additional sensitivity, not numeric RIN.',
 'Full mapped-gene background fitted before target extraction. The small Braak 0 group limits precision.',
 'Missing targets retained in planned BH denominator. No favorable-result filtering.'),file.path(out,'methods.txt'))
