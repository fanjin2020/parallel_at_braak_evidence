#!/usr/bin/env Rscript
# Re-test a current GO family and describe gene-overlap redundancy.
options(stringsAsFactors=FALSE)
if(.Platform$OS.type=='windows' && !nzchar(Sys.setlocale('LC_CTYPE','.UTF8')))stop('A UTF-8-capable R runtime is required for non-ASCII paths')
args<-commandArgs(trailingOnly=TRUE)
i_args<-match('--args-file',args)
if(!is.na(i_args))args<-readLines(args[i_args+1L],encoding='UTF-8',warn=FALSE)
arg<-function(k,default=NULL){i<-match(k,args);if(is.na(i))default else args[[i+1L]]}
if(!requireNamespace("limma",quietly=TRUE))stop("limma required")
input<-arg("--transcriptomic-results");go<-arg("--go-dir");out<-arg("--output-dir")
root<-arg("--project-root")
if(is.null(root)||!dir.exists(file.path(root,"data/public/gse131617")))
  stop("--project-root with rebuilt discovery inputs is required for pathway-region analyses")
if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))stop("Choose empty output directory")
dir.create(out,recursive=TRUE)
ranks<-read.csv(file.path(input,"gene_level_ranks.csv"),colClasses=c(entrez_id="character"))
cell<-arg('--cell-results')
if(!is.null(cell)) {
 cf<-file.path(cell,'gene_level_ranks_cell_adjusted.csv')
 if(!file.exists(cf))stop('Cell-adjusted gene ranks unavailable; do not assume a completed adjustment')
 extra<-read.csv(cf,colClasses=c(entrez_id='character'))
 ranks<-rbind(ranks,extra[,names(ranks),drop=FALSE])
}
members<-read.csv(file.path(go,"go_gene_membership.csv"),colClasses="character")
names_go<-read.csv(file.path(go,"go_terms.csv"),colClasses=c(go_id="character"))
sets<-split(members$entrez_id,members$go_id)
rows<-list()
for(model in unique(ranks$model)) {
  r<-ranks[ranks$model==model,];if(anyDuplicated(r$entrez_id))stop("Duplicate gene ranks")
  z<-sign(r$t)*qnorm(pmax(r$p,.Machine$double.xmin)/2,lower.tail=FALSE);names(z)<-r$entrez_id
  index<-lapply(sets,function(g)match(g,names(z)))
  if(any(vapply(index,anyNA,logical(1))))stop("Rank universe differs from GO preparation universe")
  c<-limma::cameraPR(z,index,inter.gene.cor=.01,use.ranks=FALSE,sort=FALSE)
  rows[[model]]<-data.frame(go_id=rownames(c),model=model,direction=c$Direction,p=c$PValue,q=p.adjust(c$PValue,"BH"))
}
all<-merge(do.call(rbind,rows),names_go,by="go_id",all.x=TRUE)
write.csv(all,file.path(out,"current_GO_all_tests.csv"),row.names=FALSE)
full<-all[all$model=="full",];full<-full[order(full$q,full$p,full$go_id),]
significant<-full[full$q<.05,]
# Greedy representative groups by Jaccard >=0.5; a display reduction, not a new hypothesis test.
group<-integer(nrow(significant));representatives<-integer();links<-list()
if(nrow(significant)) for(i in seq_len(nrow(significant))) {
  a<-sets[[significant$go_id[i]]];assigned<-FALSE
  for(j in representatives) {
    b<-sets[[significant$go_id[j]]];jac<-length(intersect(a,b))/length(union(a,b))
    if(jac>=.5) {
      group[i]<-group[j];assigned<-TRUE
      links[[length(links)+1L]]<-data.frame(member=significant$go_id[i],representative=significant$go_id[j],jaccard=jac)
      break
    }
  }
  if(!assigned){group[i]<-length(representatives)+1L;representatives<-c(representatives,i)}
}
significant$overlap_group<-group
significant$display_representative<-seq_len(nrow(significant))%in%representatives
write.csv(significant,file.path(out,"current_GO_significant_grouped.csv"),row.names=FALSE)
write.csv(if(length(links))do.call(rbind,links) else data.frame(member=character(),representative=character(),jaccard=numeric()),file.path(out,"GO_overlap_links.csv"),row.names=FALSE)
write.csv(data.frame(tested_sets=length(sets),significant_full=nrow(significant),display_groups=length(representatives),correlation=.01,
                    correction="BH over all current eligible BP sets, separately by sensitivity model",reduction="Jaccard >=0.5 to first lowest-q representative; significance and complete table unchanged"),file.path(out,"pathway_specification.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,"R_session_info.txt"))
cat("Current GO analysis complete:",length(sets),"tested sets\n")

if(!is.null(root)) {
  base<-file.path(root,"data/public/gse131617")
  m<-read.csv(file.path(base,"sample_manifest.csv"),check.names=FALSE)
  ex<-read.csv(gzfile(file.path(base,"expression_feature.csv.gz")),check.names=FALSE)
  representatives<-read.csv(file.path(input,"gene_representatives.csv"),colClasses=c(feature_id="character",entrez_id="character"))
  e<-as.matrix(ex[match(representatives$feature_id,ex$feature_id),-1,drop=FALSE]);storage.mode(e)<-"double";rownames(e)<-representatives$entrez_id
  m<-m[match(colnames(e),m$sample_id),]
  sd_gene<-apply(e,1,sd)
  if(any(!is.finite(sd_gene))||any(sd_gene<=0))stop("Invalid gene SD for pathway-score sensitivity")
  z<-(e-rowMeans(e))/sd_gene
  scores<-t(vapply(sets,function(g)colMeans(z[g,,drop=FALSE]),numeric(ncol(e))))
  rownames(scores)<-names(sets)
  m$cell<-factor(paste(m$brain_region,m$braak_stage,sep="__"))
  x<-model.matrix(~0+cell+age_at_death+sex_male+pmi_hours+apoe4_carrier,m)
  colnames(x)<-make.names(colnames(x))
  if(qr(x)$rank!=ncol(x))stop("Pathway interaction design rank deficient")
  rho<-limma::duplicateCorrelation(scores,x,block=m$donor_key)$consensus
  if(!is.finite(rho))stop('Nonfinite donor correlation')
  f<-limma::lmFit(scores,x,block=m$donor_key,correlation=rho)
  regions<-c("EC","FC","TC")
  c<-sapply(regions,function(r){v<-setNames(rep(0,ncol(x)),colnames(x));v[make.names(paste0("cell",r,"__V-VI"))]<-1;v[make.names(paste0("cell",r,"__0"))]<- -1;v})
  regional<-limma::eBayes(limma::contrasts.fit(f,c),robust=TRUE)
  output<-list()
  for(j in seq_along(regions)) {
    t<-limma::topTable(regional,coef=j,number=Inf,sort.by="none",confint=TRUE)
    output[[j]]<-data.frame(go_id=rownames(t),region=regions[j],score_difference=t$logFC,ci_low=t$CI.L,ci_high=t$CI.R,p=t$P.Value,q=p.adjust(t$P.Value,"BH"))
  }
  h<-limma::eBayes(limma::contrasts.fit(f,cbind(c[,1]-c[,2],c[,1]-c[,3])),robust=TRUE)
  ht<-limma::topTable(h,coef=1:2,number=Inf,sort.by="none")
  write.csv(merge(do.call(rbind,output),names_go,by="go_id"),file.path(out,"pathway_score_region_effects.csv"),row.names=FALSE)
  write.csv(data.frame(go_id=rownames(ht),F=ht$F,p=ht$P.Value,q=p.adjust(ht$P.Value,"BH")),file.path(out,"pathway_score_region_interactions.csv"),row.names=FALSE)
  writeLines("Exploratory mean standardized-gene pathway scores, not camera enrichment effect sizes or activation states. Gene standardization uses the complete discovery expression matrix. BH across all eligible pathways, separately per region and interaction family. Donor-blocked limma.",file.path(out,"pathway_score_interpretation.txt"))
}
