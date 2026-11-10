#!/usr/bin/env Rscript
# Independently supplied, gene-level normalized expression and donor metadata.
options(stringsAsFactors=FALSE)
if(.Platform$OS.type=='windows' && !nzchar(Sys.setlocale('LC_CTYPE','.UTF8')))stop('A UTF-8-capable R runtime is required for non-ASCII paths')
args<-commandArgs(trailingOnly=TRUE)
i_args<-match('--args-file',args)
if(!is.na(i_args))args<-readLines(args[i_args+1L],encoding='UTF-8',warn=FALSE)
arg<-function(k,default=NULL){i<-match(k,args);if(is.na(i))default else args[[i+1L]]}
need<-function(k){v<-arg(k);if(is.null(v)||!nzchar(v))stop("Missing ",k);v}
for(p in c("limma","statmod"))if(!requireNamespace(p,quietly=TRUE))stop("Missing ",p)
kind<-need("--validation-kind")
if(!kind%in%c("late_Braak_vs_zero","related_AD_vs_control"))stop("Explicit validation-kind required")
source<-need("--source");independence<-need("--independence-evidence")
out<-need("--output-dir");if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))stop("Choose empty output")
m<-read.csv(need("--manifest"),check.names=FALSE)
tab<-read.csv(need("--expression"),check.names=FALSE)
if(names(tab)[1]!="entrez_id"||anyDuplicated(tab$entrez_id))stop("Expression first column must be unique Entrez IDs; combine probes before validation without using disease effects")
e<-as.matrix(tab[,-1]);storage.mode(e)<-"double";rownames(e)<-as.character(tab$entrez_id)
if(!all(c("sample_id","donor_key","group")%in%names(m))||anyDuplicated(m$sample_id)||!setequal(m$sample_id,colnames(e))||any(!is.finite(e)))stop("Invalid expression/manifest")
m<-m[match(colnames(e),m$sample_id),]
covariates<-strsplit(need("--covariates"),",",fixed=TRUE)[[1]]
if(!all(covariates%in%names(m)))stop("Missing requested confounders; do not silently fit reduced model")
for(v in c("group",covariates))if(any(vapply(split(m[[v]],m$donor_key),function(x)length(unique(x[!is.na(x)]))>1L,logical(1)))&&v%in%c("group","age_at_death","sex_male"))stop("Inconsistent donor metadata")
complete<-complete.cases(m[,c("group",covariates),drop=FALSE]);excluded<-sum(!complete)
m<-m[complete,,drop=FALSE];e<-e[,complete,drop=FALSE]
m$group<-factor(m$group,levels=c("reference","case"))
if(anyNA(m$group)||length(unique(m$group))!=2L)stop("group must be reference/case")
nd<-table(m$group[!duplicated(m$donor_key)])
if(any(nd<5L))stop("At least five independent donors per contrast group required")
covariates<-covariates[vapply(covariates,function(x)length(unique(m[[x]]))>1L,logical(1))]
x<-model.matrix(as.formula(paste("~0+group",if(length(covariates))paste0("+",paste(covariates,collapse="+"))else"")),m)
if(qr(x)$rank!=ncol(x)||nrow(x)<=ncol(x)+5L)stop("Nonestimable external design")
repeated<-anyDuplicated(m$donor_key)>0L
rho<-if(repeated)limma::duplicateCorrelation(e,x,block=m$donor_key)$consensus else NA_real_
f<-if(repeated)limma::lmFit(e,x,block=m$donor_key,correlation=rho)else limma::lmFit(e,x)
c<-matrix(0,ncol(x),1,dimnames=list(colnames(x),"case_vs_reference"));c["groupcase",1]<-1;c["groupreference",1]<- -1
f<-limma::eBayes(limma::contrasts.fit(f,c),robust=TRUE)
t<-limma::topTable(f,number=Inf,sort.by="none",confint=TRUE)
all<-data.frame(entrez_id=rownames(t),log2FC=t$logFC,ci_low=t$CI.L,ci_high=t$CI.R,moderated_t=t$t,p=t$P.Value,q_genome=p.adjust(t$P.Value,"BH"))
targets<-read.csv(need("--targets"),colClasses="character")
if(!all(c("entrez_id","gene_symbol","discovery_direction")%in%names(targets))||anyDuplicated(targets$entrez_id))stop("Targets need unique entrez_id, gene_symbol, discovery_direction (up/down)")
if(!all(targets$discovery_direction%in%c("up","down")))stop("Invalid discovery direction")
result<-merge(targets,all,by="entrez_id",all.x=TRUE)
result$available<-!is.na(result$p)
result$q_target_family<-p.adjust(result$p,"BH",n=nrow(targets))
result$direction_agrees<-ifelse(result$available,sign(result$log2FC)==ifelse(result$discovery_direction=="up",1,-1),NA)
result$replication_label<-if(kind=="late_Braak_vs_zero")"Independent late-Braak contrast, subject to supplied cohort-independence evidence"else"Related AD/control association; NOT replication of Braak V-VI versus 0"
dir.create(out,recursive=TRUE)
write.csv(result,file.path(out,"external_target_results.csv"),row.names=FALSE,na="")
write.csv(all,file.path(out,"external_all_gene_results.csv"),row.names=FALSE)
write.csv(data.frame(planned_targets=nrow(result),available=sum(result$available),unavailable=sum(!result$available),
                    same_direction=sum(result$direction_agrees,na.rm=TRUE),opposite_direction=sum(!result$direction_agrees,na.rm=TRUE),
                    same_direction_q05=sum(result$direction_agrees & result$q_target_family<.05,na.rm=TRUE),
                    opposite_direction_q05=sum(!result$direction_agrees & result$q_target_family<.05,na.rm=TRUE),
                    validation_kind=kind,correction="BH over all prespecified supplied targets, including unmeasured targets in family size"),
          file.path(out,"external_validation_summary.csv"),row.names=FALSE)
go<-arg('--go-dir')
if(!is.null(go)) {
  memberships<-read.csv(file.path(go,'go_gene_membership.csv'),colClasses='character')
  sets<-lapply(split(memberships$entrez_id,memberships$go_id),function(g)intersect(g,all$entrez_id))
  sets<-sets[lengths(sets)>=10L & lengths(sets)<=300L]
  if(!length(sets))stop('No GO sets supported by external measured genes')
  z<-sign(all$moderated_t)*qnorm(pmax(all$p,.Machine$double.xmin)/2,lower.tail=FALSE);names(z)<-all$entrez_id
  index<-lapply(sets,function(g)match(g,names(z)))
  pathways<-limma::cameraPR(z,index,inter.gene.cor=.01,use.ranks=FALSE,sort=FALSE)
  pathway_table<-data.frame(go_id=rownames(pathways),n_genes=pathways$NGenes,direction=pathways$Direction,p=pathways$PValue,q=p.adjust(pathways$PValue,'BH'),validation_kind=kind)
  discovery_file<-arg('--discovery-pathways')
  if(!is.null(discovery_file)) {
    discovery<-read.csv(discovery_file)
    discovery<-discovery[discovery$model=='full',c('go_id','direction','q')]
    names(discovery)[2:3]<-c('discovery_direction','discovery_q')
    pathway_table<-merge(pathway_table,discovery,by='go_id',all.x=TRUE)
    pathway_table$direction_agrees<-pathway_table$direction==pathway_table$discovery_direction
  }
  write.csv(pathway_table,file.path(out,'external_pathway_results.csv'),row.names=FALSE)
  writeLines('Current discovery GO sets intersected with externally measured genes, 10-300 genes; BH across all retained sets. Disease contrast, tissue and platform may differ. Direction agreement is not proof of exact late-Braak replication.',file.path(out,'external_pathway_interpretation.txt'))
}
write.csv(data.frame(source=source,independence_evidence=independence,validation_kind=kind,n_reference=unname(nd["reference"]),n_case=unname(nd["case"]),n_samples=nrow(m),excluded_missing_covariates=excluded,covariates=paste(covariates,collapse=";")),file.path(out,"external_analysis_specification.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,"R_session_info.txt"))
cat("External analysis completed; check comparison compatibility before claiming replication.\n")
