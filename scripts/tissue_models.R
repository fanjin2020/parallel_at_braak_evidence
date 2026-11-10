#!/usr/bin/env Rscript
# Donor-aware sensitivity models and regional late-Braak contrasts.
options(stringsAsFactors=FALSE)
if(.Platform$OS.type=='windows' && !nzchar(Sys.setlocale('LC_CTYPE','.UTF8')))stop('A UTF-8-capable R runtime is required for non-ASCII paths')
args <- commandArgs(trailingOnly=TRUE)
i_args<-match('--args-file',args)
if(!is.na(i_args))args<-readLines(args[i_args+1L],encoding='UTF-8',warn=FALSE)
arg <- function(k, default=NULL) { i<-match(k,args); if(is.na(i)) default else args[[i+1L]] }
for(p in c("limma","statmod")) if(!requireNamespace(p,quietly=TRUE)) stop("Missing R package: ",p)
root <- normalizePath(arg("--project-root"),mustWork=TRUE)
out <- arg("--output-dir")
if(is.null(out) || (dir.exists(out) && length(list.files(out,all.files=TRUE,no..=TRUE)))) stop("Choose a new empty output directory")
dir.create(out,recursive=TRUE)
base <- file.path(root,"data/public/gse131617")
manifest <- read.csv(file.path(base,"sample_manifest.csv"),check.names=FALSE)
tab <- read.csv(gzfile(file.path(base,"expression_feature.csv.gz")),check.names=FALSE)
ann <- read.csv(file.path(base,"feature_annotation.csv"),colClasses="character")
selected <- read.csv(file.path(root,"data/contracts/gene_symbols.csv"),colClasses="character")
if(anyDuplicated(tab$feature_id)||anyDuplicated(manifest$sample_id)||anyDuplicated(ann$feature_id)) stop("Duplicate input key")
e <- as.matrix(tab[,-1L]); storage.mode(e)<-"double";rownames(e)<-tab$feature_id
if(!setequal(colnames(e),manifest$sample_id)||any(!is.finite(e))) stop("Invalid expression or sample alignment")
manifest <- manifest[match(colnames(e),manifest$sample_id),]
required <- c("donor_key","braak_stage","brain_region","age_at_death","sex_male","pmi_hours","apoe4_carrier")
if(!all(required %in% names(manifest))||anyNA(manifest[,required])) stop("Incomplete required manifest")
for(v in c("braak_stage","age_at_death","sex_male","pmi_hours","apoe4_carrier")) {
  if(any(vapply(split(manifest[[v]],manifest$donor_key),function(x) length(unique(x))!=1L,logical(1)))) stop("Donor metadata inconsistent: ",v)
}
manifest$braak_stage<-factor(manifest$braak_stage,levels=c("0","I-II","III-IV","V-VI"))
manifest$brain_region<-factor(manifest$brain_region,levels=c("EC","FC","TC"))
if(anyNA(manifest[,c('braak_stage','brain_region')]))stop('Unknown Braak stage or brain-region label')
status<-list();results<-list();diagnostics<-list()
record <- function(component,state,reason="") status[[length(status)+1L]] <<- data.frame(component=component,status=state,reason=reason)
write <- function(x,name) write.csv(x,file.path(out,paste0(name,".csv")),row.names=FALSE,na="")
donors<-manifest[!duplicated(manifest$donor_key),]
summary_rows<-list()
for(s in levels(manifest$braak_stage)) {
  d<-donors[donors$braak_stage==s,]
  for(v in c("age_at_death","sex_male","pmi_hours","apoe4_carrier")) {
    x<-d[[v]]
    summary_rows[[length(summary_rows)+1L]]<-data.frame(braak_stage=s,variable=v,n_donors=nrow(d),n_available=sum(!is.na(x)),mean=mean(x),sd=sd(x),q25=unname(quantile(x,.25)),median=median(x),q75=unname(quantile(x,.75)),minimum=min(x),maximum=max(x))
  }
}
write(do.call(rbind,summary_rows),"donor_characteristics")
donor_values<-list()
for(s in levels(manifest$braak_stage))for(v in c('age_at_death','sex_male','pmi_hours','apoe4_carrier')) {
  counts<-as.data.frame(table(donors[donors$braak_stage==s,v],useNA='ifany'))
  donor_values[[length(donor_values)+1L]]<-data.frame(braak_stage=s,variable=v,value=as.character(counts[[1]]),n_donors=counts[[2]])
}
write(do.call(rbind,donor_values),'donor_value_counts')
write(as.data.frame(table(manifest$braak_stage,manifest$brain_region)),"region_sample_counts")
metadata_file<-arg("--sample-metadata")
extra<-character()
availability<-data.frame(variable=c("rin","batch","tissue_ph","cause_of_death","array_processing_date","array_scan_date","scanner_id","cell_proportions"),available=FALSE,n_available=0L)
if(!is.null(metadata_file)) {
  meta<-read.csv(metadata_file,check.names=FALSE)
  if(!"sample_id" %in% names(meta)||anyDuplicated(meta$sample_id)||!all(manifest$sample_id %in% meta$sample_id)) stop("Supplementary metadata requires unique matching sample_id")
  meta<-meta[match(manifest$sample_id,meta$sample_id),]
  for(v in c("rin","batch","tissue_ph","cause_of_death","array_processing_date","array_scan_date","scanner_id")) if(v %in% names(meta)) {
    manifest[[v]]<-meta[[v]]
    if(v %in% c("rin","tissue_ph")) manifest[[v]]<-as.numeric(manifest[[v]]) else manifest[[v]]<-factor(manifest[[v]])
    availability[availability$variable==v,c("available","n_available")]<-list(TRUE,sum(!is.na(meta[[v]])))
    extra<-c(extra,v)
  }
  fractions<-grep("^cell_prop_",names(meta),value=TRUE)
  if(length(fractions)>=2L) {
    f<-as.matrix(meta[,fractions]);storage.mode(f)<-"double"
    if(anyNA(f)||any(f<0|f>1)||any(abs(rowSums(f)-1)>.01)) stop("Cell fractions must be complete, nonnegative and sum to 1")
    manifest[,fractions]<-f
    extra<-c(extra,sort(fractions)[-1L])
    availability[availability$variable=="cell_proportions",c("available","n_available")]<-list(TRUE,nrow(f))
    record("cell_composition","supplied_fractions_only","Externally estimated, validated fractions required; no deconvolution or functional-state inference performed here")
  }
}
write(availability,"technical_metadata_availability")
for(v in availability$variable[!availability$available]) record(v,"requires_source_metadata_or_reference","Not in supplied manifest; not imputed from Braak labels or gene effects")
if(length(extra)) {
  technical_values<-list()
  for(s in levels(manifest$braak_stage))for(v in extra) {
    counts<-as.data.frame(table(manifest[manifest$braak_stage==s,v],useNA='ifany'))
    technical_values[[length(technical_values)+1L]]<-data.frame(braak_stage=s,variable=v,value=as.character(counts[[1]]),n_tissue_samples=counts[[2]])
  }
  write(do.call(rbind,technical_values),'technical_value_counts')
}
if('--descriptions-only'%in%args) {
  write(do.call(rbind,status),'analysis_status')
  capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'))
  quit(save='no',status=0)
}

fit_model <- function(y,m,interaction=FALSE,extras=character(),age_spline=FALSE) {
  m<-droplevels(m)
  if(!all(c("0","V-VI") %in% levels(m$braak_stage))) stop("Both contrast groups required")
  nd<-table(m$braak_stage[!duplicated(m$donor_key)])
  if(any(nd[c("0","V-VI")]<3L)) stop("Fewer than 3 donors in a contrast group")
  cov<-c(if(age_spline) "splines::ns(age_at_death,df=3)" else "age_at_death","sex_male","pmi_hours","apoe4_carrier",extras)
  cov<-cov[vapply(cov,function(v) if(grepl("splines::",v)) TRUE else length(unique(m[[v]]))>1L,logical(1))]
  if(interaction) { m$cell<-factor(paste(m$brain_region,m$braak_stage,sep="__")); lead<-"cell" }
  else lead<-c("braak_stage",if(nlevels(m$brain_region)>1L) "brain_region")
  design<-model.matrix(as.formula(paste("~0+",paste(c(lead,cov),collapse="+"))),m)
  colnames(design)<-make.names(colnames(design))
  if(nrow(design)!=ncol(y)||qr(design)$rank!=ncol(design)||nrow(design)<=ncol(design)+5L) stop("Incomplete or nonestimable design")
  repeated<-anyDuplicated(m$donor_key)>0L
  rho<-if(repeated) limma::duplicateCorrelation(y,design,block=m$donor_key)$consensus else NA_real_
  f<-if(repeated) limma::lmFit(y,design,block=m$donor_key,correlation=rho) else limma::lmFit(y,design)
  if(!is.finite(rho)&&repeated) stop("Nonfinite donor correlation")
  if(interaction) {
    contrasts<-sapply(levels(m$brain_region),function(r) {
      z<-rep(0,ncol(design));names(z)<-colnames(design)
      hi<-make.names(paste0("cell",r,"__V-VI"));lo<-make.names(paste0("cell",r,"__0"))
      if(!all(c(hi,lo)%in%names(z))) stop("Missing region-stage cell")
      z[hi]<-1;z[lo]<- -1;z
    })
  } else {
    contrasts<-matrix(0,ncol(design),1,dimnames=list(colnames(design),"late_vs_zero"))
    contrasts["braak_stageV.VI",1]<-1;contrasts["braak_stage0",1]<- -1
  }
  moderated<-limma::eBayes(limma::contrasts.fit(f,contrasts),robust=TRUE)
  list(fit=moderated,base=f,contrasts=contrasts,rho=rho,n_donors=length(unique(m$donor_key)),n_samples=nrow(m),n_zero=unname(nd["0"]),n_late=unname(nd["V-VI"]))
}

extract <- function(f,label) {
  rows<-lapply(seq_len(ncol(f$fit$coefficients)),function(j) {
    tt<-limma::topTable(f$fit,coef=j,number=Inf,sort.by="none",confint=TRUE)
    data.frame(feature_id=rownames(tt),model=label,region=colnames(f$fit$coefficients)[j],log2FC=tt$logFC,ci_low=tt$CI.L,ci_high=tt$CI.R,
               moderated_t=tt$t,p=tt$P.Value,q_all_features=p.adjust(tt$P.Value,"BH"),n_donors=f$n_donors,n_samples=f$n_samples,n_zero=f$n_zero,n_late=f$n_late)
  })
  do.call(rbind,rows)
}

run_fit <- function(label,keep=rep(TRUE,nrow(manifest)),extras=character(),interaction=FALSE,age_spline=FALSE) {
  cat('Fitting',label,'\n');flush.console()
  needed<-c(required,extras);keep<-keep&complete.cases(manifest[,needed,drop=FALSE])
  tryCatch({
    f<-fit_model(e[,keep,drop=FALSE],manifest[keep,,drop=FALSE],interaction,extras,age_spline)
    r<-extract(f,label);results[[label]]<<-r
    diagnostics[[label]]<<-data.frame(model=label,n_donors=f$n_donors,n_samples=f$n_samples,n_zero=f$n_zero,n_late=f$n_late,donor_correlation=f$rho)
    write(r,paste0("features_",label));record(label,"ok");invisible(f)
  },error=function(x) {record(label,"not_estimable",conditionMessage(x));NULL})
}
models<-strsplit(arg('--feature-models','original_covariates,male_only,female_only,age_spline,age_overlap,technical,region_interaction'),',',fixed=TRUE)[[1]]
if(any(!models%in%c('original_covariates','male_only','female_only','age_spline','age_overlap','technical','region_interaction','none')))stop('Unknown feature model')
full<-if('original_covariates'%in%models)run_fit("original_covariates") else NULL
if('male_only'%in%models)run_fit("male_only",manifest$sex_male==1)
if('female_only'%in%models)run_fit("female_only",manifest$sex_male==0)
if('age_spline'%in%models)run_fit("age_spline",age_spline=TRUE)
age_ranges<-lapply(c("0","V-VI"),function(s) range(donors$age_at_death[donors$braak_stage==s]))
overlap<-c(max(vapply(age_ranges,`[`,numeric(1),1)),min(vapply(age_ranges,`[`,numeric(1),2)))
if('age_overlap'%in%models) {
  if(overlap[1]<overlap[2]) run_fit("age_overlap",manifest$age_at_death>=overlap[1]&manifest$age_at_death<=overlap[2]) else record("age_overlap","no_overlap")
}
if(length(extra)&&'technical'%in%models) {
  run_fit("technical_complete_case_reference",complete.cases(manifest[,extra,drop=FALSE]))
  run_fit("technical_adjusted",extras=extra)
}
regional<-if('region_interaction'%in%models)run_fit("region_interaction",interaction=TRUE) else NULL
if(!is.null(regional)) {
  c<-regional$contrasts
  hetero<-cbind(c[,1]-c[,2],c[,1]-c[,3])
  h<-limma::eBayes(limma::contrasts.fit(regional$base,hetero),robust=TRUE)
  tt<-limma::topTable(h,coef=1:2,number=Inf,sort.by="none")
  write(data.frame(feature_id=rownames(tt),F=tt$F,p=tt$P.Value,q_all_features=p.adjust(tt$P.Value,"BH")),"region_heterogeneity")
}
if(length(results)) {
all<-do.call(rbind,results)
sel<-merge(all,selected,by="feature_id",all=FALSE)
sel$q_selected_exploratory<-ave(sel$p,interaction(sel$model,sel$region),FUN=function(p)p.adjust(p,"BH"))
sel$inference<-"Original-data-selected targets; sensitivity estimates are not independent replication"
write(sel,"selected_gene_sensitivities")
write(sel[sel$gene_symbol!="UTY",],"non_UTY_gene_sensitivities")
write(sel[sel$gene_symbol=="UTY",],"UTY_separate_exploratory")
targets<-sel[sel$model=="original_covariates" & sel$gene_symbol!="UTY",c("entrez_id","gene_symbol","log2FC")]
targets$discovery_direction<-ifelse(targets$log2FC>0,"up","down")
write(targets[,c("entrez_id","gene_symbol","discovery_direction")],"external_validation_targets")
write(do.call(rbind,diagnostics),"model_sample_sizes")
}

# Select one uniquely mapped feature per Entrez ID by mean expression, not disease effect.
a<-ann[ann$mapping_status=="unique_entrez_id" & ann$feature_id %in% rownames(e),]
a$average<-rowMeans(e[a$feature_id,,drop=FALSE]);a<-a[order(a$entrez_id,-a$average,a$feature_id),];a<-a[!duplicated(a$entrez_id),]
gene<-e[a$feature_id,,drop=FALSE];rownames(gene)<-a$entrez_id
rank_rows<-list()
rank_models<-strsplit(arg('--rank-models','full,male,female'),',',fixed=TRUE)[[1]]
if(any(!rank_models%in%c('full','male','female','none')))stop('Unknown gene rank model')
for(label in setdiff(rank_models,'none')) {
  cat('Fitting gene-level ranks:',label,'\n');flush.console()
  keep<-if(label=="male") manifest$sex_male==1 else if(label=="female") manifest$sex_male==0 else rep(TRUE,nrow(manifest))
  tryCatch({
    f<-fit_model(gene[,keep,drop=FALSE],manifest[keep,,drop=FALSE]);tt<-limma::topTable(f$fit,number=Inf,sort.by="none")
    rank_rows[[label]]<-data.frame(entrez_id=rownames(tt),model=label,log2FC=tt$logFC,t=tt$t,p=tt$P.Value)
  },error=function(x)record(paste0("gene_rank_",label),"not_estimable",conditionMessage(x)))
}
if(length(rank_rows))write(do.call(rbind,rank_rows),"gene_level_ranks")
write(a[,c("feature_id","entrez_id","average")],"gene_representatives")
write(do.call(rbind,status),"analysis_status")
capture.output(sessionInfo(),file=file.path(out,"R_session_info.txt"))
cat("Transcriptomic sensitivities completed. Inspect analysis_status.csv.\n")
