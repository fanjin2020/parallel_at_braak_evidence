#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
if(.Platform$OS.type=='windows') Sys.setlocale('LC_CTYPE','.UTF8')
args <- commandArgs(TRUE)
arg <- function(k,default=NULL){i<-match(k,args);if(is.na(i))default else args[i+1L]}
filearg <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
scriptdir <- dirname(normalizePath(filearg,mustWork=TRUE))
source(file.path(scriptdir,'tissue_helpers.R'))
source(file.path(scriptdir,'sex_interaction_helpers.R'))
for(p in c('limma','statmod'))if(!requireNamespace(p,quietly=TRUE))stop('Missing ',p)
root <- arg('--project-root');out <- arg('--output-dir');diagnostic <- FALSE
if(is.null(root)||is.null(out))stop('Specify --project-root and --output-dir')
root <- normalizePath(root,mustWork=TRUE,winslash='/')
qc <- file.path(root,'cel_qc');joint <- file.path(root,'targeted/qc_joint')
inputs <- c(expression=file.path(qc,'rma_core_expression.csv.gz'),
 manifest=file.path(qc,'retained_sample_manifest.csv'),exclusions=file.path(qc,'technical_exclusions.csv'),
 representatives=file.path(qc,'fixed_gene_representatives.csv'),targets=file.path(root,'data/contracts/gene_symbols.csv'),
 technical=file.path(root,'data/public/cel_metadata/sample_technical_metadata.csv'),
 membership=file.path(root,'revision/go_annotation/go_gene_membership.csv'),dimensions=file.path(joint,'dimensions.csv'),
 prior=file.path(joint,'features_baseline.csv'),CA=file.path(joint,'estimated_fractions_CA_joint_quantile.csv'),
 NG=file.path(joint,'estimated_fractions_NG_joint_quantile.csv'),
 annotation=file.path(root,'data/public/gse131617/feature_annotation.csv'))
if(any(!file.exists(inputs)))stop('Required cached inputs missing: ',paste(names(inputs)[!file.exists(inputs)],collapse=', '))
signature <- list(inputs=setNames(as.character(tools::md5sum(inputs)),names(inputs)),
 code=unname(tools::md5sum(c(filearg,file.path(scriptdir,c('sex_interaction_helpers.R','tissue_helpers.R'))))),diagnostic=diagnostic)
dir.create(out,recursive=TRUE,showWarnings=FALSE);cp<-file.path(out,'checkpoints');dir.create(cp,showWarnings=FALSE)
sf <- file.path(out,'resume_signature.rds')
if(file.exists(sf)){if(!identical(readRDS(sf),signature))stop('Inputs/settings/code changed; use a new output directory')
}else saveRDS(signature,sf)
writecsv <- function(x,name)write.csv(x,file.path(out,paste0(name,'.csv')),row.names=FALSE,na='')
e <- read.csv(gzfile(inputs['expression']),check.names=FALSE)
Y <- as.matrix(e[,-1]);storage.mode(Y)<-'double';rownames(Y)<-as.character(e[[1]]);rm(e)
prior <- read.csv(inputs['prior'],colClasses=c(id='character'))
if(anyDuplicated(prior$id)||!all(prior$id%in%rownames(Y)))stop('Invalid background identifiers')
Y <- Y[match(prior$id,rownames(Y)),,drop=FALSE]
m <- read.csv(inputs['manifest']);if(anyDuplicated(m$sample_id)||!setequal(m$sample_id,colnames(Y)))stop('Retained sample mismatch')
m <- m[match(colnames(Y),m$sample_id),]
dim_expected <- read.csv(inputs['dimensions']);exclusions <- read.csv(inputs['exclusions'])
stopifnot(ncol(Y)==dim_expected$samples,nrow(Y)==dim_expected$features,
 length(unique(m$donor_key))==dim_expected$donors,!any(m$sample_id%in%exclusions$sample_id),
 ncol(Y)==212L,length(unique(m$donor_key))==71L)
meta <- read.csv(inputs['technical']);if(anyDuplicated(meta$sample_id))stop('Duplicate technical metadata')
m$array_scan_date <- factor(meta$array_scan_date[match(m$sample_id,meta$sample_id)])
models <- list(baseline=character(),scan_date='array_scan_date')
for(ref in c('CA','NG')){
 f <- read.csv(inputs[ref],check.names=FALSE)
 if(anyDuplicated(f$sample_id)||!setequal(f$sample_id,m$sample_id))stop('Fraction/sample mismatch')
 f <- f[match(m$sample_id,f$sample_id),];p <- as.matrix(f[,-1]);storage.mode(p)<-'double'
 stopifnot(all(is.finite(p)),all(p>=0),all(p<=1),max(abs(rowSums(p)-1))<1e-6,'Neurons'%in%names(f))
 keys <- setdiff(names(f),c('sample_id','Neurons'));cols <- paste0(ref,'_',make.names(keys));m[,cols] <- f[,keys]
 models[[paste0('joint_',ref,'_joint_quantile')]] <- c('array_scan_date',cols)
}
a <- read.csv(inputs['representatives'],colClasses='character');targets <- read.csv(inputs['targets'],colClasses='character')
stopifnot(!anyDuplicated(a$entrez_id),!anyDuplicated(a$feature_id),all(a$feature_id%in%rownames(Y)),
 !anyDuplicated(targets$feature_id),all(targets$feature_id%in%rownames(Y)),sum(targets$gene_symbol!='UTY')==16L,
 sum(targets$gene_symbol=='UTY')==1L)
G <- Y[a$feature_id,,drop=FALSE];rownames(G) <- a$entrez_id
membership <- read.csv(inputs['membership'],colClasses='character')
sets <- lapply(split(membership$entrez_id,membership$go_id),function(x)intersect(unique(x),rownames(G)))
sets <- sets[lengths(sets)>=10&lengths(sets)<=300]
stopifnot(length(sets)==dim_expected$GO_sets,nrow(G)==dim_expected$genes)
S <- review_scores(G,sets)
donors <- m[!duplicated(m$donor_key),]
coverage <- aggregate(rep(1L,nrow(donors)),donors[c('sex_male','braak_stage')],sum);names(coverage)[3]<-'donors'
arrays <- aggregate(rep(1L,nrow(m)),m[c('sex_male','braak_stage')],sum);names(arrays)[3]<-'arrays'
writecsv(merge(coverage,arrays,by=c('sex_male','braak_stage')),'cohort_coverage')
annotation <- read.csv(inputs['annotation'],colClasses='character')
if(anyDuplicated(annotation$feature_id))stop('Duplicate feature annotation')
features <- paths <- diagnostics <- statuses <- list();start <- Sys.time();job <- 0L
for(model in names(models))for(background in c('feature','GO_score')){
 job <- job+1L;key <- paste(model,background,sep='__');file <- file.path(cp,paste0(key,'.rds'))
 cat(sprintf('[Braak x sex %d/%d] %s %s\n',job,2L*length(models),model,background));flush.console()
 if(file.exists(file)){ans<-readRDS(file);cat('  cached\n')}else{
  ans <- tryCatch({z<-sex_interaction_fit(if(background=='feature')Y else S,m,models[[model]])
   list(status='completed',reason='',result=z)},
   not_estimable=function(e)list(status='not_estimable',reason=conditionMessage(e),result=NULL))
  tmp<-paste0(file,'.tmp');saveRDS(ans,tmp);if(!file.rename(tmp,file))stop('Checkpoint rename failed')
 }
 statuses[[key]] <- data.frame(model=model,background=background,status=ans$status,reason=ans$reason)
 if(ans$status=='completed'){
  z <- cbind(model=model,ans$result$effects)
  if(background=='feature'){
   z$entrez_id <- annotation$entrez_id[match(z$id,annotation$feature_id)]
   z$mapping_status <- annotation$mapping_status[match(z$id,annotation$feature_id)]
   z$target_symbol <- targets$gene_symbol[match(z$id,targets$feature_id)]
   z$interpretation_scope <- ifelse(!is.na(z$target_symbol)&z$target_symbol=='UTY',
     'Y-linked target: not interpreted as a sex interaction','Genome-wide exploratory interaction')
   features[[model]] <- z
  }else paths[[model]] <- z
  diagnostics[[key]] <- cbind(model=model,background=background,ans$result$diagnostic)
 }
 writecsv(do.call(rbind,statuses),'analysis_status')
 writecsv(data.frame(completed_jobs=job,planned_jobs=2L*length(models),diagnostic_only=diagnostic,
  elapsed_seconds=as.numeric(difftime(Sys.time(),start,units='secs'))),'progress')
}
bind <- function(x)if(length(x))do.call(rbind,x)else data.frame(status=character())
f <- bind(features);g <- bind(paths)
writecsv(f,'feature_sex_interactions');writecsv(g,'GO_sex_interactions');writecsv(bind(diagnostics),'model_diagnostics')
t <- if(nrow(f))f[!is.na(f$target_symbol)&f$target_symbol!='UTY',]else f
writecsv(t,'target_sex_interactions')
summaries <- list()
for(model in names(models))for(background in c('all_features','non_Y_targets','GO_scores')){
 z <- switch(background,all_features=f,non_Y_targets=t,GO_scores=g)
 if(!'model'%in%names(z))next
 z <- z[z$model==model,,drop=FALSE];if(!nrow(z))next
 summaries[[length(summaries)+1L]] <- data.frame(model=model,background=background,n=nrow(z),
  late_interaction_q05=sum(z$interaction_late_q<.05),omnibus_interaction_q05=sum(z$omnibus_q<.05),
  min_late_q=min(z$interaction_late_q),min_omnibus_q=min(z$omnibus_q))
}
writecsv(bind(summaries),'interaction_summary')
capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'))
writeLines(c('QC-retained RMA arrays only: 212 arrays, 71 donors.',
 'The sex-interaction analyses were exploratory and were specified after the initial analyses.',
 'Four-level Braak by binary recorded sex cell means, region, age, PMI and APOE4; optional scan date and fixed CA/NG fractions.',
 'Primary interaction contrast: (V-VI minus 0 in males) minus (V-VI minus 0 in females), one degree of freedom.',
 'Secondary omnibus test: three sex-difference contrasts for I-II, III-IV and V-VI versus 0, moderated F with three degrees of freedom.',
 'Donor-blocked duplicateCorrelation, lmFit and robust eBayes fitted to the full background before target extraction.',
 'Separate BH families per model and test: all 22,011 features or all 4,826 GO scores. Model sensitivities are not independent confirmations.',
 'Sex-specific contrasts share the interaction model and its covariate slopes; they need not equal separately fitted sex-subgroup estimates.',
 'GO scores average standardized expression in the fixed 14,950-gene universe; this is not CAMERA enrichment or cellular activation.',
 'The 16 non-Y target results retain full-feature q values. UTY remains in the full-background family but is excluded from target sex-interaction interpretation.',
 'Nonsignificant interaction tests do not establish equivalence. Sparse sex-by-stage cells limit precision.',
 'Missing discovery RIN, true experimental batch and cause of death are not replaced by scan date.'),file.path(out,'methods.txt'))
cat('SEX_INTERACTION_COMPLETE; inspect analysis_status.csv before interpretation\n')
