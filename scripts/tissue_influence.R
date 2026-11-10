#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
if(.Platform$OS.type=='windows') Sys.setlocale('LC_CTYPE','.UTF8')
args<-commandArgs(TRUE)
arg<-function(k,default=NULL){i<-match(k,args);if(is.na(i))default else args[i+1L]}
filearg<-grep('^--file=',commandArgs(FALSE),value=TRUE)[1]
scriptdir<-dirname(normalizePath(sub('^--file=','',filearg),mustWork=TRUE))
source(file.path(scriptdir,'tissue_helpers.R'))
for(p in c('limma','statmod'))if(!requireNamespace(p,quietly=TRUE))stop('Missing ',p)
if('--self-test'%in%args){review_self_test();quit(status=0)}
root<-arg('--project-root');out<-arg('--output-dir');diagnostic<-'--diagnostic'%in%args
if(is.null(root)||is.null(out))stop('Specify --project-root and --output-dir')
qc<-file.path(root,'cel_qc');joint<-file.path(root,'targeted/qc_joint');go<-file.path(root,'revision/go_annotation')
inputs<-c(file.path(qc,c('rma_core_expression.csv.gz','retained_sample_manifest.csv','technical_exclusions.csv','fixed_gene_representatives.csv')),
 file.path(root,'data/contracts/gene_symbols.csv'),file.path(root,'data/public/cel_metadata/sample_technical_metadata.csv'),
 file.path(go,'go_gene_membership.csv'),file.path(joint,'dimensions.csv'),
 file.path(joint,paste0('features_',c('baseline','scan_date','joint_CA_joint_quantile','joint_NG_joint_quantile'),'.csv')),
 file.path(joint,paste0('estimated_fractions_',c('CA','NG'),'_joint_quantile.csv')))
if(any(!file.exists(inputs)))stop('Required rebuilt inputs missing: ',paste(basename(inputs[!file.exists(inputs)]),collapse=', '))
signature<-list(inputs=setNames(as.character(tools::md5sum(inputs)),sub(paste0(root,'/'),'',inputs,fixed=TRUE)),
 code=unname(tools::md5sum(c(sub('^--file=','',filearg),file.path(scriptdir,'tissue_helpers.R')))),diagnostic=diagnostic)
dir.create(out,recursive=TRUE,showWarnings=FALSE);cp<-file.path(out,'checkpoints');dir.create(cp,showWarnings=FALSE)
signature_file<-file.path(out,'resume_signature.rds')
if(file.exists(signature_file)){
 if(!identical(readRDS(signature_file),signature))stop('Inputs/settings/code changed; use a new output directory')
}else saveRDS(signature,signature_file)
writecsv<-function(x,name)write.csv(x,file.path(out,paste0(name,'.csv')),row.names=FALSE,na='')
e<-read.csv(gzfile(inputs[1]),check.names=FALSE);Y<-as.matrix(e[,-1]);storage.mode(Y)<-'double';rownames(Y)<-as.character(e[[1]]);rm(e)
m<-read.csv(inputs[2]);m<-m[match(colnames(Y),m$sample_id),]
a<-read.csv(inputs[4],colClasses='character');targets<-read.csv(inputs[5],colClasses='character')
prior<-read.csv(file.path(joint,'features_baseline.csv'),colClasses=c(id='character'))
Y<-Y[match(prior$id,rownames(Y)),,drop=FALSE]
dim_expected<-read.csv(file.path(joint,'dimensions.csv'))
stopifnot(!anyNA(m$sample_id),ncol(Y)==dim_expected$samples,nrow(Y)==dim_expected$features,
 !any(m$sample_id%in%read.csv(inputs[3])$sample_id),!anyDuplicated(a$entrez_id),!anyDuplicated(targets$feature_id))
meta<-read.csv(inputs[6]);if(anyDuplicated(meta$sample_id))stop('Duplicated technical metadata')
m$array_scan_date<-factor(meta$array_scan_date[match(m$sample_id,meta$sample_id)])
models<-list(baseline=character(),scan_date='array_scan_date')
for(ref in c('CA','NG')){
 f<-read.csv(file.path(joint,paste0('estimated_fractions_',ref,'_joint_quantile.csv')),check.names=FALSE)
 if(anyDuplicated(f$sample_id)||!setequal(m$sample_id,f$sample_id))stop('Fraction/sample mismatch')
 f<-f[match(m$sample_id,f$sample_id),];p<-as.matrix(f[,-1]);stopifnot(all(is.finite(p)),all(p>=0),all(p<=1),max(abs(rowSums(p)-1))<1e-6)
 keys<-setdiff(names(f),c('sample_id','Neurons'));cols<-paste0(ref,'_',make.names(keys));m[,cols]<-f[,keys]
 models[[paste0('joint_',ref,'_joint_quantile')]]<-c('array_scan_date',cols)
}
if(diagnostic)models<-models['baseline']
G<-Y[a$feature_id,,drop=FALSE];rownames(G)<-a$entrez_id
membership<-read.csv(inputs[7],colClasses='character')
sets<-lapply(split(membership$entrez_id,membership$go_id),function(x)intersect(unique(x),rownames(G)))
sets<-sets[lengths(sets)>=10&lengths(sets)<=300]
stopifnot(length(sets)==dim_expected$GO_sets,nrow(G)==dim_expected$genes)
S<-review_scores(G,sets)
donors<-sort(unique(m$donor_key));dates<-sort(unique(as.character(m$array_scan_date)))
if(diagnostic){donors<-head(donors,2);dates<-head(dates,2)}
age0<-range(m$age_at_death[m$braak_stage=='0']);age1<-range(m$age_at_death[m$braak_stage=='V-VI'])
overlap<-c(max(age0[1],age1[1]),min(age0[2],age1[2]))
writecsv(data.frame(group=c('0','V-VI','intersection'),lower=c(age0[1],age1[1],overlap[1]),upper=c(age0[2],age1[2],overlap[2])),'age_support')
coverage<-aggregate(rep(1,nrow(m)),m[c('braak_stage','brain_region','sex_male','array_scan_date')],sum);names(coverage)[5]<-'arrays';writecsv(coverage,'design_coverage')
jobs<-list()
addjob<-function(model,kind,value=''){
 key<-sprintf('%04d',length(jobs)+1L);jobs[[key]]<<-list(key=key,model=model,kind=kind,value=value)
}
for(model in names(models)){
 for(kind in c('reference','feature_region','GO_region','male','female','age_overlap','age_spline'))addjob(model,kind)
 for(d in donors)addjob(model,'omit_donor',d)
 for(d in dates)addjob(model,'omit_scan_date',d)
}
states<-list();results<-list();started<-Sys.time()
for(job in jobs){
 key<-job$key;file<-file.path(cp,paste0(key,'.rds'))
 cat(sprintf('[QC review %s/%d] %s %s %s\n',key,length(jobs),job$model,job$kind,job$value));flush.console()
 if(file.exists(file)){ans<-readRDS(file);cat('  cached\n')}else{
  ans<-tryCatch({
   keep<-switch(job$kind,male=m$sex_male==1,female=m$sex_male==0,
    age_overlap=m$age_at_death>=overlap[1]&m$age_at_death<=overlap[2],
    omit_donor=m$donor_key!=job$value,omit_scan_date=as.character(m$array_scan_date)!=job$value,rep(TRUE,nrow(m)))
   is_go<-job$kind=='GO_region';regional<-job$kind%in%c('feature_region','GO_region')
   yy<-if(is_go)S else Y
   z<-review_fit(yy[,keep,drop=FALSE],m[keep,,drop=FALSE],models[[job$model]],regional,job$kind=='age_spline')
   reproduction<-NULL
   if(job$kind=='reference'){
    old<-read.csv(file.path(joint,paste0('features_',job$model,'.csv')),colClasses=c(id='character'))
    old<-old[match(z$effects$id,old$id),];diff<-c(effect=max(abs(z$effects$effect-old$log2FC)),p=max(abs(z$effects$p-old$p)),q=max(abs(z$effects$q-old$q)))
    if(any(!is.finite(diff))||any(diff>1e-6))stop('QC reference reproduction mismatch')
    reproduction<-data.frame(statistic=names(diff),max_absolute_difference=unname(diff),passed=TRUE)
   }
   if(!is_go)z$effects<-review_targets(z$effects,targets)
   list(status='completed',reason='',result=z,reproduction=reproduction)
  },not_estimable=function(e)list(status='not_estimable',reason=conditionMessage(e),result=NULL))
  tmp<-paste0(file,'.tmp');saveRDS(ans,tmp);if(!file.rename(tmp,file))stop('Checkpoint rename failed')
 }
 states[[key]]<-data.frame(job=key,model=job$model,kind=job$kind,omitted_unit=job$value,status=ans$status,reason=ans$reason)
 results[[key]]<-ans
 writecsv(do.call(rbind,states),'analysis_status')
 writecsv(data.frame(completed_jobs=length(states),planned_jobs=length(jobs),elapsed_seconds=as.numeric(difftime(Sys.time(),started,units='secs')),diagnostic_only=diagnostic),'progress')
}
bind<-function(x)if(length(x))do.call(rbind,x)else data.frame(status=character())
effects<-interactions<-go_effects<-go_interactions<-diags<-checks<-list()
for(job in jobs){
 z<-results[[job$key]];if(z$status!='completed')next
 context<-data.frame(model=job$model,kind=job$kind,omitted_unit=job$value)
 diags[[job$key]]<-cbind(context,z$result$diagnostic)
 if(!is.null(z$reproduction))checks[[job$key]]<-cbind(model=job$model,z$reproduction)
 if(job$kind=='GO_region'){
   go_effects[[job$key]]<-cbind(context,z$result$effects)
   go_interactions[[job$key]]<-cbind(context,z$result$interaction)
 }else{
   effects[[job$key]]<-cbind(context,z$result$effects)
   if(!is.null(z$result$interaction))interactions[[job$key]]<-cbind(context,z$result$interaction)
 }
}
writecsv(bind(effects),'target_effects');writecsv(bind(interactions),'feature_interactions')
writecsv(bind(go_effects),'GO_region_effects');writecsv(bind(go_interactions),'GO_interactions')
writecsv(bind(diags),'model_diagnostics');writecsv(bind(checks),'reference_checks')
z<-bind(effects);summaries<-list()
for(model in names(models))for(unit in c('omit_donor','omit_scan_date'))for(id in targets$feature_id){
 ref<-z[z$model==model&z$kind=='reference'&z$feature_id==id,];v<-z[z$model==model&z$kind==unit&z$feature_id==id,]
 planned<-if(unit=='omit_donor')length(donors)else length(dates)
 if(nrow(ref)!=1)stop('Missing QC reference fit')
 summaries[[length(summaries)+1L]]<-data.frame(model=model,omission=unit,feature_id=id,gene_symbol=ref$gene_symbol,
  planned=planned,estimable=nrow(v),reference_effect=ref$effect,reference_q=ref$q,
  min_effect=if(nrow(v))min(v$effect)else NA_real_,max_effect=if(nrow(v))max(v$effect)else NA_real_,
  max_absolute_change=if(nrow(v))max(abs(v$effect-ref$effect))else NA_real_,
  direction_reversals=sum(sign(v$effect)!=sign(ref$effect)),q05_count=sum(v$q<.05),
  interpretation='Descriptive influence checks; no donor/date removed from the primary analysis')
}
writecsv(bind(summaries),'influence_summary')
capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'))
writeLines(c('All models use QC-retained RMA arrays; full feature or GO-score background fitted before target extraction.',
 'Moderated F tests compare the V-VI minus 0 effect across EC, FC and TC (two independent contrasts).',
 'BH over all features or all GO scores separately per model; regional effects have separate full-background BH per region.',
 'GO scores are averages of standardized expression, not CAMERA enrichment effects or cellular activation.',
 'Age overlap uses the intersection of the observed Braak 0/V-VI age ranges; spline age uses three degrees of freedom.',
 'Estimated cell fractions are held fixed during omissions; uncertainty from fraction estimation is not propagated.',
 'Nonestimable fits remain explicit. Absence of significance is not evidence of equal regional effects.',
 'Partial resumption is permitted only with identical input and code fingerprints.'),file.path(out,'methods.txt'))
cat('QC_REVIEW_COMPLETE; inspect analysis_status.csv before interpretation\n')
