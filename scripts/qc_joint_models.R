#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
args<-commandArgs(TRUE)
arg<-function(k,default=NULL){i<-match(k,args);if(is.na(i))default else args[i+1L]}
expected<-c(samples=as.integer(arg('--expected-samples',212)),features=as.integer(arg('--expected-features',22011)),genes=as.integer(arg('--expected-genes',14950)),GO=as.integer(arg('--expected-go',4826)))
for(p in c('limma','statmod','dtangle'))if(!requireNamespace(p,quietly=TRUE))stop('Missing package: ',p)

fit_background<-function(Y,d,extra=character()){
 if(ncol(Y)!=nrow(d)||any(!is.finite(Y)))stop('Expression/design mismatch or nonfinite expression')
 d<-droplevels(d)
 cov<-c('brain_region','age_at_death','sex_male','pmi_hours','apoe4_carrier',extra)
 if(anyNA(d[,c('braak_stage','donor_key',cov)]))stop('Missing model covariates')
 cov<-cov[vapply(cov,function(k)length(unique(d[[k]]))>1L,logical(1))]
 if(!all(c('0','V-VI')%in%d$braak_stage))stop('Braak comparison groups absent')
 if(length(unique(d$donor_key))<3L||!anyDuplicated(d$donor_key))stop('Insufficient repeated donor support')
 X<-model.matrix(reformulate(c('braak_stage',cov),intercept=FALSE),d)
 if(qr(X)$rank<ncol(X)||nrow(X)<=ncol(X)+5L)stop('Rank-deficient or under-supported design')
 rho<-limma::duplicateCorrelation(Y,X,block=d$donor_key)$consensus
 if(!is.finite(rho))stop('Nonfinite donor correlation')
 f<-limma::lmFit(Y,X,block=d$donor_key,correlation=rho)
 contrast<-matrix(0,ncol(X),1)
 contrast[match('braak_stageV-VI',colnames(X)),1]<-1
 contrast[match('braak_stage0',colnames(X)),1]<- -1
 tt<-limma::topTable(limma::eBayes(limma::contrasts.fit(f,contrast),robust=TRUE),number=Inf,sort.by='none',confint=TRUE)
 if(any(!is.finite(tt$P.Value)))stop('Nonfinite full-background P values')
 list(table=data.frame(id=rownames(tt),log2FC=tt$logFC,ci_low=tt$CI.L,ci_high=tt$CI.R,
      t=tt$t,p=tt$P.Value,q=p.adjust(tt$P.Value,'BH')),
      diagnostic=data.frame(n_samples=ncol(Y),n_donors=length(unique(d$donor_key)),
      n_features=nrow(Y),n_parameters=ncol(X),residual_df=nrow(X)-ncol(X),
      donor_correlation=rho,condition_number=kappa(X)))
}


root<-arg('--project-root');qc<-arg('--qc-results');refdir<-arg('--reference-dir')
go<-arg('--go-dir');out<-arg('--output-dir');metadata<-arg('--sample-metadata')
if(is.null(out)||dir.exists(out)&&length(list.files(out)))stop('Choose an empty output directory')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write<-function(x,name)write.csv(x,file.path(out,paste0(name,'.csv')),row.names=FALSE,na='')
status<-list();progress<-function(component,state,reason=''){
 status[[component]]<<-data.frame(component=component,status=state,reason=reason,updated=format(Sys.time(),'%Y-%m-%d %H:%M:%S'))
 write(do.call(rbind,status),'analysis_status')
 cat(format(Sys.time(),'%H:%M:%S'),component,state,reason,'\n');flush.console()
}
progress('inputs','running')
e<-read.csv(gzfile(file.path(qc,'rma_core_expression.csv.gz')),check.names=FALSE)
Y<-as.matrix(e[,-1]);storage.mode(Y)<-'double';rownames(Y)<-as.character(e[[1]])
rm(e);gc()
m<-read.csv(file.path(qc,'retained_sample_manifest.csv'))
exclusions<-read.csv(file.path(qc,'technical_exclusions.csv'))
original<-read.csv(file.path(root,'data/public/gse131617/feature_annotation.csv'),colClasses='character')
original_manifest<-read.csv(file.path(root,'data/public/gse131617/sample_manifest.csv'))
stopifnot(!anyDuplicated(rownames(Y)),!anyDuplicated(colnames(Y)),!anyDuplicated(m$sample_id),
 setequal(colnames(Y),m$sample_id),all(m$sample_id%in%original_manifest$sample_id),all(is.finite(Y)))
m<-m[match(colnames(Y),m$sample_id),]
if(nrow(m)!=expected['samples'])stop('Unexpected retained array count: ',nrow(m),' expected ',expected['samples'])
if(any(m$sample_id%in%exclusions$sample_id))stop('An excluded technical-failure array is present')
if(!setequal(setdiff(original_manifest$sample_id,exclusions$sample_id),m$sample_id))stop('Retained samples do not equal original minus documented technical exclusions')
Y<-Y[rownames(Y)%in%original$feature_id,,drop=FALSE]
if(nrow(Y)!=expected['features']||!setequal(rownames(Y),original$feature_id))stop('Full original feature background was not preserved')
a<-read.csv(file.path(qc,'fixed_gene_representatives.csv'),colClasses='character')
stopifnot(!anyDuplicated(a$entrez_id),!anyDuplicated(a$feature_id),all(a$feature_id%in%rownames(Y)))
G<-Y[a$feature_id,,drop=FALSE];rownames(G)<-a$entrez_id
if(nrow(G)!=expected['genes'])stop('Unexpected fixed gene universe: ',nrow(G))
targets<-read.csv(file.path(root,'data/contracts/gene_symbols.csv'),colClasses='character')
if(!all(targets$feature_id%in%rownames(Y)))stop('A fixed target feature is absent')
meta<-read.csv(metadata);stopifnot(!anyDuplicated(meta$sample_id),all(m$sample_id%in%meta$sample_id))
meta<-meta[match(m$sample_id,meta$sample_id),]
if(!'array_scan_date'%in%names(meta)||anyNA(meta$array_scan_date)||any(trimws(meta$array_scan_date)==''))stop('Complete scan dates required')
m$array_scan_date<-factor(meta$array_scan_date)
m$braak_stage<-factor(m$braak_stage,levels=c('0','I-II','III-IV','V-VI'));m$brain_region<-factor(m$brain_region)
members<-read.csv(file.path(go,'go_gene_membership.csv'),colClasses='character')
sets0<-lapply(split(members$entrez_id,members$go_id),unique)
sets<-lapply(sets0,intersect,y=rownames(G));size<-lengths(sets)
write(data.frame(go_id=names(sets),original_members=lengths(sets0),measurable_members=size,tested=size>=10&size<=300),'GO_coverage')
sets<-sets[size>=10&size<=300];if(!length(sets))stop('No GO sets with 10-300 measurable genes')
prior_go<-read.csv(file.path(qc,'GO_coverage.csv'),colClasses=c(go_id='character'))
tested_flag<-tolower(trimws(as.character(prior_go$tested)))
if(!all(tested_flag%in%c('true','false','1','0')))stop('Invalid prior GO tested flag')
if(length(sets)!=expected['GO']||!setequal(names(sets),prior_go$go_id[tested_flag%in%c('true','1')]))stop('GO test family differs from the earlier QC run')
write(a,'fixed_gene_representatives');write(m,'model_sample_manifest')
write(data.frame(features=nrow(Y),genes=nrow(G),samples=nrow(m),donors=length(unique(m$donor_key)),GO_sets=length(sets)),'dimensions')
pathways<-function(tt,model){
 z<-sign(tt$t)*qnorm(pmax(tt$p,.Machine$double.xmin)/2,lower.tail=FALSE);names(z)<-tt$id
 idx<-lapply(sets,function(s)match(s,names(z)))
 if(any(vapply(idx,anyNA,logical(1))))stop('GO background mismatch')
 ct<-limma::cameraPR(z,idx,inter.gene.cor=.01,use.ranks=FALSE,sort=FALSE)
 data.frame(go_id=rownames(ct),model=model,n_genes=ct$NGenes,direction=ct$Direction,p=ct$PValue,q=p.adjust(ct$PValue,'BH'))
}
progress('inputs','completed')
progress('baseline_fit','running')
baseline_feature<-fit_background(Y,m)
progress('baseline_fit','completed')
models<-list(baseline=character(),scan_date='array_scan_date')
markers_out<-list();fractions_diag<-list()
for(reference in c('CA','NG')){
 ref<-read.csv(file.path(refdir,paste0('reference_',reference,'.csv')),check.names=FALSE,colClasses=c(entrez_id='character'))
 stopifnot(!anyDuplicated(ref$entrez_id));rownames(ref)<-ref$entrez_id;ref<-as.matrix(ref[,-1]);storage.mode(ref)<-'double'
 common<-setdiff(intersect(rownames(G),rownames(ref)),targets$entrez_id)
 if(length(common)<1000||any(!is.finite(ref))||any(ref<0))stop('Invalid cell reference')
 R0<-t(log2(ref[common,,drop=FALSE]+.5));B0<-t(G[common,,drop=FALSE])
 markers<-dtangle::find_markers(Y=R0,references=R0,marker_method='diff',data_type='rna-seq')
 selected<-lapply(markers$L,function(x)head(x,100L))
 if(any(lengths(selected)<10))stop('Insufficient reference markers')
 for(j in seq_along(selected))markers_out[[length(markers_out)+1]]<-data.frame(reference=reference,cell_type=rownames(R0)[j],entrez_id=colnames(R0)[selected[[j]]])
 for(mode in c('joint_quantile','native_logscale')){
  label<-paste(reference,mode,sep='_');progress(paste0('fractions_',label),'running')
  if(mode=='joint_quantile'){
   joined<-t(limma::normalizeBetweenArrays(t(rbind(B0,R0)),method='quantile'))
   B<-joined[seq_len(nrow(B0)),,drop=FALSE];R<-joined[nrow(B0)+seq_len(nrow(R0)),,drop=FALSE]
  }else{B<-B0;R<-R0}
  ans<-dtangle::dtangle(Y=B,references=R,markers=selected,n_markers=lengths(selected),gamma=.7)
  props<-ans$estimates;colnames(props)<-rownames(R0)
  stopifnot(all(is.finite(props)),all(props>=0),all(props<=1),max(abs(rowSums(props)-1))<1e-6,'Neurons'%in%colnames(props))
  f<-data.frame(sample_id=m$sample_id,props,check.names=FALSE);write(f,paste0('estimated_fractions_',label))
  keep<-setdiff(colnames(props),'Neurons');new<-paste0(label,'_',make.names(keep));m[,new]<-props[,keep,drop=FALSE]
  models[[paste0('composition_',label)]]<-new;models[[paste0('joint_',label)]]<-c('array_scan_date',new)
  fractions_diag[[label]]<-data.frame(reference=reference,mode=mode,n_reference_genes=length(common),gamma=.7,max_fraction=max(props),min_fraction=min(props))
  progress(paste0('fractions_',label),'completed')
 }
}
write(do.call(rbind,markers_out),'reference_markers');write(do.call(rbind,fractions_diag),'deconvolution_diagnostics')
selected_results<-list();go_results<-list();diagnostics<-list()
for(label in names(models)){
 progress(label,'running')
 tryCatch({
  f<-if(label=='baseline')baseline_feature else fit_background(Y,m,models[[label]]);write(f$table,paste0('features_',label))
  g<-fit_background(G,m,models[[label]]);write(g$table,paste0('genes_',label))
  selected_results[[label]]<-merge(transform(f$table,feature_id=id,model=label),targets,by='feature_id')
  go_results[[label]]<-pathways(g$table,label);write(go_results[[label]],paste0('GO_',label))
  diagnostics[[paste0(label,'_features')]]<-cbind(model=label,background='features',f$diagnostic)
  diagnostics[[paste0(label,'_genes')]]<-cbind(model=label,background='genes',g$diagnostic)
  write(do.call(rbind,selected_results),'target_effects_full_background');write(do.call(rbind,go_results),'GO_all_models')
  progress(label,'completed')
 },error=function(e)progress(label,'not_estimable',conditionMessage(e)))
}
loro<-list()
for(label in c('baseline','scan_date','joint_CA_joint_quantile','joint_NG_joint_quantile')){
 for(region in levels(m$brain_region)){
  key<-paste('LORO',label,region,sep='_');progress(key,'running')
  tryCatch({
   keep<-m$brain_region!=region
   g<-fit_background(G[,keep,drop=FALSE],m[keep,],models[[label]])
   loro[[key]]<-transform(pathways(g$table,label),excluded_region=region)
   write(loro[[key]],key);write(do.call(rbind,loro),'GO_LORO_all_models')
   progress(key,'completed')
  },error=function(e)progress(key,'not_estimable',conditionMessage(e)))
 }
}
if(length(diagnostics))write(do.call(rbind,diagnostics),'model_diagnostics')
if(length(go_results)){
 all<-do.call(rbind,go_results);baseline<-all[all$model=='baseline',]
 summaries<-lapply(names(go_results),function(k){z<-merge(baseline,go_results[[k]],by='go_id',suffixes=c('_baseline','_model'));data.frame(model=k,tested=nrow(z),significant=sum(z$q_model<.05),baseline_significant=sum(z$q_baseline<.05),baseline_same_direction_significant=sum(z$q_baseline<.05&z$q_model<.05&z$direction_baseline==z$direction_model),baseline_opposite_direction_significant=sum(z$q_baseline<.05&z$q_model<.05&z$direction_baseline!=z$direction_model))})
 write(do.call(rbind,summaries),'GO_model_summary')
}
write(as.data.frame(table(m$array_scan_date,m$braak_stage,m$brain_region)),'scan_date_braak_region_counts')
capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'))
writeLines(c('QC-retained RMA arrays only; same sample and fixed original gene representatives across models.',
 'CA/NG fractions re-estimated on the QC-passed matrix; original targets excluded from reference marker selection.',
 'Correlation and empirical Bayes moderation estimated on full feature/gene backgrounds, before target extraction.',
 'Feature BH and gene BH families are separate. Every eligible GO set is tested per model and region-omission fit.',
 'LORO retains full-sample estimated fractions and the fixed gene/GO universe; fraction-estimation uncertainty is not propagated.',
 'Scan date is a technical proxy, not RIN or laboratory batch. Estimated fractions are not measured cell counts.',
 'No donor/scan-date leave-one-out expansion or new CEL download is performed.'),file.path(out,'methods.txt'))
st<-do.call(rbind,status);progress('workflow',if(any(st$status!='completed'))'completed_with_nonestimable_models' else 'completed')
if(any(st$status!='completed'))quit(status=3)
