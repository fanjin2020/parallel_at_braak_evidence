#!/usr/bin/env Rscript
# Reference-based composition sensitivity; estimated fractions are not measured cell counts.
options(stringsAsFactors=FALSE)
if(.Platform$OS.type=='windows')Sys.setlocale('LC_CTYPE','.UTF8')
args<-commandArgs(trailingOnly=TRUE);i<-match('--args-file',args)
if(!is.na(i))args<-readLines(args[i+1L],encoding='UTF-8',warn=FALSE)
arg<-function(k,default=NULL){i<-match(k,args);if(is.na(i))default else args[i+1L]}
for(p in c('dtangle','limma','statmod'))if(!requireNamespace(p,quietly=TRUE))stop('Missing package: ',p)
root<-arg('--project-root');refdir<-arg('--reference-dir');out<-arg('--output-dir')
if(dir.exists(out)&&length(list.files(out)))stop('Choose an empty output directory')
dir.create(out,recursive=TRUE)
base<-file.path(root,'data/public/gse131617')
m<-read.csv(file.path(base,'sample_manifest.csv'),check.names=FALSE)
e<-read.csv(gzfile(file.path(base,'expression_feature.csv.gz')),check.names=FALSE)
ann<-read.csv(file.path(base,'feature_annotation.csv'),colClasses='character')
targets<-read.csv(file.path(root,'data/contracts/gene_symbols.csv'),colClasses='character')
y<-as.matrix(e[,-1]);rownames(y)<-e$feature_id;storage.mode(y)<-'double'
m<-m[match(colnames(y),m$sample_id),]
if(anyNA(m$sample_id)||any(!is.finite(y)))stop('Invalid expression or alignment')
m$braak_stage<-factor(m$braak_stage,levels=c('0','I-II','III-IV','V-VI'))
m$brain_region<-factor(m$brain_region)
a<-ann[ann$mapping_status=='unique_entrez_id' & ann$feature_id%in%rownames(y),]
a$average<-rowMeans(y[a$feature_id,,drop=FALSE]);a<-a[order(a$entrez_id,-a$average,a$feature_id),];a<-a[!duplicated(a$entrez_id),]
gene<-y[a$feature_id,,drop=FALSE];rownames(gene)<-a$entrez_id
write<-function(x,n)write.csv(x,file.path(out,paste0(n,'.csv')),row.names=FALSE,na='')
fit<-function(Y,extra=NULL){
 d<-if(is.null(extra))m else cbind(m,extra);cols<-if(is.null(extra))character() else names(extra)
 X<-model.matrix(as.formula(paste('~0+braak_stage+brain_region+age_at_death+sex_male+pmi_hours+apoe4_carrier',paste(cols,collapse='+'),sep=if(length(cols))'+' else '')),d)
 if(qr(X)$rank!=ncol(X))stop('Composition-adjusted design is rank deficient')
 rho<-limma::duplicateCorrelation(Y,X,block=d$donor_key)$consensus
 f<-limma::lmFit(Y,X,block=d$donor_key,correlation=rho)
 z<-matrix(0,ncol(X),1);z[which(colnames(X)=='braak_stageV-VI'),1]<-1;z[which(colnames(X)=='braak_stage0'),1]<- -1
 tt<-limma::topTable(limma::eBayes(limma::contrasts.fit(f,z),robust=TRUE),number=Inf,sort.by='none',confint=TRUE)
 data.frame(feature_id=rownames(tt),log2FC=tt$logFC,ci_low=tt$CI.L,ci_high=tt$CI.R,t=tt$t,p=tt$P.Value,q=p.adjust(tt$P.Value,'BH'),donor_correlation=rho,n_donors=length(unique(d$donor_key)),n_samples=nrow(d),design_condition_number=kappa(X))
}
target_ids<-targets$feature_id
effects<-list();summaries<-list();markers_out<-list();diagnostics<-list();status<-list()
full_reference<-fit(y)
write(full_reference,'features_unadjusted_composition')
effects[['unadjusted_composition']]<-transform(full_reference[match(target_ids,full_reference$feature_id),],model='unadjusted_composition',family='BH over all features; targets extracted after full-background fit')
for(reference in c('CA','NG')) {
 ref<-read.csv(file.path(refdir,paste0('reference_',reference,'.csv')),check.names=FALSE,colClasses=c(entrez_id='character'))
 rownames(ref)<-ref$entrez_id;ref<-as.matrix(ref[,-1]);storage.mode(ref)<-'double'
 common<-intersect(rownames(gene),rownames(ref));common<-setdiff(common,targets$entrez_id)
 if(length(common)<1000)stop('Insufficient reference coverage')
 r<-t(log2(ref[common,,drop=FALSE]+.5));bulk<-t(gene[common,,drop=FALSE])
 # Reference-only marker ranking excludes the original target genes to reduce direct circularity.
 markers<-dtangle::find_markers(Y=r,references=r,marker_method='diff',data_type='rna-seq')
 selected<-lapply(markers$L,function(x)head(x,100L))
 if(any(lengths(selected)<10))stop('Insufficient markers for a reference cell type')
 for(j in seq_along(selected))markers_out[[length(markers_out)+1L]]<-data.frame(reference=reference,cell_type=rownames(r)[j],entrez_id=colnames(r)[selected[[j]]])
 for(mode in c('joint_quantile','native_logscale')) {
  label<-paste(reference,mode,sep='_');cat('Deconvolving',label,'\n');flush.console()
  tryCatch({
   if(mode=='joint_quantile'){
    joined<-t(limma::normalizeBetweenArrays(t(rbind(bulk,r)),method='quantile'))
    B<-joined[seq_len(nrow(bulk)),,drop=FALSE];R<-joined[nrow(bulk)+seq_len(nrow(r)),,drop=FALSE]
   } else {B<-bulk;R<-r}
   ans<-dtangle::dtangle(Y=B,references=R,markers=selected,n_markers=lengths(selected),gamma=.7)
   props<-ans$estimates;colnames(props)<-rownames(r)
   if(any(!is.finite(props))||any(props<0)||max(abs(rowSums(props)-1))>1e-6)stop('Invalid fraction estimates')
   fractions<-data.frame(sample_id=m$sample_id,props,check.names=FALSE);names(fractions)[-1]<-paste0('cell_prop_',make.names(colnames(props)))
   write(fractions,paste0('estimated_fractions_',label))
   for(s in levels(m$braak_stage))for(cell in colnames(props)){
    v<-props[m$braak_stage==s,cell]
    summaries[[length(summaries)+1L]]<-data.frame(model=label,braak_stage=s,cell_type=cell,n=length(v),mean=mean(v),sd=sd(v),minimum=min(v),maximum=max(v))
   }
   diagnostics[[label]]<-data.frame(model=label,n_common_genes=length(common),n_target_genes_excluded=length(intersect(targets$entrez_id,intersect(rownames(gene),rownames(ref)))),min_fraction=min(props),max_fraction=max(props),samples_fraction_above_095=sum(apply(props,1,max)>.95),reference_condition_number=kappa(t(R)),gamma=.7)
   extra<-as.data.frame(props[,setdiff(colnames(props),'Neurons'),drop=FALSE]);names(extra)<-paste0('cell_prop_',make.names(names(extra)))
   full_adjusted<-fit(y,extra);write(full_adjusted,paste0('features_',label))
   effects[[label]]<-transform(full_adjusted[match(target_ids,full_adjusted$feature_id),],model=label,family='BH over all features; targets extracted after full-background fit')
   cat('Fitting genome-wide composition-adjusted model\n');flush.console()
   genome<-fit(gene,extra);names(genome)[1]<-'entrez_id';genome$model<-paste0('cell_composition_',label)
   write(genome,paste0('gene_level_ranks_',label))
   if(label=='CA_joint_quantile'){
    genome$model<-'cell_composition_CA';write(genome,'gene_level_ranks_cell_adjusted')
   }
   status[[label]]<-data.frame(component=label,status='ok',reason='Cross-platform exploratory sensitivity; fractions and their uncertainty are model-dependent')
  },error=function(e)status[[label]]<<-data.frame(component=label,status='not_estimable',reason=conditionMessage(e)))
 }
}
all<-merge(do.call(rbind,effects),targets,by='feature_id')
write(all,'selected_gene_composition_sensitivity');write(all[all$gene_symbol!='UTY',],'non_UTY_composition_sensitivity');write(all[all$gene_symbol=='UTY',],'UTY_composition_sensitivity')
write(do.call(rbind,summaries),'fraction_summary');write(do.call(rbind,markers_out),'reference_markers');write(do.call(rbind,diagnostics),'deconvolution_diagnostics');write(do.call(rbind,status),'analysis_status')
file.copy(file.path(refdir,'reference_provenance.json'),out)
writeLines(c('Exploratory snRNA-reference deconvolution of RMA microarray data.','Six broad cell classes; Neurons is the omitted composition covariate.','References: Sutton 2022 Data 5 CA (Hodge 2019) and NG (Nagy 2020).','Top 100 reference-only diff markers per class; original 17 targets excluded from markers.','Primary sensitivity: joint quantile normalization of log2 bulk and log2(reference + 0.5), dtangle gamma 0.7.','Native-logscale and alternate-reference results are reported even when inconsistent.','These are relative expression-based composition estimates, not validated histological cell fractions.','Nominal intervals condition on estimated composition and original gene selection; cross-platform, reference-region, RNA-content and nuclear-capture biases remain.','Composition can be a disease consequence; attenuation after adjustment does not establish confounding or mediation.'),file.path(out,'analysis_specification.txt'))
capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'))
