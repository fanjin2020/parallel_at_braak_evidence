#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
if(.Platform$OS.type=='windows' && !nzchar(Sys.setlocale('LC_CTYPE','.UTF8')))stop('UTF-8 locale required')
args<-commandArgs(trailingOnly=TRUE)
i<-match('--args-file',args)
if(!is.na(i))args<-readLines(args[i+1L],encoding='UTF-8',warn=FALSE)
arg<-function(k){i<-match(k,args);if(is.na(i))stop('Missing ',k);args[i+1L]}
for(p in c('limma','statmod'))if(!requireNamespace(p,quietly=TRUE))stop('Missing ',p)
root<-arg('--project-root');input<-arg('--transcriptomic-results');go<-arg('--go-dir');out<-arg('--output-dir')
if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))stop('Choose an empty output directory')
dir.create(out,recursive=TRUE)
base<-file.path(root,'data/public/gse131617')
m<-read.csv(file.path(base,'sample_manifest.csv'),check.names=FALSE)
tab<-read.csv(gzfile(file.path(base,'expression_feature.csv.gz')),check.names=FALSE)
rep<-read.csv(file.path(input,'gene_representatives.csv'),colClasses=c(feature_id='character',entrez_id='character'))
if(anyDuplicated(rep$entrez_id)||anyDuplicated(rep$feature_id)||anyDuplicated(m$sample_id))stop('Duplicate key')
e<-as.matrix(tab[match(rep$feature_id,tab$feature_id),-1,drop=FALSE]);storage.mode(e)<-'double';rownames(e)<-rep$entrez_id
m<-m[match(colnames(e),m$sample_id),]
if(anyNA(m)||any(!is.finite(e)))stop('Incomplete expression or manifest')
m$braak_stage<-factor(m$braak_stage,levels=c('0','I-II','III-IV','V-VI'))
m$brain_region<-factor(m$brain_region,levels=c('EC','FC','TC'))
membership<-read.csv(file.path(go,'go_gene_membership.csv'),colClasses='character')
sets<-split(membership$entrez_id,membership$go_id)
index<-lapply(sets,function(x)match(x,rownames(e)))
if(any(vapply(index,anyNA,logical(1))))stop('GO universe mismatch')
results<-list();ranks<-list();diagnostics<-list()
for(region in c('EC','FC','TC')) {
 cat('Fitting current-GO leave-out:',region,'\n');flush.console()
 keep<-m$brain_region!=region;d<-droplevels(m[keep,]);y<-e[,keep,drop=FALSE]
 x<-model.matrix(~0+braak_stage+brain_region+age_at_death+sex_male+pmi_hours+apoe4_carrier,d)
 colnames(x)<-make.names(colnames(x))
 if(qr(x)$rank!=ncol(x))stop('Rank-deficient leave-out model')
 rho<-limma::duplicateCorrelation(y,x,block=d$donor_key)$consensus
 if(!is.finite(rho))stop('Nonfinite donor correlation')
 fit<-limma::lmFit(y,x,block=d$donor_key,correlation=rho)
 contrast<-matrix(0,ncol(x),1,dimnames=list(colnames(x),'late_vs_zero'))
 contrast['braak_stageV.VI',1]<-1;contrast['braak_stage0',1]<- -1
 fit<-limma::eBayes(limma::contrasts.fit(fit,contrast),robust=TRUE)
 tt<-limma::topTable(fit,number=Inf,sort.by='none')
 z<-sign(tt$t)*qnorm(pmax(tt$P.Value,.Machine$double.xmin)/2,lower.tail=FALSE)
 ca<-limma::cameraPR(z,index,inter.gene.cor=.01,use.ranks=FALSE,sort=FALSE)
 results[[region]]<-data.frame(go_id=rownames(ca),omitted_region=region,direction=ca$Direction,p=ca$PValue,q=p.adjust(ca$PValue,'BH'))
 ranks[[region]]<-data.frame(entrez_id=rownames(tt),omitted_region=region,log2FC=tt$logFC,t=tt$t,p=tt$P.Value)
 diagnostics[[region]]<-data.frame(omitted_region=region,n_donors=length(unique(d$donor_key)),n_samples=nrow(d),donor_correlation=rho,tested_sets=length(sets))
 write.csv(results[[region]],file.path(out,paste0('GO_leave_',region,'.csv')),row.names=FALSE)
}
full<-read.csv(file.path(arg('--pathway-results'),'current_GO_all_tests.csv'))
summary<-full[full$model=='full',c('go_id','name','n_genes','direction','p','q')]
for(region in names(results)) {
 r<-results[[region]];j<-match(summary$go_id,r$go_id)
 if(anyNA(j))stop('Incomplete GO leave-out family')
 summary[[paste0('direction_without_',region)]]<-r$direction[j]
 summary[[paste0('q_without_',region)]]<-r$q[j]
}
summary$all_leaveout_same_direction<-Reduce(`&`,lapply(names(results),function(r)summary$direction==summary[[paste0('direction_without_',r)]]))
summary$full_q05_direction_consistent<-summary$q<.05 & summary$all_leaveout_same_direction
summary$interpretation<-'Within-cohort directional sensitivity, not independent replication; each leave-out BH family contains all current eligible GO sets'
write.csv(summary,file.path(out,'current_GO_direction_consistency.csv'),row.names=FALSE)
write.csv(do.call(rbind,results),file.path(out,'current_GO_LORO_tests.csv'),row.names=FALSE)
write.csv(do.call(rbind,ranks),file.path(out,'gene_ranks_leave_region_out.csv'),row.names=FALSE)
write.csv(do.call(rbind,diagnostics),file.path(out,'model_diagnostics.csv'),row.names=FALSE)
write.csv(data.frame(component='current_GO_LORO',status='completed',reason='All three donor-blocked two-region refits; no requirement for leave-out significance'),file.path(out,'analysis_status.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'))
cat('Current-GO leave-out checks complete\n')
