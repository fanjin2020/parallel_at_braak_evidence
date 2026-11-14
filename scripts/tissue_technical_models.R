#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
args<-commandArgs(TRUE)
j<-match('--args-file',args)
if(!is.na(j))args<-readLines(args[j+1],encoding='UTF-8')
arg<-function(k,default=NULL){j<-match(k,args);if(is.na(j))default else args[j+1]}
for(p in c('limma','statmod'))if(!requireNamespace(p,quietly=TRUE))stop(p,' required')
root<-arg('--project-root');out<-arg('--output-dir');go<-arg('--go-dir')
meta_file<-arg('--sample-metadata');cell_dir<-arg('--cell-results')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write<-function(x,n)write.csv(x,file.path(out,paste0(n,'.csv')),row.names=FALSE,na='')
base<-file.path(root,'data/public/gse131617')
m<-read.csv(file.path(base,'sample_manifest.csv'))
e<-read.csv(gzfile(file.path(base,'expression_feature.csv.gz')),check.names=FALSE)
a<-read.csv(file.path(base,'feature_annotation.csv'),colClasses='character')
targets<-read.csv(file.path(root,'data/contracts/gene_symbols.csv'),colClasses='character')
Y<-as.matrix(e[,-1]);storage.mode(Y)<-'double';rownames(Y)<-e[[1]]
m<-m[match(colnames(Y),m$sample_id),]
stopifnot(!anyNA(m$sample_id),!anyDuplicated(m$sample_id),all(is.finite(Y)))
meta<-read.csv(meta_file)
stopifnot(!anyDuplicated(meta$sample_id),all(m$sample_id%in%meta$sample_id))
meta<-meta[match(m$sample_id,meta$sample_id),]
if(!'array_scan_date'%in%names(meta)||anyNA(meta$array_scan_date))stop('Complete scan dates required')
m$array_scan_date<-factor(meta$array_scan_date)
m$braak_stage<-factor(m$braak_stage,levels=c('0','I-II','III-IV','V-VI'))
m$brain_region<-factor(m$brain_region)
a<-a[a$mapping_status=='unique_entrez_id'&a$feature_id%in%rownames(Y),]
a$average<-rowMeans(Y[a$feature_id,,drop=FALSE])
a<-a[order(a$entrez_id,-a$average,a$feature_id),];a<-a[!duplicated(a$entrez_id),]
G<-Y[a$feature_id,,drop=FALSE];rownames(G)<-a$entrez_id
write(a,'gene_representatives')
members<-read.csv(file.path(go,'go_gene_membership.csv'),colClasses='character')
terms<-read.csv(file.path(go,'go_terms.csv'))
sets<-split(members$entrez_id,members$go_id)
if(any(!unique(members$entrez_id)%in%rownames(G)))stop('GO universe mismatch')
write(data.frame(n_features=nrow(Y),n_genes=nrow(G),n_GO=length(sets),
                 n_donors=length(unique(m$donor_key)),n_samples=nrow(m)),'dimensions')
models<-list(original=character(),scan_date='array_scan_date')
for(ref in c('CA','NG'))for(mode in c('joint_quantile','native_logscale')){
 label<-paste(ref,mode,sep='_');file<-file.path(cell_dir,paste0('estimated_fractions_',label,'.csv'))
 if(!file.exists(file))stop('Missing fractions: ',basename(file))
 f<-read.csv(file,check.names=FALSE);stopifnot(!anyDuplicated(f$sample_id))
 f<-f[match(m$sample_id,f$sample_id),];cols<-grep('^cell_prop_',names(f),value=TRUE)
 v<-as.matrix(f[,cols]);stopifnot(!anyNA(v),all(v>=0),all(v<=1),all(abs(rowSums(v)-1)<1e-5))
 omit<-grep('Neurons$',cols,value=TRUE)
 if(length(omit)!=1)stop('One Neurons reference category required')
 keep<-setdiff(cols,omit);new<-paste(label,keep,sep='_');m[,new]<-f[,keep,drop=FALSE]
 models[[paste0('composition_',label)]]<-new
 models[[paste0('joint_',label)]]<-c('array_scan_date',new)
}
status<-list();diagnostics<-list();selected<-list();go_rows<-list();rank_rows<-list()
fit<-function(mat,d,extra){
 d<-droplevels(d)
 if(!all(c('0','V-VI')%in%d$braak_stage))stop('Contrast groups absent')
 base_cov<-c('brain_region','age_at_death','sex_male','pmi_hours','apoe4_carrier',extra)
 base_cov<-base_cov[vapply(base_cov,function(k)length(unique(d[[k]]))>1,logical(1))]
 if(anyNA(d[,c('braak_stage',base_cov)]))stop('Missing design covariates')
 X<-model.matrix(reformulate(c('braak_stage',base_cov),intercept=FALSE),d)
 if(qr(X)$rank<ncol(X)||nrow(X)<=ncol(X)+5)stop('Nonestimable design')
 rho<-limma::duplicateCorrelation(mat,X,block=d$donor_key)$consensus
 if(!is.finite(rho))stop('Nonfinite donor correlation')
 f<-limma::lmFit(mat,X,block=d$donor_key,correlation=rho)
 c<-matrix(0,ncol(X),1)
 c[match('braak_stageV-VI',colnames(X)),1]<-1
 c[match('braak_stage0',colnames(X)),1]<- -1
 tt<-limma::topTable(limma::eBayes(limma::contrasts.fit(f,c),robust=TRUE),
                     number=Inf,sort.by='none',confint=TRUE)
 tab<-data.frame(id=rownames(tt),log2FC=tt$logFC,ci_low=tt$CI.L,ci_high=tt$CI.R,
                 t=tt$t,p=tt$P.Value,q=p.adjust(tt$P.Value,'BH'))
 z<-X;numeric_cols<-which(apply(X,2,sd)>0)
 z[,numeric_cols]<-scale(X[,numeric_cols,drop=FALSE],center=FALSE,scale=apply(X[,numeric_cols,drop=FALSE],2,sd))
 list(table=tab,diagnostic=data.frame(n_samples=ncol(mat),n_donors=length(unique(d$donor_key)),
      n_features=nrow(mat),n_parameters=ncol(X),donor_correlation=rho,
      condition_number=kappa(X),scaled_condition_number=kappa(z)))
}
pathways<-function(tt,label){
 z<-sign(tt$t)*qnorm(pmax(tt$p,.Machine$double.xmin)/2,lower.tail=FALSE);names(z)<-tt$id
 idx<-lapply(sets,function(s)match(s,names(z)))
 if(any(vapply(idx,anyNA,logical(1))))stop('Missing genes in GO family')
 ct<-limma::cameraPR(z,idx,inter.gene.cor=.01,use.ranks=FALSE,sort=FALSE)
 data.frame(go_id=rownames(ct),model=label,direction=ct$Direction,p=ct$PValue,
            q=p.adjust(ct$PValue,'BH'))
}
for(label in names(models)){
 cat('Full-background model:',label,'\n');flush.console()
 tryCatch({
  f<-fit(Y,m,models[[label]]);g<-fit(G,m,models[[label]])
  write(f$table,paste0('features_',label));write(g$table,paste0('genes_',label))
  diagnostics[[paste0(label,'_features')]]<-cbind(model=label,background='features',f$diagnostic)
  diagnostics[[paste0(label,'_genes')]]<-cbind(model=label,background='genes',g$diagnostic)
  selected[[label]]<-merge(transform(f$table,feature_id=id,model=label),targets,by='feature_id')
  rank_rows[[label]]<-transform(g$table,model=label,entrez_id=id)
  go_rows[[label]]<-pathways(g$table,label)
  status[[label]]<-data.frame(component=label,status='completed',reason='')
 },error=function(e)status[[label]]<<-data.frame(component=label,status='not_estimable',reason=conditionMessage(e)))
}
if(length(selected))write(do.call(rbind,selected),'target_effects_full_background')
if(length(rank_rows))write(do.call(rbind,rank_rows),'gene_ranks')
if(length(go_rows)){
 all<-merge(do.call(rbind,go_rows),terms,by='go_id',all.x=TRUE)
 write(all,'GO_all_models')
 s<-all[all$model=='original'&all$q<.05,];s<-s[order(s$q,s$p,s$go_id),]
 reps<-character();groups<-integer(nrow(s))
 for(i in seq_len(nrow(s))){
  jac<-vapply(reps,function(k)length(intersect(sets[[k]],sets[[s$go_id[i]]]))/length(union(sets[[k]],sets[[s$go_id[i]]])),numeric(1))
  if(any(jac>=.5))groups[i]<-which(jac>=.5)[1] else {reps<-c(reps,s$go_id[i]);groups[i]<-length(reps)}
 }
 s$display_group<-groups;s$display_representative<-s$go_id%in%reps;write(s,'GO_grouped')
 write(head(s[s$display_representative,],12),'GO_display_12')
}
write(as.data.frame(table(m$array_scan_date,m$braak_stage,m$brain_region)), 'scan_date_braak_region_counts')
pc<-prcomp(t(G),center=TRUE,scale.=FALSE)
variance<-pc$sdev^2/sum(pc$sdev^2)
pc_tests<-list()
for(i in seq_len(min(10,ncol(pc$x)))){
 d<-m;d$score<-pc$x[,i]
 for(k in c('array_scan_date','braak_stage','brain_region')){
  f<-lm(reformulate(k,'score'),d);null<-lm(score~1,d)
  v<-anova(null,f)
  pc_tests[[length(pc_tests)+1]]<-data.frame(pc=i,variance_fraction=variance[i],variable=k,
      descriptive_r2=summary(f)$r.squared,nominal_p=v$`Pr(>F)`[2],
      note='Sample-level descriptive association; repeated donor tissues not independent')
 }
}
write(do.call(rbind,pc_tests),'PCA_metadata_associations')
influence<-list();loro<-list()
for(label in intersect(c('original','scan_date','composition_CA_joint_quantile','joint_CA_joint_quantile','joint_NG_joint_quantile'),names(models))){
 for(region in levels(m$brain_region)){
  keep<-m$brain_region!=region
  tryCatch({
   g<-fit(G[,keep,drop=FALSE],m[keep,],models[[label]])
   loro[[length(loro)+1]]<-transform(pathways(g$table,label),excluded_region=region)
  },error=function(e)status[[paste(label,region)]]<<-data.frame(component=paste(label,'without',region),status='not_estimable',reason=conditionMessage(e)))
 }
 if(label%in%c('original','scan_date','joint_CA_joint_quantile')){
  for(kind in c('donor_key','array_scan_date'))for(id in unique(as.character(m[[kind]]))){
   keep<-as.character(m[[kind]])!=id
   tryCatch({
    f<-fit(Y[,keep,drop=FALSE],m[keep,],models[[label]])
    t<-f$table[f$table$id%in%targets$feature_id,]
    t$model<-label;t$exclusion_type<-kind;t$exclusion_index<-match(id,unique(as.character(m[[kind]])))
    influence[[length(influence)+1]]<-t
   },error=function(e)status[[paste(label,kind,id)]]<<-data.frame(component=paste(label,kind,match(id,unique(as.character(m[[kind]])))),status='not_estimable',reason=conditionMessage(e)))
  }
 }
}
if(length(loro))write(do.call(rbind,loro),'GO_LORO_all_models')
if(length(influence))write(do.call(rbind,influence),'target_leave_one_out')
write(do.call(rbind,diagnostics),'model_diagnostics');write(do.call(rbind,status),'analysis_status')
capture.output(sessionInfo(),file=file.path(out,'R_session_info.txt'))
writeLines(c('Each model estimates correlation and moderation on its full background.',
 'Probe and gene BH families remain distinct. Composition fractions are estimated, not measured.',
 'Scan date is an acquisition proxy, not RNA extraction batch or RIN.',
 'All GO sets are tested separately in each model and LORO refit.',
 'Selected target intervals remain conditional on prior target selection.',
 'PCA metadata tests are descriptive; no claim of donor-independent significance.'),
 file.path(out,'methods.txt'))
