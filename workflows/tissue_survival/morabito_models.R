args <- commandArgs(TRUE)
if (length(args)!=3) stop('Expected pseudobulk_directory go_annotation_directory output_directory')
for (p in c('edgeR','limma','statmod')) if(!requireNamespace(p,quietly=TRUE)) stop('Missing package: ',p)
src<-args[1]; go<-args[2]; out<-args[3]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write<-function(x,name)write.csv(x,file.path(out,paste0(name,'.csv')),row.names=FALSE,na='')
donors<-read.csv(file.path(src,'donors.csv'),stringsAsFactors=FALSE)
genes<-read.csv(file.path(src,'genes.csv'),stringsAsFactors=FALSE)
cells<-read.csv(file.path(src,'donor_cell_counts.csv'),check.names=FALSE)
types<-c('MG','ASC','ODC','OPC','EX','INH','PER.END')
targets<-c('MTHFD2','OXTR','HOPX','CXCL10','IRF1','CD14','TAP1','PTPN6','TNFAIP2','ICAM1','HAMP','ITPKC','JAK3','GMFG','HCK','CEBPB')
if(anyDuplicated(donors$donor_id)||anyDuplicated(genes$gene_id))stop('Duplicate donor or gene IDs')
cells<-cells[match(donors$donor_id,cells$donor_id),,drop=FALSE]
if(anyNA(cells$donor_id))stop('Cell counts not aligned')
members<-read.csv(file.path(go,'go_gene_membership.csv'),colClasses='character')
dictionary<-read.csv(file.path(go,'human_gene_dictionary.csv'),colClasses='character')
terms<-sort(unique(members$go_id))
if(length(terms)!=4826L)stop('Expected the retained 4826-set GO family; do not replace it with a significant-only subset')
duplicate_symbols<-duplicated(dictionary$gene_symbol)|duplicated(dictionary$gene_symbol,fromLast=TRUE)
dictionary<-dictionary[!duplicate_symbols,,drop=FALSE]
gene_symbol_unique<-!(duplicated(genes$gene_symbol)|duplicated(genes$gene_symbol,fromLast=TRUE))
mapping<-dictionary$entrez_id[match(genes$gene_symbol,dictionary$gene_symbol)]
mapping[!gene_symbol_unique]<-NA_character_
mapping[duplicated(mapping)|duplicated(mapping,fromLast=TRUE)]<-NA_character_
names(mapping)<-genes$gene_id
sets<-split(members$entrez_id,members$go_id)
status<-list(); results<-list(); pathways<-list(); summaries<-list()
counter<-0L
for(kind in types){
  counts<-as.matrix(read.csv(gzfile(file.path(src,paste0('counts_',kind,'.csv.gz'))),row.names=1,check.names=FALSE))
  storage.mode(counts)<-'double'
  if(!identical(colnames(counts),donors$donor_id)||!identical(rownames(counts),genes$gene_id))stop('Counts not aligned')
  if(anyNA(counts)||any(counts<0)||any(counts!=round(counts)))stop('Invalid pseudobulk counts')
  for(model in c('age_sex_batch_RIN','plus_PMI')){
    counter<-counter+1L; key<-paste(kind,model,sep='__')
    keep<-cells[[kind]]>=20 & donors$braak %in% c(1,2,5,6)
    if(model=='plus_PMI')keep<-keep & !is.na(donors$pmi)
    m<-donors[keep,,drop=FALSE]; m$late<-as.integer(m$braak>=5); m$batch<-factor(m$batch)
    low<-sum(m$late==0); high<-sum(m$late==1); reason<-NULL; fit<-NULL; table<-NULL
    if(min(low,high)<6)reason<-'Fewer than six independent donors in a contrast group'
    formula<-if(model=='plus_PMI')~late+age+sex_male+batch+rin+pmi else ~late+age+sex_male+batch+rin
    X<-NULL
    if(is.null(reason)){
      design<-tryCatch(model.matrix(formula,m),error=function(e)e)
      if(inherits(design,'error'))reason<-paste('Design not estimable:',conditionMessage(design)) else {
        X<-design
        if(nrow(X)!=nrow(m)||qr(X)$rank!=ncol(X)||nrow(X)-ncol(X)<5)reason<-'Incomplete/rank-deficient design or residual df below five; no covariates removed'
      }
    }
    if(is.null(reason)){
      attempt<-tryCatch({
        d<-edgeR::DGEList(counts=counts[,keep,drop=FALSE]); expressed<-edgeR::filterByExpr(d,design=X)
        if(sum(expressed)<100)stop('Fewer than 100 expressed features')
        d<-edgeR::calcNormFactors(d[expressed,,keep.lib.sizes=FALSE],method='TMM')
        v<-limma::voom(d,X,plot=FALSE)
        f<-limma::eBayes(limma::lmFit(v,X),robust=TRUE)
        t<-limma::topTable(f,coef='late',number=Inf,sort.by='none',confint=TRUE)
        list(fit=f,table=t)
      },error=function(e)e)
      if(inherits(attempt,'error'))reason<-conditionMessage(attempt) else {fit<-attempt$fit;table<-attempt$table}
    }
    status[[key]]<-data.frame(cell_type=kind,model=model,n_low=low,n_high=high,n_donors=nrow(m),
      status=if(is.null(reason))'estimated' else 'not_estimable',reason=if(is.null(reason))'' else reason,
      genes_tested=if(is.null(reason))nrow(table) else NA_integer_,contrast='Braak V-VI versus I-II',
      formula=paste(deparse(formula),collapse=''),primary=kind=='MG' & model=='age_sex_batch_RIN')
    r<-data.frame(cell_type=kind,model=model,gene_symbol=targets,status='not_estimable',
      logFC=NA_real_,CI_L=NA_real_,CI_R=NA_real_,PValue=NA_real_,FDR_genomewide=NA_real_)
    p<-data.frame(cell_type=kind,model=model,go_id=terms,n_genes=0L,status='model_not_estimable',
      direction=NA_character_,PValue=NA_real_)
    if(is.null(reason)){
      target_ids<-genes$gene_id[match(targets,genes$gene_symbol)]
      ambiguous<-!gene_symbol_unique[match(targets,genes$gene_symbol)]
      idx<-match(target_ids,rownames(table));idx[which(ambiguous)]<-NA_integer_
      r$status<-ifelse(is.na(idx),'filtered_low_expression','tested')
      r$status[is.na(match(targets,genes$gene_symbol))]<-'absent_from_input'
      r$status[which(ambiguous)]<-'ambiguous_gene_symbol'
      r$logFC<-table$logFC[idx];r$CI_L<-table$CI.L[idx];r$CI_R<-table$CI.R[idx]
      r$PValue<-table$P.Value[idx];r$FDR_genomewide<-table$adj.P.Val[idx]
      z<-fit$t[,'late'];names(z)<-rownames(fit$t)
      gid<-mapping[names(z)]; ok<-!is.na(gid); z<-z[ok];names(z)<-gid[ok]
      index<-lapply(sets,function(x)which(names(z)%in%x));sizes<-lengths(index)
      valid<-sizes>=10 & sizes<=300
      p$n_genes<-sizes[match(p$go_id,names(sizes))]
      p$status<-'insufficient_measured_members'
      if(any(valid)){
        cam<-limma::cameraPR(z,index[valid],inter.gene.cor=.01,use.ranks=FALSE,sort=FALSE)
        ix<-match(p$go_id,rownames(cam)); good<-!is.na(ix)
        p$status[good]<-'tested';p$direction[good]<-cam$Direction[ix[good]];p$PValue[good]<-cam$PValue[ix[good]]
      }
      summaries[[key]]<-data.frame(cell_type=kind,model=model,genomewide_features=nrow(table),
        unique_GO_mapped_features=length(z),GO_tested=sum(p$status=='tested'))
    }
    p$q_fixed_GO_family<-p.adjust(ifelse(is.na(p$PValue),1,p$PValue),method='BH',n=4826)
    results[[key]]<-r;pathways[[key]]<-p
    write(do.call(rbind,status),'GSE174367_model_status')
    cat(counter,'/14',kind,model,status[[key]]$status,'n=',nrow(m),'\n')
  }
}
res<-do.call(rbind,results)
res$q_all_planned_targets<-p.adjust(ifelse(is.na(res$PValue),1,res$PValue),method='BH')
res$q_primary_MG<-NA_real_;primary<-res$cell_type=='MG' & res$model=='age_sex_batch_RIN'
res$q_primary_MG[primary]<-p.adjust(ifelse(is.na(res$PValue[primary]),1,res$PValue[primary]),method='BH',n=16)
res$discovery_direction<-'down';res$direction_match<-ifelse(res$status=='tested',res$logFC<0,NA)
gores<-do.call(rbind,pathways)
gores$q_all_cell_models<-p.adjust(ifelse(is.na(gores$PValue),1,gores$PValue),method='BH')
write(res,'GSE174367_target_results');write(gores,'GSE174367_GO_all_models')
n0<-sum(donors$braak==0)
write(data.frame(status=if(n0==0)'not_estimable' else 'not_analyzed',contrast='Braak V-VI versus 0',n_stage_0=n0,
                 reason=if(n0==0)'No Braak-0 donors in this independent cohort' else 'Metadata differs from the inspected release; review the contrast plan before additional fitting'),'GSE174367_exact_contrast_status')
if(length(summaries))write(do.call(rbind,summaries),'GSE174367_model_dimensions')
writeLines(capture.output(sessionInfo()),file.path(out,'GSE174367_R_sessionInfo.txt'))
