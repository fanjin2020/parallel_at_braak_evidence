#!/usr/bin/env Rscript
# Extract published group summaries; these are not sample-level RIN measurements.
args<-commandArgs(trailingOnly=TRUE)
if('--args-file'%in%args)args<-readLines(args[match('--args-file',args)+1L],encoding='UTF-8',warn=FALSE)
arg<-function(k){i<-match(k,args);if(is.na(i))stop('Missing argument: ',k);args[i+1L]}
d<-readxl::read_excel(arg('--source-xls'),sheet='Supple_Table_S1',col_names=FALSE,.name_repair='minimal')
out<-arg('--output-dir');dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(as.character(d[13,3])!='RIN')stop('Unexpected source sheet structure')
rows<-list()
for(pair in list(c('EC',14),c('TC',16),c('FC',18))) {
  region<-pair[1];i<-as.integer(pair[2])
  if(as.character(d[i,3])!=region)stop('Region row mismatch')
  for(j in 1:4)rows[[length(rows)+1L]]<-data.frame(brain_region=region,braak_stage=c('0','I-II','III-IV','V-VI')[j],
    n_donors=c(13,20,19,19)[j],rin_mean=as.numeric(d[i,j+4]),rin_sd=as.numeric(d[i+1,j+4]),
    published_four_group_kw_p=as.numeric(sub('[a-z]$','',as.character(d[i,10]))),
    source='Miyashita et al. 2014, Supplementary Table S1; DOI:10.1038/tp.2014.35',
    interpretation='Published aggregate summary; not individual metadata, not a newly calculated test, not evidence that RNA-quality confounding is absent')
}
write.csv(do.call(rbind,rows),file.path(out,'published_RIN_by_braak_region.csv'),row.names=FALSE)
