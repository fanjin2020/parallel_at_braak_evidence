#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
args<-commandArgs(TRUE)
arg<-function(k){i<-match(k,args);if(is.na(i))stop('Missing ',k);args[i+1L]}
src<-normalizePath(arg('--package-dir'),mustWork=TRUE)
out<-arg('--output-dir')
if(dir.exists(out)&&length(list.files(out)))stop('Choose an empty export directory')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
inventory<-list()
for(n in c('GDSCALE','MEDHIST','RECCMEDS','MODHACH','BACKMEDS','INITHEALTH')){
 p<-file.path(src,'data',paste0(n,'.rda'))
 if(!file.exists(p))stop('Missing local ADNIMERGE2 table: ',n)
 e<-new.env();load(p,envir=e);d<-e[[n]]
 if(!is.data.frame(d)||!nrow(d))stop('Invalid table: ',n)
 write.csv(d,file.path(out,paste0(n,'.csv')),row.names=FALSE,na='',fileEncoding='UTF-8')
 inventory[[n]]<-data.frame(table=n,rows=nrow(d),package_version=read.dcf(file.path(src,'DESCRIPTION'))[1,'Version'],
                           md5=unname(tools::md5sum(p)))
 cat(n,nrow(d),'rows exported locally\n')
}
write.csv(do.call(rbind,inventory),file.path(out,'source_inventory.csv'),row.names=FALSE)
