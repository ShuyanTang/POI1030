#!/bin/Rscript
library(ggplot2)


data<-read.table('collapsing.tsv',header = T, sep = "\t")
data$Gene<-as.factor(data$Gene)


burden_test<-function(dat,type,gp){

  ctrl_name=paste0(type,'_',gp)
  case_name=paste0(type,'_Case')
  dat$p<-c()
  dat$OR<-c()
	side = 'greater'

  for (i in 1:length(dat$Gene)) {
    ctrl=dat[,c(ctrl_name)][i]
    case=dat[,c(case_name)][i]
    
    ctrl_an=dat[,c(paste0('AN_',gp))][i]
    case_an=dat[,c('AN_Case')][i]

    m<-matrix(c(case,ctrl,case_an-case,ctrl_an-ctrl),nrow=2,ncol=2)
    f<-fisher.test(m,alternative = side )
    dat$p[i]<-f$p.value
    dat$OR[i]<-f$estimate
  }
  dat2<-dat[dat[,case_name]>0,c('Gene','p','OR')]
  dat2$p.adj<-p.adjust(dat2$p,method = "BH")
  names(dat2)[names(dat2) == 'p.adj'] <- paste0(type,'.',gp,'.P_adj')
  names(dat2)[names(dat2) == 'p'] <- paste0(type,'.',gp,'.P')
  names(dat2)[names(dat2) == 'OR'] <- paste0(type,'.',gp,'.OR')

  data<<-merge(data,dat2,all.x=T,by.x='Gene')
}


## LoF model
burden_test(data,'LoF','Control')

## missense model
burden_test(data,'Cadd_20','Control')
burden_test(data,'Cadd_10','Control')
burden_test(data,'SPM','Control')
burden_test(data,'REVEL','Control')
burden_test(data,'M.CAP','Control')
burden_test(data,'MetaSVM','Control')

write.table(data,file='genelist.burden_result.tsv',sep = "\t",quote = F,row.names = F)
