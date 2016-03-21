##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file makes a table with the concentration
# of coliphage by beach and indicator
##########################################


rm(list=ls())
library(foreign)

# --------------------------------------
# paste n followed by percent in parentheses
# input numeric variables
# output string
# --------------------------------------
n.p.paren=function(n,p){
  paste(n," (",as.numeric(sprintf("%0.1f",p*100)),")",sep="")
}

# --------------------------------------
# load the and pre-preprocess the 
# analysis dataset
# (refer to the base functions script
# for details on the pre-processing)
# --------------------------------------
wq=read.csv("~/Documents/CRG/coliphage/13beaches-data/temp/beaches-coli-ent-wq.csv",stringsAsFactors=TRUE)


#-------------------------------------------------------
# Summary of log10 concentration
#-------------------------------------------------------
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

wq.table=function(ind,cond){
  data=wq
  data=subset(data,data$risk==cond)
  x=data[[ind]]
  #min max
  min=sprintf("%0.1f",quantile(x[x>=0],probs=c(0,1),na.rm=TRUE)[1])
  if(quantile(x[x>=0],probs=c(0,1),na.rm=TRUE)[2]<1){
    max=sprintf("%0.1f",quantile(x[x>=0],probs=c(0,1),na.rm=TRUE)[2])
  }else{
    max=sprintf("%0.0f",quantile(x[x>=0],probs=c(0,1),na.rm=TRUE)[2])
  }
  #geometric mean
  gm=sprintf("%0.2f",gm_mean(x[x>=0]))
  #non-detects
  x_nd=data[[paste(ind,"_nd",sep="")]]
  nd=sprintf("%0.0f",sum(x_nd == 'Below detection', na.rm=TRUE))
  #n
  n=sum(table(x))
  
  out=data.frame(n=paste(n),min=min,max=max,gm=gm,nd=nd)
  rownames(out)=NULL
  return(out)
}

fmc1601tab.low=wq.table(ind="fmc1601",cond="Low")
fmc1601tab.high=wq.table(ind="fmc1601",cond="High")
fmc1602tab.low=wq.table(ind="fmc1602",cond="Low")
fmc1602tab.high=wq.table(ind="fmc1602",cond="High")
fpc1601tab.low=wq.table(ind="fpc1601",cond="Low")
fpc1601tab.high=wq.table(ind="fpc1601",cond="High")
fpc1602tab.low=wq.table(ind="fpc1602",cond="Low")
fpc1602tab.high=wq.table(ind="fpc1602",cond="High")

wq.table=rbind(fmc1601tab.low,fmc1601tab.high,
               fmc1602tab.low,fmc1602tab.high,
               fpc1601tab.low,fpc1601tab.high,
               fpc1602tab.low,fpc1602tab.high)
wq.table$lab=c(rep(c("~~~~Low risk","~~~~High risk"),4))
wq.table=rbind(rep(NA,ncol(wq.table)),wq.table[1:2,],
               rep(NA,ncol(wq.table)),wq.table[3:4,], 
               rep(NA,ncol(wq.table)),wq.table[5:6,],
               rep(NA,ncol(wq.table)),wq.table[7:8,])
wq.table$lab[1]="\\textbf{Somatic coliphage (EPA 1601)}"
wq.table$lab[4]="\\textbf{Somatic coliphage (EPA 1602)}"
wq.table$lab[7]="\\textbf{Male-specific coliphage (EPA 1601)}"
wq.table$lab[10]="\\textbf{Male-specific coliphage (EPA 1602)}"
wq.table=wq.table[,c(6,1:5)]

rownames(wq.table)=NULL

save(wq.table,file="~/Documents/CRG/coliphage/13beaches-data/Temp/wqtable_risk.RData")

write.csv(wq.table,file="~/Documents/CRG/coliphage/Results/Tables/Table2.csv",na="",row.names=FALSE)







