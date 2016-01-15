##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file makes a table with the concentration
# of coliphage by beach and indicator
##########################################


rm(list=ls())
library(foreign)

setwd("~/Dropbox/Coliphage/")

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
load("~/dropbox/coliphage/data/temp/beaches-coli-ent-wq.RData")


#-------------------------------------------------------
# Summary of log10 concentration
#-------------------------------------------------------
gm_mean = function(a){
  if(length(a)==0){
    NA
  }else{
    prod(a)^(1/length(a))
  }
}

wq.table=function(data,ind,beach){
  x.beach=data[data$beach==beach & !is.na(data$beach),]
  x=x.beach[[ind]]
  
  #min max
  min=sprintf("%0.1f",quantile(x[x>=0],probs=c(0,1),na.rm=TRUE)[1])
  if(quantile(x[x>=0],probs=c(0,1),na.rm=TRUE)[2]<1){
  max=sprintf("%0.1f",quantile(x[x>=0],probs=c(0,1),na.rm=TRUE)[2])
  }else{
  max=sprintf("%0.0f",quantile(x[x>=0],probs=c(0,1),na.rm=TRUE)[2])
  }
  #geometric mean
  x=x[!is.na(x)]
  gm=sprintf("%0.1f",gm_mean(x[x>=0]))
  #non-detects
  x_nd=x.beach[[paste(ind,"_nd",sep="")]]
  nd=sprintf("%0.0f",sum(x_nd == 'Below detection', na.rm=TRUE))
  #n
  n=length(x)

  out=data.frame(n=paste(n),min=min,max=max,gm=gm,nd=nd)
  rownames(out)=NULL
  return(out)
}


# avalon
a.1601fmc.tab=wq.table(data=wq, beach="Avalon", ind="fmc1601")
a.1602fmc.tab=wq.table(data=wq, beach="Avalon", ind="fmc1602")
a.1601fpc.tab=wq.table(data=wq, beach="Avalon", ind="fpc1601")
a.1602fpc.tab=wq.table(data=wq, beach="Avalon", ind="fpc1602")

# doheny
d.1601fmc.tab=wq.table(data=wq, beach="Doheny", ind="fmc1601")
d.1602fmc.tab=wq.table(data=wq, beach="Doheny", ind="fmc1602")
d.1601fpc.tab=wq.table(data=wq, beach="Doheny", ind="fpc1601")
d.1602fpc.tab=wq.table(data=wq, beach="Doheny", ind="fpc1602")

# malibu
m.1601fpc.tab=wq.table(data=wq, beach="Malibu", ind="fpc1601")
m.1602fpc.tab=wq.table(data=wq, beach="Malibu", ind="fpc1602")

# mission bay
mb.1601fmc.tab=wq.table(data=wq, beach="Mission Bay", ind="fmc1601")
mb.1601fpc.tab=wq.table(data=wq, beach="Mission Bay", ind="fpc1601")

# fairhope
f.1601fpc.tab=wq.table(data=wq, beach="Fairhope", ind="fpc1601")

# goddard
g.1601fpc.tab=wq.table(data=wq, beach="Goddard", ind="fpc1601")

wq.table=rbind(a.1601fmc.tab,d.1601fmc.tab,mb.1601fmc.tab,
                a.1602fmc.tab,d.1602fmc.tab,
                a.1601fpc.tab,d.1601fpc.tab,m.1601fpc.tab,mb.1601fpc.tab,f.1601fpc.tab,g.1601fpc.tab,
                a.1602fpc.tab,d.1602fpc.tab)
wq.table$lab=c("~~~Avalon","~~~Doheny","~~~Mission Bay",
                "~~~Avalon","~~~Doheny",
                "~~~Avalon","~~~Doheny","~~~Malibu","~~~Mission Bay","~~~Fairhope","~~~Goddard",
                "~~~Avalon","~~~Doheny")
wq.table=rbind(rep(NA,ncol(wq.table)),wq.table[1:3,],
                rep(NA,ncol(wq.table)),wq.table[4:5,], 
                rep(NA,ncol(wq.table)),wq.table[6:11,],
                rep(NA,ncol(wq.table)),wq.table[12:13,])
wq.table$lab[1]="Somatic coliphage (EPA 1601)"
wq.table$lab[5]="Somatic coliphage (EPA 1602)"
wq.table$lab[8]="Male-specific coliphage (EPA 1601)"
wq.table$lab[15]="Male-specific coliphage (EPA 1602)"
wq.table=wq.table[,c(6,1:5)]

rownames(wq.table)=NULL

save(wq.table,file="~/Dropbox/Coliphage/Data/Temp/wqtable.RData")

wq.table[which(wq.table$lab=="~~~Avalon"),"lab"]="Avalon"
wq.table[which(wq.table$lab=="~~~Doheny"),"lab"]="Doheny"
wq.table[which(wq.table$lab=="~~~Malibu"),"lab"]="Malibu"
wq.table[which(wq.table$lab=="~~~Mission Bay"),"lab"]="Mission Bay"
wq.table[which(wq.table$lab=="~~~Fairhope"),"lab"]="Fairhope"
wq.table[which(wq.table$lab=="~~~Goddard"),"lab"]="Goddard"

write.csv(wq.table,file="~/Dropbox/Coliphage/Results/Tables/Table1.csv",na="",row.names=FALSE)






