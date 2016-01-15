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
wq=read.csv("~/Dropbox/13beaches/data/final/13beaches-wq-samples.csv",stringsAsFactors=FALSE)

# restrict to 6 beaches with coliphage data
adm.list=c("Avalon","Doheny","Malibu")
other.list=c("Mission Bay","Fairhope","Goddard")

wq.adm=wq[wq$beach %in% adm.list,]
wq.oth=wq[wq$beach %in% other.list,]

# subset to relevant indicators
ind.list=c("14FPC051AV071","14FPC051DO071","14FPC051AV071","16FPC081AV071",
 "14FMC051AV071","16FMC081AV071","12ENT041AV071","12ENT041AV081","14FPC051DO071",
 "16FPC081DO071","14FMC051DO071","16FMC081DO071","15ENT041DO071","15ENT041DO081",
 "16FMC081MA091","14FPC051MA091","16FPC081MA091","12ENT041MA091")

wq.adm=wq.adm[wq.adm$groupindex %in% ind.list,]

all.inds=rbind(wq.adm,wq.oth)

all.inds$label=NULL
all.inds$label[all.inds$groupindex=="14FMC051AV071"]="F- Coliphage (EPA 1601)"
all.inds$label[all.inds$groupindex=="14FMC051DO071"]="F- Coliphage (EPA 1601)"

all.inds$label[all.inds$groupindex=="16FMC081AV071"]="F- Coliphage (EPA 1602)"
all.inds$label[all.inds$groupindex=="16FMC081DO071"]="F- Coliphage (EPA 1602)"
all.inds$label[all.inds$groupindex=="16FMC081MA091"]="F- Coliphage (EPA 1602)"

all.inds$label[all.inds$groupindex=="14FPC051AV071"]="F+ Coliphage (EPA 1601)"
all.inds$label[all.inds$groupindex=="14FPC051DO071"]="F+ Coliphage (EPA 1601)"
all.inds$label[all.inds$groupindex=="14FPC051MA091"]="F+ Coliphage (EPA 1601)"

all.inds$label[all.inds$groupindex=="16FPC081AV071"]="F+ Coliphage (EPA 1602)"
all.inds$label[all.inds$groupindex=="16FPC081DO071"]="F+ Coliphage (EPA 1602)"
all.inds$label[all.inds$groupindex=="16FPC081MA091"]="F+ Coliphage (EPA 1602)"

all.inds$label[all.inds$groupindex=="12ENT041AV071"]="EPA 1600 Enterococcus 1"
all.inds$label[all.inds$groupindex=="12ENT041AV081"]="EPA 1600 Enterococcus 2"
all.inds$label[all.inds$groupindex=="15ENT041DO071"]="EPA 1600 Enterococcus 1"
all.inds$label[all.inds$groupindex=="15ENT041DO081"]="EPA 1600 Enterococcus 2"
all.inds$label[all.inds$groupindex=="12ENT041MA091"]="EPA 1600 Enterococcus 1"

all.inds$label[all.inds$analysismethod=="EPA 1601" & 
   all.inds$indicator=="Fminus coliphage"]="F- Coliphage (EPA 1601)"
all.inds$label[all.inds$analysismethod=="EPA 1601" & 
   all.inds$indicator=="Fplus coliphage"]="F+ Coliphage (EPA 1601)"
all.inds$label[all.inds$analysismethod=="EPA 1602" & 
   all.inds$indicator=="Fminus coliphage"]="F- Coliphage (EPA 1602)"
all.inds$label[all.inds$analysismethod=="EPA 1602" & 
   all.inds$indicator=="Fplus coliphage"]="F+ Coliphage (EPA 1602)"


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

wq.table=function(data,label,beach){
  x=data[data$label==label,]
  x=x[x$beach==beach & !is.na(x$beach),]
  #min max
  min=sprintf("%0.1f",quantile(x$result[x$result>=0],probs=c(0,1),na.rm=TRUE)[1])
  if(quantile(x$result[x$result>=0],probs=c(0,1),na.rm=TRUE)[2]<1){
  max=sprintf("%0.1f",quantile(x$result[x$result>=0],probs=c(0,1),na.rm=TRUE)[2])
  }else{
  max=sprintf("%0.0f",quantile(x$result[x$result>=0],probs=c(0,1),na.rm=TRUE)[2])
  }
  #geometric mean
  result=x$result[!is.na(x$result)]
  gm=sprintf("%0.1f",gm_mean(result[result>=0]))
  #non-detects
  qual=x$qualifier
  nd=sprintf("%0.0f",sum(x$qualifier == '<', na.rm=TRUE))
  #n
  n=nrow(x)
  #if(n==nd){
   # min="--"
   # max="--"
   # gm="--"
  #}
  out=data.frame(n=paste(n),min=min,max=max,gm=gm,nd=nd)
  rownames(out)=NULL
  return(out)
}

# avalon
a.1601fmc.tab=wq.table(data=all.inds, beach="Avalon", label="F- Coliphage (EPA 1601)")
a.1602fmc.tab=wq.table(data=all.inds, beach="Avalon", label="F- Coliphage (EPA 1602)")
a.1601fpc.tab=wq.table(data=all.inds, beach="Avalon", label="F+ Coliphage (EPA 1601)")
a.1602fpc.tab=wq.table(data=all.inds, beach="Avalon", label="F+ Coliphage (EPA 1602)")

# doheny
d.1601fmc.tab=wq.table(data=all.inds, beach="Doheny", label="F- Coliphage (EPA 1601)")
d.1602fmc.tab=wq.table(data=all.inds, beach="Doheny", label="F- Coliphage (EPA 1602)")
d.1601fpc.tab=wq.table(data=all.inds, beach="Doheny", label="F+ Coliphage (EPA 1601)")
d.1602fpc.tab=wq.table(data=all.inds, beach="Doheny", label="F+ Coliphage (EPA 1602)")

# malibu
#1601fmc not measured in malibu
m.1602fmc.tab=wq.table(data=all.inds, beach="Malibu", label="F- Coliphage (EPA 1602)")
m.1601fpc.tab=wq.table(data=all.inds, beach="Malibu", label="F+ Coliphage (EPA 1601)")
m.1602fpc.tab=wq.table(data=all.inds, beach="Malibu", label="F+ Coliphage (EPA 1602)")

# mission bay
mb.1601fmc.tab=wq.table(data=all.inds, beach="Mission Bay", label="F- Coliphage (EPA 1601)")
mb.1601fpc.tab=wq.table(data=all.inds, beach="Mission Bay", label="F+ Coliphage (EPA 1601)")

# fairhope
f.1601fpc.tab=wq.table(data=all.inds, beach="Fairhope", label="F+ Coliphage (EPA 1601)")

# goddard
g.1601fpc.tab=wq.table(data=all.inds, beach="Goddard", label="F+ Coliphage (EPA 1601)")

wq.table=rbind(a.1601fmc.tab,d.1601fmc.tab,mb.1601fmc.tab,
                a.1602fmc.tab,d.1602fmc.tab,m.1602fmc.tab,
                a.1601fpc.tab,d.1601fpc.tab,m.1601fpc.tab,mb.1601fpc.tab,f.1601fpc.tab,g.1601fpc.tab,
                a.1602fpc.tab,d.1602fpc.tab,m.1602fpc.tab)
wq.table$lab=c("~~~Avalon","~~~Doheny","~~~Mission Bay",
                "~~~Avalon","~~~Doheny","~~~Malibu",
                "~~~Avalon","~~~Doheny","~~~Malibu","~~~Mission Bay","~~~Fairhope","~~~Goddard",
                "~~~Avalon","~~~Doheny","~~~Malibu")
wq.table=rbind(rep(NA,ncol(wq.table)),wq.table[1:3,],
                rep(NA,ncol(wq.table)),wq.table[4:6,], 
                rep(NA,ncol(wq.table)),wq.table[7:12,],
                rep(NA,ncol(wq.table)),wq.table[13:15,])
wq.table$lab[1]="F- Coliphage (EPA 1601)"
wq.table$lab[5]="F- Coliphage (EPA 1602)"
wq.table$lab[9]="F+ Coliphage (EPA 1601)"
wq.table$lab[16]="F+ Coliphage (EPA 1602)"
wq.table=wq.table[,c(6,1:5)]

save(wq.table,file="~/Dropbox/Coliphage/Data/Temp/wqtable.RData")





