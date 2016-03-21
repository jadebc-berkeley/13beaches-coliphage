##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file makes a table with the concentration
# of coliphage by beach and indicator
##########################################


rm(list=ls())
library(foreign)

source("~/Documents/CRG/coliphage/13beaches-coliphage/src/figures/theme_complete_bw.R")


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
wq=read.csv("~/Documents/CRG/coliphage/13beaches-data/temp/beaches-coli-ent-wq.csv")


#-------------------------------------------------------
# Summary of log10 concentration
#-------------------------------------------------------
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
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
wq.table$lab=c("Avalon","Doheny","Mission Bay",
               "Avalon","Doheny",
               "Avalon","Doheny","Malibu","Mission Bay","Fairhope","Goddard",
               "Avalon","Doheny")
wq.table=wq.table[,c(6,1:5)]

rownames(wq.table)=NULL


# drop min column since it's not informative
wq.table=wq.table[,-which(colnames(wq.table)=="min")]

wq.table$n=as.numeric(as.character(wq.table$n))
wq.table$nd=as.numeric(as.character(wq.table$nd))
wq.table$pdetect=(wq.table$n - wq.table$nd)/wq.table$n

wq.table$lab=as.factor(wq.table$lab)
wq.table$pdetect.char=paste(as.numeric(sprintf("%3.0f",wq.table$pdetect*100)),"%",sep="")

# rename coliphage
wq.table$ind=c(rep("Somatic Coliphage (EPA 1601)",3),
               rep("Somatic Coliphage (EPA 1602)",2),
               rep("Male-Specific Coliphage (EPA 1601)",6),
               rep("Male-Specific Coliphage (EPA 1602)",2))

# order ind
wq.table$ind.f=factor(wq.table$ind, levels=c("Male-Specific Coliphage (EPA 1602)",
                                           "Male-Specific Coliphage (EPA 1601)",
                                           "Somatic Coliphage (EPA 1602)",
                                           "Somatic Coliphage (EPA 1601)"))

green="#529C7A"
blue="#5CAFDB"
yellow="#C4E051"
grey="#B4B6B8"
mycolors=c(green,blue,yellow,grey)
mycolors=c(grey,yellow,blue,green)


pdf("~/Documents/CRG/coliphage/results/figures/coliphage_detected.pdf",height=5,width=9)
ggplot(wq.table,aes(x=ind.f,y=pdetect))+geom_bar(aes(fill=ind.f),color="black",stat="identity")+
  geom_text(aes(label=pdetect.char),hjust=-0.15,size=4)+
  facet_wrap(~lab)+coord_flip()+theme_complete_bw()+xlab("")+
  ylab("Percent of samples that detected coliphage")+
  scale_fill_manual(values=mycolors,guide=FALSE)+
  scale_y_continuous(limits=c(0,1),breaks=c(0,.25,.5,.75,1),
                     labels=c("0","25%","50%","75%","100%"))
dev.off()

