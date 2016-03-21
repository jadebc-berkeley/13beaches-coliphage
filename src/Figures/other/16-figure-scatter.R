##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 8/17/15

# Scatter plot of enterococcus by coliphage
# concentration
##########################################


rm(list=ls())
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)
source("~/Documents/CRG/coliphage/13beaches-coliphage/src/figures/theme_complete_bw.R")


# --------------------------------------
# load the and pre-preprocess the 
# analysis dataset
# (refer to the base functions script
# for details on the pre-processing)
# --------------------------------------
wq=read.csv("~/Documents/CRG/coliphage/13beaches-data/temp/beaches-coli-ent-wq.csv",stringsAsFactors=TRUE)

wq$fmc1601[wq$fmc1601_nd=="Below detection"]=0.1
wq$fmc1602[wq$fmc1602_nd=="Below detection"]=0.1
wq$fpc1601[wq$fpc1601_nd=="Below detection"]=0.1
wq$fpc1602[wq$fpc1602_nd=="Below detection"]=0.1
wq$entero[wq$entero_nd=="Below detection"]=0.1

wq$logfmc1601=log(wq$fmc1601,base=10)
wq$logfmc1602=log(wq$fmc1602,base=10)
wq$logfpc1601=log(wq$fpc1601,base=10)
wq$logfpc1602=log(wq$fpc1602,base=10)
wq$logentero=log(wq$entero,base=10)

wq=wq[,c("beach","logentero","logfmc1601","logfmc1602",
         "logfpc1601","logfpc1602")]
wq.long=melt(wq,ids=c("beach"))
wq.long$logentero=rep(wq.long[wq.long$variable=="logentero","value"],5)
wq.long=subset(wq.long,wq.long$variable!="logentero")

wq.long$variable=as.character(wq.long$variable)
wq.long$variable[wq.long$variable=="logfmc1601"]="Somatic Coliphage (EPA 1601)"
wq.long$variable[wq.long$variable=="logfmc1602"]="Somatic Coliphage (EPA 1602)"
wq.long$variable[wq.long$variable=="logfpc1601"]="Male-Specific Coliphage (EPA 1601)"
wq.long$variable[wq.long$variable=="logfpc1602"]="Male-Specific Coliphage (EPA 1602)"

# reorder coliphage type
wq.long$variable.f=factor(wq.long$variable, levels=c("Somatic Coliphage (EPA 1601)",
                                                     "Somatic Coliphage (EPA 1602)",
                                                     "Male-Specific Coliphage (EPA 1601)",
                                                     "Male-Specific Coliphage (EPA 1602)"))


pdf("~/Documents/CRG/coliphage/results/figures/scatter.pdf",height=5,width=7)
ggplot(wq.long,aes(x=logentero,y=value))+geom_point(alpha=0.3)+
  facet_wrap(~variable.f)+theme_complete_bw()+
  scale_x_continuous(limits=c(-1.1,5),breaks=c(-1,0,1,2,3,4))+
  scale_y_continuous(limits=c(-1.3,4),breaks=c(-1,0,1,2,3))+
  xlab("Log10 Enterococcus Concentration")+
  ylab("Log10 Coliphage Concentration")
dev.off()





