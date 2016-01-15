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
source("~/dropbox/coliphage/programs/figures/theme_complete_bw.R")


setwd("~/Dropbox/Coliphage/")

# --------------------------------------
# load the and pre-preprocess the 
# analysis dataset
# (refer to the base functions script
# for details on the pre-processing)
# --------------------------------------
load("~/dropbox/coliphage/data/temp/beaches-coli-ent-wq.RData")

wq$logfmc1601=log(wq$fmc1601,base=10)
wq$logfmc1602=log(wq$fmc1602,base=10)
wq$logfpc1601=log(wq$fpc1601,base=10)
wq$logfpc1602=log(wq$fpc1602,base=10)
wq$logentero=log(wq$entero,base=10)



fmc1601=ggplot(wq,aes(y=logfmc1601,x=logentero))+geom_point(alpha=0.3)+theme_complete_bw()+
  xlab("")+ylab("")+
  ggtitle("F- Coliphage (EPA 1601)")

fmc1602=ggplot(wq,aes(y=logfmc1602,x=logentero))+geom_point(alpha=0.3)+theme_complete_bw()+
  xlab("")+ylab("")+
  ggtitle("F- Coliphage (EPA 1602)")

fpc1601=ggplot(wq,aes(y=logfpc1601,x=logentero))+geom_point(alpha=0.3)+theme_complete_bw()+
  xlab("")+ylab("")+
  ggtitle("F+ Coliphage (EPA 1601)")

fpc1602=ggplot(wq,aes(y=logfpc1602,x=logentero))+geom_point(alpha=0.3)+theme_complete_bw()+
  xlab("")+ylab("")+
  ggtitle("F+ Coliphage (EPA 1602)")

pdf("~/dropbox/coliphage/results/figures/scatter.pdf",height=9,width=11)
grid.arrange(fmc1601,fmc1602,fpc1601,fpc1602,nrow=2,ncol=2)
dev.off()

wq=wq[,c("beach","logentero","logfmc1601","logfmc1602",
         "logfpc1601","logfpc1602")]
wq.long=melt(wq,ids=c("beach"))
wq.long$logentero=rep(wq.long[wq.long$variable=="logentero","value"],5)
wq.long=subset(wq.long,wq.long$variable!="logentero")

wq.long$variable=as.character(wq.long$variable)
wq.long$variable[wq.long$variable=="logfmc1601"]="F- Coliphage (EPA 1601)"
wq.long$variable[wq.long$variable=="logfmc1602"]="F- Coliphage (EPA 1602)"
wq.long$variable[wq.long$variable=="logfpc1601"]="F+ Coliphage (EPA 1601)"
wq.long$variable[wq.long$variable=="logfpc1602"]="F+ Coliphage (EPA 1602)"



ggplot(wq.long,aes(x=logentero,y=value))+geom_point(alpha=0.3)+
  facet_wrap(~variable)+theme_complete_bw()+
  xlab("Log10 Enterococcus Concentration")+
  ylab("Log10 Coliphage Concentration")



