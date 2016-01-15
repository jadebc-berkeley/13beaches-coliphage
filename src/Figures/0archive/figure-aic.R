##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file creates a figure comparing AIC
# values for different models
##########################################

library(ggplot2)
library(grid)

rm(list=ls())

load("~/dropbox/coliphage/results/rawoutput/regress-10day-body.Rdata")
load("~/dropbox/coliphage/results/rawoutput/regress-10day-body-entero.Rdata")
load("~/dropbox/coliphage/results/rawoutput/regress-10day-body-joint.Rdata")

source("~/dropbox/coliphage/programs/figures/theme_complete_bw.R")

# ------------------------------------------------
# combine output
# ------------------------------------------------
aic=c(aic.fmc1601,aic.fmc1602,aic.fpc1601,aic.fpc1602,
      aic.entero.fmc1601,aic.entero.fmc1602,aic.entero.fpc1601,aic.entero.fpc1602,
      aic.fmc1601.int,aic.fmc1602.int,aic.fpc1601.int,aic.fpc1602.int)
aic.low=c(aic.fmc1602.low,aic.fpc1601.low,aic.fpc1602.low,
          aic.entero.low.fmc1602,aic.entero.low.fpc1601,aic.entero.low.fpc1602,
          aic.fmc1602.low.int,aic.fpc1601.low.int,aic.fpc1602.low.int)
aic.high=c(aic.fmc1602.high,aic.fpc1601.high,aic.fpc1602.high,
           aic.entero.high.fmc1602,aic.entero.high.fpc1601,aic.entero.high.fpc1602,           
           aic.fmc1602.high.int,aic.fpc1601.high.int,aic.fpc1602.high.int)

aic.tab=data.frame(c(aic,aic.low,aic.high))

aic.tab$ind=c(rep(c("F- Coliphage (EPA 1601)","F- Coliphage (EPA 1602)",
      "F+ Coliphage (EPA 1601)", "F+ Coliphage (EPA 1602)"),3),
      rep(c("F- Coliphage (EPA 1602)",
            "F+ Coliphage (EPA 1601)", "F+ Coliphage (EPA 1602)"),6))
aic.tab$ind.f=factor(aic.tab$ind,levels=c("F+ Coliphage (EPA 1602)",
          "F+ Coliphage (EPA 1601)","F- Coliphage (EPA 1602)",
          "F- Coliphage (EPA 1601)"))

aic.tab$type.f=as.factor(c(rep("Single indicator",4),rep("Enterococcus Only",4),rep("Joint indicator",4),
                         rep("Single indicator",3),rep("Enterococcus Only",3),rep("Joint indicator",3),
                         rep("Single indicator",3),rep("Enterococcus Only",3),rep("Joint indicator",3)))

# order panels
aic.tab$cond=as.factor(c(rep("All conditions",12),rep("Low risk conditions",9),
               rep("High risk conditions",9)))
aic.tab$cond.f=factor(aic.tab$cond,levels=c("All conditions","Low risk conditions",
                                            "High risk conditions"))

# create indicator for min within each panel
colnames(aic.tab)[1]=c("aic")
min.all=min(aic.tab$aic[aic.tab$cond=="All conditions"])
min.low=min(aic.tab$aic[aic.tab$cond=="Low risk conditions"])
min.high=min(aic.tab$aic[aic.tab$cond=="High risk conditions"])

aic.tab$min=0
aic.tab$min[aic.tab$aic==min.all & aic.tab$cond=="All conditions"]=1
aic.tab$min[aic.tab$aic==min.low & aic.tab$cond=="Low risk conditions"]=1
aic.tab$min[aic.tab$aic==min.high & aic.tab$cond=="High risk conditions"]=1

aic.tab$min.f=factor(aic.tab$min)
levels(aic.tab$min.f)=c("AIC value","Minimum AIC value")

#Color vector; minimum is blue
mycol=c("black","#0066FF")
names(mycol)=levels(aic.tab$min.f)
colScale <- scale_colour_manual(name = "minimum",values = mycol)


pdf("~/dropbox/coliphage/results/figures/aicplot.pdf",height=4,width=11)
ggplot(aic.tab,aes(x=ind.f,y=aic,group=type.f,col=min.f))+geom_point(aes(shape=type.f),
  position=position_dodge(0.4))+
  facet_wrap(~cond.f,ncol=4)+coord_flip()+
scale_color_manual("",breaks="Minimum AIC value",values=mycol)+
scale_shape_manual("",breaks=rev(levels(aic.tab$type.f)),values=c(17,19,1))+  
  ylab("Akaike Information Criterion")+theme_complete_bw()+xlab(" ")+
  ylim(1000,9000)+theme(legend.title=element_blank())
dev.off()





