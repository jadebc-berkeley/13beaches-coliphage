##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file creates a figure comparing AIC
# values for different models
##########################################

library(ggplot2)
library(grid)

rm(list=ls())

load("~/dropbox/coliphage/results/rawoutput/regress-10day-body-pool.Rdata")
load("~/dropbox/coliphage/results/rawoutput/regress-10day-body-entero-pool.Rdata")
load("~/dropbox/coliphage/results/rawoutput/regress-10day-body-joint-pool.Rdata")

source("~/dropbox/coliphage/programs/figures/theme_complete_bw.R")

# ------------------------------------------------
# combine output
# ------------------------------------------------
aic=c(aic.fmc,aic.fpc,aic.entero.fmc,aic.entero.fpc,
      aic.fmc.int,aic.fpc.int)
aic.low=c(aic.fmc.low,aic.fpc.high,aic.entero.low.fmc,aic.entero.low.fpc,
          aic.fmc.low.int,aic.fpc.low.int)
aic.high=c(aic.fmc.high,aic.fpc.high,aic.entero.high.fmc,aic.entero.high.fpc,
           aic.fmc.high.int,aic.fpc.high.int)

aic.tab=data.frame(c(aic,aic.low,aic.high))

aic.tab$ind=c(rep(c("F- Coliphage","F+ Coliphage"),9))

aic.tab$type=as.factor(c(rep("Coliphage detected",2),rep("Enterococcus > 35 CFU/100 ml",2),
                         rep("Coliphage detected $\\&$ 35 CFU/100 ml",2),
                         rep("Coliphage detected",2),rep("Enterococcus  > 35 CFU/100 ml",2),
                         rep("Coliphage detected $\\&$ 35 CFU/100 ml",2),
                         rep("Coliphage detected",2),rep("Enterococcus  > 35 CFU/100 ml",2),
                         rep("Coliphage detected $\\&$ 35 CFU/100 ml",2)))
aic.tab$type.f=factor(aic.tab$type,levels=c("Coliphage detected","Coliphage detected $\\&$ 35 CFU/100 ml",
                                            "Enterococcus > 35 CFU/100 ml"))

# order panels
aic.tab$cond=as.factor(c(rep("All conditions",6),rep("Low risk conditions",6),
               rep("High risk conditions",6)))
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


pdf("~/dropbox/coliphage/results/figures/aicplot_pool.pdf",height=4,width=11)
ggplot(aic.tab,aes(x=ind,y=aic,group=type.f,col=min.f))+geom_point(aes(shape=type.f),
  position=position_dodge(0.4))+
  facet_wrap(~cond.f,ncol=4)+coord_flip()+
scale_color_manual("",breaks="Minimum AIC value",values=mycol)+
scale_shape_manual("",breaks=rev(levels(aic.tab$type.f)),values=c(17,19,1))+  
  ylab("Akaike Information Criterion")+theme_complete_bw()+xlab(" ")+
  ylim(1000,9000)+theme(legend.title=element_blank())
dev.off()





