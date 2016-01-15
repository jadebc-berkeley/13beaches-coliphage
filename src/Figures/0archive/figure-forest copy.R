##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# Plots of the RR across the range of concentration
##########################################


rm(list=ls())
library(ggplot2)
library(plyr)
library(scales)
library(grid)

load("~/dropbox/coliphage/results/rawoutput/regress-10day-body.Rdata")
load("~/dropbox/coliphage/results/rawoutput/regress-10day-body-entero.Rdata")

source("~/dropbox/coliphage/programs/figures/theme_complete_bw.R")
# ------------------------------------------------
# function to make table row with exponentiated point estimate
# and 95% ci in parentheses
# ------------------------------------------------
mkdf=function(out){
  pt.est=out[2,1]
  lb=pt.est-qnorm(.975)*out[2,2]
  ub=pt.est+qnorm(.975)*out[2,2]
  data.frame(pt.est=exp(pt.est),lb=exp(lb),ub=exp(ub))
}

# ------------------------------------------------
# combine output
# ------------------------------------------------
# estimates pooled across beach  --------------------------------------------
gici10.body.pool=list(overall.fit10.fmc1601, overall.fit10.fmc1602,
                      overall.fit10.fpc1601, overall.fit10.fpc1602)

gici10.body.high=list(overall.fit10.fmc1602.high,
                           overall.fit10.fpc1601.high, overall.fit10.fpc1602.high)

gici10.body.low=list(overall.fit10.fmc1602.low,
                          overall.fit10.fpc1601.low, overall.fit10.fpc1602.low)

# results pooled across beach
gici10.plot.pool=ldply(lapply(gici10.body.pool,mkdf),data.frame)
gici10.plot.high=ldply(lapply(gici10.body.high,mkdf),data.frame)
gici10.plot.low=ldply(lapply(gici10.body.low,mkdf),data.frame)

# entero results
gici10.comb.entero=list(overall.fit10.entero, overall.fit10.entero.high,overall.fit10.entero.low)

gici10.plot.entero=ldply(lapply(gici10.comb.entero,mkdf),data.frame)

# combine results into one data frame
gici10.plot=rbind(gici10.plot.pool,gici10.plot.high,gici10.plot.low,gici10.plot.entero)

gici10.plot$lab=as.factor(c(rep("All conditions",4),rep("High risk conditions",3),
                  rep("Low risk conditions",3),"All conditions","High risk conditions",
                  "Low risk conditions"))

gici10.plot$ind=as.factor(c("F- Coliphage (EPA 1601)","F- Coliphage (EPA 1602)",
                                "F+ Coliphage (EPA 1601)","F+ Coliphage (EPA 1602)", 
                                "F- Coliphage (EPA 1602)",
                                "F+ Coliphage (EPA 1601)","F+ Coliphage (EPA 1602)",
                                "F- Coliphage (EPA 1602)",
                                "F+ Coliphage (EPA 1601)","F+ Coliphage (EPA 1602)",
                            rep("Enterococcus Only",3)))

gici10.plot$type=as.factor(c(rep("coliphage",10),rep("enterococcus",3)))

# order beach factor
gici10.plot$ind.f=factor(gici10.plot$ind,
            levels=c("Enterococcus Only","F+ Coliphage (EPA 1602)","F+ Coliphage (EPA 1601)",
                     "F- Coliphage (EPA 1602)","F- Coliphage (EPA 1601)"))

# order panels
gici10.plot$lab.f=factor(gici10.plot$lab,levels=c("All conditions",
        "Low risk conditions","High risk conditions"))

pdf("~/dropbox/coliphage/results/figures/forestplot.pdf",height=4,width=11)
ggplot(gici10.plot,aes(x=ind.f,y=pt.est,group=lab.f))+
  geom_errorbar(aes(ymin=lb,ymax=ub,color=type),width=0.3,position=position_dodge(width=0.4))+facet_wrap(~lab.f)+
  geom_point(aes(color=type,shape=type))+
  coord_flip()+geom_hline(yintercept=1,linetype="dotted")+
  scale_color_manual(breaks=gici10.plot$type,values=c("black","#0066FF"),guide=FALSE)+
  scale_shape_manual(breaks=gici10.plot$type,values=c(19,17),guide=FALSE)+  
  scale_y_log10(breaks=c(0.1,0.2,0.5,1,2,5),limits=c(0.2,5))+theme_complete_bw()+
  ylab("Cumulative incidence ratio")+xlab("")
dev.off()

save(gici10.plot,file="~/dropbox/coliphage/results/figures/forestplot.RData")




