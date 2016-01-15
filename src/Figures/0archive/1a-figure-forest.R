##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# Forest plot of CIRs for single indicators
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
gici10.body.pool=list(overall.fit10.fmc1601, overall.fit10.fmc1602,
                      overall.fit10.fpc1601, overall.fit10.fpc1602)

gici10.body.high=list(overall.fit10.fmc1602.high,
                           overall.fit10.fpc1601.high, overall.fit10.fpc1602.high)

gici10.body.low=list(overall.fit10.fmc1602.low,
                          overall.fit10.fpc1601.low, overall.fit10.fpc1602.low)

# entero results
gici10.entero.pool=list(overall.fit10.entero.fmc1601,overall.fit10.entero.fmc1602,
                        overall.fit10.entero.fpc1601,overall.fit10.entero.fpc1602)
gici10.entero.high=list(overall.fit10.entero.high.fmc1602,overall.fit10.entero.high.fpc1601,
                        overall.fit10.entero.high.fpc1602)
gici10.entero.low=list(overall.fit10.entero.low.fmc1602,overall.fit10.entero.low.fpc1601,
                       overall.fit10.entero.low.fpc1602)
  
# ------------------------------------------------
# format output
# ------------------------------------------------
gici10.plot.pool=ldply(lapply(gici10.body.pool,mkdf),data.frame)
gici10.plot.high=ldply(lapply(gici10.body.high,mkdf),data.frame)
gici10.plot.low=ldply(lapply(gici10.body.low,mkdf),data.frame)

gici10.plot.entero.pool=ldply(lapply(gici10.entero.pool,mkdf),data.frame)
gici10.plot.entero.high=ldply(lapply(gici10.entero.high,mkdf),data.frame)
gici10.plot.entero.low=ldply(lapply(gici10.entero.low,mkdf),data.frame)

# combine results into one data frame
gici10.plot=rbind(gici10.plot.pool,gici10.plot.high,gici10.plot.low,
                  gici10.plot.entero.pool,gici10.plot.entero.high,gici10.plot.entero.low)

gici10.plot$lab=as.factor(rep(c(rep("All conditions",4),rep("High risk conditions",3),
                  rep("Low risk conditions",3)),2))

gici10.plot$ind=as.factor(rep(c("F- Coliphage (EPA 1601)","F- Coliphage (EPA 1602)",
                                "F+ Coliphage (EPA 1601)","F+ Coliphage (EPA 1602)", 
                                "F- Coliphage (EPA 1602)",
                                "F+ Coliphage (EPA 1601)","F+ Coliphage (EPA 1602)",
                                "F- Coliphage (EPA 1602)",
                                "F+ Coliphage (EPA 1601)","F+ Coliphage (EPA 1602)"),2))

gici10.plot$type=as.factor(c(rep("Coliphage detected",10),rep("Enterococcus > 35 CFU/100 ml",10)))
gici10.plot$type.f <- factor(gici10.plot$type, levels=rev(levels(gici10.plot$type)) )

# order beach factor
gici10.plot$ind.f=factor(gici10.plot$ind,
            levels=c("F- Coliphage (EPA 1601)","F- Coliphage (EPA 1602)",
                     "F+ Coliphage (EPA 1601)","F+ Coliphage (EPA 1602)"))

# order panels
gici10.plot$lab.f=factor(gici10.plot$lab,levels=c("High risk conditions",
        "Low risk conditions","All conditions"))

forest=ggplot(gici10.plot,aes(x=lab.f,y=pt.est))+
  geom_errorbar(aes(ymin=lb,ymax=ub,color=type.f),width=0.3,
  position=position_dodge(width=0.4))+
  facet_wrap(~ind.f,ncol=4)+
  geom_point(aes(color=type.f,shape=type.f),position=position_dodge(width=0.4))+
  coord_flip()+ geom_hline(yintercept=1,linetype="dotted")+
  scale_color_manual("",breaks=rev(levels(gici10.plot$type.f)),values=c("#0066FF","black"))+
  scale_shape_manual("",breaks=rev(levels(gici10.plot$type.f)),values=c(17,19))+  
  scale_y_log10(breaks=c(0.1,0.2,0.5,1,2,5),limits=c(0.2,5))+theme_complete_bw()+
  ylab("Cumulative incidence ratio  (95% CI)")+xlab("")+
  ggtitle("Single indicators")

save(forest,file="~/dropbox/coliphage/results/rawoutput/gg_forest.RData")
