##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# Forest plot of CIRs for single indicators

# Results pooled across assay and beach
##########################################


rm(list=ls())
library(ggplot2)
library(plyr)
library(scales)
library(grid)

load("~/Documents/CRG/coliphage/results/rawoutput/regress-10day-body-pool-both.Rdata")
load("~/Documents/CRG/coliphage/results/rawoutput/regress-10day-body-entero-pool-both.Rdata")
load("~/Documents/CRG/coliphage/results/rawoutput/regress-10day-body-joint-pool-both.Rdata")  

source("~/Documents/CRG/coliphage/13beaches-coliphage/src/figures/theme_complete_bw.R")
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

mkdf.joint=function(out){
  row=grep("ent4",rownames(out))
  pt.est=out[row,1]
  lb=pt.est-qnorm(.975)*out[4,2]
  ub=pt.est+qnorm(.975)*out[4,2]
  data.frame(pt.est=exp(pt.est),lb=exp(lb),ub=exp(ub))
}

# ------------------------------------------------
# combine output
# ------------------------------------------------
# single coliphage results
gici10.body.pool=overall.fit10
gici10.body.high=overall.fit10.high
gici10.body.low=overall.fit10.low

# entero results
gici10.entero.pool=overall.fit10.entero
gici10.entero.high=overall.fit10.entero.high
gici10.entero.low=overall.fit10.entero.low

# joint indicator results
gici10.body.int=overall.fit10.int
gici10.body.high.int=overall.fit10.high.int
gici10.body.low.int=overall.fit10.low.int

# ------------------------------------------------
# format output
# ------------------------------------------------
gici10.plot.pool=mkdf(gici10.body.pool)
gici10.plot.high=mkdf(gici10.body.high)
gici10.plot.low=mkdf(gici10.body.low)

gici10.plot.entero.pool=mkdf(gici10.entero.pool)
gici10.plot.entero.high=mkdf(gici10.entero.high)
gici10.plot.entero.low=mkdf(gici10.entero.low)

gici10.plot.int=mkdf.joint(gici10.body.int)
gici10.plot.int.high=mkdf.joint(gici10.body.high.int)
gici10.plot.int.low=mkdf.joint(gici10.body.low.int)

# combine results into one data frame
gici10.plot=rbind(gici10.plot.pool,gici10.plot.high,gici10.plot.low,
                  gici10.plot.entero.pool,gici10.plot.entero.high,gici10.plot.entero.low,
                  gici10.plot.int,gici10.plot.int.high,gici10.plot.int.low)

gici10.plot$lab=as.factor(rep(c("All conditions","High risk conditions",
                                "Low risk conditions"),3))

gici10.plot$type=as.factor(c(rep("Coliphage detected",3),rep("Enterococcus > 35 CFU/100 ml",3),
                             rep("Coliphage detected &\nEnterococcus > 35 CFU/100 ml",3)))

# order panels
gici10.plot$lab.f=factor(gici10.plot$lab,levels=c("High risk conditions",
                                                  "Low risk conditions","All conditions"))

# order indicators
gici10.plot$type.f=factor(gici10.plot$type,
                          levels=c("Coliphage detected &\nEnterococcus > 35 CFU/100 ml",
                                                    "Enterococcus > 35 CFU/100 ml",
                                                    "Coliphage detected"))


pdf("~/Documents/CRG/coliphage/results/figures/forestplots_pool_both.pdf",height=5,width=9)
ggplot(gici10.plot,aes(x=lab.f,y=pt.est))+
  geom_errorbar(aes(ymin=lb,ymax=ub,color=type.f),width=0.3,
                position=position_dodge(width=0.4))+
  geom_point(aes(color=type.f,shape=type.f),position=position_dodge(width=0.4))+
  coord_flip()+ geom_hline(yintercept=1,linetype="dotted")+
  scale_color_manual("",breaks=rev(levels(gici10.plot$type.f)),values=c("#048F2E","#0066FF","black"))+
  scale_shape_manual("",breaks=rev(levels(gici10.plot$type.f)),values=c(15,17,19))+  
  scale_y_log10(breaks=c(0.2,0.5,1,2),limits=c(0.3,3))+theme_complete_bw()+
  ylab("Cumulative incidence ratio  (95% CI)")+xlab("")
dev.off()
