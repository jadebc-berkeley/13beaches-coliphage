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

load("~/Documents/CRG/coliphage/results/rawoutput/regress-10day-body-pool.Rdata")
load("~/Documents/CRG/coliphage/results/rawoutput/regress-10day-body-entero-pool.Rdata")
load("~/Documents/CRG/coliphage/results/rawoutput/regress-10day-body-joint-pool.Rdata")  

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
gici10.body.pool=list(overall.fit10.fmc, overall.fit10.fpc)
gici10.body.high=list(overall.fit10.fmc.high, overall.fit10.fpc.high)
gici10.body.low=list(overall.fit10.fmc.low, overall.fit10.fpc.low)

# entero results
gici10.entero.pool=list(overall.fit10.entero.fmc,overall.fit10.entero.fpc)
gici10.entero.high=list(overall.fit10.entero.high.fmc,overall.fit10.entero.high.fpc)
gici10.entero.low=list(overall.fit10.entero.low.fmc,overall.fit10.entero.low.fpc)

# joint indicator results
gici10.body.int=list(overall.fit10.fmc.int, overall.fit10.fpc.int)
gici10.body.high.int=list(overall.fit10.fmc.high.int, overall.fit10.fpc.high.int)
gici10.body.low.int=list(overall.fit10.fmc.low.int, overall.fit10.fpc.low.int)

# ------------------------------------------------
# format output
# ------------------------------------------------
gici10.plot.pool=ldply(lapply(gici10.body.pool,mkdf),data.frame)
gici10.plot.high=ldply(lapply(gici10.body.high,mkdf),data.frame)
gici10.plot.low=ldply(lapply(gici10.body.low,mkdf),data.frame)

gici10.plot.entero.pool=ldply(lapply(gici10.entero.pool,mkdf),data.frame)
gici10.plot.entero.high=ldply(lapply(gici10.entero.high,mkdf),data.frame)
gici10.plot.entero.low=ldply(lapply(gici10.entero.low,mkdf),data.frame)

gici10.plot.int=ldply(lapply(gici10.body.int,mkdf.joint),data.frame)
gici10.plot.int.high=ldply(lapply(gici10.body.high.int,mkdf.joint),data.frame)
gici10.plot.int.low=ldply(lapply(gici10.body.low.int,mkdf.joint),data.frame)

# combine results into one data frame
gici10.plot=rbind(gici10.plot.pool,gici10.plot.high,gici10.plot.low,
                  gici10.plot.entero.pool,gici10.plot.entero.high,gici10.plot.entero.low,
                  gici10.plot.int,gici10.plot.int.high,gici10.plot.int.low)

gici10.plot$lab=as.factor(rep(c(rep("All conditions",2),rep("High risk conditions",2),
                                rep("Low risk conditions",2)),3))

gici10.plot$ind=as.factor(rep(c("Somatic Coliphage","Male-Specific Coliphage"),9))

gici10.plot$type=as.factor(c(rep("Coliphage detected",6),rep("Enterococci > 35 CFU/100 ml",6),
                             rep("Coliphage detected &\nEnterococci > 35 CFU/100 ml",6)))

# order beach factor
gici10.plot$ind.f=factor(gici10.plot$ind,
                         levels=c("Somatic Coliphage","Male-Specific Coliphage"))

# order panels
gici10.plot$lab.f=factor(gici10.plot$lab,levels=c("High risk conditions",
                                                  "Low risk conditions","All conditions"))

# order indicators
gici10.plot$type.f=factor(gici10.plot$type,
                          levels=c("Coliphage detected &\nEnterococci > 35 CFU/100 ml",
                                                    "Enterococci > 35 CFU/100 ml",
                                                    "Coliphage detected"))


pdf("~/Documents/CRG/coliphage/results/figures/forestplots_pool.pdf",height=5,width=9)
ggplot(gici10.plot,aes(x=lab.f,y=pt.est))+
  geom_errorbar(aes(ymin=lb,ymax=ub,color=type.f),width=0.3,
                position=position_dodge(width=0.4))+
  facet_wrap(~ind.f,ncol=4)+
  geom_point(aes(color=type.f,shape=type.f),position=position_dodge(width=0.4))+
  coord_flip()+ geom_hline(yintercept=1,linetype="dotted")+
  scale_color_manual("",breaks=rev(levels(gici10.plot$type.f)),values=c("#048F2E","#0066FF","black"))+
  scale_shape_manual("",breaks=rev(levels(gici10.plot$type.f)),values=c(15,17,19))+  
  scale_y_log10(breaks=c(0.2,0.5,1,2),limits=c(0.3,3))+theme_complete_bw()+
  ylab("Cumulative incidence ratio  (95% CI)")+xlab("")
dev.off()



pdf("~/Documents/CRG/coliphage/results/figures/forestplots_pool_present.pdf",height=3,width=8)
ggplot(gici10.plot,aes(x=lab.f,y=pt.est))+
  geom_errorbar(aes(ymin=lb,ymax=ub,color=type.f),width=0.3,
                position=position_dodge(width=0.4))+
  facet_wrap(~ind.f,ncol=4)+
  geom_point(aes(color=type.f,shape=type.f),position=position_dodge(width=0.4))+
  coord_flip()+ geom_hline(yintercept=1,linetype="dotted")+
  scale_color_manual("",breaks=rev(levels(gici10.plot$type.f)),values=c("#048F2E","#0066FF","black"))+
  scale_shape_manual("",breaks=rev(levels(gici10.plot$type.f)),values=c(15,17,19))+  
  scale_y_log10(breaks=c(0.2,0.5,1,2),limits=c(0.3,3))+theme_complete_bw()+
  ylab("Cumulative incidence ratio  (95% CI)")+xlab("")
dev.off()
