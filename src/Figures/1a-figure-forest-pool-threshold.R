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

load("~/dropbox/coliphage/results/rawoutput/regress-10day-body-pool-threshold.Rdata")

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
# single coliphage results
gici10.body.pool.25=list(overall.fit10.fmc.25, overall.fit10.fpc.25)
gici10.body.pool.50=list(overall.fit10.fmc.50, overall.fit10.fpc.50)
gici10.body.pool.75=list(overall.fit10.fmc.75, overall.fit10.fpc.75)

gici10.body.high.25=list(overall.fit10.fmc.high.25, overall.fit10.fpc.high.25)
gici10.body.high.50=list(overall.fit10.fmc.high.50, overall.fit10.fpc.high.50)
gici10.body.high.75=list(overall.fit10.fmc.high.75, overall.fit10.fpc.high.75)

gici10.body.low.25=list(overall.fit10.fmc.low.25, overall.fit10.fpc.low.25)
gici10.body.low.50=list(overall.fit10.fmc.low.50, overall.fit10.fpc.low.50)
gici10.body.low.75=list(overall.fit10.fmc.low.75, overall.fit10.fpc.low.75)

# ------------------------------------------------
# format output
# ------------------------------------------------
gici10.plot.pool.25=ldply(lapply(gici10.body.pool.25,mkdf),data.frame)
gici10.plot.pool.50=ldply(lapply(gici10.body.pool.50,mkdf),data.frame)
gici10.plot.pool.75=ldply(lapply(gici10.body.pool.75,mkdf),data.frame)

gici10.plot.high.25=ldply(lapply(gici10.body.high.25,mkdf),data.frame)
gici10.plot.high.50=ldply(lapply(gici10.body.high.50,mkdf),data.frame)
gici10.plot.high.75=ldply(lapply(gici10.body.high.75,mkdf),data.frame)

gici10.plot.low.25=ldply(lapply(gici10.body.low.25,mkdf),data.frame)
gici10.plot.low.50=ldply(lapply(gici10.body.low.50,mkdf),data.frame)
gici10.plot.low.75=ldply(lapply(gici10.body.low.75,mkdf),data.frame)

# combine results into one data frame
gici10.plot=rbind(gici10.plot.pool.25,gici10.plot.high.25,gici10.plot.low.25,
                  gici10.plot.pool.50,gici10.plot.high.50,gici10.plot.low.50,
                  gici10.plot.pool.75,gici10.plot.high.75,gici10.plot.low.75)

gici10.plot$lab=as.factor(rep(c(rep("All conditions",2),rep("High risk conditions",2),
                                rep("Low risk conditions",2)),3))

gici10.plot$ind=as.factor(rep(c("Somatic Coliphage","Male-Specific Coliphage"),9))

gici10.plot$type=as.factor(c(rep(">25% of samples detected coliphage",6),
                             rep(">50% of samples detected coliphage",6),
                             rep(">75% of samples detected coliphage",6)))

# order beach factor
gici10.plot$ind.f=factor(gici10.plot$ind,
                         levels=c("Somatic Coliphage","Male-Specific Coliphage"))

# order panels
gici10.plot$lab.f=factor(gici10.plot$lab,levels=c("High risk conditions",
                                                  "Low risk conditions","All conditions"))

# order indicators
gici10.plot$type.f=factor(gici10.plot$type,
                          levels=c(">75% of samples detected coliphage",
                                   ">50% of samples detected coliphage",
                                   ">25% of samples detected coliphage"))


pdf("~/dropbox/coliphage/results/figures/forestplots_pool_threshold.pdf",height=5,width=9)
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
