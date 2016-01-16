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

load("~/Documents/CRG/coliphage/results/rawoutput/regress-10day-body-threshold.Rdata")

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

# ------------------------------------------------
# combine output
# ------------------------------------------------
# single coliphage results
gici10.body.pool.25=list(overall.fit10.fmc1601.25, overall.fit10.fmc1602.25,
                      overall.fit10.fpc1601.25, overall.fit10.fpc1602.25)
gici10.body.pool.50=list(overall.fit10.fmc1601.50, overall.fit10.fmc1602.50,
                      overall.fit10.fpc1601.50, overall.fit10.fpc1602.50)
gici10.body.pool.75=list(overall.fit10.fmc1601.75, overall.fit10.fmc1602.75,
                      overall.fit10.fpc1601.75, overall.fit10.fpc1602.75)

gici10.body.high.25=list(overall.fit10.fmc1602.high.25,
                      overall.fit10.fpc1601.high.25, overall.fit10.fpc1602.high.25)
gici10.body.high.50=list(overall.fit10.fmc1602.high.50,
                      overall.fit10.fpc1601.high.50, overall.fit10.fpc1602.high.50)
gici10.body.high.75=list(overall.fit10.fmc1602.high.75,
                      overall.fit10.fpc1601.high.75, overall.fit10.fpc1602.high.75)

gici10.body.low.25=list(overall.fit10.fmc1602.low.25,
                     overall.fit10.fpc1601.low.25, overall.fit10.fpc1602.low.25)
gici10.body.low.50=list(overall.fit10.fmc1602.low.50,
                     overall.fit10.fpc1601.low.50, overall.fit10.fpc1602.low.50)
gici10.body.low.75=list(overall.fit10.fmc1602.low.75,
                     overall.fit10.fpc1601.low.75, overall.fit10.fpc1602.low.75)



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

gici10.plot$lab=as.factor(rep(c(rep("All conditions",4),rep("High risk conditions",3),
                                rep("Low risk conditions",3)),3))

gici10.plot$ind=as.factor(rep(c("F- Coliphage (EPA 1601)","F- Coliphage (EPA 1602)",
                                "F+ Coliphage (EPA 1601)","F+ Coliphage (EPA 1602)", 
                                "F- Coliphage (EPA 1602)",
                                "F+ Coliphage (EPA 1601)","F+ Coliphage (EPA 1602)",
                                "F- Coliphage (EPA 1602)",
                                "F+ Coliphage (EPA 1601)","F+ Coliphage (EPA 1602)"),3))


gici10.plot$type=c(rep(">25% of samples detected coliphage",10),
                             rep(">50% of samples detected coliphage",10),
                             rep(">75% of samples detected coliphage",10))


# order beach factor
gici10.plot$ind.f=factor(gici10.plot$ind,
                         levels=c("F- Coliphage (EPA 1601)","F- Coliphage (EPA 1602)",
                                  "F+ Coliphage (EPA 1601)","F+ Coliphage (EPA 1602)"))

# order panels
gici10.plot$lab.f=factor(gici10.plot$lab,levels=c("High risk conditions",
                                                  "Low risk conditions","All conditions"))

# order indicators
# order indicators
gici10.plot$type.f=factor(gici10.plot$type,
                          levels=c(">75% of samples detected coliphage",
                                   ">50% of samples detected coliphage",
                                   ">25% of samples detected coliphage"))



pdf("~/Documents/CRG/coliphage/results/figures/forestplots-threshold.pdf",height=5,width=14)
ggplot(gici10.plot,aes(x=lab.f,y=pt.est))+
  geom_errorbar(aes(ymin=lb,ymax=ub,color=type.f),width=0.3,
                position=position_dodge(width=0.4))+
  facet_wrap(~ind.f,ncol=4)+
  geom_point(aes(color=type.f,shape=type.f),position=position_dodge(width=0.4))+
  coord_flip()+ geom_hline(yintercept=1,linetype="dotted")+
  scale_color_manual("",breaks=rev(levels(gici10.plot$type.f)),values=c("#048F2E","#0066FF","black"))+
  scale_shape_manual("",breaks=rev(levels(gici10.plot$type.f)),values=c(15,17,19))+  
  scale_y_log10(breaks=c(0.1,0.2,0.5,1,2,5,10),limits=c(0.08,13))+theme_complete_bw()+
  ylab("Cumulative incidence ratio  (95% CI)")+xlab("")
dev.off()



