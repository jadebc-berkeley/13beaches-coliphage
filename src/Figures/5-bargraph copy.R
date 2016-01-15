##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 11/10/15

# Bar graph comparing risk of illness
# under different conditions
##########################################


rm(list=ls())
library(ggplot2)
library(plyr)
library(scales)
library(grid)

source("~/dropbox/coliphage/programs/figures/theme_complete_bw.R")

load("~/dropbox/coliphage/results/rawoutput/bargraph-output.Rdata")

bargraph=data.frame(pt.est=c(ns.bs$bootest, s.bs$bootest, fmc.bs$bootest, 
    fpc.bs$bootest, fmc.ent.bs$bootest, fpc.ent.bs$bootest))

bargraph$lb=c(ns.bs$boot95lb, s.bs$boot95lb, fmc.bs$boot95lb, 
    fpc.bs$boot95lb, fmc.ent.bs$boot95lb, fpc.ent.bs$boot95lb)

bargraph$ub=c(ns.bs$boot95ub, s.bs$boot95ub, fmc.bs$boot95ub, 
    fpc.bs$boot95ub, fmc.ent.bs$boot95ub, fpc.ent.bs$boot95ub)

bargraph$lab=c("Non-swimmers","Swimmers not \nexposed to coliphage",
               "Swimmers exposed to \nsomatic coliphage",
               "Swimmers exposed to \nmale-specific coliphage",
               "Swimmers exposed to \nsomatic coliphage and \nEnterococcus > 35 CFU/100 ml",
               "Swimmers exposed to \nmale-specific coliphage \nand Enterococcus > \n35 CFU/100 ml")

# order label
bargraph$lab.f=factor(bargraph$lab,
        levels=c("Non-swimmers","Swimmers not \nexposed to coliphage",
                 "Swimmers exposed to \nsomatic coliphage",
                 "Swimmers exposed to \nmale-specific coliphage",
                 "Swimmers exposed to \nsomatic coliphage and \nEnterococcus > 35 CFU/100 ml",
                 "Swimmers exposed to \nmale-specific coliphage \nand Enterococcus > \n35 CFU/100 ml"))
        
        

bargraph$lab=as.factor(bargraph$lab)

bargraph$pt.est=bargraph$pt.est*100
bargraph$lb=bargraph$lb*100
bargraph$ub=bargraph$ub*100

bargraph$pt.est.lab=paste(as.numeric(sprintf("%0.1f",bargraph$pt.est)),"%",sep="")

somatic=bargraph[c(1:3,5),]
male=bargraph[c(1:2,4,6),]


somatic.plot=ggplot(somatic,aes(x=lab.f,y=pt.est))+geom_bar(fill="grey",stat="identity",width=.6)+
  ylab("Probability of Gastrointestinal Illness (%)")+theme_complete_bw()+
  geom_errorbar(aes(ymin=lb,ymax=ub),width=0.25)+
  geom_text(data=somatic, mapping=aes(x=lab,y=pt.est,label=pt.est.lab),vjust=2.7)+
  ylim(c(0,8.5))+xlab("")+ggtitle("Somatic Coliphage")

male.plot=ggplot(male,aes(x=lab.f,y=pt.est))+geom_bar(fill="grey",stat="identity",width=.6)+
  ylab("Probability of Gastrointestinal Illness (%)")+theme_complete_bw()+
  geom_errorbar(aes(ymin=lb,ymax=ub),width=0.25)+
  geom_text(data=male, mapping=aes(x=lab,y=pt.est,label=pt.est.lab),vjust=2.7)+
  ylim(c(0,8.5))+xlab("")+ggtitle("Male-Specific Coliphage")

pdf("~/dropbox/coliphage/results/figures/bargraph.pdf",height=9,width=8)
grid.arrange(somatic.plot, male.plot)
dev.off()

