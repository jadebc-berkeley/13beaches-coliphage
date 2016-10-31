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
library(gridExtra)

source("~/Documents/CRG/coliphage/13beaches-coliphage/src/figures/theme_complete_bw.R")

load("~/Documents/CRG/coliphage/results/rawoutput/bargraph-output.Rdata")

# bargraph=data.frame(pt.est=c(ns, swim, s.som, s.male, s.som.ent, s.male.ent))
# 
# bargraph$lb=c(lb.ns, lb.swim, lb.som, lb.male, lb.som.ent, lb.male.ent)
# 
# bargraph$ub=c(ub.ns, ub.swim, ub.som, ub.male, ub.som.ent, ub.male.ent)

bargraph=data.frame(pt.est=c(ns[1,1], swim[1,1], s.som[1,1], s.male[1,1], s.som.ent[1,1], s.male.ent[1,1]))
bargraph$lb=c(ns[1,2], swim[1,2], s.som[1,2], s.male[1,2], s.som.ent[1,2], s.male.ent[1,2])
bargraph$ub=c(ns[1,3], swim[1,3], s.som[1,3], s.male[1,3], s.som.ent[1,3], s.male.ent[1,3])

bargraph$lab=c("Non-swimmers","Swimmers not \nexposed to coliphage",
               "Swimmers exposed to \nsomatic coliphage",
               "Swimmers exposed to \nmale-specific coliphage",
               "Swimmers exposed to \nsomatic coliphage and \nenterococci > 35 CFU/100 ml",
               "Swimmers exposed to \nmale-specific coliphage \nand enterococci > \n35 CFU/100 ml")

# order label
bargraph$lab.f=factor(bargraph$lab,
        levels=c("Non-swimmers","Swimmers not \nexposed to coliphage",
                 "Swimmers exposed to \nsomatic coliphage",
                 "Swimmers exposed to \nmale-specific coliphage",
                 "Swimmers exposed to \nsomatic coliphage and \nenterococci > 35 CFU/100 ml",
                 "Swimmers exposed to \nmale-specific coliphage \nand enterococci > \n35 CFU/100 ml"))
        
        

bargraph$lab=as.factor(bargraph$lab)

bargraph$pt.est=bargraph$pt.est*100
bargraph$lb=bargraph$lb*100
bargraph$ub=bargraph$ub*100

bargraph$pt.est.lab=paste(round(bargraph$pt.est,digits=1),"%",sep="")


somatic=bargraph[c(1:3,5),]
male=bargraph[c(1:2,4,6),]


somatic.plot=ggplot(somatic,aes(x=lab.f,y=pt.est))+geom_bar(fill="grey",stat="identity",width=.6)+
  ylab("Probability of Gastrointestinal Illness (%)")+theme_complete_bw()+
  geom_errorbar(aes(ymin=lb,ymax=ub),width=0.25)+
  geom_text(data=somatic, mapping=aes(x=lab,y=pt.est,label=pt.est.lab),vjust=c(3,5.4,3.5,5.5))+
  ylim(c(0,10))+xlab("")+ggtitle("Somatic Coliphage")

male.plot=ggplot(male,aes(x=lab.f,y=pt.est))+geom_bar(fill="grey",stat="identity",width=.6)+
  ylab("Probability of Gastrointestinal Illness (%)")+theme_complete_bw()+
  geom_errorbar(aes(ymin=lb,ymax=ub),width=0.25)+
  geom_text(data=male, mapping=aes(x=lab,y=pt.est,label=pt.est.lab),vjust=c(3,5.4,3.5,5.5))+
  ylim(c(0,10))+xlab("")+ggtitle("Male-Specific Coliphage")

pdf("~/Documents/CRG/coliphage/results/figures/bargraph.pdf",height=9,width=8)
grid.arrange(somatic.plot, male.plot)
dev.off()

