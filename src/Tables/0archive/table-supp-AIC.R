##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file creates a table comparing AIC
# values for different models
##########################################

rm(list=ls())

load("~/dropbox/coliphage/results/rawoutput/regress-10day-body.Rdata")
load("~/dropbox/coliphage/results/rawoutput/regress-10day-body-entero.Rdata")
load("~/dropbox/coliphage/results/rawoutput/regress-10day-body-joint.Rdata")

# ------------------------------------------------
# combine output
# ------------------------------------------------
aic=c(aic.fmc1601,aic.fmc1602,aic.fpc1601,aic.fpc1602,
      aic.entero.fmc1601,aic.entero.fmc1602,aic.entero.fpc1601,aic.entero.fpc1602,
      aic.fmc1601.int,aic.fmc1602.int,aic.fpc1601.int,aic.fpc1602.int)
aic.low=c(NA,aic.fmc1602.low,aic.fpc1601.low,aic.fpc1602.low,
          NA,aic.entero.low.fmc1602,aic.entero.low.fpc1601,aic.entero.low.fpc1602,
          NA,aic.fmc1602.low.int,aic.fpc1601.low.int,aic.fpc1602.low.int)
aic.high=c(NA,aic.fmc1602.high,aic.fpc1601.high,aic.fpc1602.high,
           NA,aic.entero.high.fmc1602,aic.entero.high.fpc1601,aic.entero.high.fpc1602,           
           NA,aic.fmc1602.high.int,aic.fpc1601.high.int,aic.fpc1602.high.int)

aic.tab=data.frame(c(aic,aic.low,aic.high))

aic.tab$ind=c(rep(c("F- Coliphage (EPA 1601)","F- Coliphage (EPA 1602)",
                    "F+ Coliphage (EPA 1601)", "F+ Coliphage (EPA 1602)"),9))

aic.tab$type=as.factor(c(rep("Coliphage only",4),rep("Enterococcus Only",4),rep("Joint indicator",4),
                           rep("Coliphage only",4),rep("Enterococcus Only",4),rep("Joint indicator",4),
                           rep("Coliphage only",4),rep("Enterococcus Only",4),rep("Joint indicator",4)))
aic.tab$type.f=factor(aic.tab$type,levels=c("Coliphage only","Joint indicator",
                                            "Enterococcus Only"))

aic.tab$type=NULL

aic.tab$cond=as.factor(c(rep("All conditions",12),rep("Low risk conditions",12),
                         rep("High risk conditions",12)))

aic.tab.out=cbind(aic.tab[aic.tab$cond=="All conditions",1:3],
        aic.tab[aic.tab$cond=="Low risk conditions",1],
        aic.tab[aic.tab$cond=="High risk conditions",1])

colnames(aic.tab.out)=c("all","ind","type.f","low","high")
aic.tab.out=aic.tab.out[,c("ind","type.f","all","low","high")]

aic.tab.out=aic.tab.out[order(aic.tab.out$ind,aic.tab.out$type.f),]

save(aic.tab.out,file="~/dropbox/coliphage/results/tables/aic.RData")


