##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file creates a table comparing AIC
# values for different models
##########################################

rm(list=ls())

load("~/Documents/CRG/coliphage/results/rawoutput/regress-10day-body-pool.Rdata")
load("~/Documents/CRG/coliphage/results/rawoutput/regress-10day-body-entero-pool.Rdata")
load("~/Documents/CRG/coliphage/results/rawoutput/regress-10day-body-joint-pool.Rdata")

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

aic.tab$cond=as.factor(c(rep("All conditions",6),rep("Low risk conditions",6),
                         rep("High risk conditions",6)))

aic.tab.out=cbind(aic.tab[aic.tab$cond=="All conditions",1:3],
        aic.tab[aic.tab$cond=="Low risk conditions",1],
        aic.tab[aic.tab$cond=="High risk conditions",1])

colnames(aic.tab.out)=c("all","ind","type","low","high")
aic.tab.out=aic.tab.out[,c("ind","type","all","low","high")]

aic.tab.out=aic.tab.out[order(aic.tab.out$ind,aic.tab.out$type),]

save(aic.tab.out,file="~/Documents/CRG/coliphage/results/tables/aic-pool.RData")


