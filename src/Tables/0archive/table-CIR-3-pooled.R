##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file makes a table with regression output
# For a table with the combined and interaction
# results for body immersion
##########################################


rm(list=ls())

setwd("~/Dropbox/Coliphage/")  

load("results/rawoutput/regress-3day-body.Rdata")
load("results/rawoutput/regress-3day-body-entero.Rdata")


# ------------------------------------------------
# combine output
# ------------------------------------------------

# n's --------------------------------------------
# pooled n's
gici3.n.body.pool=c(all.n3.fmc1601,all.n3.fmc1602,
                     all.n3.fpc1601,all.n3.fpc1602)

gici3.n.body.high=c(all.n3.fmc1602.high,
                     all.n3.fpc1601.high,all.n3.fpc1602.high)
  
gici3.n.body.low=c(all.n3.fmc1602.low,
                  all.n3.fpc1601.low,all.n3.fpc1602.low)

entero.n=c(all.n3.entero35,all.n3.entero35.high,all.n3.entero35.low)


# estimates pooled across beach  --------------------------------------------
gici3.body.pool=list(overall.fit3.fmc1601, overall.fit3.fmc1602,
                      overall.fit3.fpc1601, overall.fit3.fpc1602)

gici3.body.pool.high=list(overall.fit3.fmc1602.high,
                      overall.fit3.fpc1601.high, overall.fit3.fpc1602.high)

gici3.body.pool.low=list(overall.fit3.fmc1602.low,
                      overall.fit3.fpc1601.low, overall.fit3.fpc1602.low)

entero.pool=list(overall.fit3.entero,overall.fit3.entero.low,
                 overall.fit3.entero.high)
  
# ------------------------------------------------
# function to make table row with exponentiated point estimate
# and 95% ci in parentheses
# ------------------------------------------------
mkrow.pool=function(out){
  pt.est=out[2,1]
  lb=pt.est-qnorm(.975)*out[2,2]
  ub=pt.est+qnorm(.975)*out[2,2]
  paste(sprintf("%0.2f",exp(pt.est))," (",
        sprintf("%0.2f",exp(lb)),",",
        sprintf("%0.2f",exp(ub)), ")",sep="")
}

# ------------------------------------------------
# convert results into table format
# pooled results
# ------------------------------------------------
# results pooled across beach
gici3.tab.body.pool=data.frame(combined=unlist(lapply(gici3.body.pool,mkrow.pool)))
gici3.tab.body.pool.high=data.frame(combined=unlist(lapply(gici3.body.pool.high,mkrow.pool)))
gici3.tab.body.pool.low=data.frame(combined=unlist(lapply(gici3.body.pool.low,mkrow.pool)))
entero.tab.pool=data.frame(combined=unlist(lapply(entero.pool,mkrow.pool)))
entero.tab.pool$combined=as.character(entero.tab.pool$combined)

gici3.tab.body.pool.high$combined=as.character(gici3.tab.body.pool.high$combined)
gici3.tab.body.pool.high=rbind("",gici3.tab.body.pool.high)
gici3.tab.body.pool.low$combined=as.character(gici3.tab.body.pool.low$combined)
gici3.tab.body.pool.low=rbind("",gici3.tab.body.pool.low)

gici3.pool=data.frame(cbind(gici3.n.body.pool,gici3.tab.body.pool,
                      c("",gici3.n.body.high),gici3.tab.body.pool.high,
                      c("",gici3.n.body.low),gici3.tab.body.pool.low))
colnames(gici3.pool)=NULL
gici3.pool[,1]=as.character(gici3.pool[,1])
gici3.pool[,2]=as.character(gici3.pool[,2])
gici3.pool[,3]=as.character(gici3.pool[,3])
gici3.pool[,4]=as.character(gici3.pool[,4])
gici3.pool[,5]=as.character(gici3.pool[,5])
gici3.pool[,6]=as.character(gici3.pool[,6])

entero=c(entero.n[1],entero.tab.pool[1,],entero.n[2],entero.tab.pool[2,],
    entero.n[3],entero.tab.pool[3,])

gici3.pool=rbind(gici3.pool,entero)
lab=c("F- Coliphage (EPA 1601)","F- Coliphage (EPA 1602)",
      "F+ Coliphage (EPA 1601)","F+ Coliphage (EPA 1602)",
      "Enterococcus")
gici3.pool=cbind(lab,gici3.pool)
save(gici3.pool,file="~/Dropbox/Coliphage/Results/Tables/CIR-pool-3.RData")






