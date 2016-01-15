##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file makes a table with regression output
# For a table with the combined and interaction
# results for body immersion
##########################################


rm(list=ls())

setwd("~/Dropbox/Coliphage/")  

load("results/rawoutput/regress-10day-body.Rdata")
load("results/rawoutput/regress-10day-head.Rdata")
load("results/rawoutput/regress-10day-swall.Rdata")

# ------------------------------------------------
# combine output
# ------------------------------------------------

# n's --------------------------------------------
# pooled n's
gici10.n.body.pool=c(all.n10.fmc1601,all.n10.fmc1602,
                     all.n10.fpc1601,all.n10.fpc1602)

gici10.n.head.pool=c(all.n10.fmc1601head,all.n10.fmc1602head,
                     all.n10.fpc1601head,all.n10.fpc1602head)

gici10.n.swall.pool=c(all.n10.fmc1601swall,all.n10.fmc1602swall,
                      all.n10.fpc1601swall,all.n10.fpc1602swall)


# estimates pooled across beach  --------------------------------------------
gici10.body.pool=list(overall.fit10.fmc1601, overall.fit10.fmc1602,
                      overall.fit10.fpc1601, overall.fit10.fpc1602)

gici10.head.pool=list(overall.fit10.fmc1601head, overall.fit10.fmc1602head,
                      overall.fit10.fpc1601head, overall.fit10.fpc1602head)

gici10.swall.pool=list(overall.fit10.fmc1601swall, overall.fit10.fmc1602swall,
                       overall.fit10.fpc1601swall, overall.fit10.fpc1602swall)

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
gici10.tab.body.pool=data.frame(combined=unlist(lapply(gici10.body.pool,mkrow.pool)))
gici10.tab.body.pool.high=data.frame(combined=unlist(lapply(gici10.body.pool.high,mkrow.pool)))
gici10.tab.body.pool.low=data.frame(combined=unlist(lapply(gici10.body.pool.low,mkrow.pool)))

gici10.tab.head.pool=data.frame(combined=unlist(lapply(gici10.head.pool,mkrow.pool)))
gici10.tab.swall.pool=data.frame(combined=unlist(lapply(gici10.swall.pool,mkrow.pool)))

gici10.exp.pool=cbind(gici10.n.body.pool,gici10.tab.body.pool,
                      ,gici10.tab.body.pool.high,
                      gici10.n.swall.pool,gici10.tab.body.pool.low)


#save(gici10.exp.pool,file="~/Dropbox/Coliphage/Results/Tables/reg-exposure-10.RData")







