##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file makes a table with regression output
# For a table with the combined and interaction
# results for body immersion

# Joint indicators

# Version with:
# 1 = P(Y|coliphage absent, entero<35)
# 2 = P(Y|coliphage present, entero<35)
# 3 = P(Y|coliphage absent, entero>35)
# 4 = P(Y|coliphage present, entero>35)
##########################################


rm(list=ls())

setwd("~/Dropbox/Coliphage/")  

load("~/dropbox/coliphage/results/rawoutput/regress-10day-body-joint-unadj.Rdata")


# ------------------------------------------------
# combine output
# ------------------------------------------------

# n's --------------------------------------------
gici10.n.joint=c(all.n10.fmc1601,all.n10.fmc1602,
                     all.n10.fpc1601,all.n10.fpc1602)

gici10.n.joint.high=c(all.n10.fmc1602.high,
                     all.n10.fpc1601.high,all.n10.fpc1602.high)

gici10.n.joint.low=c(all.n10.fmc1602.low,
                    all.n10.fpc1601.low,all.n10.fpc1602.low)


# estimates  --------------------------------------------
gici10.joint.int=list(overall.fit10.fmc1601.int, overall.fit10.fmc1602.int,
    overall.fit10.fpc1601.int, overall.fit10.fpc1602.int)

gici10.joint.high.int=list(overall.fit10.fmc1602.high.int,
    overall.fit10.fpc1601.high.int, overall.fit10.fpc1602.high.int)

gici10.joint.low.int=list(overall.fit10.fmc1602.low.int,
    overall.fit10.fpc1601.low.int, overall.fit10.fpc1602.low.int)


# ------------------------------------------------
# function to make table row with exponentiated point estimate
# and 95% ci in parentheses

# first row is RR if entero<35 & coli present
# second row is RR if enter>35 & coli not present
# third row is RR if entero>35 & coli present
# ------------------------------------------------
mkjoint=function(out){
  pt.est.1.0=out[2,1]
  lb.1.0=pt.est.1.0-qnorm(.975)*out[2,2]
  ub.1.0=pt.est.1.0+qnorm(.975)*out[2,2]
  row.1.0=paste(sprintf("%0.2f",exp(pt.est.1.0))," (",
        sprintf("%0.2f",exp(lb.1.0)),",",
        sprintf("%0.2f",exp(ub.1.0)), ")",sep="")
  
  pt.est.0.1=out[3,1]
  lb.0.1=pt.est.0.1-qnorm(.975)*out[3,2]
  ub.0.1=pt.est.0.1+qnorm(.975)*out[3,2]
  row.0.1=paste(sprintf("%0.2f",exp(pt.est.0.1))," (",
        sprintf("%0.2f",exp(lb.0.1)),",",
        sprintf("%0.2f",exp(ub.0.1)), ")",sep="")
  
  pt.est.1.1=out[4,1]
  lb.1.1=pt.est.1.1-qnorm(.975)*out[4,2]
  ub.1.1=pt.est.1.1+qnorm(.975)*out[4,2]
  row.1.1=paste(sprintf("%0.2f",exp(pt.est.1.1))," (",
        sprintf("%0.2f",exp(lb.1.1)),",",
        sprintf("%0.2f",exp(ub.1.1)), ")",sep="")  
  
  return(c("ref",row.1.0,row.0.1,row.1.1))
}

# ------------------------------------------------
# convert results into table format
# pooled results
# ------------------------------------------------
# results pooled across beach
gici10.joint.tab=mapply(mkjoint,gici10.joint.int)
gici10.joint.tab=c(gici10.joint.tab[,1],gici10.joint.tab[,2],
                   gici10.joint.tab[,3],gici10.joint.tab[,4])
gici10.joint.tab.high=mapply(mkjoint,gici10.joint.high.int)
gici10.joint.tab.high=c(gici10.joint.tab.high[,1],gici10.joint.tab.high[,2],
                   gici10.joint.tab.high[,3])
gici10.joint.tab.low=mapply(mkjoint,gici10.joint.low.int)
gici10.joint.tab.low=c(gici10.joint.tab.low[,1],gici10.joint.tab.low[,2],
                   gici10.joint.tab.low[,3])

gici10.joint.tab.high=c(rep("",4),gici10.joint.tab.high)
gici10.joint.tab.low=c(rep("",4),gici10.joint.tab.low)


gici10.n.joint.tab=format(gici10.n.joint,scientific=FALSE,big.mark=",")
gici10.n.joint.low.tab=format(gici10.n.joint.low,scientific=FALSE,big.mark=",")
gici10.n.joint.low.tab=c(rep("",4),gici10.n.joint.low.tab)
gici10.n.joint.high.tab=format(gici10.n.joint.high,scientific=FALSE,big.mark=",")
gici10.n.joint.high.tab=c(rep("",4),gici10.n.joint.high.tab)

gici10.joint.tab.all=data.frame(cbind(gici10.n.joint.tab,gici10.joint.tab,
                                      gici10.n.joint.low.tab,gici10.joint.tab.low,
                                      gici10.n.joint.high.tab,gici10.joint.tab.high))

colnames(gici10.joint.tab.all)=c("n","All conditions","n","Low risk conditions",
                                 "n","High risk conditions")

ind=c("F- Coliphage (EPA 1601)",rep("",3),
      "F- Coliphage (EPA 1602)",rep("",3),
      "F+ Coliphage (EPA 1601)",rep("",3),
      "F+ Coliphage (EPA 1602)",rep("",3))

labC=rep(c("Not detected", "Detected","Not detected","Detected"),4)

labE=rep(c("<35 ml/100 CFU","<35 ml/100 CFU",">35 ml/100 CFU",">35 ml/100 CFU"),4)

gici10.joint.tab.all=cbind(ind,labC,labE,gici10.joint.tab.all)

# removing result with sparse data
gici10.joint.tab.all[["High risk conditions"]]=as.character(gici10.joint.tab.all[["High risk conditions"]])
gici10.joint.tab.all[11,9]="(not estimated)"
  
#save(gici10.joint.tab.all,file="~/Dropbox/Coliphage/Results/Tables/CIR-joint-10.RData")







