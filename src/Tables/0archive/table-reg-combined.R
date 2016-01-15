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
load("results/rawoutput/regress-10day-body.Rdata")

# ------------------------------------------------
# combine output
# ------------------------------------------------
# n's
gici3.n.comb=c(do.n3.fmc1601, mb.n3.fmc1601,
				av.n3.fmc1602, do.n3.fmc1602,
				av.n3.fpc1601, do.n3.fpc1601, fa.n3.fpc1601, ma.n3.fpc1601, mb.n3.fpc1601,
                av.n3.fpc1602, do.n3.fpc1602)

gici10.n.comb=c(do.n10.fmc1601, mb.n10.fmc1601,
				av.n10.fmc1602, do.n10.fmc1602,
				av.n10.fpc1601, do.n10.fpc1601, fa.n10.fpc1601, ma.n10.fpc1601, mb.n10.fpc1601,
                av.n10.fpc1602, do.n10.fpc1602)

# pooled across berm/groundwater
gici3.comb=list(dofit3.fmc1601, mbfit3.fmc1601,
				avfit3.fmc1602, dofit3.fmc1602,
				avfit3.fpc1601, dofit3.fpc1601, fafit3.fpc1601, mafit3.fpc1601, mbfit3.fpc1601,
                avfit3.fpc1602, dofit3.fpc1602)

gici10.comb=list(dofit10.fmc1601, mbfit10.fmc1601,
        avfit10.fmc1602, dofit10.fmc1602,
        avfit10.fpc1601, dofit10.fpc1601, fafit10.fpc1601, mafit10.fpc1601, mbfit10.fpc1601,
                avfit10.fpc1602, dofit10.fpc1602)

# berm open/ groundwater above median
gici3.int1=list(avfit3.fmc1602.gw.1, dofit3.fmc1602.berm.1,avfit3.fpc1601.gw.1,
        dofit3.fpc1601.berm.1,avfit3.fpc1602.gw.1,  dofit3.fpc1602.berm.1,mafit3.fpc1601.berm.1)

gici10.int1=list(avfit10.fmc1602.gw.1, dofit10.fmc1602.berm.1,avfit10.fpc1601.gw.1,
        dofit10.fpc1601.berm.1,avfit10.fpc1602.gw.1,  dofit10.fpc1602.berm.1,mafit10.fpc1601.berm.1)

# berm closed / groundwater below median
gici3.int0=list(avfit3.fmc1602.gw.0,dofit3.fmc1602.berm.0, avfit3.fpc1601.gw.0, 
        dofit3.fpc1601.berm.0, avfit3.fpc1602.gw.0, dofit3.fpc1602.berm.0,mafit3.fpc1601.berm.0)

gici10.int0=list(avfit10.fmc1602.gw.0,dofit10.fmc1602.berm.0, avfit10.fpc1601.gw.0, 
        dofit10.fpc1601.berm.0, avfit10.fpc1602.gw.0, dofit10.fpc1602.berm.0,mafit10.fpc1601.berm.0)

# ------------------------------------------------
# function to make table row with exponentiated point estimate
# and 95% ci in parentheses
# ------------------------------------------------
mkrow=function(out){
  pt.est=out$fit[2,1]
  lb=pt.est-qnorm(.975)*out$fit[2,2]
  ub=pt.est+qnorm(.975)*out$fit[2,2]
  paste(sprintf("%0.2f",exp(pt.est))," (",
        sprintf("%0.2f",exp(lb)),",",
        sprintf("%0.2f",exp(ub)), ")",sep="")
}

# ------------------------------------------------
# convert results into table format
# ------------------------------------------------

# 3-day gici
gici3.tab.comb=data.frame(combined=unlist(lapply(gici3.comb,mkrow)))

label=c("Doheny","Mission Bay",
        "Avalon","Doheny",
        "Avalon","Doheny","Fairhope","Malibu","Mission Bay",
        "Avalon","Doheny")

gici3.comb.print=cbind(label,gici3.n.comb,gici3.tab.comb)
gici3.comb.print$label=as.character(gici3.comb.print$label)

colnames(gici3.comb.print)=c("label","N","comb")

# stratify by berm/groundwater
gici3.tab.int1=data.frame(int1=unlist(lapply(gici3.int1,mkrow)))
gici3.tab.int0=data.frame(int1=unlist(lapply(gici3.int0,mkrow)))

# add blank rows for other beaches
gici3.tab.int1.print=c("--","--",as.character(gici3.tab.int1[c(1:4),1]),
      "--",as.character(gici3.tab.int1[5,1]),"--",as.character(gici3.tab.int1[c(6:7),1]))
# add blank rows for other beaches
gici3.tab.int0.print=c("--","--",as.character(gici3.tab.int0[c(1:4),1]),
      "--",as.character(gici3.tab.int0[5,1]),"--",as.character(gici3.tab.int0[c(6:7),1]))


gici3.tab.all=cbind(gici3.comb.print,gici3.tab.int0.print,gici3.tab.int1.print)

# 10-day gici
gici10.tab.comb=data.frame(combined=unlist(lapply(gici10.comb,mkrow)))


gici10.comb.print=cbind(label,gici10.n.comb,gici10.tab.comb)
gici10.comb.print$label=as.character(gici10.comb.print$label)

colnames(gici10.comb.print)=c("label","N","comb")

# stratify by berm/groundwater
gici10.tab.int1=data.frame(int1=unlist(lapply(gici10.int1,mkrow)))
gici10.tab.int0=data.frame(int1=unlist(lapply(gici10.int0,mkrow)))

# add blank rows for other beaches
gici10.tab.int1.print=c("--","--",as.character(gici10.tab.int1[c(1:4),1]),
   "--",as.character(gici10.tab.int1[5,1]),"--",as.character(gici10.tab.int1[c(6:7),1]))

# add blank rows for other beaches
gici10.tab.int0.print=c("--","--",as.character(gici10.tab.int0[c(1:4),1]),
   "--",as.character(gici10.tab.int0[5,1]),"--",as.character(gici10.tab.int0[c(6:7),1]))

gici10.tab.all=cbind(gici10.comb.print,gici10.tab.int0.print,gici10.tab.int1.print)

save(gici3.tab.all,gici10.tab.all,file="~/Dropbox/Coliphage/Results/Tables/reg-combined.RData")






