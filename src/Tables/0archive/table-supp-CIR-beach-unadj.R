##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file makes a table with regression output
# For a table with the combined and interaction
# results for body immersion

# unadjusted analyses
##########################################

rm(list=ls())
setwd("~/Dropbox/Coliphage/")  

load("results/rawoutput/regress-3day-body-beach.Rdata")
load("results/rawoutput/regress-10day-body-beach.Rdata")

# ------------------------------------------------
# combine output
# ------------------------------------------------
# n's
gici3.n.comb=c(do.n3.fmc1601, mb.n3.fmc1601,
				av.n3.fmc1602, do.n3.fmc1602,
				av.n3.fpc1601, do.n3.fpc1601, fa.n3.fpc1601, ma.n3.fpc1601, mb.n3.fpc1601,
        av.n3.fpc1602, do.n3.fpc1602)

gici3.n.int0=c(av.n3.fmc1602.int0, do.n3.fmc1602.int0,
        av.n3.fpc1601.int0, do.n3.fpc1601.int0, ma.n3.fpc1601.int0,
        av.n3.fpc1602.int0, do.n3.fpc1602.int0)

gici3.n.int1=c(av.n3.fmc1602.int1, do.n3.fmc1602.int1,
        av.n3.fpc1601.int1, do.n3.fpc1601.int1, ma.n3.fpc1601.int1,
        av.n3.fpc1602.int1, do.n3.fpc1602.int1)

gici10.n.comb=c(do.n10.fmc1601, mb.n10.fmc1601,
				av.n10.fmc1602, do.n10.fmc1602,
				av.n10.fpc1601, do.n10.fpc1601, fa.n10.fpc1601, ma.n10.fpc1601, mb.n10.fpc1601,
        av.n10.fpc1602, do.n10.fpc1602)

gici10.n.int0=c(av.n10.fmc1602.int0, do.n10.fmc1602.int0,
        av.n10.fpc1601.int0, do.n10.fpc1601.int0,  ma.n10.fpc1601.int0, 
        av.n10.fpc1602.int0, do.n10.fpc1602.int0)

gici10.n.int1=c(av.n10.fmc1602.int1, do.n10.fmc1602.int1,
        av.n10.fpc1601.int1, do.n10.fpc1601.int1, ma.n10.fpc1601.int1,
        av.n10.fpc1602.int1, do.n10.fpc1602.int1)

# estimates pooled across conditions
gici3.comb=list(dofit3.fmc1601, mbfit3.fmc1601,
				avfit3.fmc1602, dofit3.fmc1602,
				avfit3.fpc1601, dofit3.fpc1601, fafit3.fpc1601, mafit3.fpc1601, mbfit3.fpc1601,
                avfit3.fpc1602, dofit3.fpc1602)

gici10.comb=list(dofit10.fmc1601, mbfit10.fmc1601,
        avfit10.fmc1602, dofit10.fmc1602,
        avfit10.fpc1601, dofit10.fpc1601, fafit10.fpc1601, mafit10.fpc1601, mbfit10.fpc1601,
                avfit10.fpc1602, dofit10.fpc1602)

# estimates when berm open/ groundwater above median
gici3.int1=list(avfit3.fmc1602.gw.1, dofit3.fmc1602.berm.1,avfit3.fpc1601.gw.1,
        dofit3.fpc1601.berm.1,avfit3.fpc1602.gw.1,  dofit3.fpc1602.berm.1,mafit3.fpc1601.berm.1)

gici10.int1=list(avfit10.fmc1602.gw.1, dofit10.fmc1602.berm.1,avfit10.fpc1601.gw.1,
        dofit10.fpc1601.berm.1,avfit10.fpc1602.gw.1,  dofit10.fpc1602.berm.1,mafit10.fpc1601.berm.1)

# estimates when berm closed / groundwater below median
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

# 3-day gici ----------------------------------------------------
# pooled results
gici3.tab.comb=data.frame(combined=unlist(lapply(gici3.comb,mkrow)))
gici3.tab.comb$combined=as.character(gici3.tab.comb$combined)

# stratified results
gici3.tab.int1=data.frame(int1=unlist(lapply(gici3.int1,mkrow)))
gici3.tab.int0=data.frame(int0=unlist(lapply(gici3.int0,mkrow)))
gici3.tab.int1$int1=as.character(gici3.tab.int1$int1)
gici3.tab.int0$int0=as.character(gici3.tab.int0$int0)

# combine the pooled and high/low risk estimates into one object
gici3.tab.all=c(gici3.tab.comb$combined[1:3],gici3.tab.int0$int0[1],gici3.tab.int1$int1[1],
                gici3.tab.comb$combined[4],gici3.tab.int0$int0[2],gici3.tab.int1$int1[2],
                gici3.tab.comb$combined[5],gici3.tab.int0$int0[3],gici3.tab.int1$int1[3],
                gici3.tab.comb$combined[6],gici3.tab.int0$int0[4],gici3.tab.int1$int1[4],
                gici3.tab.comb$combined[7:8],gici3.tab.int0$int0[5],gici3.tab.int1$int1[5],
                gici3.tab.comb$combined[9:10],gici3.tab.int0$int0[6],gici3.tab.int1$int1[6],
                gici3.tab.comb$combined[11],gici3.tab.int0$int0[7],gici3.tab.int1$int1[7])

gici3.n.all=c(gici3.n.comb[1:3],gici3.n.int0[1],gici3.n.int1[1],
                gici3.n.comb[4],gici3.n.int0[2],gici3.n.int1[2],
                gici3.n.comb[5],gici3.n.int0[3],gici3.n.int1[3],
                gici3.n.comb[6],gici3.n.int0[4],gici3.n.int1[4],
                gici3.n.comb[7:8],gici3.n.int0[5],gici3.n.int1[5],
                gici3.n.comb[9:10],gici3.n.int0[6],gici3.n.int1[6],
                gici3.n.comb[11],gici3.n.int0[7],gici3.n.int1[7])
  
label=c("Doheny","Mission Bay",
        "Avalon","","","Doheny","","",
        "Avalon","","","Doheny","","","Fairhope","Malibu","","","Mission Bay",
        "Avalon","","","Doheny","","")

gici3.comb.print=data.frame(cbind(label,gici3.n.all,gici3.tab.all))
gici3.comb.print$label=as.character(gici3.comb.print$label)
gici3.comb.print$cond=c(rep("All conditions",3),"Groundwater below median flow",
      "Groundwater above median flow","All conditions","Berm closed","Berm open",
      "All conditions","Groundwater below median flow","Groundwater above median flow",
      "All conditions","Berm closed","Berm open","All conditions","All conditions",
      "Berm closed","Berm open","All conditions","All conditions","Groundwater below median flow",
      "Groundwater above median flow","All conditions","Berm closed","Berm open")

gici3.comb.print=gici3.comb.print[,c(1,4,2:3)]

colnames(gici3.comb.print)=c("label","cond","N","CIR")

# 10-day gici ----------------------------------------------------
gici10.tab.comb=data.frame(combined=unlist(lapply(gici10.comb,mkrow)))
gici10.tab.comb$combined=as.character(gici10.tab.comb$combined)

# stratify by berm/groundwater
gici10.tab.int1=data.frame(int1=unlist(lapply(gici10.int1,mkrow)))
gici10.tab.int0=data.frame(int0=unlist(lapply(gici10.int0,mkrow)))
gici10.tab.int1$int1=as.character(gici10.tab.int1$int1)
gici10.tab.int0$int0=as.character(gici10.tab.int0$int0)

# combine the pooled and high/low risk estimates into one object
gici10.tab.all=c(gici10.tab.comb$combined[1:3],gici10.tab.int0$int0[1],gici10.tab.int1$int1[1],
                gici10.tab.comb$combined[4],gici10.tab.int0$int0[2],gici10.tab.int1$int1[2],
                gici10.tab.comb$combined[5],gici10.tab.int0$int0[3],gici10.tab.int1$int1[3],
                gici10.tab.comb$combined[6],gici10.tab.int0$int0[4],gici10.tab.int1$int1[4],
                gici10.tab.comb$combined[7:8],gici10.tab.int0$int0[5],gici10.tab.int1$int1[5],
                gici10.tab.comb$combined[9:10],gici10.tab.int0$int0[6],gici10.tab.int1$int1[6],
                gici10.tab.comb$combined[11],gici10.tab.int0$int0[7],gici10.tab.int1$int1[7])

gici10.n.all=c(gici10.n.comb[1:3],gici10.n.int0[1],gici10.n.int1[1],
                gici10.n.comb[4],gici10.n.int0[2],gici10.n.int1[2],
                gici10.n.comb[5],gici10.n.int0[3],gici10.n.int1[3],
                gici10.n.comb[6],gici10.n.int0[4],gici10.n.int1[4],
                gici10.n.comb[7:8],gici10.n.int0[5],gici10.n.int1[5],
                gici10.n.comb[9:10],gici10.n.int0[6],gici10.n.int1[6],
                gici10.n.comb[11],gici10.n.int0[7],gici10.n.int1[7])

gici10.comb.print=data.frame(cbind(label,gici10.n.all,gici10.tab.all))
gici10.comb.print$label=as.character(gici10.comb.print$label)

gici10.comb.print$cond=c(rep("All conditions",3),"Groundwater below median flow",
      "Groundwater above median flow","All conditions","Berm closed","Berm open",
      "All conditions","Groundwater below median flow","Groundwater above median flow",
      "All conditions","Berm closed","Berm open","All conditions","All conditions",
      "Berm closed","Berm open","All conditions","All conditions","Groundwater below median flow",
      "Groundwater above median flow","All conditions","Berm closed","Berm open")

gici10.comb.print=gici10.comb.print[,c(1,4,2:3)]

colnames(gici10.comb.print)=c("label","cond","N","CIR")

save(gici3.comb.print,gici10.comb.print,file="~/Dropbox/Coliphage/Results/Tables/reg-combined-unadj.RData")






