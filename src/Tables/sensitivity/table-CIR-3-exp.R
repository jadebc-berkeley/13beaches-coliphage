##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file makes a table with regression output
# For a table with the combined and interaction
# results for body immersion

# 3-day gastrointestinal illness
##########################################


rm(list=ls())

setwd("~/Dropbox/Coliphage/")  

load("results/rawoutput/regress-3day-body-beach.Rdata")
load("results/rawoutput/regress-3day-head.Rdata")
load("results/rawoutput/regress-3day-swall.Rdata")

# ------------------------------------------------
# combine output
# ------------------------------------------------

# n's --------------------------------------------
# combined n's
gici3.n.body=c(do.n3.fmc1601, mb.n3.fmc1601,
				av.n3.fmc1602, do.n3.fmc1602,
				av.n3.fpc1601, do.n3.fpc1601, fa.n3.fpc1601, ma.n3.fpc1601, mb.n3.fpc1601,
        av.n3.fpc1602, do.n3.fpc1602)

gici3.n.head=c(do.n3.fmc1601head, mb.n3.fmc1601head,
				av.n3.fmc1602head, do.n3.fmc1602head,
				av.n3.fpc1601head, do.n3.fpc1601head, fa.n3.fpc1601head, ma.n3.fpc1601head, mb.n3.fpc1601head,
        av.n3.fpc1602head, do.n3.fpc1602head)

gici3.n.swall=c(do.n3.fmc1601swall, mb.n3.fmc1601swall,
				av.n3.fmc1602swall, do.n3.fmc1602swall,
				av.n3.fpc1601swall, do.n3.fpc1601swall, fa.n3.fpc1601swall, ma.n3.fpc1601swall, mb.n3.fpc1601swall,
        av.n3.fpc1602swall, do.n3.fpc1602swall)

# n if berm closed/groundwater below median
gici3.n.int0.body=c(av.n3.fmc1602.int0, do.n3.fmc1602.int0,
        av.n3.fpc1601.int0, do.n3.fpc1601.int0, ma.n3.fpc1601.int0,
        av.n3.fpc1602.int0, do.n3.fpc1602.int0)

gici3.n.int0.head=c(av.n3.fmc1602.int0head, do.n3.fmc1602.int0head,
        av.n3.fpc1601.int0head, do.n3.fpc1601.int0head, ma.n3.fpc1601.int0head,
        av.n3.fpc1602.int0head, do.n3.fpc1602.int0)

gici3.n.int0.swall=c(av.n3.fmc1602.int0swall, do.n3.fmc1602.int0swall,
        av.n3.fpc1601.int0swall, do.n3.fpc1601.int0swall, ma.n3.fpc1601.int0swall,
        av.n3.fpc1602.int0swall, do.n3.fpc1602.int0)

# n if berm open/groundwater above median
gici3.n.int1.body=c(av.n3.fmc1602.int1, do.n3.fmc1602.int1,
        av.n3.fpc1601.int1, do.n3.fpc1601.int1, ma.n3.fpc1601.int1, 
        av.n3.fpc1602.int1, do.n3.fpc1602.int1)

gici3.n.int1.head=c(av.n3.fmc1602.int1head, do.n3.fmc1602.int1head,
        av.n3.fpc1601.int1head, do.n3.fpc1601.int1head, ma.n3.fpc1601.int1head,
        av.n3.fpc1602.int1head, do.n3.fpc1602.int1)

gici3.n.int1.swall=c(av.n3.fmc1602.int1swall, do.n3.fmc1602.int1swall,
        av.n3.fpc1601.int1swall, do.n3.fpc1601.int1swall, ma.n3.fpc1601.int1swall,
        av.n3.fpc1602.int1swall, do.n3.fpc1602.int1)

# estimates pooled across conditions
gici3.body=list(dofit3.fmc1601, mbfit3.fmc1601,
				avfit3.fmc1602, dofit3.fmc1602,
				avfit3.fpc1601, dofit3.fpc1601, fafit3.fpc1601, mafit3.fpc1601, mbfit3.fpc1601,
                avfit3.fpc1602, dofit3.fpc1602)

gici3.head=list(dofit3.fmc1601head, mbfit3.fmc1601head,
				avfit3.fmc1602head, dofit3.fmc1602head,
				avfit3.fpc1601head, dofit3.fpc1601head, fafit3.fpc1601head, mafit3.fpc1601head, mbfit3.fpc1601head,
                avfit3.fpc1602head, dofit3.fpc1602head)

gici3.swall=list(dofit3.fmc1601swall, mbfit3.fmc1601swall,
				avfit3.fmc1602swall, dofit3.fmc1602swall,
				avfit3.fpc1601swall, dofit3.fpc1601swall, fafit3.fpc1601swall, mafit3.fpc1601swall, mbfit3.fpc1601swall,
                avfit3.fpc1602swall, dofit3.fpc1602swall)

# estimates when berm open/ groundwater above median
gici3.int1.body=list(avfit3.fmc1602.gw.1, dofit3.fmc1602.berm.1,avfit3.fpc1601.gw.1,
        dofit3.fpc1601.berm.1,avfit3.fpc1602.gw.1,  dofit3.fpc1602.berm.1,mafit3.fpc1601.berm.1)

gici3.int1.head=list(avfit3.fmc1602.gw.1.head, dofit3.fmc1602.berm.1.head,avfit3.fpc1601.gw.1.head,
        dofit3.fpc1601.berm.1.head,avfit3.fpc1602.gw.1.head,  dofit3.fpc1602.berm.1.head,mafit3.fpc1601.berm.1.head)

gici3.int1.swall=list(avfit3.fmc1602.gw.1.swall, dofit3.fmc1602.berm.1.swall,avfit3.fpc1601.gw.1.swall,
        dofit3.fpc1601.berm.1.swall,avfit3.fpc1602.gw.1.swall,  dofit3.fpc1602.berm.1.swall,mafit3.fpc1601.berm.1.swall)

# estimates when berm closed / groundwater below median
gici3.int0.body=list(avfit3.fmc1602.gw.0,dofit3.fmc1602.berm.0, avfit3.fpc1601.gw.0, 
        dofit3.fpc1601.berm.0, avfit3.fpc1602.gw.0, dofit3.fpc1602.berm.0,mafit3.fpc1601.berm.0)

gici3.int0.head=list(avfit3.fmc1602.gw.0.head, dofit3.fmc1602.berm.0.head,avfit3.fpc1601.gw.0.head,
        dofit3.fpc1601.berm.0.head,avfit3.fpc1602.gw.0.head,  dofit3.fpc1602.berm.0.head,mafit3.fpc1601.berm.1.head)

gici3.int0.swall=list(avfit3.fmc1602.gw.0.swall, dofit3.fmc1602.berm.0.swall,avfit3.fpc1601.gw.0.swall,
        dofit3.fpc1601.berm.0.swall,avfit3.fpc1602.gw.0.swall,  dofit3.fpc1602.berm.0.swall,mafit3.fpc1601.berm.1.swall)

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
# pooled results
gici3.tab.body=data.frame(combined=unlist(lapply(gici3.body,mkrow)))
gici3.tab.head=data.frame(combined=unlist(lapply(gici3.head,mkrow)))
gici3.tab.swall=data.frame(combined=unlist(lapply(gici3.swall,mkrow)))

gici3.tab.body$combined=as.character(gici3.tab.body$combined)
gici3.tab.head$combined=as.character(gici3.tab.head$combined)
gici3.tab.swall$combined=as.character(gici3.tab.swall$combined)

# stratified results
gici3.tab.int1.body=data.frame(int1=unlist(lapply(gici3.int1.body,mkrow)))
gici3.tab.int1.head=data.frame(int1=unlist(lapply(gici3.int1.head,mkrow)))
gici3.tab.int1.swall=data.frame(int1=unlist(lapply(gici3.int1.swall,mkrow)))

gici3.tab.int0.body=data.frame(int0=unlist(lapply(gici3.int0.body,mkrow)))
gici3.tab.int0.head=data.frame(int0=unlist(lapply(gici3.int0.head,mkrow)))
gici3.tab.int0.swall=data.frame(int0=unlist(lapply(gici3.int0.swall,mkrow)))

gici3.tab.int1.body$int1=as.character(gici3.tab.int1.body$int1)
gici3.tab.int1.head$int1=as.character(gici3.tab.int1.head$int1)
gici3.tab.int1.swall$int1=as.character(gici3.tab.int1.swall$int1)

gici3.tab.int0.body$int0=as.character(gici3.tab.int0.body$int0)
gici3.tab.int0.head$int0=as.character(gici3.tab.int0.head$int0)
gici3.tab.int0.swall$int0=as.character(gici3.tab.int0.swall$int0)

# combine the pooled and high/low risk estimates into one object
gici3.tab.body=c(gici3.tab.body$combined[1:3],gici3.tab.int0.body$int0[1],gici3.tab.int1.body$int1[1],
                gici3.tab.body$combined[4],gici3.tab.int0.body$int0[2],gici3.tab.int1.body$int1[2],
                gici3.tab.body$combined[5],gici3.tab.int0.body$int0[3],gici3.tab.int1.body$int1[3],
                gici3.tab.body$combined[6],gici3.tab.int0.body$int0[4],gici3.tab.int1.body$int1[4],
                gici3.tab.body$combined[7:8],gici3.tab.int0.body$int0[5],gici3.tab.int1.body$int1[5],
                gici3.tab.body$combined[9:10],gici3.tab.int0.body$int0[6],gici3.tab.int1.body$int1[6],
                gici3.tab.body$combined[11],gici3.tab.int0.body$int0[7],gici3.tab.int1.body$int1[7])

gici3.tab.head=c(gici3.tab.head$combined[1:3],gici3.tab.int0.head$int0[1],gici3.tab.int1.head$int1[1],
                gici3.tab.head$combined[4],gici3.tab.int0.head$int0[2],gici3.tab.int1.head$int1[2],
                gici3.tab.head$combined[5],gici3.tab.int0.head$int0[3],gici3.tab.int1.head$int1[3],
                gici3.tab.head$combined[6],gici3.tab.int0.head$int0[4],gici3.tab.int1.head$int1[4],
                gici3.tab.head$combined[7:8],gici3.tab.int0.head$int0[5],gici3.tab.int1.head$int1[5],
                gici3.tab.head$combined[9:10],gici3.tab.int0.head$int0[6],gici3.tab.int1.head$int1[6],
                gici3.tab.head$combined[11],gici3.tab.int0.head$int0[7],gici3.tab.int1.head$int1[7])

gici3.tab.swall=c(gici3.tab.swall$combined[1:3],gici3.tab.int0.swall$int0[1],gici3.tab.int1.swall$int1[1],
                gici3.tab.swall$combined[4],gici3.tab.int0.swall$int0[2],gici3.tab.int1.swall$int1[2],
                gici3.tab.swall$combined[5],gici3.tab.int0.swall$int0[3],gici3.tab.int1.swall$int1[3],
                gici3.tab.swall$combined[6],gici3.tab.int0.swall$int0[4],gici3.tab.int1.swall$int1[4],
                gici3.tab.swall$combined[7:8],gici3.tab.int0.swall$int0[5],gici3.tab.int1.swall$int1[5],
                gici3.tab.swall$combined[9:10],gici3.tab.int0.swall$int0[6],gici3.tab.int1.swall$int1[6],
                gici3.tab.swall$combined[11],gici3.tab.int0.swall$int0[7],gici3.tab.int1.swall$int1[7])

gici3.n.body=c(gici3.n.body[1:3],gici3.n.int0.body[1],gici3.n.int1.body[1],
                gici3.n.body[4],gici3.n.int0.body[2],gici3.n.int1.body[2],
                gici3.n.body[5],gici3.n.int0.body[3],gici3.n.int1.body[3],
                gici3.n.body[6],gici3.n.int0.body[4],gici3.n.int1.body[4],
                gici3.n.body[7:8],gici3.n.int0.body[5],gici3.n.int1.body[5],
                gici3.n.body[9:10],gici3.n.int0.body[6],gici3.n.int1.body[6],
                gici3.n.body[11],gici3.n.int0.body[7],gici3.n.int1.body[7])

gici3.n.head=c(gici3.n.head[1:3],gici3.n.int0.head[1],gici3.n.int1.head[1],
                gici3.n.head[4],gici3.n.int0.head[2],gici3.n.int1.head[2],
                gici3.n.head[5],gici3.n.int0.head[3],gici3.n.int1.head[3],
                gici3.n.head[6],gici3.n.int0.head[4],gici3.n.int1.head[4],
                gici3.n.head[7:8],gici3.n.int0.head[5],gici3.n.int1.head[5],
                gici3.n.head[9:10],gici3.n.int0.head[6],gici3.n.int1.head[6],
                gici3.n.head[11],gici3.n.int0.head[7],gici3.n.int1.head[7])

gici3.n.swall=c(gici3.n.swall[1:3],gici3.n.int0.swall[1],gici3.n.int1.swall[1],
                gici3.n.swall[4],gici3.n.int0.swall[2],gici3.n.int1.swall[2],
                gici3.n.swall[5],gici3.n.int0.swall[3],gici3.n.int1.swall[3],
                gici3.n.swall[6],gici3.n.int0.swall[4],gici3.n.int1.swall[4],
                gici3.n.swall[7:8],gici3.n.int0.swall[5],gici3.n.int1.swall[5],
                gici3.n.swall[9:10],gici3.n.int0.swall[6],gici3.n.int1.swall[6],
                gici3.n.swall[11],gici3.n.int0.swall[7],gici3.n.int1.swall[7])

label=c("Doheny","Mission Bay",
        "Avalon","","","Doheny","","",
        "Avalon","","","Doheny","","","Fairhope","Malibu","","","Mission Bay",
        "Avalon","","","Doheny","","")

gici3.exp=data.frame(cbind(label,gici3.n.body,gici3.tab.body,
    gici3.n.head,gici3.tab.head,gici3.n.swall,gici3.tab.swall))

gici3.exp$label=as.character(gici3.exp$label)
gici3.exp$cond=c(rep("All conditions",3),"Groundwater below median flow",
      "Groundwater above median flow","All conditions","Berm closed","Berm open",
      "All conditions","Groundwater below median flow","Groundwater above median flow",
      "All conditions","Berm closed","Berm open","All conditions","All conditions",
      "Berm closed","Berm open","All conditions","All conditions","Groundwater below median flow",
      "Groundwater above median flow","All conditions","Berm closed","Berm open")

gici3.exp=gici3.exp[,c(1,8,2:7)]

colnames(gici3.exp)=c("label","conditions","n-body","body","n-head","head","n-swall","swall")

save(gici3.exp,file="~/Dropbox/Coliphage/Results/Tables/reg-exposure-3.RData")





