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
# combined n's
gici10.n.body=c(do.n10.fmc1601, mb.n10.fmc1601,
				av.n10.fmc1602, do.n10.fmc1602,
				av.n10.fpc1601, do.n10.fpc1601, fa.n10.fpc1601, ma.n10.fpc1601, mb.n10.fpc1601,
        av.n10.fpc1602, do.n10.fpc1602)

gici10.n.head=c(do.n10.fmc1601head, mb.n10.fmc1601head,
				av.n10.fmc1602head, do.n10.fmc1602head,
				av.n10.fpc1601head, do.n10.fpc1601head, fa.n10.fpc1601head, ma.n10.fpc1601head, mb.n10.fpc1601head,
        av.n10.fpc1602head, do.n10.fpc1602head)

gici10.n.swall=c(do.n10.fmc1601swall, mb.n10.fmc1601swall,
				av.n10.fmc1602swall, do.n10.fmc1602swall,
				av.n10.fpc1601swall, do.n10.fpc1601swall, fa.n10.fpc1601swall, ma.n10.fpc1601swall, mb.n10.fpc1601swall,
        av.n10.fpc1602swall, do.n10.fpc1602swall)

# n if berm closed/groundwater below median
gici10.n.int0.body=c(av.n10.fmc1602.int0, do.n10.fmc1602.int0,
        av.n10.fpc1601.int0, do.n10.fpc1601.int0, ma.n10.fpc1601.int0, 
        av.n10.fpc1602.int0, do.n10.fpc1602.int0)

gici10.n.int0.head=c(av.n10.fmc1602.int0head, do.n10.fmc1602.int0head,
        av.n10.fpc1601.int0head, do.n10.fpc1601.int0head, ma.n10.fpc1601.int0head, 
        av.n10.fpc1602.int0head, do.n10.fpc1602.int0)

gici10.n.int0.swall=c(av.n10.fmc1602.int0swall, do.n10.fmc1602.int0swall,
        av.n10.fpc1601.int0swall, do.n10.fpc1601.int0swall,ma.n10.fpc1601.int0swall,
        av.n10.fpc1602.int0swall, do.n10.fpc1602.int0)

# n if berm open/groundwater above median
gici10.n.int1.body=c(av.n10.fmc1602.int1, do.n10.fmc1602.int1,
        av.n10.fpc1601.int1, do.n10.fpc1601.int1, ma.n10.fpc1601.int1, 
        av.n10.fpc1602.int1, do.n10.fpc1602.int1)

gici10.n.int1.head=c(av.n10.fmc1602.int1head, do.n10.fmc1602.int1head,
        av.n10.fpc1601.int1head, do.n10.fpc1601.int1head, ma.n10.fpc1601.int1head, 
        av.n10.fpc1602.int1head, do.n10.fpc1602.int1)

gici10.n.int1.swall=c(av.n10.fmc1602.int1swall, do.n10.fmc1602.int1swall,
        av.n10.fpc1601.int1swall, do.n10.fpc1601.int1swall, ma.n10.fpc1601.int1swall,
        av.n10.fpc1602.int1swall, do.n10.fpc1602.int1)

# estimates pooled across conditions  --------------------------------------------
gici10.body=list(dofit10.fmc1601, mbfit10.fmc1601,
				avfit10.fmc1602, dofit10.fmc1602,
				avfit10.fpc1601, dofit10.fpc1601, fafit10.fpc1601, mafit10.fpc1601, mbfit10.fpc1601,
                avfit10.fpc1602, dofit10.fpc1602)

gici10.head=list(dofit10.fmc1601head, mbfit10.fmc1601head,
				avfit10.fmc1602head, dofit10.fmc1602head,
				avfit10.fpc1601head, dofit10.fpc1601head, fafit10.fpc1601head, mafit10.fpc1601head, mbfit10.fpc1601head,
                avfit10.fpc1602head, dofit10.fpc1602head)

gici10.swall=list(dofit10.fmc1601swall, mbfit10.fmc1601swall,
				avfit10.fmc1602swall, dofit10.fmc1602swall,
				avfit10.fpc1601swall, dofit10.fpc1601swall, fafit10.fpc1601swall, mafit10.fpc1601swall, mbfit10.fpc1601swall,
                avfit10.fpc1602swall, dofit10.fpc1602swall)

# estimates when berm open/ groundwater above median
gici10.int1.body=list(avfit10.fmc1602.gw.1, dofit10.fmc1602.berm.1,avfit10.fpc1601.gw.1,
        dofit10.fpc1601.berm.1,avfit10.fpc1602.gw.1,  dofit10.fpc1602.berm.1,mafit10.fpc1601.berm.1)

gici10.int1.head=list(avfit10.fmc1602.gw.1.head, dofit10.fmc1602.berm.1.head,avfit10.fpc1601.gw.1.head,
        dofit10.fpc1601.berm.1.head,avfit10.fpc1602.gw.1.head,  dofit10.fpc1602.berm.1.head,mafit10.fpc1601.berm.1.head)

gici10.int1.swall=list(avfit10.fmc1602.gw.1.swall, dofit10.fmc1602.berm.1.swall,avfit10.fpc1601.gw.1.swall,
        dofit10.fpc1601.berm.1.swall,avfit10.fpc1602.gw.1.swall,  dofit10.fpc1602.berm.1.swall,mafit10.fpc1601.berm.1.swall)

# estimates when berm closed / groundwater below median
gici10.int0.body=list(avfit10.fmc1602.gw.0,dofit10.fmc1602.berm.0, avfit10.fpc1601.gw.0, 
        dofit10.fpc1601.berm.0, avfit10.fpc1602.gw.0, dofit10.fpc1602.berm.0,mafit10.fpc1601.berm.0)

gici10.int0.head=list(avfit10.fmc1602.gw.0.head, dofit10.fmc1602.berm.0.head,avfit10.fpc1601.gw.0.head,
        dofit10.fpc1601.berm.0.head,avfit10.fpc1602.gw.0.head,  dofit10.fpc1602.berm.0.head,mafit10.fpc1601.berm.1.head)

gici10.int0.swall=list(avfit10.fmc1602.gw.0.swall, dofit10.fmc1602.berm.0.swall,avfit10.fpc1601.gw.0.swall,
        dofit10.fpc1601.berm.0.swall,avfit10.fpc1602.gw.0.swall,  dofit10.fpc1602.berm.0.swall,mafit10.fpc1601.berm.1.swall)

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
# results stratified by beach
# ------------------------------------------------
# combined results
gici10.tab.body=data.frame(combined=unlist(lapply(gici10.body,mkrow)))
gici10.tab.head=data.frame(combined=unlist(lapply(gici10.head,mkrow)))
gici10.tab.swall=data.frame(combined=unlist(lapply(gici10.swall,mkrow)))

# combined results
gici10.tab.body=data.frame(combined=unlist(lapply(gici10.body,mkrow)))
gici10.tab.head=data.frame(combined=unlist(lapply(gici10.head,mkrow)))
gici10.tab.swall=data.frame(combined=unlist(lapply(gici10.swall,mkrow)))

gici10.tab.body$combined=as.character(gici10.tab.body$combined)
gici10.tab.head$combined=as.character(gici10.tab.head$combined)
gici10.tab.swall$combined=as.character(gici10.tab.swall$combined)

# stratified results
gici10.tab.int1.body=data.frame(int1=unlist(lapply(gici10.int1.body,mkrow)))
gici10.tab.int1.head=data.frame(int1=unlist(lapply(gici10.int1.head,mkrow)))
gici10.tab.int1.swall=data.frame(int1=unlist(lapply(gici10.int1.swall,mkrow)))

gici10.tab.int0.body=data.frame(int0=unlist(lapply(gici10.int0.body,mkrow)))
gici10.tab.int0.head=data.frame(int0=unlist(lapply(gici10.int0.head,mkrow)))
gici10.tab.int0.swall=data.frame(int0=unlist(lapply(gici10.int0.swall,mkrow)))

gici10.tab.int1.body$int1=as.character(gici10.tab.int1.body$int1)
gici10.tab.int1.head$int1=as.character(gici10.tab.int1.head$int1)
gici10.tab.int1.swall$int1=as.character(gici10.tab.int1.swall$int1)

gici10.tab.int0.body$int0=as.character(gici10.tab.int0.body$int0)
gici10.tab.int0.head$int0=as.character(gici10.tab.int0.head$int0)
gici10.tab.int0.swall$int0=as.character(gici10.tab.int0.swall$int0)

# combine the pooled and high/low risk estimates into one object
gici10.tab.body=c(gici10.tab.body$combined[1:3],gici10.tab.int0.body$int0[1],gici10.tab.int1.body$int1[1],
                gici10.tab.body$combined[4],gici10.tab.int0.body$int0[2],gici10.tab.int1.body$int1[2],
                gici10.tab.body$combined[5],gici10.tab.int0.body$int0[3],gici10.tab.int1.body$int1[3],
                gici10.tab.body$combined[6],gici10.tab.int0.body$int0[4],gici10.tab.int1.body$int1[4],
                gici10.tab.body$combined[7:8],gici10.tab.int0.body$int0[5],gici10.tab.int1.body$int1[5],
                gici10.tab.body$combined[9:10],gici10.tab.int0.body$int0[6],gici10.tab.int1.body$int1[6],
                gici10.tab.body$combined[11],gici10.tab.int0.body$int0[7],gici10.tab.int1.body$int1[7])

gici10.tab.head=c(gici10.tab.head$combined[1:3],gici10.tab.int0.head$int0[1],gici10.tab.int1.head$int1[1],
                gici10.tab.head$combined[4],gici10.tab.int0.head$int0[2],gici10.tab.int1.head$int1[2],
                gici10.tab.head$combined[5],gici10.tab.int0.head$int0[3],gici10.tab.int1.head$int1[3],
                gici10.tab.head$combined[6],gici10.tab.int0.head$int0[4],gici10.tab.int1.head$int1[4],
                gici10.tab.head$combined[7:8],gici10.tab.int0.head$int0[5],gici10.tab.int1.head$int1[5],
                gici10.tab.head$combined[9:10],gici10.tab.int0.head$int0[6],gici10.tab.int1.head$int1[6],
                gici10.tab.head$combined[11],gici10.tab.int0.head$int0[7],gici10.tab.int1.head$int1[7])

gici10.tab.swall=c(gici10.tab.swall$combined[1:3],gici10.tab.int0.swall$int0[1],gici10.tab.int1.swall$int1[1],
                gici10.tab.swall$combined[4],gici10.tab.int0.swall$int0[2],gici10.tab.int1.swall$int1[2],
                gici10.tab.swall$combined[5],gici10.tab.int0.swall$int0[3],gici10.tab.int1.swall$int1[3],
                gici10.tab.swall$combined[6],gici10.tab.int0.swall$int0[4],gici10.tab.int1.swall$int1[4],
                gici10.tab.swall$combined[7:8],gici10.tab.int0.swall$int0[5],gici10.tab.int1.swall$int1[5],
                gici10.tab.swall$combined[9:10],gici10.tab.int0.swall$int0[6],gici10.tab.int1.swall$int1[6],
                gici10.tab.swall$combined[11],gici10.tab.int0.swall$int0[7],gici10.tab.int1.swall$int1[7])

gici10.n.body=c(gici10.n.body[1:3],gici10.n.int0.body[1],gici10.n.int1.body[1],
                gici10.n.body[4],gici10.n.int0.body[2],gici10.n.int1.body[2],
                gici10.n.body[5],gici10.n.int0.body[3],gici10.n.int1.body[3],
                gici10.n.body[6],gici10.n.int0.body[4],gici10.n.int1.body[4],
                gici10.n.body[7:8],gici10.n.int0.body[5],gici10.n.int1.body[5],
                gici10.n.body[9:10],gici10.n.int0.body[6],gici10.n.int1.body[6],
                gici10.n.body[11],gici10.n.int0.body[7],gici10.n.int1.body[7])

gici10.n.head=c(gici10.n.head[1:3],gici10.n.int0.head[1],gici10.n.int1.head[1],
                gici10.n.head[4],gici10.n.int0.head[2],gici10.n.int1.head[2],
                gici10.n.head[5],gici10.n.int0.head[3],gici10.n.int1.head[3],
                gici10.n.head[6],gici10.n.int0.head[4],gici10.n.int1.head[4],
                gici10.n.head[7:8],gici10.n.int0.head[5],gici10.n.int1.head[5],
                gici10.n.head[9:10],gici10.n.int0.head[6],gici10.n.int1.head[6],
                gici10.n.head[11],gici10.n.int0.head[7],gici10.n.int1.head[7])

gici10.n.swall=c(gici10.n.swall[1:3],gici10.n.int0.swall[1],gici10.n.int1.swall[1],
                gici10.n.swall[4],gici10.n.int0.swall[2],gici10.n.int1.swall[2],
                gici10.n.swall[5],gici10.n.int0.swall[3],gici10.n.int1.swall[3],
                gici10.n.swall[6],gici10.n.int0.swall[4],gici10.n.int1.swall[4],
                gici10.n.swall[7:8],gici10.n.int0.swall[5],gici10.n.int1.swall[5],
                gici10.n.swall[9:10],gici10.n.int0.swall[6],gici10.n.int1.swall[6],
                gici10.n.swall[11],gici10.n.int0.swall[7],gici10.n.int1.swall[7])

label=c("Doheny","Mission Bay",
        "Avalon","","","Doheny","","",
        "Avalon","","","Doheny","","","Fairhope","Malibu","","","Mission Bay",
        "Avalon","","","Doheny","","")

gici10.exp=data.frame(cbind(label,gici10.n.body,gici10.tab.body,
    gici10.n.head,gici10.tab.head,gici10.n.swall,gici10.tab.swall))

gici10.exp$label=as.character(gici10.exp$label)
gici10.exp$cond=c(rep("All conditions",3),"Groundwater below median flow",
      "Groundwater above median flow","All conditions","Berm closed","Berm open",
      "All conditions","Groundwater below median flow","Groundwater above median flow",
      "All conditions","Berm closed","Berm open","All conditions","All conditions",
      "Berm closed","Berm open","All conditions","All conditions","Groundwater below median flow",
      "Groundwater above median flow","All conditions","Berm closed","Berm open")

gici10.exp=gici10.exp[,c(1,8,2:7)]

colnames(gici10.exp)=c("label","conditions","n-body","body","n-head","head","n-swall","swall")

save(gici10.exp,file="~/Dropbox/Coliphage/Results/Tables/reg-exposure-10.RData")





