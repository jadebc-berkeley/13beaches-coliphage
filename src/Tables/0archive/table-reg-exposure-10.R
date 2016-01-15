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

# combine output
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

gici10.tab.body=data.frame(combined=unlist(lapply(gici10.body,mkrow)))
gici10.tab.head=data.frame(combined=unlist(lapply(gici10.head,mkrow)))
gici10.tab.swall=data.frame(combined=unlist(lapply(gici10.swall,mkrow)))

label=c("Doheny","Mission Bay",
		"Avalon","Doheny",
		"Avalon","Doheny","Fairhope","Malibu","Mission Bay",
		"Avalon","Doheny")

gici10.exp=cbind(label,gici10.n.body,gici10.tab.body,
                gici10.n.head,gici10.tab.head,
                gici10.n.swall,gici10.tab.swall)
gici10.exp$label=as.character(gici10.exp$label)
colnames(gici10.exp)=c("label","n-body","body","n-head","head","n-swall","swall")

save(gici10.exp,file="~/Dropbox/Coliphage/Results/Tables/reg-exposure-10.RData")





