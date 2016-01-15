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
load("results/rawoutput/regress-3day-head.Rdata")
load("results/rawoutput/regress-3day-swall.Rdata")

# combine output
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

gici3.tab.body=data.frame(combined=unlist(lapply(gici3.body,mkrow)))
gici3.tab.head=data.frame(combined=unlist(lapply(gici3.head,mkrow)))
gici3.tab.swall=data.frame(combined=unlist(lapply(gici3.swall,mkrow)))

label=c("Doheny","Mission Bay",
		"Avalon","Doheny",
		"Avalon","Doheny","Fairhope","Malibu","Mission Bay",
		"Avalon","Doheny")

gici3.exp=cbind(label,gici3.n.body,gici3.tab.body,
                gici3.n.head,gici3.tab.head,
                gici3.n.swall,gici3.tab.swall)
gici3.exp$label=as.character(gici3.exp$label)
#gici3.exp.print=rbind(c("EPA 1601 F- Coliphage",rep(NA,3)),gici3.exp[c(1:2),],
#                      c("EPA 1602 F- Coliphage",rep(NA,3)),gici3.exp[c(3:4),],
#                      c("EPA 1601 F+ Coliphage",rep(NA,3)),gici3.exp[c(5:9),],
#                      c("EPA 1602 F+ Coliphage",rep(NA,3)),gici3.exp[c(10:11),])
colnames(gici3.exp)=c("label","n-body","body","n-head","head","n-swall","swall")

save(gici3.exp,file="~/Dropbox/Coliphage/Results/Tables/reg-exposure-3.RData")





