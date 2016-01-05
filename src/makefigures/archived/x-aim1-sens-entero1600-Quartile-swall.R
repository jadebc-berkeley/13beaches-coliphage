# --------------------------------------
# x-aim1-sens-entero1600-quartiles-swall.R
# ben arnold (benarnold@berkeley.edu)
#
# description:
# plot CIRs associated with Enterococcus
# EPA 1600 quartiles
# using summary regression output
#
# Sensitivity analysis looking at only
# head immersion swimmers
#
#
# version 2 (24 mar 2015)
# re-designed the main CIR figure to 
# report incidence in the figure with
# adjusted CIRs in the portion
#
# version 1 (20 feb 2015)
#
# --------------------------------------

# --------------------------------------
# input files:
#	13beaches-wq.csv
#	aim1-entero1600-Quartile-regs-swall.Rdata
#
# output files:
#	aim1-sens-entero1600-Quartile-CI-swall.pdf
# --------------------------------------

# --------------------------------------
# preamble
# --------------------------------------

rm(list=ls())
library(RColorBrewer)


# --------------------------------------
# Load the water quality dataset to
# get the mid-points of the Entero
# concentrations by each quartile
# --------------------------------------
wq <- read.csv("~/dropbox/13beaches/data/final/13beaches-wq.csv")

minQs <- tapply(wq$avgdyentero1600,wq$qavgdyentero1600,function(x) min(x,na.rm=T))
maxQs <- tapply(wq$avgdyentero1600,wq$qavgdyentero1600,function(x) max(x,na.rm=T))
labQs <- paste(sprintf("%1.0f",10^minQs)," to ",sprintf("%1.0f",10^maxQs),sep="")
labQs <- paste("Q",1:4,"\n(",labQs,")",sep="")

rngQs <- paste(sprintf("%1.0f",round(10^minQs)),"-",sprintf("%1.0f",floor(10^maxQs)),sep="")
rngQs[4] <-paste(">",sprintf("%1.0f",floor(10^maxQs))[3],sep="")

midQs <- minQs + (maxQs-minQs)/2


# --------------------------------------
# Plot Enterococcus 1600 Quartile Results
# --------------------------------------

load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-entero1600-Quartile-regs-swall.Rdata")


pdf("~/dropbox/13beaches/aim1-results/figs/aim1-sens-entero1600-Quartile-CI-swall.pdf",width=16,height=5)
op <- par(mar=c(4,7,8,0)+0.1,xpd=TRUE)
cols <- brewer.pal(9,"YlGnBu")[8:5]
ytics <- seq(0,150,by=30)
# set up an empty plot
MidPts <- barplot(1:4,names.arg=NA,border=NA,col=NA,
	ylim=range(ytics),ylab="",yaxt="n",
	las=1,bty="n"
	)
	segments(x0=0,x1=max(MidPts+0.5),y0=ytics,lty=2,lwd=1,col="gray80")
	segments(x0=mean(MidPts[3:4]),y0=min(ytics)-8,y1=max(ytics)+23,lwd=2,col="gray80")
	axis(2,at=ytics,las=1)
	mtext("Diarrhea\nIncidence\nper 1000",side=2,line=3,las=1)
	
	# calculate X coordinates relative to the mid points for each group
	xspan <- 0.37
	xplus <- c(-xspan, -xspan/3, xspan/3, xspan)  # evenly distribute 4 datapoints around each midpoint
	
	
	x0to4   <- xplus+MidPts[1]
	x5to10  <- xplus+MidPts[2]
	x11plus <- xplus+MidPts[3]
	xall    <- xplus+MidPts[4]
	allxs   <- c(x0to4,x5to10,x11plus,xall)  # for table and quartile labels in header/footer
	labx <- MidPts[1]-xspan*1.5  # for left-hand labels in the header/footer

	# plot age 0 to 4 estimates
	segments(x0=x0to4,y0=ci.0to4[,"CIlb"]*1000,y1=ci.0to4[,"CIub"]*1000,lwd=2,col=cols)
	points(x0to4,ci.0to4[,"CI"]*1000,pch=16,cex=1.75,lwd=2,col=cols)
	text(x=x0to4,y=ci.0to4[,"CI"]*1000,sprintf("%1.0f",ci.0to4[,"CI"]*1000),pos=4,cex=0.7,col=cols,font=2)
	
	# plot age 5 to 10 estimates
	segments(x0=x5to10,y0=ci.5to10[,"CIlb"]*1000,y1=ci.5to10[,"CIub"]*1000,lwd=2,col=cols)
	points(x5to10,ci.5to10[,"CI"]*1000,pch=16,cex=1.75,lwd=2,col=cols)
	text(x=x5to10,y=ci.5to10[,"CI"]*1000,sprintf("%1.0f",ci.5to10[,"CI"]*1000),pos=4,cex=0.7,col=cols,font=2)
	
	# plot age > 10 estimates
	segments(x0=x11plus,y0=ci.11plus[,"CIlb"]*1000,y1=ci.11plus[,"CIub"]*1000,lwd=2,col=cols)
	points(x11plus,ci.11plus[,"CI"]*1000,pch=16,cex=1.75,lwd=2,col=cols)
	text(x=x11plus,y=ci.11plus[,"CI"]*1000,sprintf("%1.0f",ci.11plus[,"CI"]*1000),pos=4,cex=0.7,col=cols,font=2)
	
	# plot all age estimates
	segments(x0=xall,y0=ci.all[,"CIlb"]*1000,y1=ci.all[,"CIub"]*1000,lwd=2,col=cols)
	points(xall,ci.all[,"CI"]*1000,pch=16,cex=1.75,lwd=2,col=cols)
	text(x=xall,y=ci.all[,"CI"]*1000,sprintf("%1.0f",ci.all[,"CI"]*1000),pos=4,cex=0.7,col=cols,font=2)
	
	
	# print header labels
	mtext(c("Ages\n0 to 4","Ages\n5 to 10","Ages\n>10","All\nAges"),at=MidPts,side=3,line=5.5  )
	
	mtext(expression(paste(italic("Enterococcus")," 1600 Quartile")),side=3,line=4,at=labx,adj=1,col=cols[3],cex=0.8)
	mtext(c("Q1","Q2","Q3","Q4"),side=3,line=4,at=allxs,col=cols,cex=0.9,font=2)
	
	mtext("Range (CFU/100ml)",side=3,line=3,at=labx,adj=1,col=cols[3],cex=0.8)
	mtext(rngQs,side=3,line=3,at=allxs,col=cols,cex=0.8)
	
	# Print adjusted CIRs and 95% CIs (formatted)
	cirform <- function(cirs) {
		paste(sprintf("%1.2f",cirs),sep="")
	}
	circiform <- function(circi) {
		apply(circi,1,function(x) paste("(",sprintf("%1.2f",x["CIRlb"]),", ",sprintf("%1.2f",x["CIRub"]),")",sep="") )
	}
	mtext("Adjusted CIR",side=3,line=2,at=labx,adj=1,cex=0.8,col="gray30")
	mtext(c("ref",cirform(cir.age0to4[,1]),"ref",cirform(cir.age5to10[,1]),"ref",cirform(cir.age11plus[,1]),"ref",cirform(cir.all[,1])),side=3,line=2,at=allxs,cex=0.75)
	
	mtext("(95% CI)",side=3,line=1,at=labx,adj=1,cex=0.8,col="gray30")
	mtext(c("",circiform(cir.age0to4[,2:3]),"",circiform(cir.age5to10[,2:3]),"",circiform(cir.age11plus[,2:3]),"",circiform(cir.all[,2:3])),side=3,line=1,at=allxs,cex=0.7)
	

	# print footer labels
	mtext(c("Q1","Q2","Q3","Q4"),side=1,line=0.5,at=allxs,col=cols,cex=0.9,font=2)
	
	# print table with Ns
	mtext("Incident Diarrhea Cases",side=1,line=2,at=labx,adj=1,cex=0.8,col="gray30")
	ns <- c(N.age0to4[,1],N.age5to10[,1],N.age11plus[,1],N.all[,1])
	mtext(  format(ns,big.mark=","),side=1,line=2,at=allxs+0.03,adj=1,cex=0.75    )
	
	mtext("Population At Risk",side=1,line=3,at=labx,adj=1,cex=0.8,col="gray30")
	Ns <- c(N.age0to4[,2],N.age5to10[,2],N.age11plus[,2],N.all[,2])
	mtext(  format(Ns,big.mark=","),side=1,line=3,at=allxs+0.03,adj=1,cex=0.75    )
	
par(op)
dev.off()





