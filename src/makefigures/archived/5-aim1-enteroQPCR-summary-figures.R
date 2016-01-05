# --------------------------------------
# 5-aim1-enteroQPCR-summary-figures.R
# ben arnold (benarnold@berkeley.edu)
#
# description:
# plot CIRs associated with Enterococcus
# EPA 1600 output
# using summary regression output
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
#	aim1-enteroQPCR-Quartile-regs-body.Rdata
#	aim1-enteroQPCR-470cce-regs-body.Rdata
#
# output files:
#	aim1-enteroQPCR-Quartile-CI-byage.pdf
#	aim1-enteroQPCR-470cce-CI-bypointsource.pdf
#	aim1-enteroQPCR-470cce-CI-byage.pdf
# # # # # # #	aim1-enteroQPCR-470cce-CI-byage-pointsource.pdf
# # # # # # #	aim1-enteroQPCR-470cce-CI-byage-nonpointsource.pdf
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

minQs <- tapply(wq$avgdyenteropcr,wq$qavgdyenteropcr,function(x) min(x,na.rm=T))
maxQs <- tapply(wq$avgdyenteropcr,wq$qavgdyenteropcr,function(x) max(x,na.rm=T))
labQs <- paste(sprintf("%1.0f",10^minQs)," to ",sprintf("%1.0f",10^maxQs),sep="")
labQs <- paste("Q",1:4,"\n(",labQs,")",sep="")

rngQs <- paste(sprintf("%1.0f",round(10^minQs)),"-",sprintf("%1.0f",floor(10^maxQs)),sep="")
rngQs[4] <-paste(">",sprintf("%1.0f",floor(10^maxQs))[3],sep="")

midQs <- minQs + (maxQs-minQs)/2


# --------------------------------------
# Plot Enterococcus 1600 Quartile Results
# --------------------------------------

load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-enteroQPCR-Quartile-regs-body.Rdata")


pdf("~/dropbox/13beaches/aim1-results/figs/aim1-enteroQPCR-Quartile-CI-byage.pdf",width=16,height=5)
op <- par(mar=c(4,7,8,0)+0.1,xpd=TRUE)
cols <- brewer.pal(9,"YlOrBr")[9:6]
ytics <- seq(0,100,by=20)
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
	
	mtext(expression(paste(italic("Enterococcus")," qPCR Quartile")),side=3,line=4,at=labx,adj=1,col=cols[3],cex=0.8)
	mtext(c("Q1","Q2","Q3","Q4"),side=3,line=4,at=allxs,col=cols,cex=0.9,font=2)
	
	mtext("Range (CCE/100ml)",side=3,line=3,at=labx,adj=1,col=cols[3],cex=0.8)
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

# --------------------------------------
# Plot Enterococcus qPCR <=470 / >470 Results
# All ages, stratified by beach type
# --------------------------------------


load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-enteroQPCR-470cce-regs-body.Rdata")


pdf("~/dropbox/13beaches/aim1-results/figs/aim1-enteroQPCR-470cce-CI-bypointsource.pdf",width=11,height=5)
op <- par(mar=c(5,9,8,0)+0.1,xpd=TRUE)
cols <- brewer.pal(9,"YlOrBr")[c(8,6)]
ytics <- seq(0,100,by=20)
# set up an empty plot
MidPts <- barplot(1:3,names.arg=NA,border=NA,col=NA,
	ylim=range(ytics),ylab="",yaxt="n",
	las=1,bty="n"
)
segments(x0=0.1,x1=max(MidPts+0.5),y0=ytics,lty=2,lwd=1,col="gray80")
segments(x0=mean(MidPts[2:3]),y0=min(ytics)-8,y1=max(ytics)+18,lwd=2,col="gray80")
axis(2,at=ytics,las=1)
mtext("Diarrhea\nIncidence\nper 1000",side=2,line=3,las=1)

# calculate X coordinates relative to the mid points for each group
xspan <- 0.2
xplus <- c(-xspan, xspan)  # evenly distribute 2 datapoints around each midpoint

xps   <- xplus+MidPts[1]
xnps  <- xplus+MidPts[2]
xall    <- xplus+MidPts[3]
allxs   <- c(xps,xnps,xall)  # for table and category labels in header/footer
labx <- MidPts[1]-xspan*3  # for left-hand labels in the header/footer

# plot point source estimates
segments(x0=xps,y0=ci.ps[,"CIlb"]*1000,y1=ci.ps[,"CIub"]*1000,lwd=2,col=cols)
points(xps,ci.ps[,"CI"]*1000,pch=16,cex=1.5,lwd=2,col=cols)
text(x=xps,y=ci.ps[,"CI"]*1000,sprintf("%1.0f",ci.ps[,"CI"]*1000),pos=4,cex=0.7,col=cols,font=2)

# plot non-point source estimates
segments(x0=xnps,y0=ci.nps[,"CIlb"]*1000,y1=ci.nps[,"CIub"]*1000,lwd=2,col=cols)
points(xnps,ci.nps[,"CI"]*1000,pch=16,cex=1.5,lwd=2,col=cols)
text(x=xnps,y=ci.nps[,"CI"]*1000,sprintf("%1.0f",ci.nps[,"CI"]*1000),pos=4,cex=0.7,col=cols,font=2)

# plot all conditions estimates
segments(x0=xall,y0=ci.all[,"CIlb"]*1000,y1=ci.all[,"CIub"]*1000,lwd=2,col=cols)
points(xall,ci.all[,"CI"]*1000,pch=16,cex=1.5,lwd=2,col=cols)
text(x=xall,y=ci.all[,"CI"]*1000,sprintf("%1.0f",ci.all[,"CI"]*1000),pos=4,cex=0.7,col=cols,font=2)

# print header labels
mtext(c("Identified\nPoint Source","No\nPoint Source","All\nConditions"),at=MidPts,side=3,line=5.5  )

mtext(expression(paste(italic("Enterococcus")," qPCR")),side=3,line=4,at=labx,adj=1,col=cols[2],cex=0.8)
mtext("Category",side=3,line=3,at=labx,adj=1,col=cols[2],cex=0.8)
mtext(expression(""<= 470,""> 470),side=3,line=4,at=allxs,col=cols,cex=1,font=2)
mtext("CCE/100ml",side=3,line=3,at=allxs,col=cols,cex=0.75)


# Print adjusted CIRs and 95% CIs (formatted)
cirform <- function(cirs) {
	paste(sprintf("%1.2f",cirs),sep="")
}
circiform <- function(circi) {
	paste("(",sprintf("%1.2f",circi["CIRlb"]),", ",sprintf("%1.2f",circi["CIRub"]),")",sep="")
}
mtext("Adjusted CIR",side=3,line=2,at=labx,adj=1,cex=0.8,col="gray30")
mtext(c("ref",cirform(cir.all[2,1]),"ref",cirform(cir.all[1,1]),"ref",cirform(cir.all[3,1])),side=3,line=2,at=allxs,cex=0.75)

mtext("(95% CI)",side=3,line=1,at=labx,adj=1,cex=0.8,col="gray30")
mtext(c("",circiform(cir.all[2,2:3]),"",circiform(cir.all[1,2:3]),"",circiform(cir.all[3,2:3])),side=3,line=1,at=allxs,cex=0.7)


# print footer labels
mtext(expression(""<= 470,""> 470),side=1,line=0.25,at=allxs,col=cols,cex=1,font=2)
mtext("CCE/100ml",side=1,line=1.25,at=allxs,col=cols,cex=0.75)

# print table with Ns
mtext("Incident Diarrhea Cases",side=1,line=2.25,at=labx,adj=1,cex=0.8,col="gray30")
ns <- c(N.ps[,1],N.nps[,1],N.all[,1])
mtext(  format(ns,big.mark=","),side=1,line=2.25,at=allxs+0.06,adj=1,cex=0.75    )

mtext("Population At Risk",side=1,line=3.25,at=labx,adj=1,cex=0.8,col="gray30")
Ns <- c(N.ps[,2],N.nps[,2],N.all[,2])
mtext(  format(Ns,big.mark=","),side=1,line=3.25,at=allxs+0.06,adj=1,cex=0.75    )

par(op)
dev.off()



# --------------------------------------
# Plot Enterococcus qPCR <=470 / >470 Results
# Age-stratified:
# Overall
# Non-Point Source
# Point Source
# --------------------------------------


### Basic Figure Shell for the Entero qPCR 470 CCE plots
### called repeatedly for different stratifications (point source, non-point source, overall)
### the function is pretty idiosyncratic -- just trying to avoid repeating all 50+ lines of code three times

CIplot470cce <- function(pdata,cols) {
	
	# pdata is a list of 12 data plotting objects in this order:
		# ci.all    : Matrix of cumulative incidence estimates w/ 95% CIs -- All ages
		# ci.0to4   : Matrix of cumulative incidence estimates w/ 95% CIs -- Age 0 to 4
		# ci.5to10  : Matrix of cumulative incidence estimates w/ 95% CIs -- Age 5 to 10
		# ci.11plus : Matrix of cumulative incidence estimates w/ 95% CIs -- Age 11 plus
	ci.all    <- pdata[[1]]
	ci.0to4   <- pdata[[2]]
	ci.5to10  <- pdata[[3]]
	ci.11plus <- pdata[[4]]
		# cir.all       : Vector of CIR estimate w/ 95% CIs -- All ages
		# cir.age0to4   : Vector of CIR estimate w/ 95% CIs -- Age 0 to 4
		# cir.age5to10  : Vector of CIR estimate w/ 95% CIs -- Age 5 to 10
		# cir.age11plus : Vector of CIR estimate w/ 95% CIs -- Age 11 plus
	cir.all       <- pdata[[5]]
	cir.age0to4   <- pdata[[6]]
	cir.age5to10  <- pdata[[7]]
	cir.age11plus <- pdata[[8]]
		# N.all       : Matrix of N incident cases and N exposed -- All ages
		# N.age0to4   : Matrix of N incident cases and N exposed -- Age 0 to 4
		# N.age5to10  : Matrix of N incident cases and N exposed -- Age 5 to 10
		# N.age11plus : Matrix of N incident cases and N exposed -- Age 11 plus
	N.all       <- pdata[[9]]
	N.age0to4   <- pdata[[10]]
	N.age5to10  <- pdata[[11]]
	N.age11plus <- pdata[[12]]
	
	op <- par(mar=c(5,8,8,0)+0.1,xpd=TRUE)
	ytics <- seq(0,100,by=20)
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
	xspan <- 0.2
	xplus <- c(-xspan, xspan)  # evenly distribute 2 datapoints around each midpoint
	
	x0to4   <- xplus+MidPts[1]
	x5to10  <- xplus+MidPts[2]
	x11plus <- xplus+MidPts[3]
	xall    <- xplus+MidPts[4]
	allxs   <- c(x0to4,x5to10,x11plus,xall)  # for table and category labels in header/footer
	labx <- MidPts[1]-xspan*3  # for left-hand labels in the header/footer

	
	# plot age 0 to 4 estimates
	segments(x0=x0to4,y0=ci.0to4[,"CIlb"]*1000,y1=ci.0to4[,"CIub"]*1000,lwd=2,col=cols)
	points(x0to4,ci.0to4[,"CI"]*1000,pch=16,cex=1.5,lwd=2,col=cols)
	text(x=x0to4,y=ci.0to4[,"CI"]*1000,sprintf("%1.0f",ci.0to4[,"CI"]*1000),pos=4,cex=0.7,col=cols,font=2)
	
	# plot age 5 to 10 estimates
	segments(x0=x5to10,y0=ci.5to10[,"CIlb"]*1000,y1=ci.5to10[,"CIub"]*1000,lwd=2,col=cols)
	points(x5to10,ci.5to10[,"CI"]*1000,pch=16,cex=1.5,lwd=2,col=cols)
	text(x=x5to10,y=ci.5to10[,"CI"]*1000,sprintf("%1.0f",ci.5to10[,"CI"]*1000),pos=4,cex=0.7,col=cols,font=2)
	
	# plot age > 10 estimates
	segments(x0=x11plus,y0=ci.11plus[,"CIlb"]*1000,y1=ci.11plus[,"CIub"]*1000,lwd=2,col=cols)
	points(x11plus,ci.11plus[,"CI"]*1000,pch=16,cex=1.5,lwd=2,col=cols)
	text(x=x11plus,y=ci.11plus[,"CI"]*1000,sprintf("%1.0f",ci.11plus[,"CI"]*1000),pos=4,cex=0.7,col=cols,font=2)
	
	# plot all age estimates
	segments(x0=xall,y0=ci.all[,"CIlb"]*1000,y1=ci.all[,"CIub"]*1000,lwd=2,col=cols)
	points(xall,ci.all[,"CI"]*1000,pch=16,cex=1.5,lwd=2,col=cols)
	text(x=xall,y=ci.all[,"CI"]*1000,sprintf("%1.0f",ci.all[,"CI"]*1000),pos=4,cex=0.7,col=cols,font=2)
	
	
	# print header labels
	mtext(c("Ages\n0 to 4","Ages\n5 to 10","Ages\n>10","All\nAges"),at=MidPts,side=3,line=5.5  )
	
	mtext(expression(paste(italic("Enterococcus")," qPCR")),side=3,line=4,at=labx,adj=1,col=cols[2],cex=0.8)
	mtext("Category",side=3,line=3,at=labx,adj=1,col=cols[2],cex=0.8)
	mtext(expression(""<= 470,""> 470),side=3,line=4,at=allxs,col=cols,cex=1,font=2)
	mtext("CCE/100ml",side=3,line=3,at=allxs,col=cols,cex=0.75)
	
	
	# Print adjusted CIRs and 95% CIs (formatted)
	cirform <- function(cirs) {
		paste(sprintf("%1.2f",cirs),sep="")
	}
	circiform <- function(circi) {
		paste("(",sprintf("%1.2f",circi["CIRlb"]),", ",sprintf("%1.2f",circi["CIRub"]),")",sep="")
	}
	mtext("Adjusted CIR",side=3,line=2,at=labx,adj=1,cex=0.8,col="gray30")
	mtext(c("ref",cirform(cir.age0to4[1]),"ref",cirform(cir.age5to10[1]),"ref",cirform(cir.age11plus[1]),"ref",cirform(cir.all[1])),side=3,line=2,at=allxs,cex=0.75)
	
	mtext("(95% CI)",side=3,line=1,at=labx,adj=1,cex=0.8,col="gray30")
	mtext(c("",circiform(cir.age0to4[2:3]),"",circiform(cir.age5to10[2:3]),"",circiform(cir.age11plus[2:3]),"",circiform(cir.all[2:3])),side=3,line=1,at=allxs,cex=0.7)
	
	# print footer labels
	mtext(expression(""<= 470,""> 470),side=1,line=0.25,at=allxs,col=cols,cex=1,font=2)
	mtext("CCE/100ml",side=1,line=1.25,at=allxs,col=cols,cex=0.75)
	
	# print table with Ns
	mtext("Incident Diarrhea Cases",side=1,line=2.25,at=labx,adj=1,cex=0.8,col="gray30")
	ns <- c(N.age0to4[,1],N.age5to10[,1],N.age11plus[,1],N.all[,1])
	mtext(  format(ns,big.mark=","),side=1,line=2.25,at=allxs+0.06,adj=1,cex=0.75    )
	
	mtext("Population At Risk",side=1,line=3.25,at=labx,adj=1,cex=0.8,col="gray30")
	Ns <- c(N.age0to4[,2],N.age5to10[,2],N.age11plus[,2],N.all[,2])
	mtext(  format(Ns,big.mark=","),side=1,line=3.25,at=allxs+0.06,adj=1,cex=0.75    )
	
	par(op)
	
}  # end of CIplot35cfu function


# Overall Results
all.pdata <- list(
	ci.all, ci.0to4, ci.5to10, ci.11plus,
	cir.all["Overall",], cir.age0to4["Overall",], cir.age5to10["Overall",], cir.age11plus["Overall",],
	N.all, N.age0to4, N.age5to10, N.age11plus
)
# Point source results
ps.pdata <- list(
	ci.ps, ci.ps0to4, ci.ps5to10, ci.ps11plus,
	cir.all["Point source",], cir.age0to4["Point source",], cir.age5to10["Point source",], cir.age11plus["Point source",],
	N.ps, N.ps0to4, N.ps5to10, N.ps11plus
)
# Non-point source results
nps.pdata <- list(
	ci.nps, ci.nps0to4, ci.nps5to10, ci.nps11plus,
	cir.all["Non-point source",], cir.age0to4["Non-point source",], cir.age5to10["Non-point source",], cir.age11plus["Non-point source",],
	N.nps, N.nps0to4, N.nps5to10, N.nps11plus
)


## Overall Plot
pdf("~/dropbox/13beaches/aim1-results/figs/aim1-enteroQPCR-470cce-CI-byage.pdf",width=16,height=5)
cols <- brewer.pal(9,"YlOrBr")[c(8,6)]
CIplot470cce(all.pdata,cols=cols)
dev.off()

## Point source 35 CFU Plot
pdf("~/dropbox/13beaches/aim1-results/figs/aim1-enteroQPCR-470cce-CI-byage-pointsource.pdf",width=16,height=5)
cols <- brewer.pal(9,"YlOrBr")[c(8,6)]
CIplot470cce(ps.pdata,cols=cols)
dev.off()

## Non-Point source 35 CFU Plot
pdf("~/dropbox/13beaches/aim1-results/figs/aim1-enteroQPCR-470cce-CI-byage-nonpointsource.pdf",width=16,height=5)
cols <- brewer.pal(9,"YlOrBr")[c(8,6)]
CIplot470cce(nps.pdata,cols=cols)
dev.off()






# --------------------------------------
# BEACH SPECIFIC RESULTS BELOW
# --------------------------------------


# # 
# # --------------------------------------
# # Calculate CIRs from regression output
# # --------------------------------------


# # function to get Estimates and SEs from a linear combination of regression coefficients
# lccalc <- function(lc,x,vcv) {
	# # lc : linear combination of coefficients
	# # x : log-linear model object returned from coeftest (class=coeftest)
	# # vcv : variance-covariance matrix of coefficients for robust SEs
	# est <- exp(t(lc)%*%x[,1])
	# se  <- sqrt( t(lc)%*%vcv%*%lc )
	# lb <- exp(log(est)-1.96*se)
	# ub <- exp(log(est)+1.96*se)
	# list(CIR=est,CIRlb=lb,CIRub=ub)
# }

# # function to bind together CIRs for Q2 - Q5
# cirbind <- function(fit) {
	# # fit : results returned from mpreg with fit and vcovCL objects
	# lc2 <- c(0,1,rep(0,nrow(fit$fit)-2))
	# lc3 <- c(0,0,1,rep(0,nrow(fit$fit)-3))
	# lc4 <- c(0,0,0,1,rep(0,nrow(fit$fit)-4))
	# lc5 <- c(0,0,0,0,1,rep(0,nrow(fit$fit)-5))
	# res <- rbind(
	# lccalc(lc2,fit$fit,fit$vcovCL),
	# lccalc(lc3,fit$fit,fit$vcovCL),
	# lccalc(lc4,fit$fit,fit$vcovCL),
	# lccalc(lc5,fit$fit,fit$vcovCL))
	# rownames(res) <- paste("Entero Q",2:5,sep="")
	# return(res)
# }

# # grab beach-specific CIRs and 95% CIs
# fitlist <- list(hufit,sifit,wpfit,wefit,avfit,bofit,dhfit,edfit,fafit,gdfit,mafit,mbfit,sufit)
# beach.cirs <- lapply(fitlist,cirbind)
# # grab overall CIRs and 95% CIs
# all.cirs <- cirbind(list(fit=all.head,vcovCL=all.VC))

# # collate the beach CIRs into 4 matrices, corresponding to Q2 - Q5
# q2.cirs <- q3.cirs <- q4.cirs <- q5.cirs <- matrix(NA,nrow=length(beach.cirs),ncol=3)
# for (i in 1:length(beach.cirs)) {
	# q2.cirs[i,] <- unlist(beach.cirs[[i]][1,])
	# q3.cirs[i,] <- unlist(beach.cirs[[i]][2,])
	# q4.cirs[i,] <- unlist(beach.cirs[[i]][3,])
	# q5.cirs[i,] <- unlist(beach.cirs[[i]][4,])
# }

# CIRlab <- c("Huntington","Silver","Washington\nPark","West","Avalon","Boqueron","Edgewater","Doheny","Fairhope","Goddard","Malibu","Mission\nBay","Surfside")


# # --------------------------------------
# # ALL AGES
# # plot CIRs in a large panel to show
# # dose-response by Enterococcus quintile
# # --------------------------------------

# pdf("~/dropbox/13beaches/aim1-results/figs/entero1600-head-CIR-dose-response-panel.pdf",height=8*2,width=5*2)
# lo <- layout(mat=matrix(1:14,nrow=7,ncol=2,byrow=TRUE))
# ytics <- c(0.2,0.5,1,2,4,8)
# cols <- brewer.pal(9,"YlGnBu")[5:8]

# for(i in 1:14) {
	
	# if(i<=2) {
		# op <- par(mar=c(1,12,4,2)+0.1)
	# } else {
		# op <- par(mar=c(1,12,2,2)+0.1)
	# }
	# if(i<=13){
		# plotmat <- matrix(unlist(beach.cirs[[i]]),nrow=4,ncol=3)
	# } else{
		# plotmat <- matrix(unlist(all.cirs),nrow=4,ncol=3) 
	# }
	
	
	# bhs <- barplot(1:4,ylim=range(ytics),log="y",yaxt="n",col=NA,border=NA)
	# axis(side=2,at=ytics,las=1)
	# segments(x0=0,x1=max(bhs+1),y0=1,y1=1,lty=2)
	# mtext(c(CIRlab,"Combined")[i],side=2,line=3,las=1,font=2)
	# mtext("aCIR",side=2,line=3,las=1,at=max(ytics))
	# # plot estimates
	# segments(x0=bhs,y0=plotmat[,2],y1=plotmat[,3],lwd=2,col=cols)
	# points(bhs,plotmat[,1],pch=21,bg="white",cex=1.75,lwd=2,col=cols)
	# # label estimates
	# text(bhs,plotmat[,1],sprintf("%1.2f",plotmat[,1]),cex=0.9,pos=4)
	
	# if(i<=2) {
		# mtext(labQs[2:5],side=3,line=0,at=bhs,col=cols,cex=0.8)
	# }
	
# }

# par(op)
# dev.off()






