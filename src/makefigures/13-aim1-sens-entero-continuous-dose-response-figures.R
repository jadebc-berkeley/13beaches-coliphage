# --------------------------------------
# 13-aim1-sens-entero-continuous-dose-response-figures.R
# ben arnold (benarnold@berkeley.edu)
#
# description:
# plot model-predicted
# exposure-response curves for Enterococcus
# EPA 1600  and EPA qPCR 1611 
#
# --------------------------------------

# --------------------------------------
# input files:
#	13beaches-analysis.csv
#	aim1-entero1600-continuous-regs-body.Rdata
#	aim1-enteroQPCR-continuous-regs-body.Rdata
#
# output files:
#	aim1-sens-entero1600-dose-response-curve.pdf
# aim1-sens-entero1600-dose-response-curve-age0to4.pdf
# aim1-sens-entero1600-dose-response-curve-age5to10.pdf
# aim1-sens-entero1600-dose-response-curve-age11plus.pdf
# aim1-sens-entero1600-dose-response-curve-ps.pdf
# aim1-sens-entero1600-dose-response-curve-nps.pdf
#
# aim1-sens-enteroQPCR-dose-response-curve.pdf
# aim1-sens-enteroQPCR-dose-response-curve-age0to4.pdf
# aim1-sens-enteroQPCR-dose-response-curve-age5to10.pdf
# aim1-sens-enteroQPCR-dose-response-curve-age11plus.pdf
# aim1-sens-enteroQPCR-dose-response-curve-ps.pdf
# aim1-sens-enteroQPCR-dose-response-curve-nps.pdf
# --------------------------------------

# --------------------------------------
# preamble
# --------------------------------------

rm(list=ls())
library(RColorBrewer)

source("~/dropbox/13beaches/src/aim1/0-aim1-base-functions.R")


# --------------------------------------
# general plotting function for the
# dose-response curves
# --------------------------------------

plotPy <- function(pYcurve,xtics=c(0.1,1,10,100,1000,10000),xlab,ytics,ytics2,ytics2units,main,CIRres,Entero){
	# Plotting function for an Enterococcus dose-response curve from a log-linear model
	#
	# arguments:
	# pYcurve : pY object returned from pY.boot (see base functions for details)
	# xtics   : location of X-axis tics in the plot
	# xlab    : X-axis label
	# ytics1  : location and range of Y-axis tics in the dose-response plot
	# ytics2  : location and range of the Y-axis tics in the histogram plot
	# ytics2 units : scaling factor for Y-axis on the histogram
	# main    : Title of the plot (e.g., "Total Population")
	# CIRres  : text string of CIR for a log10 increase:  "CIR (CIRlb, CIRub)"
	# Entero  : Enterococcus exposure to plot in the histogram
	
	lo <- layout(mat=matrix(1:2,nrow=2,ncol=1),heights=c(1,0.4))
	op <- par(mar=c(2,4,4,2)+0.1)
	plot(pYcurve$pX,pYcurve$bootest*100,type="n",
		ylim=range(ytics),yaxt="n",ylab="Probability of Diarrhea (%)",
		xlim=range(log10(xtics)),xaxt="n",xlab="",
		bty="n",
		)
		axis(1,at=log10(xtics),labels=xtics,las=1)
		axis(2,at=ytics,las=1)
		segments(x0=log10(xtics),y0=rep(0,length(xtics)),y1=rep(max(ytics),length(xtics)),col="gray90")
		mtext(main,side=3,line=2,font=2,cex=1.5)
		mtext(paste("CIR for a log10 increase:",CIRres),side=3,line=0)
	lines(pYcurve$pX,pYcurve$bootest*100,lwd=1.2)
	lines(pYcurve$pX,pYcurve$boot95lb*100,lty=5)
	lines(pYcurve$pX,pYcurve$boot95ub*100,lty=5)
	
	op <- par(mar=c(5,4,1,2)+0.1)
	
	hist(Entero,breaks=50,
		main="",
		xlim=range(log10(xtics)),xaxt="n",xlab=xlab,
		ylim=range(ytics2),yaxt="n",ylab="",
		)
		axis(1,at=log10(xtics),labels=xtics,las=1)
		axis(2,at=ytics2,labels=ytics2/ytics2units,las=1,cex.axis=0.75)
		mtext(paste("N Exposed\n(",ytics2units,"s)",sep=""),side=2,line=2)
		
	par(op)
	
	
}

# --------------------------------------
# formatting function for CIR results
# --------------------------------------
CIRformat <- function(x) {
	# x : vector of length 3 with CIR, CIRlb, Cub
	paste(sprintf("%1.2f",x[1])," (",sprintf("%1.2f",x[2]),", ",sprintf("%1.2f",x[3]),")",sep="")
} 



# --------------------------------------
# Load and pre-process analysis dataset to 
# enable plots of the distribution of 
# Enterococcus exposure
# (refer to the base functions script
# for details on the pre-processing)
# --------------------------------------
ad <- preprocess.13beaches("~/dropbox/13beaches/data/final/13beaches-analysis.csv")

# drop individuals with no water quality information
table(ad$nowq)
ad <- subset(ad,nowq==0)
dim(ad)

# subset to body immersion swimmers
table(ad$bodycontact)
ad <- subset(ad,ad$bodycontact=="Yes")
	dim(ad)

# --------------------------------------
# EPA 1600 results
# --------------------------------------

load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-entero1600-continuous-regs-body.Rdata")


# Total Population
pdf("~/dropbox/13beaches/aim1-results/figs/aim1-sens-entero1600-dose-response-curve.pdf",height=7,width=5)
plotPy(pYcurve.all,
	xlab="Enterococcus EPA 1600 Concentration CFU / 100ml",
	ytics=seq(0,12,by=1),
	ytics2=c(0:4)*1000,
	ytics2units=1000,
	main="Total Population",
	CIRres=CIRformat(cir.age[4,]),
	Entero=ad$entero1600
)
dev.off()


# Ages 0 to 4
pdf("~/dropbox/13beaches/aim1-results/figs/aim1-sens-entero1600-dose-response-curve-age0to4.pdf",height=7,width=5)
plotPy(pYcurve.age0to4,
	xlab="Enterococcus EPA 1600 Concentration CFU / 100ml",
	ytics=seq(0,12,by=1),
	ytics2=c(0:3)*100,
	ytics2units=100,
	main="Ages 0 to 4",
	CIRres=CIRformat(cir.age[1,]),
	Entero=ad$entero1600[ad$agestrat=="(0, 4]"]
)
dev.off()

# Ages 5 to 10
pdf("~/dropbox/13beaches/aim1-results/figs/aim1-sens-entero1600-dose-response-curve-age5to10.pdf",height=7,width=5)
plotPy(pYcurve.age5to10,
	xlab="Enterococcus EPA 1600 Concentration CFU / 100ml",
	ytics=seq(0,12,by=1),
	ytics2=seq(0,6,by=2)*100,
	ytics2units=100,
	main="Ages 5 to 10",
	CIRres=CIRformat(cir.age[2,]),
	Entero=ad$entero1600[ad$agestrat=="(4, 10]"]
)
dev.off()

# Ages > 10
pdf("~/dropbox/13beaches/aim1-results/figs/aim1-sens-entero1600-dose-response-curve-age11plus.pdf",height=7,width=5)
plotPy(pYcurve.age11plus,
	xlab="Enterococcus EPA 1600 Concentration CFU / 100ml",
	ytics=seq(0,12,by=1),
	ytics2=c(0:3)*1000,
	ytics2units=1000,
	main="Ages >10",
	CIRres=CIRformat(cir.age[3,]),
	Entero=ad$entero1600[ad$agestrat==">10"]
)
dev.off()

# Point Source Beaches
pdf("~/dropbox/13beaches/aim1-results/figs/aim1-sens-entero1600-dose-response-curve-ps.pdf",height=7,width=5)
plotPy(pYcurve.ps,
	xlab="Enterococcus EPA 1600 Concentration CFU / 100ml",
	ytics=seq(0,12,by=1),
	ytics2=seq(0,3,by=1)*1000,
	ytics2units=1000,
	main="Point Source Beaches",
	CIRres=CIRformat(cir.ps[1,]),
	Entero=ad$entero1600[ad$pointsource=="Yes"]
)
dev.off()


# Non-Point Source Beaches
pdf("~/dropbox/13beaches/aim1-results/figs/aim1-sens-entero1600-dose-response-curve-nps.pdf",height=7,width=5)
plotPy(pYcurve.nps,
	xlab="Enterococcus EPA 1600 Concentration CFU / 100ml",
	ytics=seq(0,12,by=1),
	ytics2=seq(0,3,by=1)*1000,
	ytics2units=1000,
	main="Non-Point Source Beaches",
	CIRres=CIRformat(cir.ps[2,]),
	Entero=ad$entero1600[ad$pointsource=="No"]
)
dev.off()

# --------------------------------------
# EPA 1611 qPCR results
# --------------------------------------

load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-enteroQPCR-continuous-regs-body.Rdata")


# Total Population
pdf("~/dropbox/13beaches/aim1-results/figs/aim1-sens-enteroQPCR-dose-response-curve.pdf",height=7,width=5)
plotPy(pYcurve.all,
	xlab="Enterococcus EPA 1611 qPCR Concentration CCE / 100ml",
	ytics=seq(0,12,by=1),
	ytics2=seq(0,3,by=1)*1000,
	ytics2units=1000,
	main="Total Population",
	CIRres=CIRformat(cir.age[4,]),
	Entero=ad$enteroQPCR
)
dev.off()


# Ages 0 to 4
pdf("~/dropbox/13beaches/aim1-results/figs/aim1-sens-enteroQPCR-dose-response-curve-age0to4.pdf",height=7,width=5)
plotPy(pYcurve.age0to4,
	xlab="Enterococcus EPA 1611 qPCR Concentration CCE / 100ml",
	ytics=seq(0,12,by=1),
	ytics2=seq(0,3,by=1)*100,
	ytics2units=100,
	main="Ages 0 to 4",
	CIRres=CIRformat(cir.age[1,]),
	Entero=ad$enteroQPCR[ad$agestrat=="(0, 4]"]
)
dev.off()

# Ages 5 to 10
pdf("~/dropbox/13beaches/aim1-results/figs/aim1-sens-enteroQPCR-dose-response-curve-age5to10.pdf",height=7,width=5)
plotPy(pYcurve.age5to10,
	xlab="Enterococcus EPA 1611 qPCR Concentration CCE / 100ml",
	ytics=seq(0,12,by=1),
	ytics2=seq(0,6,by=2)*100,
	ytics2units=100,
	main="Ages 5 to 10",
	CIRres=CIRformat(cir.age[2,]),
	Entero=ad$enteroQPCR[ad$agestrat=="(4, 10]"]
)
dev.off()

# Ages > 10
pdf("~/dropbox/13beaches/aim1-results/figs/aim1-sens-enteroQPCR-dose-response-curve-age11plus.pdf",height=7,width=5)
plotPy(pYcurve.age11plus,
	xlab="Enterococcus EPA 1611 qPCR Concentration CCE / 100ml",
	ytics=seq(0,12,by=1),
	ytics2=seq(0,2,by=1)*1000,
	ytics2units=1000,
	main="Ages >10",
	CIRres=CIRformat(cir.age[3,]),
	Entero=ad$enteroQPCR[ad$agestrat==">10"]
)
dev.off()

# Point Source Beaches
pdf("~/dropbox/13beaches/aim1-results/figs/aim1-sens-enteroQPCR-dose-response-curve-ps.pdf",height=7,width=5)
plotPy(pYcurve.ps,
	xlab="Enterococcus EPA 1611 qPCR Concentration CCE / 100ml",
	ytics=seq(0,12,by=1),
	ytics2=seq(0,3,by=1)*1000,
	ytics2units=1000,
	main="Point Source Beaches",
	CIRres=CIRformat(cir.ps[1,]),
	Entero=ad$enteroQPCR[ad$pointsource=="Yes"]
)
dev.off()


# Non-Point Source Beaches
pdf("~/dropbox/13beaches/aim1-results/figs/aim1-sens-enteroQPCR-dose-response-curve-nps.pdf",height=7,width=5)
plotPy(pYcurve.nps,
	xlab="Enterococcus EPA 1611 qPCR Concentration CCE / 100ml",
	ytics=seq(0,12,by=1),
	ytics2=seq(0,2,by=1)*1000,
	ytics2units=1000,
	main="Non-Point Source Beaches",
	CIRres=CIRformat(cir.ps[2,]),
	Entero=ad$enteroQPCR[ad$pointsource=="No"]
)
dev.off()


