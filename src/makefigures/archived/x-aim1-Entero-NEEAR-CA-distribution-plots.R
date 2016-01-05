



# --------------------------------------
# x-aim1-Entero-NEEAR-CA-distribution-plots.R
# ben arnold (benarnold@berkeley.edu)
#
# description:
# exploratory / confirmatory analysis:
# plot the distributions and summary 
# statistics for the Entero qPCR
# results in the NEEAR versus CA studies
#
# rationale: confirm that the distributions
# look similar and reasonable to ensure that
# pooling the data makes sense.
#
# version 2 (10 apr 2015)
# added the 2-way scatterplots
#
# version 1 (8 apr 2015)
#
# --------------------------------------

# --------------------------------------
# input files:
#    13beaches-wq.csv
#
# output files:
#    Entero-distributions-NEEAR-v-CA.pdf
# --------------------------------------


# --------------------------------------
# preamble
# --------------------------------------

rm(list=ls())
library(RColorBrewer)
library(scales)

# --------------------------------------
# load the water quality dataset
# --------------------------------------

# d <- read.csv("~/dropbox/13beaches/data/final/13beaches-analysis.csv")
wq <-read.csv("~/dropbox/13beaches/data/final/13beaches-wq.csv")

# drop Site C samples from Doheny and Malibu
# (Lagoon samples are not representative of the 
# main beach conditions that individuals were exposed to)
wq <- subset(wq,beachcode!="Doheny-C" & beachcode!="Malibu-C")

# create an indicator for NEEAR studies
wq$neear <- ifelse(wq$beach!="Avalon" & wq$beach!="Doheny" & wq$beach!="Malibu" & wq$beach!="Mission Bay",1,0)

# restrict to non-missing values on qPCR and 1600 to make density plotting easier
table(is.na(wq$avgdyenteropcr))
wq1 <- subset(wq,is.na(avgdyenteropcr)==FALSE)
dim(wq1)

table(is.na(wq$avgdyentero1600))
wq2 <- subset(wq,is.na(avgdyentero1600)==FALSE)
dim(wq2)


# --------------------------------------
# make the plot
# --------------------------------------

pdf("~/dropbox/13beaches/aim1-results/figs/Entero-distributions-NEEAR-v-CA.pdf",width=12,height=6)

# common plot parameters and layout
xtics <- c(-1,0,1,2,3,4)
ytics <- seq(0,1,by=0.1)
cols <- c("black",brewer.pal(8,"Dark2")[1])
lo <- layout(mat=matrix(1:2,nrow=1,ncol=2))

# Overlay the two densties of 1600 values for NEEAR and CA studies
plot(density(wq2$avgdyentero1600),type="n",bty="n",
	main="EPA 1600",
	xlab="",xaxt="n",xlim=range(xtics),
	ylab="",yaxt="n",ylim=range(ytics),
	las=1
	)
	axis(1,at=xtics,las=1)
	axis(2,at=ytics,las=1)
	mtext("Log 10 Enterococcus 1600 CFU/100ml",side=1,line=2.5)
	mtext("Kernel Density",side=3,line=0,adj=0,at=-1.8)
	lines(density(wq2$avgdyentero1600[wq2$neear==1]),lwd=2,col=cols[1])
	lines(density(wq2$avgdyentero1600[wq2$neear==0]),lwd=2,col=cols[2])

	text(0,0.6,"NEEAR Studies",font=2,col=cols[1])
	text(2.5,0.9,"CA Studies",font=2,col=cols[2])



# Overlay the two densties of qPCR values for NEEAR and CA studies

plot(density(wq1$avgdyenteropcr),type="n",bty="n",
	main="EPA qPCR 1611",
	xlab="",xaxt="n",xlim=range(xtics),
	ylab="",yaxt="n",ylim=range(ytics),
	las=1
	)
	axis(1,at=xtics,las=1)
	axis(2,at=ytics,las=1)
	mtext("Log 10 Enterococcus qPCR CCE/100ml",side=1,line=2.5)
	mtext("Kernel Density",side=3,line=0,adj=0,at=-1.8)
	lines(density(wq1$avgdyenteropcr[wq1$neear==1]),lwd=2,col=cols[1])
	lines(density(wq1$avgdyenteropcr[wq1$neear==0]),lwd=2,col=cols[2])

	text(3,0.9,"NEEAR Studies",font=2,col=cols[1])
	text(0.5,0.5,"CA Studies",font=2,col=cols[2])

dev.off()


# --------------------------------------
# create a scatter plot of 1600 vs qPCR
# --------------------------------------
wq$entero470 <- ifelse(10^(wq$avgdyenteropcr)>470,1,0)

table(wq$neear,wq$entero470)

pdf("~/dropbox/13beaches/aim1-results/figs/Entero-XY-NEEAR-v-CA.pdf",width=12,height=4)

lo <- layout(mat=matrix(1:3,nrow=1,ncol=3))
cols <- brewer.pal(8,"Dark2")[c(8,1)]
ytics <- xtics <- seq(-1,4,by=1)
# combined plot
plot(wq$avgdyentero1600,wq$avgdyenteropcr,type="n",
	xaxt="n",xlab="log10 Enterococcus EPA 1600, CFU/100 ml",xlim=range(xtics),
	yaxt="n",ylab="log10 Enterococcus EPA qPCR 1611, CCE/100 ml",ylim=range(ytics),
	bty="n",las=1,main="Combined Plot"
	)
	axis(1,at=xtics,las=1)
	axis(2,at=ytics,las=1)
	segments(x0=log10(35),y0=min(ytics),y1=max(ytics),lty=2,lwd=1.5,col="gray40")
	segments(y0=log10(470),x0=min(xtics),x1=max(xtics),lty=2,lwd=1.5,col="gray40")
	text(x=max(xtics),y=log10(470)+0.2,"470 CCE",adj=1,col="gray20")
	text(x=log10(35),y=max(ytics)-0.2,"35 CFU",pos=4,col="gray20")
	
	points(wq$avgdyentero1600[wq$neear==1],wq$avgdyenteropcr[wq$neear==1],pch=16,cex=0.7,col=alpha(cols[1],0.5))
	points(wq$avgdyentero1600[wq$neear==0],wq$avgdyenteropcr[wq$neear==0],pch=16,cex=0.7,col=alpha(cols[2],0.5))

	legend("bottomright",legend=c("NEEAR Studies","CA Studies"),pch=16,col=alpha(cols,0.5),bty="n")
	
# NEEAR plot
plot(wq$avgdyentero1600,wq$avgdyenteropcr,type="n",
	xaxt="n",xlab="log10 Enterococcus EPA 1600, CFU/100 ml",xlim=range(xtics),
	yaxt="n",ylab="log10 Enterococcus EPA qPCR 1611, CCE/100 ml",ylim=range(ytics),
	bty="n",las=1,main="NEEAR Studies"
	)
	axis(1,at=xtics,las=1)
	axis(2,at=ytics,las=1)
	segments(x0=log10(35),y0=min(ytics),y1=max(ytics),lty=2,lwd=1.5,col="gray40")
	segments(y0=log10(470),x0=min(xtics),x1=max(xtics),lty=2,lwd=1.5,col="gray40")
	text(x=max(xtics),y=log10(470)+0.2,"470 CCE",adj=1,col="gray20")
	text(x=log10(35),y=max(ytics)-0.2,"35 CFU",pos=4,col="gray20")
	
	points(wq$avgdyentero1600[wq$neear==1],wq$avgdyenteropcr[wq$neear==1],pch=16,cex=0.7,col=alpha(cols[1],0.5))

# CA plot
plot(wq$avgdyentero1600,wq$avgdyenteropcr,type="n",
	xaxt="n",xlab="log10 Enterococcus EPA 1600, CFU/100 ml",xlim=range(xtics),
	yaxt="n",ylab="log10 Enterococcus EPA qPCR 1611, CCE/100 ml",ylim=range(ytics),
	bty="n",las=1,main="CA Studies"
	)
	axis(1,at=xtics,las=1)
	axis(2,at=ytics,las=1)
	segments(x0=log10(35),y0=min(ytics),y1=max(ytics),lty=2,lwd=1.5,col="gray40")
	segments(y0=log10(470),x0=min(xtics),x1=max(xtics),lty=2,lwd=1.5,col="gray40")
	text(x=max(xtics),y=log10(470)+0.2,"470 CCE",adj=1,col="gray20")
	text(x=log10(35),y=max(ytics)-0.2,"35 CFU",pos=4,col="gray20")
	
	points(wq$avgdyentero1600[wq$neear==0],wq$avgdyenteropcr[wq$neear==0],pch=16,cex=0.7,col=alpha(cols[2],0.5))

			
dev.off()




