##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# Plots of the RR across the range of concentration

# Results pooled across beach
##########################################


rm(list=ls())
library(ggplot2)

#######################################
# Probabiliy of illness plots - code from Ben
#######################################

# --------------------------------------
# formatting function for CIR results
# --------------------------------------
CIRformat <- function(x) {
  # x : vector of length 3 with CIR, CIRlb, Cub
  paste(sprintf("%1.2f",x[1])," (",sprintf("%1.2f",x[2]),", ",sprintf("%1.2f",x[3]),")",sep="")
} 


# --------------------------------------
# get CIs in one object for graphing dose
# response figures
# --------------------------------------
# function to get CIR Estimates and CIs from simple stratified models
getCIR <- function(x) {
  # x : log-linear model object returned from coeftest (class=coeftest)
  # NOTE: assumes exposure of interest is the first covariate and there are no interactions
  est <- exp(x[2,1])
  se  <- x[2,2]  
  lb <- exp(log(est)-1.96*se)
  ub <- exp(log(est)+1.96*se)
  res <- c(est,lb,ub)
  return(res)
}



# --------------------------------------
# general plotting function for the
# dose-response curves
# --------------------------------------
plotPy <- function(pYcurve,xtics=c(0.1,1,10,100,1000),xlab,ytics,ytics2,ytics2units,breaks,main,CIRres,Exp){
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
  # Exp  : Exposure to plot in the histogram
  
  lo <- layout(mat=matrix(1:2,nrow=2,ncol=1),heights=c(1,0.4))
  op <- par(mar=c(2,4,4,2)+0.1)
  plot(pYcurve$pX,pYcurve$bootest*100,type="n",
       ylim=range(ytics),yaxt="n",ylab="Probability of Gastrointestinal Illness (%)",
       xlim=range(log10(xtics)),xaxt="n",xlab="",
       bty="n",
  )
  axis(1,at=log10(xtics),labels=xtics,las=1)
  axis(2,at=ytics,las=1)
  segments(x0=log10(xtics),y0=rep(0,length(xtics)),y1=rep(max(ytics),length(xtics)),col="gray90")
  mtext(main,side=3,line=2,font=2,cex=1.15)
  mtext(paste("CIR for a log10 increase:",CIRres),side=3,line=0)
  lines(pYcurve$pX,pYcurve$bootest*100,lwd=1.2)
  lines(pYcurve$pX,pYcurve$boot95lb*100,lty=5)
  lines(pYcurve$pX,pYcurve$boot95ub*100,lty=5)
  
  op <- par(mar=c(5,4,1,2)+0.1)
  
  hist(Exp,breaks=breaks,
       main="",
       xlim=range(log10(xtics)),xaxt="n",xlab=xlab,
       ylim=range(ytics2),yaxt="n",ylab="",
  )
  axis(1,at=log10(xtics),labels=xtics,las=1)
  axis(2,at=ytics2,labels=ytics2/ytics2units,las=1,cex.axis=0.75)
  
  if(ytics2units>1){
    mtext(paste("N Exposed\n(",ytics2units,"s)",sep=""),side=2,line=2)
  }
  
  par(op)
  
  
}


# --------------------------------------
# general plotting function for the
# dose-response curves with effect modification
# --------------------------------------
plotPy.em <- function(pYcurve1,pYcurve0,xtics=c(0.1,1,10,100,1000),xlab,ytics,ytics2a,ytics2aunits,ytics2b,ytics2bunits,breaksa,breaksb,main,CIRres1,CIRres0,Exp1,Exp0,lab1,lab0){
  # Plotting function for an Enterococcus dose-response curve from a log-linear model
  #
  # arguments:
  # pYcurve1 : pY object returned from pY.boot in level 1 of effect modifier (see base functions for details)
  # pYcurve0 : pY object returned from pY.boot in level 0 of effect modifier  (see base functions for details)
  # xtics   : location of X-axis tics in the plot
  # xlab    : X-axis label
  # ytics1  : location and range of Y-axis tics in the dose-response plot
  # ytics2  : location and range of the Y-axis tics in the histogram plot
  # ytics2 units : scaling factor for Y-axis on the histogram
  # main    : Title of the plot (e.g., "Total Population")
  # CIRres  : text string of CIR for a log10 increase:  "CIR (CIRlb, CIRub)"
  # Exp  : Exposure to plot in the histogram
  # lab1 : Label for level 1 of effect modifier
  # lab0 : Label for level 0 of effect modifier
  
  # combine x and y for two pY objects
  comb=data.frame(pX=c(pYcurve1$pX,pYcurve0$pX))
  comb$pY=c(pYcurve1$bootest*100,pYcurve0$bootest*100)
  comb$lb=c(pYcurve1$boot95lb*100,pYcurve0$boot95lb*100)
  comb$ub=c(pYcurve1$boot95ub*100,pYcurve0$boot95ub*100)
  
  comb$em=as.factor(c(rep(1,length(pYcurve1$pX)),rep(0,length(pYcurve0$pX))))
  par(oma=c(0,0,5,0))
  lo <- layout(mat=matrix(c(1,3,2,4),nrow=2,ncol=2),heights=c(1,0.5))
  op <- par(mar=c(2,4,1,2)+0.1)
  
  plot(comb$pX[comb$em==1],comb$pY[comb$em==1],type="n",
       ylim=range(ytics),yaxt="n",ylab="Probability of Gastrointestinal Illness (%)",
       xlim=range(log10(xtics)),xaxt="n",xlab="",
       bty="n",
  )
  axis(1,at=log10(xtics),labels=xtics,las=1)
  axis(2,at=ytics,las=1)
  segments(x0=log10(xtics),y0=rep(0,length(xtics)),y1=rep(max(ytics),length(xtics)),col="gray90")
  lines(comb$pX[comb$em==1],comb$pY[comb$em==1],lwd=1.2,col="#0066CC")
  lines(comb$pX[comb$em==1],comb$lb[comb$em==1],lwd=1.2,col="#0066CC",lty=5)
  lines(comb$pX[comb$em==1],comb$ub[comb$em==1],lwd=1.2,col="#0066CC",lty=5)
  mtext(paste(lab1),side=3,line=2,font=2)  
  mtext(paste("CIR for a log10 increase:",CIRres1),side=3,line=0)
  
  
  plot(comb$pX[comb$em==0],comb$pY[comb$em==0],col=comb$em,type="n",
       ylim=range(ytics),yaxt="n",ylab="Probability of Gastrointestinal Illness (%)",
       xlim=range(log10(xtics)),xaxt="n",xlab="",
       bty="n",
  )
  axis(1,at=log10(xtics),labels=xtics,las=1)
  axis(2,at=ytics,las=1)
  segments(x0=log10(xtics),y0=rep(0,length(xtics)),y1=rep(max(ytics),length(xtics)),col="gray90")
  lines(comb$pX[comb$em==0],comb$pY[comb$em==0],lwd=1.2,col="#993300")
  lines(comb$pX[comb$em==0],comb$lb[comb$em==0],lwd=1.2,col="#993300",lty=5)
  lines(comb$pX[comb$em==0],comb$ub[comb$em==0],lwd=1.2,col="#993300",lty=5)
  mtext(paste(lab0),side=3,line=2,font=2)  
  mtext(paste("CIR for a log10 increase:",CIRres0),side=3,line=0)
  
  
  mtext(main,font=2,cex=1.15,outer=TRUE,line=3)
  
  op <- par(mar=c(5,4,1,2)+0.1)
  
  hist(Exp1,breaks=breaksa,
       main="",
       xlim=range(log10(xtics)),xaxt="n",xlab="",
       ylim=range(ytics2a),yaxt="n",ylab="",
  )
  axis(1,at=log10(xtics),labels=xtics,las=1)
  axis(2,at=ytics2a,labels=ytics2a/ytics2aunits,las=1,cex.axis=0.75)
  
  hist(Exp0,breaks=breaksb,
       main="",
       xlim=range(log10(xtics)),xaxt="n",xlab="",
       ylim=range(ytics2b),yaxt="n",ylab="",
  )
  axis(1,at=log10(xtics),labels=xtics,las=1)
  axis(2,at=ytics2b,labels=ytics2b/ytics2bunits,las=1,cex.axis=0.75)
  
  mtext(xlab,outer=TRUE,side=1,line=-2)
  
  if(ytics2aunits>1){
    mtext(paste("N Exposed\n(",ytics2aunits,"s)",sep=""),side=2,line=2)
  }
  
  par(op)
  
}


# --------------------------------------
# load the and pre-preprocess the 
# analysis dataset
# (refer to the base functions script
# for details on the pre-processing)
# --------------------------------------
beaches13=read.csv("~/Documents/CRG/coliphage/13beaches-data/final/13beaches-analysis.csv")

# load base functions
source("~/Documents/CRG/coliphage/13beaches-coliphage/src/Analysis/0-base-functions.R")

data=preprocess.6beaches(beaches13)

# restrict to 6 beaches with coliphage data
beach.list=c("Avalon","Doheny","Malibu","Mission Bay",
             "Fairhope","Goddard")

all=data[data$beach %in% beach.list,]

all=subset(all,nowq==0)
 # subset to non-missing exposure categories
# to make the robust CI calcs work
all=subset(all,all$bodycontact=="Yes")

load("~/Documents/CRG/coliphage/results/rawoutput/regress-10day-continuous-body.Rdata")



# --------------------------------------
# make figures
# --------------------------------------

pdf("~/Documents/CRG/coliphage/results/figures/dose-response-pool-fmc1601.pdf",height=7,width=5)
plotPy(all.fmc1601.pY,
       xlab="Concentration PFU / 100ml",
       ytics=seq(0,15,by=1),
       ytics2=c(0:3)*250,
       ytics2units=1,
       breaks=40,
       main="Somatic Coliphage (EPA 1601)",
       CIRres=CIRformat(getCIR(overall.fit10.fmc1601)),
       Exp=all$fmc1601[all$fmc1601>-1]
)
dev.off()

pdf("~/Documents/CRG/coliphage/results/figures/dose-response-pool-fmc1602-all.pdf",height=7,width=5)
plotPy(all.fmc1602.pY,
       xlab="Concentration PFU / 100ml",
       ytics=seq(0,35,by=5),
       ytics2=seq(0,1000,by=200),
       ytics2units=1,
       breaks=40,
       main="Somatic Coliphage (EPA 1601)",
       CIRres=CIRformat(getCIR(overall.fit10.fmc1602)),
       Exp=all$fmc1602[all$fmc1602>-1]
)
dev.off()

pdf("~/Documents/CRG/coliphage/results/figures/dose-response-pool-fmc1602.pdf",height=6,width=7)
plotPy.em(all.fmc1602.pY.low,all.fmc1602.pY.high,
          xlab="Concentration PFU / 100ml",
          xtics=c(0.1,1,10,100,1000,10000),
          ytics=seq(0,50,by=10),
          ytics2a=c(0:3)*250,
          ytics2aunits=1,
          ytics2b=c(0:3)*250,
          ytics2bunits=1,
          breaksa=30,
          breaksb=20,
          main="Somatic Coliphage (EPA 1602)",
          CIRres1=CIRformat(getCIR(overall.fit10.fmc1602.low)),
          CIRres0=CIRformat(getCIR(overall.fit10.fmc1602.high)),
          Exp1=all$fmc1602[all$risk=="Low" & all$fmc1602>-1],
          Exp0=all$fmc1602[all$risk=="High" & all$fmc1602>-1],
          lab1="Low risk conditions",lab0="High risk conditions"
)
dev.off()

pdf("~/Documents/CRG/coliphage/results/figures/dose-response-pool-fpc1601-all.pdf",height=7,width=5)
plotPy(all.fpc1601.pY,
       xlab="Concentration PFU / 100ml",
       ytics=seq(0,35,by=5),
       ytics2=seq(0,1000,by=200),
       ytics2units=1,
       breaks=40,
       main="Male-Specific Coliphage (EPA 1601)",
       CIRres=CIRformat(getCIR(overall.fit10.fpc1601)),
       Exp=all$fpc1601[all$fpc1601>-1]
)
dev.off()

pdf("~/Documents/CRG/coliphage/results/figures/dose-response-pool-fpc1601.pdf",height=6,width=7)
plotPy.em(all.fpc1601.pY.low,all.fpc1601.pY.high,
          xlab="Concentration PFU / 100ml",
          xtics=c(0.1,1,10,100,1000),
          ytics=seq(0,35,by=5),
          ytics2a=c(0:5)*250,
          ytics2aunits=1,
          ytics2b=c(0:5)*250,
          ytics2bunits=1,
          breaksa=35,
          breaksb=15,
          main="Male-Specific Coliphage (EPA 1601)",
          CIRres1=CIRformat(getCIR(overall.fit10.fpc1601.low)),
          CIRres0=CIRformat(getCIR(overall.fit10.fpc1601.high)),
          Exp1=all$fmc1601[all$risk=="Low" & all$fpc1601>-1],
          Exp0=all$fmc1601[all$risk=="High" & all$fpc1601>-1],
          lab1="Low risk conditions",lab0="High risk conditions"
)
dev.off()

pdf("~/Documents/CRG/coliphage/results/figures/dose-response-pool-fpc1602-all.pdf",height=7,width=5)
plotPy(all.fpc1602.pY,
       xlab="Concentration PFU / 100ml",
       ytics=seq(0,110,by=20),
       ytics2=c(0:4)*100,
       ytics2units=1,
       breaks=40,
       main="Male-Specific Coliphage (EPA 1602)",
       CIRres=CIRformat(getCIR(overall.fit10.fpc1602)),
       Exp=all$fpc1602[all$fpc1602>-1]
)
dev.off()

pdf("~/Documents/CRG/coliphage/results/figures/dose-response-pool-fpc1602.pdf",height=6,width=7)
plotPy.em(all.fpc1602.pY.low,all.fpc1602.pY.high,
          xlab="Concentration PFU / 100ml",
          xtics=c(0.1,1,10,100,1000),
          ytics=seq(0,110,by=20),
          ytics2a=c(0:5)*250,
          ytics2aunits=1,
          ytics2b=c(0:5)*250,
          ytics2bunits=1,
          breaksa=30,
          breaksb=25,
          main="Male-Specific Coliphage (EPA 1602)",
          CIRres1=CIRformat(getCIR(overall.fit10.fpc1602.low)),
          CIRres0=CIRformat(getCIR(overall.fit10.fpc1602.high)),
          Exp1=all$fpc1602[all$risk=="Low" & all$fpc1602>-1],
          Exp0=all$fpc1602[all$risk=="High" & all$fpc1602>-1],
          lab1="Low risk conditions",lab0="High risk conditions"
)
dev.off()
