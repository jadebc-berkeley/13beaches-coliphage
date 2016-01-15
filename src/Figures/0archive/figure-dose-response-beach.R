##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# Plots of the RR across the range of concentration

# Results stratified by beach
##########################################


rm(list=ls())
library(ggplot2)
library(foreign)

setwd("~/Dropbox/Coliphage/")

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
plotPy.em <- function(pYcurve1,pYcurve0,xtics=c(0.1,1,10,100,1000),xlab,ytics,ytics2,ytics2units,breaks,main,CIRres1,CIRres0,Exp1,Exp0,lab1,lab0){
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
  lines(comb$pX[comb$em==1],comb$pY[comb$em==1],lwd=1.2,col="#993300")
  lines(comb$pX[comb$em==1],comb$lb[comb$em==1],lwd=1.2,col="#993300",lty=5)
  lines(comb$pX[comb$em==1],comb$ub[comb$em==1],lwd=1.2,col="#993300",lty=5)
  mtext(paste(lab1),side=3,line=2,font=2)  
  mtext(paste("CIR for a log10 increase:",CIRres1),side=3,line=0)
  
  
  plot(comb$pX[comb$em==0],comb$pY[comb$em==0],col=comb$em,type="n",
       ylim=range(ytics),yaxt="n",ylab="Probability of Diarrhea (%)",
       xlim=range(log10(xtics)),xaxt="n",xlab="",
       bty="n",
  )
  axis(1,at=log10(xtics),labels=xtics,las=1)
  axis(2,at=ytics,las=1)
  segments(x0=log10(xtics),y0=rep(0,length(xtics)),y1=rep(max(ytics),length(xtics)),col="gray90")
  lines(comb$pX[comb$em==0],comb$pY[comb$em==0],lwd=1.2,col="#0066CC")
  lines(comb$pX[comb$em==0],comb$lb[comb$em==0],lwd=1.2,col="#0066CC",lty=5)
  lines(comb$pX[comb$em==0],comb$ub[comb$em==0],lwd=1.2,col="#0066CC",lty=5)
  mtext(paste(lab0),side=3,line=2,font=2)  
  mtext(paste("CIR for a log10 increase:",CIRres0),side=3,line=0)
  
  
  mtext(main,font=2,cex=1.15,outer=TRUE,line=3)
  
  op <- par(mar=c(5,4,1,2)+0.1)
  
  hist(Exp1,breaks=breaks,
       main="",
       xlim=range(log10(xtics)),xaxt="n",xlab="",
       ylim=range(ytics2),yaxt="n",ylab="",
  )
  axis(1,at=log10(xtics),labels=xtics,las=1)
  axis(2,at=ytics2,labels=ytics2/ytics2units,las=1,cex.axis=0.75)
  
  hist(Exp0,breaks=breaks,
       main="",
       xlim=range(log10(xtics)),xaxt="n",xlab="",
       ylim=range(ytics2),yaxt="n",ylab="",
  )
  axis(1,at=log10(xtics),labels=xtics,las=1)
  axis(2,at=ytics2,labels=ytics2/ytics2units,las=1,cex.axis=0.75)
  
  mtext(xlab,outer=TRUE,side=1,line=-2)
  
  if(ytics2units>1){
    mtext(paste("N Exposed\n(",ytics2units,"s)",sep=""),side=2,line=2)
  }
  
  par(op)
  
}


# --------------------------------------
# load the and pre-preprocess the 
# analysis dataset
# (refer to the base functions script
# for details on the pre-processing)
# --------------------------------------
beaches13=read.dta("~/Dropbox/13beaches/data/final/13beaches-analysis.dta")

# load base functions
source("Programs/Analysis/0-base-functions.R")

data=preprocess.6beaches(beaches13)

# restrict to 6 beaches with coliphage data
beach.list=c("Avalon","Doheny","Malibu","Mission Bay",
             "Fairhope","Goddard")

all=data[data$beach %in% beach.list,]
avalon=data[data$beach %in% "Avalon",]
doheny=data[data$beach %in% "Doheny",]
malibu=data[data$beach %in% "Malibu",]
mission=data[data$beach %in% "Mission Bay",]
fairhope=data[data$beach %in% "Fairhope",]
goddard=data[data$beach %in% "Goddard",]

data.list=list(all=all,avalon=avalon,doheny=doheny,mission=mission,
               malibu=malibu,goddard=goddard,fairhope=fairhope)

data.list=lapply(data.list,function(df){
  # drop individuals with no water quality information
  df=subset(df,nowq==0)
  # subset to non-missing exposure categories
  # to make the robust CI calcs work
  df=subset(df,df$bodycontact=="Yes")
})

# convert from list back to data frames
list2env(data.list ,.GlobalEnv)

load("~/dropbox/coliphage/results/rawoutput/regress-3day-continuous-body.Rdata")


# --------------------------------------
# make figures
# --------------------------------------
  pdf("~/dropbox/coliphage/results/figures/dose-response-av-fmc1602.pdf",height=6,width=7)
  plotPy.em(av.fmc1602.pY.gw.1,av.fmc1602.pY.gw.0,
          xlab="Concentration PFU / 100ml",
          xtics=c(0.1,1,10,100,1000,10000),
          ytics=seq(0,35,by=5),
          ytics2=c(0:4)*250,
          ytics2units=1,
          breaks=45,
          main="Avalon Beach, F- Coliphage (EPA 1602)",
          CIRres1=CIRformat(getCIR(avfit10.fmc1602.pY.gw.1$fit)),
          CIRres0=CIRformat(getCIR(avfit10.fmc1602.pY.gw.0$fit)),
          Exp1=avalon$fmc1602[avalon$groundwater=="Above median flow"],
          Exp0=avalon$fmc1602[avalon$groundwater=="Below median flow"],
          lab1="Groundwater above median",lab0="Groundwater below median"
  )
  dev.off()

  pdf("~/dropbox/coliphage/results/figures/dose-response-av-fpc1601.pdf",height=6,width=7)
  plotPy.em(av.fpc1601.pY.gw.1,av.fpc1601.pY.gw.0,
         xlab="Concentration PFU / 100ml",
         ytics=seq(0,70,by=10),
         ytics2=c(0:3)*250,
         ytics2units=1,
         breaks=25,
         main="Avalon Beach, F+ Coliphage (EPA 1601)",
         CIRres1=CIRformat(getCIR(avfit10.fpc1601.pY.gw.1$fit)),
         CIRres0=CIRformat(getCIR(avfit10.fpc1601.pY.gw.0$fit)),
         Exp1=avalon$fpc1601[avalon$groundwater=="Above median flow"],
         Exp0=avalon$fpc1601[avalon$groundwater=="Below median flow"],
         lab1="Groundwater above median",lab0="Groundwater below median"
  )
  dev.off()
  
  pdf("~/dropbox/coliphage/results/figures/dose-response-av-fpc1602.pdf",height=6,width=7)
  plotPy.em(av.fpc1602.pY.gw.1,av.fpc1602.pY.gw.0,
         xlab="Concentration PFU / 100ml",
         xtics=c(0.1,1,10,100),
         ytics=seq(0,90,by=10),
         ytics2=c(0:4)*500,
         ytics2units=1,
         breaks=20,
         main="Avalon Beach, F+ Coliphage (EPA 1602)",
         CIRres1=CIRformat(getCIR(avfit10.fpc1602.pY.gw.1$fit)),
         CIRres0=CIRformat(getCIR(avfit10.fpc1602.pY.gw.0$fit)),
         Exp1=avalon$fpc1602[avalon$groundwater=="Above median flow"],
         Exp0=avalon$fpc1602[avalon$groundwater=="Below median flow"],
         lab1="Groundwater above median",lab0="Groundwater below median"
  )
  dev.off()
  
  pdf("~/dropbox/coliphage/results/figures/dose-response-do-fmc1601.pdf",height=7,width=5)
  plotPy(do.fmc1601.pY,
         xlab="Concentration PFU / 100ml",
         ytics=seq(0,15,by=1),
         ytics2=c(0:2)*250,
         ytics2units=1,
         breaks=50,
         main="Doheny Beach, F- Coliphage (EPA 1601)",
         CIRres=CIRformat(getCIR(dofit10.fmc1601.pY$fit)),
         Exp=doheny$fmc1601
  )
  dev.off()
  
  pdf("~/dropbox/coliphage/results/figures/dose-response-do-fmc1602.pdf",height=6,width=7)
  plotPy.em(do.fmc1602.pY.berm.1,do.fmc1602.pY.berm.0,
         xlab="Concentration PFU / 100ml",
         xtics=c(0.1,1,10,100,1000,10000),
         ytics=seq(0,35,by=5),
         ytics2=c(0:3)*250,
         ytics2units=1,
         breaks=40,
         main="Doheny Beach, F- Coliphage (EPA 1602)",
         CIRres1=CIRformat(getCIR(dofit10.fmc1602.pY.berm.1$fit)),
         CIRres0=CIRformat(getCIR(dofit10.fmc1602.pY.berm.0$fit)),
         Exp1=doheny$fmc1602[doheny$berm=="Open"],
         Exp0=doheny$fmc1602[doheny$berm=="Closed"],
         lab1="Berm open",lab0="Berm closed"
  )
  dev.off()
  
  pdf("~/dropbox/coliphage/results/figures/dose-response-do-fpc1601.pdf",height=6,width=7)
  plotPy.em(do.fpc1601.pY.berm.1,do.fpc1601.pY.berm.0,
         xlab="Concentration PFU / 100ml",
         ytics=seq(0,70,by=10),
         ytics2=c(0:5)*250,
         ytics2units=1,
         breaks=30,
         main="Doheny Beach, F+ Coliphage (EPA 1601)",
         CIRres1=CIRformat(getCIR(dofit10.fpc1601.pY.berm.1$fit)),
         CIRres0=CIRformat(getCIR(dofit10.fpc1601.pY.berm.0$fit)),
         Exp1=doheny$fpc1601[doheny$berm=="Open"],
         Exp0=doheny$fpc1601[doheny$berm=="Closed"],
         lab1="Berm open",lab0="Berm closed"
)
  dev.off()
  
  pdf("~/dropbox/coliphage/results/figures/dose-response-do-fpc1602.pdf",height=6,width=7)
  plotPy.em(do.fpc1602.pY.berm.1,do.fpc1602.pY.berm.0,
         xlab="Concentration PFU / 100ml",
         xtics=c(0.1,1,10,100),
         ytics=seq(0,90,by=10),
         ytics2=c(0:7)*500,
         ytics2units=1,
         breaks=50,
         main="Doheny Beach, F+ Coliphage (EPA 1602)",
         CIRres1=CIRformat(getCIR(dofit10.fpc1602.pY.berm.1$fit)),
         CIRres0=CIRformat(getCIR(dofit10.fpc1602.pY.berm.0$fit)),
         Exp1=doheny$fpc1602[doheny$berm=="Open"],
         Exp0=doheny$fpc1602[doheny$berm=="Closed"],
         lab1="Berm open",lab0="Berm closed"
  )
  dev.off()
  
  pdf("~/dropbox/coliphage/results/figures/dose-response-fa-fpc1601.pdf",height=7,width=5)
  plotPy(fa.fpc1601.pY,
         xlab="Concentration PFU / 100ml",
         ytics=seq(0,70,by=10),
         ytics2=c(0:3)*25,
         ytics2units=1,
         breaks=40,
         main="Fairhope Beach, F+ Coliphage (EPA 1601)",
         CIRres=CIRformat(getCIR(fafit10.fpc1601.pY$fit)),
         Exp=fairhope$fpc1601
  )
  dev.off()
  
  pdf("~/dropbox/coliphage/results/figures/dose-response-go-fpc1601.pdf",height=7,width=5)
  plotPy(go.fpc1601.pY,
         xlab="Concentration PFU / 100ml",
         ytics=seq(0,70,by=10),
         ytics2=c(0:4)*50,
         ytics2units=1,
         breaks=40,
         main="Goddard Beach, F+ Coliphage (EPA 1601)",
         CIRres=CIRformat(getCIR(gofit10.fpc1601.pY$fit)),
         Exp=goddard$fpc1601
  )
  dev.off()
  
  pdf("~/dropbox/coliphage/results/figures/dose-response-ma-fmc1602.pdf",height=6,width=7)
  plotPy.em(ma.fmc1602.pY.berm.1,ma.fmc1602.pY.berm.0,
         xlab="Concentration PFU / 100ml",
         ytics=seq(0,70,by=10),
         ytics2=c(0:4)*500,
         ytics2units=1,
         breaks=35,
         main="Malibu Beach, F- Coliphage (EPA 1602)",
         CIRres1=CIRformat(getCIR(mafit10.fmc1602.pY.berm.1$fit)),
         CIRres0=CIRformat(getCIR(mafit10.fmc1602.pY.berm.0$fit)),
         Exp1=malibu$fmc1602[malibu$berm=="Open"],
         Exp0=malibu$fmc1602[malibu$berm=="Closed"],
         lab1="Berm open",lab0="Berm closed"
  )
  dev.off()
    
  pdf("~/dropbox/coliphage/results/figures/dose-response-ma-fpc1601.pdf",height=6,width=7)
  plotPy.em(ma.fpc1601.pY.berm.1,ma.fpc1601.pY.berm.0,
         xlab="Concentration PFU / 100ml",
         ytics=seq(0,70,by=10),
         ytics2=c(0:6)*250,
         ytics2units=1,
         breaks=35,
         main="Malibu Beach, F+ Coliphage (EPA 1601)",
         CIRres1=CIRformat(getCIR(mafit10.fpc1601.pY.berm.1$fit)),
         CIRres0=CIRformat(getCIR(mafit10.fpc1601.pY.berm.0$fit)),
         Exp1=malibu$fpc1601[malibu$berm=="Open"],
         Exp0=malibu$fpc1601[malibu$berm=="Closed"],
         lab1="Berm open",lab0="Berm closed"
  )
  dev.off()

  pdf("~/dropbox/coliphage/results/figures/dose-response-mb-fmc1601.pdf",height=7,width=5)
  plotPy(mb.fmc1601.pY,
         xlab="Concentration PFU / 100ml",
         ytics=seq(0,15,by=1),
         ytics2=c(0:4)*500,
         ytics2units=1,
         breaks=30,
         main="Mission Bay Beach, F- Coliphage (EPA 1601)",
         CIRres=CIRformat(getCIR(mbfit10.fmc1601.pY$fit)),
         Exp=mission$fmc1601
  )
  dev.off()
  
  pdf("~/dropbox/coliphage/results/figures/dose-response-mb-fpc1601.pdf",height=7,width=5)
  plotPy(mb.fpc1601.pY,
         xlab="Concentration PFU / 100ml",
         ytics=seq(0,70,by=10),
         ytics2=c(0:4)*1000,
         ytics2units=1,
         breaks=15,
         main="Mission Bay Beach, F+ Coliphage (EPA 1601)",
         CIRres=CIRformat(getCIR(mbfit10.fpc1601.pY$fit)),
         Exp=mission$fpc1601
  )
  dev.off()

