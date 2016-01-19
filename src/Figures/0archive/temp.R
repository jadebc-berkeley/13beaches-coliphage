# load("~/dropbox/coliphage/results/rawoutput/regress-10day-continuous-body-risk.Rdata")
# 
# all.fpc1602.pY.low.risk=all.fpc1602.pY.low
# all.fpc1602.pY.high.risk=all.fpc1602.pY.high
# load("~/dropbox/coliphage/results/rawoutput/regress-10day-continuous-body.Rdata")
# 


# --------------------------------------
# general plotting function for the
# dose-response curves with effect modification
# --------------------------------------
#plotPy.em <- function(pYcurve1,pYcurve0,xtics=c(0.1,1,10,100,1000),xlab,ytics,ytics2a,ytics2aunits,ytics2b,ytics2bunits,breaksa,breaksb,main,CIRres1,CIRres0,Exp1,Exp0,lab1,lab0){
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
  comb$rd=c(rdcurve1$bootest*100,rdcurve0$bootest*100)


  comb$em=as.factor(c(rep(1,length(pYcurve1$pX)),rep(0,length(pYcurve0$pX))))
#  par(oma=c(0,0,5,0))
  par(oma=c(0,2,5,5)) 
  lo <- layout(mat=matrix(c(1,3,2,4),nrow=2,ncol=2),heights=c(1,0.5))
  
  # Low-risk plot --------------
  plot(comb$pX[comb$em==1],comb$pY[comb$em==1],type="l",
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
  
  # double-axis for risk difference
  par(new = TRUE)
  plot(comb$pX[comb$em==1],comb$rd[comb$em==1],type="n",
       ylim=range(ytics),yaxt="n",ylab="",
       xlim=range(log10(xtics)),xaxt="n",xlab="",
       bty="n",
  )
  axis(4,at=ytics,las=1)
  mtext("Risk Difference", side=4, line=2, cex=0.75)
  
  # High-risk plot --------------
  plot(comb$pX[comb$em==0],comb$rd[comb$em==0],col=comb$em,type="n",
       ylim=range(ytics),yaxt="n",ylab="",
       xlim=range(log10(xtics)),xaxt="n",xlab="",
       bty="n",
  )
  axis(1,at=log10(xtics),labels=xtics,las=1)
  axis(2,at=ytics,las=1)
  mtext("Probability of Gastrointestinal Illness (%)", side=2, line=2, cex=0.8)
  segments(x0=log10(xtics),y0=rep(0,length(xtics)),y1=rep(max(ytics),length(xtics)),col="gray90")
  lines(comb$pX[comb$em==0],comb$pY[comb$em==0],lwd=1.2,col="#993300")
  lines(comb$pX[comb$em==0],comb$lb[comb$em==0],lwd=1.2,col="#993300",lty=5)
  lines(comb$pX[comb$em==0],comb$ub[comb$em==0],lwd=1.2,col="#993300",lty=5)
  mtext(paste(lab0),side=3,line=2,font=2)  
  mtext(paste("CIR for a log10 increase:",CIRres0),side=3,line=0)
  
  # double-axis for risk difference
  par(new = TRUE)  
  plot(comb$pX[comb$em==0],comb$rd[comb$em==0],col=comb$em,type="l",
       ylim=range(ytics),yaxt="n",ylab="",
       xlim=range(log10(xtics)),xaxt="n",xlab="",
       bty="n",
  )
  axis(4,at=ytics,las=1)
  mtext("Risk Difference", side=4, line=2, cex=0.75)
  
  
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
  
# }
# 
# plotPy.em(
#   pYcurve1=all.fpc1602.pY.low
#   pYcurve0=all.fpc1602.pY.high
#rdcurve1=all.fpc1602.pY.low.risk
#rdcurve0=all.fpc1602.pY.high.risk
#           xlab="Concentration PFU / 100ml"
#           xtics=c(0.1,1,10,100)
#ytics=seq(-10,110,10)
#           ytics2a=c(0:5)*50
#           ytics2aunits=1
#           ytics2b=c(0:5)*50
#           ytics2bunits=1
#           breaksa=30
#           breaksb=30
#           main="Male-Specific Coliphage (EPA 1602)"
#           CIRres1=CIRformat(getCIR(overall.fit10.fpc1602.low))
#           CIRres0=CIRformat(getCIR(overall.fit10.fpc1602.high))
#           Exp1=all$fpc1602[all$risk=="Low" & all$fpc1602>-1]
#           Exp0=all$fpc1602[all$risk=="High" & all$fpc1602>-1]
#           lab1="Low risk conditions"
#   lab0="High risk conditions"
# )