# --------------------------------------
# aim1-12-entero1600-figs.R
# ben arnold (benarnold@berkeley.edu)
#
# description:
# plot CIRs associated with Enterococcus
# EPA 1600 output
# using summary regression output
#
#
#
# version 1 (20 feb 2015)
#
# --------------------------------------


# --------------------------------------
# preamble
# --------------------------------------

rm(list=ls())
library(RColorBrewer)


# --------------------------------------
# Load the water quality dataset to
# get the mid-points of the Entero
# concentrations by each quintile
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
# Load the regression output
# --------------------------------------

load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-entero1600-Quartile-regs-body.Rdata")


# --------------------------------------
# CIR plot of summary estimates
# --------------------------------------

pdf("~/dropbox/13beaches/aim1-results/figs/aim1-entero1600-Quartile-CIRs-agestratified.pdf",width=10,height=3)
op <- par(mar=c(6,7,3,0)+0.1)
cols <- brewer.pal(9,"YlGnBu")[5:8]
ytics <- c(0.5,1,2,4)
# set up an empty plot
MidPts <- barplot(1:4,names.arg=NA,border=NA,col=NA,
	log="y",ylog=T,ylim=c(0.5,4.5),ylab="",yaxt="n",
	las=1,bty="n"
	)
	axis(2,at=ytics,las=1)
	segments(x0=0,x1=max(MidPts+0.5),y0=1,lty=2,lwd=1.5)
	segments(x0=mean(MidPts[1:2]),y0=0.2,y1=4.5,lwd=2,col="gray80")
	# mtext(title,side=2,line=6,at=4,font=2,adj=0.5,las=1,cex=1.5)
	mtext("aCIR",side=3,line=0.5,at=-0.2)
	mtext(c("All\nAges","Ages\n0 to 4","Ages\n5 to 10","Ages\n>10"),at=MidPts,side=3,line=0  )
	
	# calculate X coordinates relative to the mid points for each group
	xspan <- 0.37
	xplus <- c(-xspan, -xspan/3, xspan/3, xspan)
	
	
	# plot all age estimates
	xall <- xplus+MidPts[1]
	segments(x0=xall,y0=c(NA,cir.all[,"CIRlb"]),y1=c(NA,cir.all[,"CIRub"]),lwd=2,col=cols)
	points(xall,c(1,cir.all[,"CIR"]),pch=16,bg="white",cex=1.75,lwd=2,col=cols)
	
	# plot age 0 to 4 estimates
	x0to4 <- xplus+MidPts[2]
	segments(x0= x0to4,y0=c(NA,cir.age0to4[,"CIRlb"]),y1=c(NA,cir.age0to4[,"CIRub"]),lwd=2,col=cols)
	points(x0to4,c(1,cir.age0to4[,"CIR"]),pch=16,bg="white",cex=1.75,lwd=2,col=cols)
	
	# plot age 5 to 10 estimates
	x5to10 <- xplus+MidPts[3]
	segments(x0=x5to10,y0=c(NA,cir.age5to10[,"CIRlb"]),y1=c(NA,cir.age5to10[,"CIRub"]),lwd=2,col=cols)
	points(x5to10,c(1,cir.age5to10[,"CIR"]),pch=16,bg="white",cex=1.75,lwd=2,col=cols)
	
	# plot age > 10 estimates
	x11plus <- xplus+MidPts[4]
	segments(x0=x11plus,y0=c(NA,cir.age11plus[,"CIRlb"]),y1=c(NA,cir.age11plus[,"CIRub"]),lwd=2,col=cols)
	points(x11plus,c(1,cir.age11plus[,"CIR"]),pch=16,bg="white",cex=1.75,lwd=2,col=cols)
	
	# print labels and Ns
	allxs <- c(xall,x0to4,x5to10,x11plus)
	labx <- MidPts[1]-xspan*1.6
	mtext(expression(paste(italic("Enterococcus")," Quartile")),side=1,line=0.5,at=labx,adj=1,col=cols[3],cex=0.8)
	mtext(c("Q1","Q2","Q3","Q4"),side=1,line=0.5,at=allxs,col=cols,cex=0.9,font=2)
	
	mtext("Range (CFU/100ml)",side=1,line=1.5,at=labx,adj=1,col=cols[3],cex=0.8)
	mtext(rngQs,side=1,line=1.5,at=allxs,col=cols,cex=0.8)
	
	mtext("Incident Diarrhea Cases",side=1,line=3,at=labx,adj=1,cex=0.8,col="gray40")
	ns <- c(N.all[,1],N.age0to4[,1],N.age5to10[,1],N.age11plus[,1])
	mtext(  format(ns,big.mark=","),side=1,line=3,at=allxs+0.05,adj=1,cex=0.75    )
	
	mtext("Population At Risk",side=1,line=4,at=labx,adj=1,cex=0.8,col="gray40")
	Ns <- c(N.all[,2],N.age0to4[,2],N.age5to10[,2],N.age11plus[,2])
	mtext(  format(Ns,big.mark=","),side=1,line=4,at=allxs+0.05,adj=1,cex=0.75    )
	
	mtext("Incidence per 1000",side=1,line=5,at=labx,adj=1,cex=0.8,col="gray40")
	mtext(  sprintf("%1.0f",(ns/Ns)*1000),side=1,line=5,at=allxs+0.05,adj=1,cex=0.75    )
	

par(op)
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






