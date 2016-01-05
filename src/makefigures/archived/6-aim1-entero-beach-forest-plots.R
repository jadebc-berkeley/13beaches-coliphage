# --------------------------------------
# 6-aim1-entero-beach-forest-plots
# ben arnold (benarnold@berkeley.edu)
#
# description:
# plot CIRs associated with Enterococcus exposure
# using summary regression output
# stratified by beach and by point versus non-point source
#
# version 1 (1 jun 2015)
#
# --------------------------------------

# --------------------------------------
# input files:
#	aim1-entero1600-Quartile-regs-body.RData
#	aim1-enteroQPCR-Quartile-regs-body.RData
#
#	aim1-entero1600-35cfu-regs-body.RData
#   aim1-enteroQPCR-470cce-regs-body.RData
#
# output files:
#	 aim1-entero1600-Quartile-beach-forest.pdf
#	 aim1-enteroQPCR-Quartile-beach-forest.pdf
#	 aim1-entero1600-QPCR-regulatory-beach-forest.pdf
# --------------------------------------

# --------------------------------------
# preamble
# --------------------------------------

rm(list=ls())
library(RColorBrewer)



# --------------------------------------
# EPA 1600 quartile CIR forest plot
# --------------------------------------

load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-entero1600-Quartile-regs-body.Rdata")

# ------------------------------------
# Calculate CIRs from the
# beach output
# ------------------------------------

# function to get Estimates and SEs from a linear combination of regression coefficients
lccalc <- function(lc,x,vcv) {
	# lc : linear combination of coefficients
	# x : log-linear model object returned from coeftest (class=coeftest)
	# vcv : variance-covariance matrix of coefficients for robust SEs
	est <- exp(t(lc)%*%x[,1])
	se  <- sqrt( t(lc)%*%vcv%*%lc )
	lb <- exp(log(est)-1.96*se)
	ub <- exp(log(est)+1.96*se)
	list(CIR=est,CIRlb=lb,CIRub=ub)
}

# function to bind together CIRs for Q2 - Q4
cirbind <- function(fit) {
	# fit : results returned from mpreg with fit and vcovCL objects
	lc2 <- c(0,1,rep(0,nrow(fit$fit)-2))
	lc3 <- c(0,0,1,rep(0,nrow(fit$fit)-3))
	lc4 <- c(0,0,0,1,rep(0,nrow(fit$fit)-4))
	res <- rbind(
	lccalc(lc2,fit$fit,fit$vcovCL),
	lccalc(lc3,fit$fit,fit$vcovCL),
	lccalc(lc4,fit$fit,fit$vcovCL))
	rownames(res) <- paste("Entero Q",2:4,sep="")
	return(res)
}

# grab beach-specific CIRs and 95% CIs
fitlist <- list(hufit,sifit,wpfit,wefit,avfit,bofit,dhfit,edfit,fafit,gdfit,mafit,mbfit,sufit)
beach.cirs <- lapply(fitlist,cirbind)

# collate the beach CIRs into 3 matrices, corresponding to Q2 - Q4
q2.cirs <- q3.cirs <- q4.cirs  <- matrix(NA,nrow=length(beach.cirs),ncol=3)
for (i in 1:length(beach.cirs)) {
	q2.cirs[i,] <- unlist(beach.cirs[[i]][1,])
	q3.cirs[i,] <- unlist(beach.cirs[[i]][2,])
	q4.cirs[i,] <- unlist(beach.cirs[[i]][3,])
}


# --------------------------------------
# Forest Plot
# Beach-specific CIR estimates 
# --------------------------------------

# Add beach labels
# CAREFUL: this label corresponds to the order of the blist.body/head/swall lists earlier in the code
fCIRlab <- c("Huntington","Silver","Washington Park","West","Avalon","Boqueron","Edgewater","Doheny","Fairhope","Goddard","Malibu","Mission Bay","Surfside")


# Order beaches by point source and non-point source, then by the body immersion CIR
nonpoint <- c(0,0,0,0,0,0,0,1,0,0,1,1,1)
bord <- rev(order(nonpoint,q2.cirs[,1]))
q2CIRs <- q2.cirs[bord,]
q3CIRs <- q3.cirs[bord,]
q4CIRs <- q4.cirs[bord,]
fCIRlab <- fCIRlab[bord]

# Add in the combined, summary estimates (subgroups + overall)
q2CIRs <- rbind(cir.all[1,],cir.nps[1,], q2CIRs[1:4,],cir.ps[1,], q2CIRs[5:13,])
q3CIRs <- rbind(cir.all[2,],cir.nps[2,], q3CIRs[1:4,],cir.ps[2,], q3CIRs[5:13,])
q4CIRs <- rbind(cir.all[3,],cir.nps[3,], q4CIRs[1:4,],cir.ps[3,], q4CIRs[5:13,])
fCIRlab <- c("All Combined","No Point Source Combined",fCIRlab[1:4],"Point Source Combined",fCIRlab[5:13])



pdf("~/dropbox/13beaches/aim1-results/figs/aim1-entero1600-Quartile-beach-forest.pdf",height=6,width=11)

lo <- layout(mat=matrix(1:4,nrow=1,ncol=4),widths=c(1.05,1,1,1))
op <- par(mar=c(5,2,2,2)+0.1)
fcol <- brewer.pal(9,"YlGnBu")[5]
mcol <- brewer.pal(9,"YlGnBu")[6]
scol <- brewer.pal(11,"Spectral")[11]
cols <- c(scol,rep(mcol,5),rep(fcol,10))
pchs <- c(5,5,rep(15,4),5,rep(15,9))

xtics <- c(0.25,0.5,1,2,4,8)

# Beach Labels
bhs <- barplot(1:16,horiz=TRUE,xlim=c(0,1),xaxt="n",col=NA,border=NA)
text(x=rep(1,length(bhs)),y=bhs,fCIRlab,cex=1.5,adj=1,col=c(scol,mcol,rep("gray40",4),fcol,rep("gray40",9)),font=c(2,2,rep(1,4),2,rep(1,9)) )

# Q2 vs Q1
bhs <- barplot(1:16,horiz=TRUE,xlim=range(xtics),log="x",xaxt="n",col=NA,border=NA)
axis(side=1,at=xtics,las=1)
segments(x0=1,x1=1,y0=0,y1=max(bhs+1),lty=2)
mtext("Adjusted CIR",side=1,line=3,las=1)
mtext("Q2 vs Q1",side=3,line=0,at=1,col="gray40",font=2)
# scale the area of the point size by 1/SE of the estimates
# calc 1/SE
oneoverSE <- 1.96/(log(q2CIRs[,3])-log(q2CIRs[,1]))
# now scale to area of a square rather than an edge
# and further scale the cex factor so the smallest = 1
cexscale <- sqrt(oneoverSE)
cexscale <- cexscale/min(cexscale)
# label the combined estimates
segments(x0=q2CIRs[1,1],x1=q2CIRs[1,1],y0=0,y1=max(bhs)+1,lty=2,col=scol)
text(x=q2CIRs[c(1,2,7),1]+0.3,y=bhs[c(1,2,7)]+0.3,paste(sprintf("%1.2f",q2CIRs[c(1,2,7),1])," (",sprintf("%1.2f",q2CIRs[c(1,2,7),2]),", ", sprintf("%1.2f",q2CIRs[c(1,2,7),3]) ,")",sep=""),pos=4,cex=0.8,col=c(scol,mcol,fcol) )
# plot the estimates
segments(y0=bhs,x0=q2CIRs[,2],x1=q2CIRs[,3],lwd=2,col=cols)
points(q2CIRs[,1],bhs,cex=1.5*cexscale,pch=pchs,lwd=2,bg="white",col=cols)


# Q3 vs Q1
bhs <- barplot(1:16,horiz=TRUE,xlim=range(xtics),log="x",xaxt="n",col=NA,border=NA)
axis(side=1,at=xtics,las=1)
segments(x0=1,x1=1,y0=0,y1=max(bhs+1),lty=2)
mtext("Adjusted CIR",side=1,line=3,las=1)
# mtext(fCIRlab,side=2,line=0,at=bhs,cex=0.9,las=1)
mtext("Q3 vs Q1",side=3,line=0,at=1,col="gray40",font=2)
# scale the area of the point size by 1/SE of the estimates
# calc 1/SE
oneoverSE <- 1.96/(log(q3CIRs[,3])-log(q3CIRs[,1]))
# now scale to area of a square rather than an edge
# and further scale the cex factor so the smallest = 1
cexscale <- sqrt(oneoverSE)
cexscale <- cexscale/min(cexscale)
# label the combined estimates
segments(x0=q3CIRs[1,1],x1=q3CIRs[1,1],y0=0,y1=max(bhs)+1,lty=2,col=scol)
text(x=q3CIRs[c(1,2,7),1]+0.3,y=bhs[c(1,2,7)]+0.3,paste(sprintf("%1.2f",q3CIRs[c(1,2,7),1])," (",sprintf("%1.2f",q3CIRs[c(1,2,7),2]),", ", sprintf("%1.2f",q3CIRs[c(1,2,7),3]) ,")",sep=""),pos=4,cex=0.8,col=c(scol,mcol,fcol) )
# plot the estimates
segments(y0=bhs,x0=q3CIRs[,2],x1=q3CIRs[,3],lwd=2,col=cols)
points(q3CIRs[,1],bhs,cex=1.5*cexscale,pch=pchs,lwd=2,bg="white",col=cols)


# Q4 versus Q1
bhs <- barplot(1:16,horiz=TRUE,xlim=range(xtics),log="x",xaxt="n",col=NA,border=NA)
axis(side=1,at=xtics,las=1)
segments(x0=1,x1=1,y0=0,y1=max(bhs+1),lty=2)
mtext("Adjusted CIR",side=1,line=3,las=1)
# mtext(fCIRlab,side=2,line=0,at=bhs,cex=0.9,las=1)
mtext("Q4 vs Q1",side=3,line=0,at=1,col="gray40",font=2)
# scale the area of the point size by 1/SE of the estimates
# calc 1/SE
oneoverSE <- 1.96/(log(q4CIRs[,3])-log(q4CIRs[,1]))
# now scale to area of a square rather than an edge
# and further scale the cex factor so the smallest = 1
cexscale <- sqrt(oneoverSE)
cexscale <- cexscale/min(cexscale)
# label the combined estimates
segments(x0=q4CIRs[1,1],x1=q4CIRs[1,1],y0=0,y1=max(bhs)+1,lty=2,col=scol)
text(x=q4CIRs[c(1,2,7),1]+0.3,y=bhs[c(1,2,7)]+0.3,paste(sprintf("%1.2f",q4CIRs[c(1,2,7),1])," (",sprintf("%1.2f",q4CIRs[c(1,2,7),2]),", ", sprintf("%1.2f",q4CIRs[c(1,2,7),3]) ,")",sep=""),pos=4,cex=0.8,col=c(scol,mcol,fcol) )
# plot the estimates
segments(y0=bhs,x0=q4CIRs[,2],x1=q4CIRs[,3],lwd=2,col=cols)
points(q4CIRs[,1],bhs,cex=1.5*cexscale,pch=pchs,lwd=2,bg="white",col=cols)

par(op)
dev.off()


# --------------------------------------
# EPA QPCR quartile CIR forest plot
# --------------------------------------

load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-enteroQPCR-Quartile-regs-body.Rdata")

# ------------------------------------
# Calculate CIRs from the
# beach output
# ------------------------------------


# grab beach-specific CIRs and 95% CIs
fitlist <- list(hufit,sifit,wpfit,wefit,avfit,bofit,dhfit,edfit,fafit,gdfit,mafit,mbfit,sufit)
beach.cirs <- lapply(fitlist,cirbind)

# collate the beach CIRs into 3 matrices, corresponding to Q2 - Q4
q2.cirs <- q3.cirs <- q4.cirs  <- matrix(NA,nrow=length(beach.cirs),ncol=3)
for (i in 1:length(beach.cirs)) {
	q2.cirs[i,] <- unlist(beach.cirs[[i]][1,])
	q3.cirs[i,] <- unlist(beach.cirs[[i]][2,])
	q4.cirs[i,] <- unlist(beach.cirs[[i]][3,])
}


# --------------------------------------
# Forest Plot
# Beach-specific CIR estimates 
# --------------------------------------

# Add beach labels
# CAREFUL: this label corresponds to the order of the blist.body/head/swall lists earlier in the code
fCIRlab <- c("Huntington","Silver","Washington Park","West","Avalon","Boqueron","Edgewater","Doheny","Fairhope","Goddard","Malibu","Mission Bay","Surfside")


# Order beaches by point source and non-point source, then by the body immersion CIR
nonpoint <- c(0,0,0,0,0,0,0,1,0,0,1,1,1)
bord <- rev(order(nonpoint,q2.cirs[,1]))
q2CIRs <- q2.cirs[bord,]
q3CIRs <- q3.cirs[bord,]
q4CIRs <- q4.cirs[bord,]
fCIRlab <- fCIRlab[bord]

# Add in the combined, summary estimates (subgroups + overall)
q2CIRs <- rbind(cir.all[1,],cir.nps[1,], q2CIRs[1:4,],cir.ps[1,], q2CIRs[5:13,])
q3CIRs <- rbind(cir.all[2,],cir.nps[2,], q3CIRs[1:4,],cir.ps[2,], q3CIRs[5:13,])
q4CIRs <- rbind(cir.all[3,],cir.nps[3,], q4CIRs[1:4,],cir.ps[3,], q4CIRs[5:13,])
fCIRlab <- c("All Combined","No Point Source Combined",fCIRlab[1:4],"Point Source Combined",fCIRlab[5:13])



pdf("~/dropbox/13beaches/aim1-results/figs/aim1-enteroQPCR-Quartile-beach-forest.pdf",height=6,width=11)

lo <- layout(mat=matrix(1:4,nrow=1,ncol=4),widths=c(1.05,1,1,1))
op <- par(mar=c(5,2,2,2)+0.1)



fcol <- brewer.pal(9,"YlOrBr")[5]
mcol <- brewer.pal(9,"YlOrBr")[7]
scol <- brewer.pal(9,"YlOrBr")[9]
cols <- c(scol,rep(mcol,5),rep(fcol,10))
pchs <- c(5,5,rep(15,4),5,rep(15,9))

xtics <- c(0.25,0.5,1,2,4,8)

# Beach Labels
bhs <- barplot(1:16,horiz=TRUE,xlim=c(0,1),xaxt="n",col=NA,border=NA)
text(x=rep(1,length(bhs)),y=bhs,fCIRlab,cex=1.5,adj=1,col=c(scol,mcol,rep("gray40",4),fcol,rep("gray40",9)),font=c(2,2,rep(1,4),2,rep(1,9)) )

# Q2 vs Q1
bhs <- barplot(1:16,horiz=TRUE,xlim=range(xtics),log="x",xaxt="n",col=NA,border=NA)
axis(side=1,at=xtics,las=1)
segments(x0=1,x1=1,y0=0,y1=max(bhs+1),lty=2)
mtext("Adjusted CIR",side=1,line=3,las=1)
mtext("Q2 vs Q1",side=3,line=0,at=1,col="gray40",font=2)
# scale the area of the point size by 1/SE of the estimates
# calc 1/SE
oneoverSE <- 1.96/(log(q2CIRs[,3])-log(q2CIRs[,1]))
# now scale to area of a square rather than an edge
# and further scale the cex factor so the smallest = 1
cexscale <- sqrt(oneoverSE)
cexscale <- cexscale/min(cexscale)
# label the combined estimates
segments(x0=q2CIRs[1,1],x1=q2CIRs[1,1],y0=0,y1=max(bhs)+1,lty=2,col=scol)
text(x=q2CIRs[c(1,2,7),1]+0.3,y=bhs[c(1,2,7)]+0.3,paste(sprintf("%1.2f",q2CIRs[c(1,2,7),1])," (",sprintf("%1.2f",q2CIRs[c(1,2,7),2]),", ", sprintf("%1.2f",q2CIRs[c(1,2,7),3]) ,")",sep=""),pos=4,cex=0.8,col=c(scol,mcol,fcol) )
# plot the estimates
segments(y0=bhs,x0=q2CIRs[,2],x1=q2CIRs[,3],lwd=2,col=cols)
points(q2CIRs[,1],bhs,cex=1.5*cexscale,pch=pchs,lwd=2,bg="white",col=cols)


# Q3 vs Q1
bhs <- barplot(1:16,horiz=TRUE,xlim=range(xtics),log="x",xaxt="n",col=NA,border=NA)
axis(side=1,at=xtics,las=1)
segments(x0=1,x1=1,y0=0,y1=max(bhs+1),lty=2)
mtext("Adjusted CIR",side=1,line=3,las=1)
# mtext(fCIRlab,side=2,line=0,at=bhs,cex=0.9,las=1)
mtext("Q3 vs Q1",side=3,line=0,at=1,col="gray40",font=2)
# scale the area of the point size by 1/SE of the estimates
# calc 1/SE
oneoverSE <- 1.96/(log(q3CIRs[,3])-log(q3CIRs[,1]))
# now scale to area of a square rather than an edge
# and further scale the cex factor so the smallest = 1
cexscale <- sqrt(oneoverSE)
cexscale <- cexscale/min(cexscale)
# label the combined estimates
segments(x0=q3CIRs[1,1],x1=q3CIRs[1,1],y0=0,y1=max(bhs)+1,lty=2,col=scol)
text(x=q3CIRs[c(1,2,7),1]+0.3,y=bhs[c(1,2,7)]+0.3,paste(sprintf("%1.2f",q3CIRs[c(1,2,7),1])," (",sprintf("%1.2f",q3CIRs[c(1,2,7),2]),", ", sprintf("%1.2f",q3CIRs[c(1,2,7),3]) ,")",sep=""),pos=4,cex=0.8,col=c(scol,mcol,fcol) )
# plot the estimates
segments(y0=bhs,x0=q3CIRs[,2],x1=q3CIRs[,3],lwd=2,col=cols)
points(q3CIRs[,1],bhs,cex=1.5*cexscale,pch=pchs,lwd=2,bg="white",col=cols)


# Q4 versus Q1
bhs <- barplot(1:16,horiz=TRUE,xlim=range(xtics),log="x",xaxt="n",col=NA,border=NA)
axis(side=1,at=xtics,las=1)
segments(x0=1,x1=1,y0=0,y1=max(bhs+1),lty=2)
mtext("Adjusted CIR",side=1,line=3,las=1)
# mtext(fCIRlab,side=2,line=0,at=bhs,cex=0.9,las=1)
mtext("Q4 vs Q1",side=3,line=0,at=1,col="gray40",font=2)
# scale the area of the point size by 1/SE of the estimates
# calc 1/SE
oneoverSE <- 1.96/(log(q4CIRs[,3])-log(q4CIRs[,1]))
# now scale to area of a square rather than an edge
# and further scale the cex factor so the smallest = 1
cexscale <- sqrt(oneoverSE)
cexscale <- cexscale/min(cexscale)
# label the combined estimates
segments(x0=q4CIRs[1,1],x1=q4CIRs[1,1],y0=0,y1=max(bhs)+1,lty=2,col=scol)
text(x=q4CIRs[c(1,2,7),1]+0.3,y=bhs[c(1,2,7)]+0.3,paste(sprintf("%1.2f",q4CIRs[c(1,2,7),1])," (",sprintf("%1.2f",q4CIRs[c(1,2,7),2]),", ", sprintf("%1.2f",q4CIRs[c(1,2,7),3]) ,")",sep=""),pos=4,cex=0.8,col=c(scol,mcol,fcol) )
# plot the estimates
segments(y0=bhs,x0=q4CIRs[,2],x1=q4CIRs[,3],lwd=2,col=cols)
points(q4CIRs[,1],bhs,cex=1.5*cexscale,pch=pchs,lwd=2,bg="white",col=cols)

par(op)
dev.off()


# --------------------------------------
# EPA regulatory guidelines forest plot
# --------------------------------------

# grab results for EPA 1600
load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-entero1600-35cfu-regs-body.Rdata")

# function to bind together CIRs for binary exposure (above/below regulator limit)
getcir <- function(fit) {
	# fit : results returned from mpreg with fit and vcovCL objects
	lc <- c(0,1,rep(0,nrow(fit$fit)-2))
	res <- lccalc(lc,fit$fit,fit$vcovCL)
	return(res)
}

# grab beach-specific CIRs and 95% CIs
fitlist <- list(hufit,sifit,wpfit,wefit,avfit,bofit,dhfit,edfit,fafit,gdfit,mafit,mbfit,sufit)
cirs.1600 <- t(sapply(fitlist,getcir))
pool.1600 <- cir.all


# grab results for EPA 1611 QPCR
load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-enteroQPCR-470cce-regs-body.Rdata")

# grab beach-specific CIRs and 95% CIs
fitlist <- list(hufit,sifit,wpfit,wefit,avfit,bofit,dhfit,edfit,fafit,gdfit,mafit,mbfit,sufit)
cirs.QPCR <- t(sapply(fitlist,getcir))
pool.QPCR <- cir.all


# Add beach labels
# CAREFUL: this label corresponds to the order of the blist.body/head/swall lists earlier in the code
fCIRlab <- c("Huntington","Silver","Washington Park","West","Avalon","Boqueron","Edgewater","Doheny","Fairhope","Goddard","Malibu","Mission Bay","Surfside")


# Order beaches by point source and non-point source, then by the EPA 1600 CIRs
nonpoint <- c(0,0,0,0,0,0,0,1,0,0,1,1,1)
bord <- rev(order(nonpoint,unlist(cirs.1600[,1])))
CIRs1600 <- cirs.1600[bord,]
CIRsQPCR <- cirs.QPCR[bord,]
fCIRlab <- fCIRlab[bord]

# Add in the combined, summary estimates (subgroups + overall)
CIRs1600 <- rbind(pool.1600[3,],pool.1600[1,], CIRs1600[1:4,],pool.1600[2,], CIRs1600[5:13,])
CIRsQPCR <- rbind(pool.QPCR[3,],pool.QPCR[1,], CIRsQPCR[1:4,],pool.QPCR[2,], CIRsQPCR[5:13,])
fCIRlab <- c("All Combined","No Point Source Combined",fCIRlab[1:4],"Point Source Combined",fCIRlab[5:13])

# set completely unstable estimates close to zero to NA to allow for log scaling
CIRs1600[CIRs1600<0.01] <- NA
CIRsQPCR[CIRsQPCR<0.01] <- NA

# ensure the results matrices are numeric for calculations
CIRs1600 <- matrix(as.numeric(CIRs1600),nrow=16,ncol=3)
CIRsQPCR <- matrix(as.numeric(CIRsQPCR),nrow=16,ncol=3)


pdf("~/dropbox/13beaches/aim1-results/figs/aim1-entero1600-QPCR-regulatory-beach-forest.pdf",height=6,width=11)

lo <- layout(mat=matrix(1:3,nrow=1,ncol=3),widths=c(0.8,1,1))
op <- par(mar=c(5,2,3,2)+0.1)

fcol <- "black"
mcol <- "black"
scol <- "black"

pchs <- c(5,5,rep(15,4),5,rep(15,9))

xtics <- c(0.25,0.5,1,2,4,8)

# Beach Labels
bhs <- barplot(1:16,horiz=TRUE,xlim=c(0,1),xaxt="n",col=NA,border=NA)
text(x=rep(1,length(bhs)),y=bhs,fCIRlab,cex=1.5,adj=1,col=c(scol,mcol,rep("gray40",4),fcol,rep("gray40",9)),font=c(2,2,rep(1,4),2,rep(1,9)) )

# Entero 1600 > 35 CFU / 100ml
fcol <- brewer.pal(9,"YlGnBu")[5]
mcol <- brewer.pal(9,"YlGnBu")[6]
scol <- brewer.pal(11,"Spectral")[11]
cols <- c(scol,rep(mcol,5),rep(fcol,10))

bhs <- barplot(1:16,horiz=TRUE,xlim=range(xtics),log="x",xaxt="n",col=NA,border=NA)
axis(side=1,at=xtics,las=1)
segments(x0=1,x1=1,y0=0,y1=max(bhs+1),lty=2)
mtext("Adjusted CIR",side=1,line=3,las=1)
mtext("EPA 1600\n>35 CFU/100ml",side=3,line=0,at=1,col="gray40",font=2)
# scale the area of the point size by 1/SE of the estimates
# calc 1/SE
oneoverSE <- 1.96/(log(CIRs1600[,3])-log(CIRs1600[,1]))
# now scale to area of a square rather than an edge
# and further scale the cex factor so the smallest = 1
cexscale <- sqrt(oneoverSE)
cexscale <- cexscale/min(cexscale,na.rm=T)
# label the combined estimates
segments(x0=CIRs1600[1,1],x1=CIRs1600[1,1],y0=0,y1=max(bhs)+1,lty=2,col=scol)
text(x=CIRs1600[c(1,2,7),1]+0.3,y=bhs[c(1,2,7)]+0.3,paste(sprintf("%1.2f",CIRs1600[c(1,2,7),1])," (",sprintf("%1.2f",CIRs1600[c(1,2,7),2]),", ", sprintf("%1.2f",CIRs1600[c(1,2,7),3]) ,")",sep=""),pos=4,cex=0.8,col=c(scol,mcol,fcol) )
# plot the estimates
segments(y0=bhs,x0=CIRs1600[,2],x1=CIRs1600[,3],lwd=2,col=cols)
points(CIRs1600[,1],bhs,cex=1.5*cexscale,pch=pchs,lwd=2,bg="white",col=cols)


# Entero QPCR > 470 CCE / 100ml
fcol <- brewer.pal(9,"YlOrBr")[5]
mcol <- brewer.pal(9,"YlOrBr")[7]
scol <- brewer.pal(9,"YlOrBr")[9]
cols <- c(scol,rep(mcol,5),rep(fcol,10))

bhs <- barplot(1:16,horiz=TRUE,xlim=range(xtics),log="x",xaxt="n",col=NA,border=NA)
axis(side=1,at=xtics,las=1)
segments(x0=1,x1=1,y0=0,y1=max(bhs+1),lty=2)
mtext("Adjusted CIR",side=1,line=3,las=1)
# mtext(fCIRlab,side=2,line=0,at=bhs,cex=0.9,las=1)
mtext("EPA 1611 qPCR\n>470 CCE/100ml",side=3,line=0,at=1,col="gray40",font=2)
# scale the area of the point size by 1/SE of the estimates
# calc 1/SE
oneoverSE <- 1.96/(log(CIRsQPCR[,3])-log(CIRsQPCR[,1]))
# now scale to area of a square rather than an edge
# and further scale the cex factor so the smallest = 1
cexscale <- sqrt(oneoverSE)
cexscale <- cexscale/min(cexscale,na.rm=T)
# label the combined estimates
segments(x0=CIRsQPCR[1,1],x1=CIRsQPCR[1,1],y0=0,y1=max(bhs)+1,lty=2,col=scol)
text(x=CIRsQPCR[c(1,2,7),1]+0.3,y=bhs[c(1,2,7)]+0.3,paste(sprintf("%1.2f",CIRsQPCR[c(1,2,7),1])," (",sprintf("%1.2f",CIRsQPCR[c(1,2,7),2]),", ", sprintf("%1.2f",CIRsQPCR[c(1,2,7),3]) ,")",sep=""),pos=4,cex=0.8,col=c(scol,mcol,fcol) )
# plot the estimates
segments(y0=bhs,x0=CIRsQPCR[,2],x1=CIRsQPCR[,3],lwd=2,col=cols)
points(CIRsQPCR[,1],bhs,cex=1.5*cexscale,pch=pchs,lwd=2,bg="white",col=cols)


par(op)
dev.off()





