# --------------------------------------
# aim1-6-swim-exposure-figs.R
# ben arnold (benarnold@berkeley.edu)
#
# description:
# plot CIRs associated with swim exposure
# using summary regression output
#
# Caution:
# loaded input for each group (all ages,
# 0-10, 11+) has identical objects, so
# loading regression output will mask
# existing output in use unless it is
# renamed. The output has the identical
# object names to save programming time.
#
#
# version 1 (10 feb 2015)
#
# --------------------------------------


# --------------------------------------
# preamble
# --------------------------------------

rm(list=ls())
library(RColorBrewer)



# --------------------------------------
# ALL AGES
# load the regression output and
# saved workspace
# --------------------------------------

load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-swim-exposure-regs-allages.Rdata")



# --------------------------------------
# Order the beach results for plotting
# --------------------------------------
CIRlist <- list(hu.cir,si.cir,wp.cir,we.cir,av.cir,bo.cir,ed.cir,dh.cir,fa.cir,gd.cir,ma.cir,mb.cir,su.cir,all.cir)
CIRlab <- c("Huntington","Silver","Washington\nPark","West","Avalon","Boqueron","Edgewater","Doheny","Fairhope","Goddard","Malibu","Mission\nBay","Surfside","Combined")



# --------------------------------------
# ALL AGES
# plot CIRs in a large panel to show
# dose-response by swim exposure
# --------------------------------------


pdf("~/dropbox/13beaches/aim1-results/figs/swim-exposure-CIR-dose-response-panel.pdf",height=8*2,width=5*2)

lo <- layout(mat=matrix(1:14,nrow=7,ncol=2,byrow=TRUE))
ytics <- c(0.5,1,1.5,2,4,8)
cols <- brewer.pal(9,"YlOrRd")[5:7]

for(i in 1:length(CIRlist)) {
	
	if(i<=2) {
		op <- par(mar=c(1,12,4,2)+0.1)
	} else {
		op <- par(mar=c(1,12,2,2)+0.1)
	}
	
	bhs <- barplot(1:3,ylim=range(ytics),log="y",yaxt="n",col=NA,border=NA)
	axis(side=2,at=ytics,las=1)
	segments(x0=0,x1=max(bhs+1),y0=1,y1=1,lty=2)
	mtext(CIRlab[i],side=2,line=3,las=1,font=2)
	mtext("aCIR",side=2,line=3,las=1,at=max(ytics))
	# plot estimates
	segments(x0=bhs,y0=CIRlist[[i]][,2],y1=CIRlist[[i]][,3],lwd=2,col=cols)
	points(bhs,CIRlist[[i]][,1],pch=21,bg="white",cex=1.75,lwd=2,col=cols)
	# label estimates
	text(bhs,CIRlist[[i]][,1],sprintf("%1.2f",CIRlist[[i]][,1]),cex=0.9,pos=4)
	
	if(i<=2) {
		mtext(c("Body\nImmersion","Head\nImmersion","Swallowed\nWater"),side=3,line=0,at=bhs,col=cols)
	}
	
}
par(op)
dev.off()

# --------------------------------------
# ALL AGES
# Forest plot of CIRs
# --------------------------------------

fCIRlab <- c("Huntington","Silver","Washington Park","West","Avalon","Boqueron","Edgewater","Doheny","Fairhope","Goddard","Malibu","Mission Bay","Surfside")

bCIRs <- rbind(hu.cir[1,],si.cir[1,],wp.cir[1,],we.cir[1,],av.cir[1,],bo.cir[1,],ed.cir[1,],dh.cir[1,],fa.cir[1,],gd.cir[1,],ma.cir[1,],mb.cir[1,],su.cir[1,])

hCIRs <- rbind(hu.cir[2,],si.cir[2,],wp.cir[2,],we.cir[2,],av.cir[2,],bo.cir[2,],ed.cir[2,],dh.cir[2,],fa.cir[2,],gd.cir[2,],ma.cir[2,],mb.cir[2,],su.cir[2,])

sCIRs <- rbind(hu.cir[3,],si.cir[3,],wp.cir[3,],we.cir[3,],av.cir[3,],bo.cir[3,],ed.cir[3,],dh.cir[3,],fa.cir[3,],gd.cir[3,],ma.cir[3,],mb.cir[3,],su.cir[3,])



# Order beaches by fresh / marine, then by the head immersion CIR
marinelab <- c(rep(0,4),rep(1,9))
bord <- rev(order(marinelab,hCIRs[,1]))
bCIRs <- bCIRs[bord,]
hCIRs <- hCIRs[bord,]
sCIRs <- sCIRs[bord,]
fCIRlab <- fCIRlab[bord]

# Add in the combined, summary estimates (subgroups + overall)
bCIRs <- rbind(all.cir[1,],fm.cir[1,4:6],bCIRs[1:9,],fm.cir[1,1:3],bCIRs[10:13,])
hCIRs <- rbind(all.cir[2,],fm.cir[2,4:6],hCIRs[1:9,],fm.cir[2,1:3],hCIRs[10:13,])
sCIRs <- rbind(all.cir[3,],fm.cir[3,4:6],sCIRs[1:9,],fm.cir[3,1:3],sCIRs[10:13,])
fCIRlab <- c("All Combined","Marine Combined",fCIRlab[1:9],"Freshwater Combined",fCIRlab[10:13])



pdf("~/dropbox/13beaches/aim1-results/figs/swim-exposure-CIR-forest-allages.pdf",height=6,width=9)

lo <- layout(mat=matrix(1:4,nrow=1,ncol=4),widths=c(1.5,1,1,1))
op <- par(mar=c(5,2,2,2)+0.1)

fcol <- brewer.pal(11,"RdYlGn")[10]
mcol <- brewer.pal(11,"Spectral")[10]
scol <- brewer.pal(11,"Spectral")[11]
cols <- c(scol,rep(mcol,10),rep(fcol,5))
pchs <- c(5,5,rep(21,9),5,rep(21,4))

xtics <- c(0.5,1,2,4,8)

# Beach Labels
bhs <- barplot(1:16,horiz=TRUE,xlim=c(0,1),xaxt="n",col=NA,border=NA)
text(x=rep(1,length(bhs)),y=bhs,fCIRlab,cex=1.5,adj=1,col=c(scol,mcol,rep("gray40",9),fcol,rep("gray40",4)),font=c(2,2,rep(1,9),2,rep(1,4)) )

# Body Immersion
bhs <- barplot(1:16,horiz=TRUE,xlim=range(xtics),log="x",xaxt="n",col=NA,border=NA)
axis(side=1,at=xtics,las=1)
segments(x0=1,x1=1,y0=0,y1=max(bhs+1),lty=2)
mtext("aCIR",side=1,line=3,las=1)
mtext("Body Immersion",side=3,line=0,at=1,col="gray40",font=2)
# scale the area of the point size by 1/SE of the estimates
# calc 1/SE
oneoverSE <- 1.96/(log(bCIRs[,3])-log(bCIRs[,1]))
# now scale to area of a circle rather than radius
# and further scale the cex factor so the smallest = 1
cexscale <- sqrt(oneoverSE/pi)
cexscale <- cexscale/min(cexscale)
# label the combined estimate
segments(x0=bCIRs[1,1],x1=bCIRs[1,1],y0=0,y1=max(bhs+1),lty=2,col=scol)
text(x=bCIRs[1,1],y=0,paste(sprintf("%1.2f",bCIRs[1,1])," (",sprintf("%1.2f",bCIRs[1,2]),", ", sprintf("%1.2f",bCIRs[1,3]) ,")",sep=""),pos=4,cex=0.8,col=scol )
# plot the estimates
segments(y0=bhs,x0=bCIRs[,2],x1=bCIRs[,3],lwd=2,col=cols)
points(bCIRs[,1],bhs,cex=1.5*cexscale,pch=pchs,lwd=2,bg="white",col=cols)


# Head Immersion
bhs <- barplot(1:16,horiz=TRUE,xlim=range(xtics),log="x",xaxt="n",col=NA,border=NA)
axis(side=1,at=xtics,las=1)
segments(x0=1,x1=1,y0=0,y1=max(bhs+1),lty=2)
mtext("aCIR",side=1,line=3,las=1)
# mtext(fCIRlab,side=2,line=0,at=bhs,cex=0.9,las=1)
mtext("Head Immersion",side=3,line=0,at=1,col="gray40",font=2)
# scale the area of the point size by 1/SE of the estimates
# calc 1/SE
oneoverSE <- 1.96/(log(hCIRs[,3])-log(hCIRs[,1]))
# now scale to area of a circle rather than radius
# and further scale the cex factor so the smallest = 1
cexscale <- sqrt(oneoverSE/pi)
cexscale <- cexscale/min(cexscale)
# label the combined estimate
segments(x0=hCIRs[1,1],x1=hCIRs[1,1],y0=0,y1=max(bhs+1),lty=2,col=scol)
text(x=hCIRs[1,1],y=0,paste(sprintf("%1.2f",hCIRs[1,1])," (",sprintf("%1.2f",hCIRs[1,2]),", ", sprintf("%1.2f",hCIRs[1,3]) ,")",sep=""),pos=4,cex=0.8,col=scol )
# plot the estimates
segments(y0=bhs,x0=hCIRs[,2],x1=hCIRs[,3],lwd=2,col=cols)
points(hCIRs[,1],bhs,cex=1.5*cexscale,pch=pchs,lwd=2,bg="white",col=cols)


# Swallowed Water
bhs <- barplot(1:16,horiz=TRUE,xlim=range(xtics),log="x",xaxt="n",col=NA,border=NA)
axis(side=1,at=xtics,las=1)
segments(x0=1,x1=1,y0=0,y1=max(bhs+1),lty=2)
mtext("aCIR",side=1,line=3,las=1)
# mtext(fCIRlab,side=2,line=0,at=bhs,cex=0.9,las=1)
mtext("Swallowed Water",side=3,line=0,at=1,col="gray40",font=2)
# scale the area of the point size by 1/SE of the estimates
# calc 1/SE
oneoverSE <- 1.96/(log(sCIRs[,3])-log(sCIRs[,1]))
# now scale to area of a circle rather than radius
# and further scale the cex factor so the smallest = 1
cexscale <- sqrt(oneoverSE/pi)
cexscale <- cexscale/min(cexscale)
# label the combined estimate
segments(x0=sCIRs[1,1],x1=sCIRs[1,1],y0=0,y1=max(bhs+1),lty=2,col=scol)
text(x=sCIRs[1,1],y=0,paste(sprintf("%1.2f",sCIRs[1,1])," (",sprintf("%1.2f",sCIRs[1,2]),", ", sprintf("%1.2f",sCIRs[1,3]) ,")",sep=""),pos=4,cex=0.8,col=scol )
# plot the estimates
segments(y0=bhs,x0=sCIRs[,2],x1=sCIRs[,3],lwd=2,col=cols)
points(sCIRs[,1],bhs,cex=1.5*cexscale,pch=pchs,lwd=2,bg="white",col=cols)

par(op)
dev.off()

# --------------------------------------
# Load the age-stratified results
# Ages 0-10 years
# --------------------------------------


load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-swim-exposure-regs-age0to10.Rdata")

# beach-specific CIR estimates
CIRlab <- c("Huntington","Silver","Washington Park","West","Avalon","Boqueron","Edgewater","Doheny","Fairhope","Goddard","Malibu","Mission Bay","Surfside")
bCIRs.0to10 <- rbind(hu.cir[1,],si.cir[1,],wp.cir[1,],we.cir[1,],av.cir[1,],bo.cir[1,],ed.cir[1,],dh.cir[1,],fa.cir[1,],gd.cir[1,],ma.cir[1,],mb.cir[1,],su.cir[1,])

hCIRs.0to10 <- rbind(hu.cir[2,],si.cir[2,],wp.cir[2,],we.cir[2,],av.cir[2,],bo.cir[2,],ed.cir[2,],dh.cir[2,],fa.cir[2,],gd.cir[2,],ma.cir[2,],mb.cir[2,],su.cir[2,])

sCIRs.0to10 <- rbind(hu.cir[3,],si.cir[3,],wp.cir[3,],we.cir[3,],av.cir[3,],bo.cir[3,],ed.cir[3,],dh.cir[3,],fa.cir[3,],gd.cir[3,],ma.cir[3,],mb.cir[3,],su.cir[3,])

# combined CIR estimates
all.cir.0to10 <- all.cir
fm.cir.0to10 <- fm.cir

# --------------------------------------
# Load the age-stratified results
# Ages 11 + years
# --------------------------------------


load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-swim-exposure-regs-age11plus.Rdata")

# beach-specific CIR estimates
CIRlab <- c("Huntington","Silver","Washington Park","West","Avalon","Boqueron","Edgewater","Doheny","Fairhope","Goddard","Malibu","Mission Bay","Surfside")

bCIRs.11plus <- rbind(hu.cir[1,],si.cir[1,],wp.cir[1,],we.cir[1,],av.cir[1,],bo.cir[1,],ed.cir[1,],dh.cir[1,],fa.cir[1,],gd.cir[1,],ma.cir[1,],mb.cir[1,],su.cir[1,])

hCIRs.11plus <- rbind(hu.cir[2,],si.cir[2,],wp.cir[2,],we.cir[2,],av.cir[2,],bo.cir[2,],ed.cir[2,],dh.cir[2,],fa.cir[2,],gd.cir[2,],ma.cir[2,],mb.cir[2,],su.cir[2,])

sCIRs.11plus <- rbind(hu.cir[3,],si.cir[3,],wp.cir[3,],we.cir[3,],av.cir[3,],bo.cir[3,],ed.cir[3,],dh.cir[3,],fa.cir[3,],gd.cir[3,],ma.cir[3,],mb.cir[3,],su.cir[3,])

# combined CIR estimates
all.cir.11plus <- all.cir
fm.cir.11plus <- fm.cir



# --------------------------------------
# AGE-STRATIFIED
# Forest plot of CIRs
# general function (repeated below)
# --------------------------------------


ageforest <- function(CIR0to10,CIR11plus,CIRlab)
{
	# CIR0to10  : matrix of CIRs, columns = CIR, lower95, upper95 for ages 0 to 10
	# CIR11plus : matrix of CIRs, columns = CIR, lower95, upper95 for ages 0 to 10
	# CIRlab    : beach and summary text labels for CIR estimates
	
	# forest plot specs (dimensions, colors, plot range)
	lo <- layout(mat=matrix(1:3,nrow=1,ncol=3),widths=c(1,1,1))
	op <- par(mar=c(5,2,2,2)+0.1)
	fcol <- brewer.pal(11,"RdYlGn")[10]
	mcol <- brewer.pal(11,"Spectral")[10]
	scol <- brewer.pal(11,"Spectral")[11]
	cols <- c(scol,rep(mcol,10),rep(fcol,5))
	pchs <- c(5,5,rep(21,9),5,rep(21,4))
	xtics <- c(0.5,1,2,4,8)
	
	# Beach Labels
	bhs <- barplot(1:16,horiz=TRUE,xlim=c(0,1),xaxt="n",col=NA,border=NA)
	text(x=rep(1,length(bhs)),y=bhs,CIRlab,cex=1.5,adj=1,col=c(scol,mcol,rep("gray40",9),fcol,rep("gray40",4)),font=c(2,2,rep(1,9),2,rep(1,4)) )
	
	
	# Ages 0 to 10 years
	bhs <- barplot(1:16,horiz=TRUE,xlim=range(xtics),log="x",xaxt="n",col=NA,border=NA)
	axis(side=1,at=xtics,las=1)
	segments(x0=1,x1=1,y0=0,y1=max(bhs+1),lty=2)
	mtext("aCIR",side=1,line=3,las=1)
	mtext("Age 0 to 10 years",side=3,line=0,at=1,col="gray40",font=2)
	# label the combined estimate
	segments(x0=CIR0to10[1,1],x1= CIR0to10[1,1],y0=0,y1=max(bhs+1),lty=2,col=scol)
	text(x=CIR0to10[1,1],y=0,paste(sprintf("%1.2f",CIR0to10[1,1])," (",sprintf("%1.2f",CIR0to10[1,2]),", ", sprintf("%1.2f",CIR0to10[1,3]) ,")",sep=""),pos=4,cex=0.8,col=scol )
	# scale the area of the point size by 1/SE of the estimates
	# calc 1/SE
	oneoverSE <- 1.96/(log(CIR0to10[,3])-log(CIR0to10[,1]))
	# now scale to area of a circle rather than radius
	# and further scale the cex factor so the smallest = 1
	cexscale <- sqrt(oneoverSE/pi)
	cexscale <- cexscale/min(cexscale)
	# plot the estimates
	segments(y0=bhs,x0=CIR0to10[,2],x1=CIR0to10[,3],lwd=2,col=cols)
	points(CIR0to10[,1],bhs,cex=1.5*cexscale,pch=pchs,lwd=2,bg="white",col=cols)
	
	
	# Ages 11 years and older
	bhs <- barplot(1:16,horiz=TRUE,xlim=range(xtics),log="x",xaxt="n",col=NA,border=NA)
	axis(side=1,at=xtics,las=1)
	segments(x0=1,x1=1,y0=0,y1=max(bhs+1),lty=2)
	mtext("aCIR",side=1,line=3,las=1)
	mtext("Age 11 years and older",side=3,line=0,at=1,col="gray40",font=2)
	# label the combined estimate
	segments(x0=CIR11plus[1,1],x1=CIR11plus[1,1],y0=0,y1=max(bhs+1),lty=2,col=scol)
	text(x=CIR11plus[1,1],y=0,paste(sprintf("%1.2f",CIR11plus[1,1])," (",sprintf("%1.2f",CIR11plus[1,2]),", ", sprintf("%1.2f",CIR11plus[1,3]) ,")",sep=""),pos=4,cex=0.8,col=scol )
	# scale the area of the point size by 1/SE of the estimates
	# calc 1/SE
	oneoverSE <- 1.96/(log(CIR11plus[,3])-log(CIR11plus[,1]))
	# now scale to area of a circle rather than radius
	# and further scale the cex factor so the smallest = 1 
	cexscale <- sqrt(oneoverSE/pi)
	cexscale <- cexscale/min(cexscale)
	# plot the estimates
	segments(y0=bhs,x0=CIR11plus[,2],x1=CIR11plus[,3],lwd=2,col=cols)
	points(CIR11plus[,1],bhs,cex=1.5*cexscale,pch=pchs,lwd=2,bg="white",col=cols)
	
	par(op)
}


# --------------------------------------
# Body Immersion
# --------------------------------------

# Order beaches by fresh / marine, then by the CIR for <=10 year olds
marinelab <- c(rep(0,4),rep(1,9))
bord <- rev(order(marinelab,bCIRs.0to10[,1]))
bCIRs.0to10 <- bCIRs.0to10[bord,]
bCIRs.11plus <- bCIRs.11plus[bord,]
bCIRlab <- CIRlab[bord]

# Add in the combined, summary estimates (subgroups + overall)
pCIRs.0to10 <- rbind(all.cir.0to10[1,],fm.cir.0to10[1,4:6],bCIRs.0to10[1:9,],fm.cir.0to10[1,1:3],bCIRs.0to10[10:13,])
pCIRs.11plus <- rbind(all.cir.11plus[1,],fm.cir.11plus[1,4:6],bCIRs.11plus[1:9,],fm.cir.11plus[1,1:3],bCIRs.11plus[10:13,])
pCIRlab <- c("All Combined","Marine Combined",bCIRlab[1:9],"Freshwater Combined",bCIRlab[10:13])

pdf("~/dropbox/13beaches/aim1-results/figs/swim-exposure-CIR-forest-agestrat-body.pdf",height=6,width=7)
	ageforest(CIR0to10=pCIRs.0to10,CIR11plus=pCIRs.11plus,CIRlab=pCIRlab)
dev.off()

# --------------------------------------
# Head Immersion
# --------------------------------------


# Order beaches by fresh / marine, then by the CIR for <=10 year olds
marinelab <- c(rep(0,4),rep(1,9))
hord <- rev(order(marinelab,hCIRs.0to10[,1]))
hCIRs.0to10 <- hCIRs.0to10[hord,]
hCIRs.11plus <- hCIRs.11plus[hord,]
hCIRlab <- CIRlab[hord]


# Add in the combined, summary estimates (subgroups + overall)
pCIRs.0to10 <- rbind(all.cir.0to10[2,],fm.cir.0to10[2,4:6],hCIRs.0to10[1:9,],fm.cir.0to10[2,1:3],hCIRs.0to10[10:13,])
pCIRs.11plus <- rbind(all.cir.11plus[2,],fm.cir.11plus[2,4:6],hCIRs.11plus[1:9,],fm.cir.11plus[2,1:3],hCIRs.11plus[10:13,])
pCIRlab <- c("All Combined","Marine Combined",hCIRlab[1:9],"Freshwater Combined",hCIRlab[10:13])


pdf("~/dropbox/13beaches/aim1-results/figs/swim-exposure-CIR-forest-agestrat-head.pdf",height=6,width=7)
	ageforest(CIR0to10=pCIRs.0to10,CIR11plus=pCIRs.11plus,CIRlab=pCIRlab)
dev.off()


# --------------------------------------
# Swallowed Water
# --------------------------------------

# Order beaches by fresh / marine, then by the CIR for <=10 year olds
marinelab <- c(rep(0,4),rep(1,9))
sord <- rev(order(marinelab,sCIRs.0to10[,1]))
sCIRs.0to10 <- sCIRs.0to10[sord,]
sCIRs.11plus <- sCIRs.11plus[sord,]
sCIRlab <- CIRlab[sord]


# Add in the combined, summary estimates (subgroups + overall)
pCIRs.0to10 <- rbind(all.cir.0to10[3,],fm.cir.0to10[3,4:6],sCIRs.0to10[1:9,],fm.cir.0to10[3,1:3],sCIRs.0to10[10:13,])
pCIRs.11plus <- rbind(all.cir.11plus[3,],fm.cir.11plus[3,4:6],sCIRs.11plus[1:9,],fm.cir.11plus[3,1:3],sCIRs.11plus[10:13,])
pCIRlab <- c("All Combined","Marine Combined",sCIRlab[1:9],"Freshwater Combined",sCIRlab[10:13])

pdf("~/dropbox/13beaches/aim1-results/figs/swim-exposure-CIR-forest-agestrat-swal.pdf",height=6,width=7)
	ageforest(CIR0to10=pCIRs.0to10,CIR11plus=pCIRs.11plus,CIRlab=pCIRlab)
dev.off()



