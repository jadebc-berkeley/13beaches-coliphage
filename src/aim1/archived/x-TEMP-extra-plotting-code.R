





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
 

midQs <- minQs + (maxQs-minQs)/2

# --------------------------------------
# Load the regression output
# --------------------------------------

load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-entero1600-Quartile-regs-body.Rdata")




# --------------------------------------
# CIR plot of summary estimates
# plots CIR by age and point / non-point
# source
# --------------------------------------

CIRplot <- function(npsCIR,psCIR,poolCIR,title,ylab=FALSE,xlab=FALSE) {
	
	cols <- brewer.pal(9,"YlGnBu")[6:8]
	
	ytics <- c(0.2,0.5,1,2,4)
	# set up an empty plot
	MidPts <- barplot(1:3,names.arg=NA,border=NA,col=NA,
		log="y",ylog=T,ylim=c(0.2,4.5),ylab="",yaxt="n",
		las=1,bty="n"
		)
		axis(2,at=ytics,las=1)
		segments(x0=0,x1=max(MidPts+0.4),y0=1,lty=2)
		mtext(title,side=2,line=6,at=4,font=2,adj=0.5,las=1,cex=1.5)
		mtext("aCIR",side=3,line=0.5,at=-0.2)
		if(ylab==TRUE) mtext(labQs[2:4],side=3,line=1,at=MidPts,col=cols)
		if(xlab==TRUE) {
			mtext(labQs[2:4],side=1,line=1,at=MidPts,col=cols)
			mtext(expression(paste(italic("Enterococcus")," EPA 1600 Quartile (range CFU/100ml)")),side=1,line=3 )
		}
		
		# plot non-point source estimates
		npx <- MidPts-0.2
		segments(x0=npx,y0=npsCIR[,"CIRlb"],y1=npsCIR[,"CIRub"],lwd=1.5,col=cols)
		points(npx,npsCIR[,"CIR"],pch=23,bg="white",cex=1.5,lwd=1.5,col=cols)
		
		# plot point source estimates
		psx <- MidPts
		segments(x0=psx,y0=psCIR[,"CIRlb"],y1=psCIR[,"CIRub"],lwd=1.5,col=cols)
		points(psx,psCIR[,"CIR"],pch=24,bg="white",cex=1.5,lwd=1.5,col=cols)
		
		# plot overall estimates
		ox <- MidPts+0.2
		segments(x0=ox,y0=poolCIR[,"CIRlb"],y1=poolCIR[,"CIRub"],lwd=1.5,col=cols)
		points(ox,poolCIR[,"CIR"],pch=16,cex=1.5,lwd=1.5,col=cols)
		
}

lo <- layout(mat=matrix(1:4,nrow=4,ncol=1))
op <- par(mar=c(2,10,4,2))
CIRplot(cir.nps,cir.ps,cir.all,title="All\nAges",ylab=TRUE)
op <- par(mar=c(2,10,2,2))
CIRplot(cir.nps0to4,cir.ps0to4,cir.age0to4,title="Ages\n0 to 4")
CIRplot(cir.nps5to10,cir.ps5to10,cir.age5to10,title="Ages\n5 to 10")
op <- par(mar=c(4,10,2,2))
CIRplot(cir.nps11plus,cir.ps11plus,cir.age11plus,title="Ages\n>10",xlab=TRUE)
par(op)

