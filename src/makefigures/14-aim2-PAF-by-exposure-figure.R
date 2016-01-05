# --------------------------------------
# 14-aim2-PAF-by-exposure-figure.R
# ben arnold (benarnold@berkeley.edu)
#
# description:
# plot the population attributable fraction (PAF)
# that corresponds with different %s of exposure
# in the population for swimming exposure and
# for swimming in contaminated water exposure
#
# --------------------------------------

# --------------------------------------
# input files:
#	aim1-entero1600-35cfu-regs-body.RData
#	aim1-swim-exposure-regs-body.Rdata
#	aim2-PARswimex-diar.Rdata
#
# output files:
# aim2-PAF-swimex-diarrhea.pdf
#	aim2-PAF-entero1600-diarrhea.pdf
# aim2-PAF-entero1600noswim-diarrhea.pdf
# --------------------------------------

# --------------------------------------
# preamble
# --------------------------------------

rm(list=ls())
library(RColorBrewer)

# source the base functions
source("~/dropbox/13beaches/src/aim1/0-aim1-base-functions.R")

# --------------------------------------
# PAF function
# calculate the PAF given an RR and
# a proportion exposed
# --------------------------------------

PAFcalc <- function(Pe,RR) {
	# calculate the population attributable fraction
	# under a counterfactual where the exposure is removed
	# see section 4.7 of Jewell 2004 Statistics for Epidemiology
	# RR : relative risk associated with exposure
	# Pe : proportion of the population exposed
	PAF <- Pe*(RR-1) / (1 + Pe*(RR-1))
	return(PAF)

}

PAFcalc2 <- function(Pe,p0,pobs,RR) {
	# calculate the population attributable fraction
	# for an exposure variable with three levels (0, 1, 2)
	# under a counterfactual where exposure level 2 is re-assigned 
	# to level 0 (leaving level 1 the same)
	# This is akin to assigning swimmers with Entero >35 CFU to
	# be non-swimmers
	# see section 4.7 of Jewell 2004 Statistics for Epidemiology
	# for the basic sketch of how to solve for this, though this
	# is a slightly more complicated case
	# note that PAF = [P(obs) - P(counterfactual)] / P(obs)
	# where P(obs) = Pe_0*p0 + Pe_1*p1 + Pe_2*p2
	# and P(count) = (Pe_0 + Pe_2)*p0 + Pe_1*p1
	# Thus, PAF = Pe_2(p2-p0)/P(obs) = Pe_2*p0(RR-1)/P(obs)
	# 
	# Pe   : proportion of the population exposed to level 2
	# p0   : risk in the unexposed
	# pobs : observed risk in the population
	# RR   : CIR for exposure level 2 versus unexposed (level 0)
	PAF <- Pe*p0*(RR-1)/pobs
	return(PAF)
}



# --------------------------------------
# load the swim exposure regression output
# and PAR output (to get total Ns)
# --------------------------------------

load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-swim-exposure-regs-body.Rdata")
load("~/dropbox/13beaches/aim2-results/rawoutput/aim2-PARswimex-diar.Rdata")





# --------------------------------------
# Calculate the PAF across a range of
# exposures for different subgroups
# --------------------------------------
Pe <- seq(0,1,by=0.01)

paf.diar.age0to4   <- sapply(Pe,PAFcalc,RR=CIRs$CIRoverall[2,1])
paf.diar.age5to10  <- sapply(Pe,PAFcalc,RR=CIRs$CIRoverall[3,1])
paf.diar.age11plus <- sapply(Pe,PAFcalc,RR=CIRs$CIRoverall[4,1])

# print model estimates
cbind(Pe,paf.diar.age0to4,paf.diar.age5to10,paf.diar.age11plus)

# --------------------------------------
# Grab the empirical exposure levels to
# get actual PAFs to plot points
# --------------------------------------
pe <- N.body[,2]/N.total[,2]

paf.diar.age0to4e   <- PAFcalc(Pe=pe[2],RR=CIRs$CIRoverall[2,1])
paf.diar.age5to10e  <- PAFcalc(Pe=pe[3],RR=CIRs$CIRoverall[3,1])
paf.diar.age11pluse <- PAFcalc(Pe=pe[4],RR=CIRs$CIRoverall[4,1])

# --------------------------------------
# Plot the curves and points for 
# diarrhea and days missed of daily
# activities
# --------------------------------------

pdf("~/dropbox/13beaches/aim2-results/figs/aim2-PAF-swimex-diarrhea.pdf")
# lo <- layout(mat=matrix(1:2,nrow=1,ncol=2))
op <- par(mar=c(4,4,4,2)+0.1)
ytics <- seq(0,0.6,by=0.1)
cols <- brewer.pal(5,"Dark2")

plot(Pe,paf.diar.age0to4,type="n",
	bty="n",
	xaxt="n",xlab="",xlim=c(0,1),
	yaxt="n",ylab="",ylim=range(ytics),
	las=1
	)
	axis(1,at=seq(0,1,by=0.1),labels=seq(0,100,by=10))
	axis(2,at=ytics,labels=sprintf("%1.0f",ytics*100),las=1)
	mtext("Proportion of Beachgoers Exposed to Swimming (%)",side=1,line=3)
	mtext("Diarrhea Episodes",side=3,line=1,at=-0.14,adj=0)
	mtext("Population Attributable Fraction (%)",side=3,line=0,at=-0.14,adj=0)
	# mtext("A)",side=3,line=2,at=-0.14,adj=0,font=1,cex=1.25)

	lines(Pe,paf.diar.age0to4,lty=1,lwd=2,col=cols[1])
	lines(Pe,paf.diar.age5to10,lty=1,lwd=2,col=cols[2])
	lines(Pe,paf.diar.age11plus,lty=1,lwd=2,col=cols[3])
	
	points(pe[2],paf.diar.age0to4e,pch=21,bg="white",lwd=2,cex=1.4,col=cols[1])
	points(pe[3],paf.diar.age5to10e,pch=21,bg="white",lwd=2,cex=1.4,col=cols[2])
	points(pe[4],paf.diar.age11pluse,pch=21,bg="white",lwd=2,cex=1.4,col=cols[3])
	
	text(1,0.35,"Ages 0 - 4",adj=1,col=cols[1])
	text(1,0.52,"Ages 5 - 10",adj=1,col=cols[2])
	text(1,0.275,"Ages >10",adj=1,col=cols[3])
	
	segments(x0=pe[3]-0.03,x1=pe[3]-0.07,y0=paf.diar.age5to10e,col="gray60")
	text(pe[3]-0.07,paf.diar.age5to10e,"Circles mark observed exposure and PAF estimates",col="gray40",cex=0.8,pos=2)


par(op)
dev.off()


# --------------------------------------
# load the Entero 1600 regression output
# --------------------------------------

load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-entero1600-35cfu-regs-body.Rdata")


# --------------------------------------
# Calculate the PAF across a range of
# exposures for different subgroups
# --------------------------------------
Pe <- seq(0,1,by=0.01)

paf.age0to4   <- sapply(Pe,PAFcalc,RR=cir.age0to4[3,1])
paf.age5to10  <- sapply(Pe,PAFcalc,RR=cir.age5to10[3,1])
paf.age11plus <- sapply(Pe,PAFcalc,RR=cir.age11plus[3,1])

paf.ps  <- sapply(Pe,PAFcalc,RR=cir.all[2,1])
paf.nps <- sapply(Pe,PAFcalc,RR=cir.all[1,1])

# print model estimates
cbind(Pe,paf.age0to4,paf.age5to10,paf.age11plus)

# print model estimates
cbind(Pe,paf.ps,paf.nps)

# --------------------------------------
# Grab the empirical exposure levels to
# get actual PAFs to plot points
# --------------------------------------
pe.0to4   <- N.age0to4[2,2]/(sum(N.age0to4[,2]))
pe.5to10  <- N.age5to10[2,2]/(sum(N.age5to10[,2]))
pe.11plus <- N.age11plus[2,2]/(sum(N.age11plus[,2]))
pe.ps     <- N.ps[2,2]/(sum(N.ps[,2]))
pe.nps    <- N.nps[2,2]/(sum(N.nps[,2]))

paf.age0to4e   <- PAFcalc(Pe=pe.0to4,RR=cir.age0to4[3,1])
paf.age5to10e  <- PAFcalc(Pe=pe.5to10,RR=cir.age5to10[3,1])
paf.age11pluse <- PAFcalc(Pe=pe.11plus,RR=cir.age11plus[3,1])
paf.pse <- PAFcalc(Pe=pe.ps,RR=cir.all[2,1])
paf.npse <- PAFcalc(Pe=pe.nps,RR=cir.all[1,1])

# --------------------------------------
# Plot the curves and points
# --------------------------------------

pdf("~/dropbox/13beaches/aim2-results/figs/aim2-PAF-entero1600-diarrhea.pdf",width=11,height=5)
lo <- layout(mat=matrix(1:2,nrow=1,ncol=2))
op <- par(mar=c(4,4,4,2)+0.1)
ytics <- seq(0,0.25,by=0.05)
cols <- brewer.pal(5,"Dark2")

plot(Pe,paf.age0to4,type="n",
	bty="n",
	xaxt="n",xlab="",xlim=c(0,1),
	yaxt="n",ylab="",ylim=range(ytics),
	las=1
	)
	axis(1,at=seq(0,1,by=0.1),labels=seq(0,100,by=10))
	axis(2,at=ytics,labels=sprintf("%1.0f",ytics*100),las=1)
	mtext("Proportion of Swimmers Exposed >35 CFU/100ml (%)",side=1,line=3)	
	mtext("Diarrhea Episodes by Age Group",side=3,line=1,at=-0.14,adj=0)
	mtext("Population Attributable Fraction (%)",side=3,line=0,at=-0.14,adj=0)
	mtext("A)",side=3,line=2,at=-0.14,adj=0,font=1,cex=1.25)
	

	lines(Pe,paf.age0to4,lty=1,lwd=2,col=cols[1])
	lines(Pe,paf.age5to10,lty=1,lwd=2,col=cols[2])
	lines(Pe,paf.age11plus,lty=1,lwd=2,col=cols[3])
	
	points(pe.0to4,paf.age0to4e,pch=21,bg="white",lwd=2,cex=1.4,col=cols[1])
	points(pe.5to10,paf.age5to10e,pch=21,bg="white",lwd=2,cex=1.4,col=cols[2])
	points(pe.11plus,paf.age11pluse,pch=21,bg="white",lwd=2,cex=1.4,col=cols[3])
	
	text(1,0.255,"Ages 0 - 4",adj=1,srt=18,col=cols[1])
	text(1,0.215,"Ages 5 - 10",adj=1,srt=18,col=cols[2])
	text(1,0.16,"Ages >10",adj=1,srt=18,col=cols[3])
	
	segments(x0=pe.0to4,y0=paf.age0to4e+0.01,y1=paf.age0to4e+0.1,col="gray60")
	text(pe.0to4,paf.age0to4e+0.1,"Circles mark\n observed exposure\n and PAF estimates",col="gray40",cex=0.8,pos=3)

plot(Pe,paf.ps,type="n",
	bty="n",
	xaxt="n",xlab="",xlim=c(0,1),
	yaxt="n",ylab="",ylim=range(ytics),
	las=1
	)
	axis(1,at=seq(0,1,by=0.1),labels=seq(0,100,by=10))
	axis(2,at=ytics,labels=sprintf("%1.0f",ytics*100),las=1)
	mtext("Proportion of Swimmers Exposed >35 CFU/100ml (%)",side=1,line=3)
	mtext("Diarrhea Episodes by Pollution Type",side=3,line=1,at=-0.14,adj=0)
	mtext("Population Attributable Fraction (%)",side=3,line=0,at=-0.14,adj=0)
	mtext("B)",side=3,line=2,at=-0.14,adj=0,font=1,cex=1.25)
	
	lines(Pe,paf.ps,lty=1,lwd=2,col=cols[4])
	lines(Pe,paf.nps,lty=1,lwd=2,col=cols[5])
	
	points(pe.ps, paf.pse,pch=21,bg="white",lwd=2,cex=1.4,col=cols[4])
	points(pe.nps,paf.npse,pch=21,bg="white",lwd=2,cex=1.4,col=cols[5])

	text(1,0.22,"Point Source",adj=1,srt=0,col=cols[4])
	text(1,0.05,"Non-Point Source",adj=1,srt=0,col=cols[5])

par(op)
dev.off()



# --------------------------------------
# load the Entero 1600 regression output
# with non-swimmers as the reference
# category
#
# The PAF calculations are slightly
# more complicated because of the 
# 3-level, categorical exposure
# --------------------------------------

# load the diarrhea CIR estimates to generate the curves
load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-entero1600-noswim35cfu-regs-body.Rdata")

# load the full analysis dataset to calculate the baseline
# risk of diarrhea and GI illness among non-swimmers
ad <- preprocess.13beaches("~/dropbox/13beaches/data/final/13beaches-analysis.csv")
	table(ad$swim35)
	table(is.na(ad$swim35))
ad <- subset(ad,is.na(ad$swim35)==FALSE)
	dim(ad)

# grab the overall prevalence
pobs.diar <- tapply(ad$diarrheaci10,ad$agestrat,mean)
pobs.gi   <- tapply(ad$gici10,ad$agestrat,mean)

# grab the prevalence in nonswimmers by age category
p0.diar <- tapply(ad$diarrheaci10[ad$anycontact=="No"],ad$agestrat[ad$anycontact=="No"],mean)
p0.gi   <- tapply(ad$gici10[ad$anycontact=="No"],ad$agestrat[ad$anycontact=="No"],mean)

# Grab empirical exposure %s to plot them
pe.tab <- table(ad$swim35,ad$agestrat)
pe35plus <- pe.tab[3,]/colSums(pe.tab)

# --------------------------------------
# Calculate the PAF across a range of
# exposures for different subgroups
# --------------------------------------
Pe <- seq(0,0.5,by=0.01)
RR.diar.0to4   <- exp(age0to4$fit[3,1])
RR.diar.5to10  <- exp(age5to10$fit[3,1])
RR.diar.11plus <- exp(age11plus$fit[3,1])

paf.age0to4    <- sapply(Pe,PAFcalc2,p0=p0.diar[3],pobs=pobs.diar[3],RR=RR.diar.0to4)
paf.age5to10   <- sapply(Pe,PAFcalc2,p0=p0.diar[2],pobs=pobs.diar[2],RR=RR.diar.5to10)
paf.age11plus  <- sapply(Pe,PAFcalc2,p0=p0.diar[1],pobs=pobs.diar[1],RR=RR.diar.11plus)

# print model estimates
cbind(Pe,paf.age0to4,paf.age5to10,paf.age11plus)

# --------------------------------------
# Grab the empirical exposure levels to
# get actual PAFs to plot points
# --------------------------------------
paf.age0to4e   <- PAFcalc2(Pe=pe35plus[3],p0=p0.diar[3],pobs=pobs.diar[3],RR=RR.diar.0to4)
paf.age5to10e  <- PAFcalc2(Pe=pe35plus[2],p0=p0.diar[2],pobs=pobs.diar[2],RR=RR.diar.5to10)
paf.age11pluse <- PAFcalc2(Pe=pe35plus[1],p0=p0.diar[1],pobs=pobs.diar[1],RR=RR.diar.11plus)


# --------------------------------------
# load the GI illness CIR estimates to 
# generate the GI curves
# --------------------------------------
load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-entero1600-noswim35cfu-regs-body-GI-illness.Rdata")

# --------------------------------------
# Calculate the PAF across a range of
# exposures for different subgroups
# --------------------------------------
Pe <- seq(0,0.5,by=0.01)
RR.gi.0to4   <- exp(age0to4$fit[3,1])
RR.gi.5to10  <- exp(age5to10$fit[3,1])
RR.gi.11plus <- exp(age11plus$fit[3,1])

paf.gi.age0to4    <- sapply(Pe,PAFcalc2,p0=p0.gi[3],pobs=pobs.gi[3],RR=RR.gi.0to4)
paf.gi.age5to10   <- sapply(Pe,PAFcalc2,p0=p0.gi[2],pobs=pobs.gi[2],RR=RR.gi.5to10)
paf.gi.age11plus  <- sapply(Pe,PAFcalc2,p0=p0.gi[1],pobs=pobs.gi[1],RR=RR.gi.11plus)

# print model estimates
cbind(Pe,paf.gi.age0to4,paf.gi.age5to10,paf.gi.age11plus)

# --------------------------------------
# Grab the empirical exposure levels to
# get actual PAFs to plot points
# --------------------------------------
paf.gi.age0to4e   <- PAFcalc2(Pe=pe35plus[3],p0=p0.gi[3],pobs=pobs.gi[3],RR=RR.gi.0to4)
paf.gi.age5to10e  <- PAFcalc2(Pe=pe35plus[2],p0=p0.gi[2],pobs=pobs.gi[2],RR=RR.gi.5to10)
paf.gi.age11pluse <- PAFcalc2(Pe=pe35plus[1],p0=p0.gi[1],pobs=pobs.gi[1],RR=RR.gi.11plus)




pdf("~/dropbox/13beaches/aim2-results/figs/aim2-PAF-entero1600noswim.pdf",width=6,height=11)
lo <- layout(mat=matrix(1:2,nrow=2,ncol=1))
op <- par(mar=c(4,4,4,2)+0.1)
ytics <- seq(0,0.4,by=0.1)
cols <- brewer.pal(5,"Dark2")

# Diarrhea
plot(Pe,paf.age0to4,type="n",
	bty="n",
	xaxt="n",xlab="",xlim=c(0,max(Pe)),
	yaxt="n",ylab="",ylim=range(ytics),
	las=1
	)
	axis(1,at=seq(0,max(Pe),by=0.1),labels=seq(0,max(Pe)*100,by=10))
	axis(2,at=ytics,labels=sprintf("%1.0f",ytics*100),las=1)
	mtext("Proportion of Beachgoers Exposed >35 CFU/100ml (%)",side=1,line=3)	
	mtext("Diarrhea Episodes by Age",side=3,line=0,at=-0.06,adj=0)
	mtext("Population Attributable Fraction (%)",side=3,line=1,at=-0.06,adj=0)
	mtext("A)",side=3,line=2,at=-0.06,adj=0,font=1,cex=1.25)
	

	lines(Pe,paf.age0to4,lty=1,lwd=2,col=cols[1])
	lines(Pe,paf.age5to10,lty=1,lwd=2,col=cols[2])
	lines(Pe,paf.age11plus,lty=1,lwd=2,col=cols[3])
	
	points(pe35plus[3],paf.age0to4e,pch=21,bg="white",lwd=2,cex=1.4,col=cols[1])
	points(pe35plus[2],paf.age5to10e,pch=21,bg="white",lwd=2,cex=1.4,col=cols[2])
	points(pe35plus[1],paf.age11pluse,pch=21,bg="white",lwd=2,cex=1.4,col=cols[3])
	
	text(max(Pe),max(paf.age0to4)+0.015,"Ages 0 - 4",adj=1,srt=20,col=cols[1])
	text(max(Pe),max(paf.age5to10)+0.02,"Ages 5 - 10",adj=1,srt=20,col=cols[2])
	text(max(Pe),max(paf.age11plus)+0.02,"Ages >10",adj=1,srt=20,col=cols[3])
	
	segments(x0=pe35plus[2],y0=paf.age5to10e+0.02,y1=paf.age5to10e+0.15,col="gray60")
	text(pe35plus[2],paf.age5to10e+0.15,"Circles mark\n observed exposure\n and PAF estimates",col="gray40",cex=0.8,pos=3)

# GI illness
plot(Pe,paf.age0to4,type="n",
	bty="n",
	xaxt="n",xlab="",xlim=c(0,max(Pe)),
	yaxt="n",ylab="",ylim=range(ytics),
	las=1
	)
	axis(1,at=seq(0,max(Pe),by=0.1),labels=seq(0,max(Pe)*100,by=10))
	axis(2,at=ytics,labels=sprintf("%1.0f",ytics*100),las=1)
	mtext("Proportion of Beachgoers Exposed >35 CFU/100ml (%)",side=1,line=3)	
	mtext("Gastrointestinal Illness Episodes by Age",side=3,line=0,at=-0.06,adj=0)
	mtext("Population Attributable Fraction (%)",side=3,line=1,at=-0.06,adj=0)
	mtext("B)",side=3,line=2,at=-0.06,adj=0,font=1,cex=1.25)
	

	lines(Pe,paf.gi.age0to4,lty=1,lwd=2,col=cols[1])
	lines(Pe,paf.gi.age5to10,lty=1,lwd=2,col=cols[2])
	lines(Pe,paf.gi.age11plus,lty=1,lwd=2,col=cols[3])
	
	points(pe35plus[3],paf.gi.age0to4e,pch=21,bg="white",lwd=2,cex=1.4,col=cols[1])
	points(pe35plus[2],paf.gi.age5to10e,pch=21,bg="white",lwd=2,cex=1.4,col=cols[2])
	points(pe35plus[1],paf.gi.age11pluse,pch=21,bg="white",lwd=2,cex=1.4,col=cols[3])
	
	text(max(Pe),max(paf.gi.age0to4)+0.02,"Ages 0 - 4",adj=1,srt=15,col=cols[1])
	text(max(Pe),max(paf.gi.age5to10)+0.02,"Ages 5 - 10",adj=1,srt=15,col=cols[2])
	text(max(Pe),max(paf.gi.age11plus)-0.03,"Ages >10",adj=1,srt=15,col=cols[3])


par(op)
dev.off()


