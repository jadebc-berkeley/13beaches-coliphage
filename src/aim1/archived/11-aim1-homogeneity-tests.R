
# ------------------------------------
# aim1-11-homogeneity-tests.R
#
# Calculate Cochran's Q statistic and
# the related Higgens' I squared for
# swim exposure analyses and for
# the Enterococcus qunitle analyses
#
# version 3 (1 jun 2015)
# updated script to handle updated raw output
#
# version 2 (21 feb 2015)
# added tests for enterococcus 1600
#
# version 1 (16 feb 2015)
# ------------------------------------

# ------------------------------------
# preamble
# ------------------------------------
rm(list=ls())


# ------------------------------------
# Swim exposure analyses
# ------------------------------------

# ------------------------------------
# Calculate CIRs from beach output
# ------------------------------------
# swim exposure regs w/o interactions
# (need to sum 2 coefficients: anycontact + higher level of contact)
CIR.swimex <- function(fo,vcv) {
	# fo : fit object returned from lbreg or mpreg
	# vcv: variance covariance matrix from lbreg or mpreg
	nr <- nrow(vcv)
	lc <- c(0,1,1,rep(0,nr-3)) # linear combination of betas for the estimate
	est <- exp(t(lc)%*%fo[,1])
	se <- sqrt( t(lc)%*%vcv%*%lc )
	lb <- exp(log(est)-1.96*se)
	ub <- exp(log(est)+1.96*se)
	res <- c(est,lb,ub)
	return(res)
}

# beach name labels that correspond to the beach extraction order used below
fCIRlab <- c("Huntington","Silver","Washington Park","West","Avalon","Boqueron","Edgewater","Doheny","Fairhope","Goddard","Malibu","Mission Bay","Surfside")


# body immersion results
load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-swim-exposure-regs-body.Rdata")
bCIRs <- rbind(
	CIR.swimex(hu.fit$fit,hu.fit$vcovCL),
	CIR.swimex(si.fit$fit,si.fit$vcovCL),
	CIR.swimex(wp.fit$fit,wp.fit$vcovCL),
	CIR.swimex(we.fit$fit,we.fit$vcovCL),
	CIR.swimex(av.fit$fit,av.fit$vcovCL),
	CIR.swimex(bo.fit$fit,bo.fit$vcovCL),
	CIR.swimex(ed.fit$fit,ed.fit$vcovCL),
	CIR.swimex(dh.fit$fit,dh.fit$vcovCL),
	CIR.swimex(fa.fit$fit,fa.fit$vcovCL),
	CIR.swimex(gd.fit$fit,gd.fit$vcovCL),
	CIR.swimex(ma.fit$fit,ma.fit$vcovCL),
	CIR.swimex(mb.fit$fit,mb.fit$vcovCL),
	CIR.swimex(su.fit$fit,su.fit$vcovCL)
	)
cir.body <- CIRs$CIRoverall[1,1]
	
# head immersion results
load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-swim-exposure-regs-head.Rdata")
hCIRs <- rbind(
	CIR.swimex(hu.fit$fit,hu.fit$vcovCL),
	CIR.swimex(si.fit$fit,si.fit$vcovCL),
	CIR.swimex(wp.fit$fit,wp.fit$vcovCL),
	CIR.swimex(we.fit$fit,we.fit$vcovCL),
	CIR.swimex(av.fit$fit,av.fit$vcovCL),
	CIR.swimex(bo.fit$fit,bo.fit$vcovCL),
	CIR.swimex(ed.fit$fit,ed.fit$vcovCL),
	CIR.swimex(dh.fit$fit,dh.fit$vcovCL),
	CIR.swimex(fa.fit$fit,fa.fit$vcovCL),
	CIR.swimex(gd.fit$fit,gd.fit$vcovCL),
	CIR.swimex(ma.fit$fit,ma.fit$vcovCL),
	CIR.swimex(mb.fit$fit,mb.fit$vcovCL),
	CIR.swimex(su.fit$fit,su.fit$vcovCL)
	)
cir.head <- CIRs$CIRoverall[1,1]
		
# swallowed water results
load("~/dropbox/13beaches/aim1-results/rawoutput/aim1-swim-exposure-regs-swall.Rdata")
sCIRs <- rbind(
	CIR.swimex(hu.fit$fit,hu.fit$vcovCL),
	CIR.swimex(si.fit$fit,si.fit$vcovCL),
	CIR.swimex(wp.fit$fit,wp.fit$vcovCL),
	CIR.swimex(we.fit$fit,we.fit$vcovCL),
	CIR.swimex(av.fit$fit,av.fit$vcovCL),
	CIR.swimex(bo.fit$fit,bo.fit$vcovCL),
	CIR.swimex(ed.fit$fit,ed.fit$vcovCL),
	CIR.swimex(dh.fit$fit,dh.fit$vcovCL),
	CIR.swimex(fa.fit$fit,fa.fit$vcovCL),
	CIR.swimex(gd.fit$fit,gd.fit$vcovCL),
	CIR.swimex(ma.fit$fit,ma.fit$vcovCL),
	CIR.swimex(mb.fit$fit,mb.fit$vcovCL),
	CIR.swimex(su.fit$fit,su.fit$vcovCL)
	)
cir.swall <- CIRs$CIRoverall[1,1]

rownames(bCIRs) <- rownames(hCIRs) <- rownames(sCIRs) <- fCIRlab

# ------------------------------------
# Function to calculate Q and I2
# ------------------------------------

I2calc <- function(studyCIRs,overallCIR) {
	# calculate Cochran's Q, P(Q), and Higgens I2
	# given study-specific CIRs and an overall CIR
	#
	# See Hedges & Vivea 1998, equation (7)
	# and Hidgens 2003
	#
	# studyCIRs : a 3 column matrix with CIR, CIRlb, CIRub as columns
	# overallCIR : the overall CIR across studies (scalar)
	
	# calculate the weights for each study
	# equal to the inverse variance of each estimate
	wi = 1/ ( ( (log(studyCIRs[,3]) - log(studyCIRs[,1]))/1.96 ) ^2 )
	
	# calculate squared deviations for each estimate 
	# from the overall estimate on the log scale
	bdiff2 <- ( log(studyCIRs[,1]) - log(overallCIR)  )^2
	
	# weight the squared differences by the inverse variance
	bwdiff2 <- bdiff2*wi
	
	# calculate the weighted mean squared deviations
	Q <- sum(bwdiff2)
	
	# calculate a P-value using the Chi-squared distribution with k-1 df
	k <- nrow(studyCIRs)
	P <- 1-pchisq(Q,df=k-1)
	
	# calculate I2
	I2 <- 100*( (Q - (k-1))/Q)
	
	cat("\n Cochran's Q = ",Q,"\n prob(Q,",k-1,")   = ",P,"\n I squared   = ",I2,"\n\n",sep="")
	
	return(rbind(Q,P,I2))

}


# ------------------------------------
# run tests for
# body immersion
# head immersion
# swallowed water
# ------------------------------------
bI <- I2calc(bCIRs,cir.body)
hI <- I2calc(hCIRs,cir.head)
sI <- I2calc(sCIRs,cir.swall)



# ------------------------------------
# Enterococcus exposure analyses
# ------------------------------------

# ------------------------------------
# load the Enterococcus EPA 1600
# quartile regression results
# ------------------------------------

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

# ------------------------------------
# run tests for
# Q2 Enterococcus
# Q3 Enterococcus
# Q4 Enterococcus
# ------------------------------------
I2 <- I2calc(q2.cirs,cir.all[1,1])
I3 <- I2calc(q3.cirs,cir.all[2,1])
I4 <- I2calc(q4.cirs,cir.all[3,1])












