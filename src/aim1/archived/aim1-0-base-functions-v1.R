

# --------------------------------------
# aim1-0-base-functions.R
# ben arnold (benarnold@berkeley.edu)
#
# description:
# base functions used in many Aim 1 analyses
#
#
# version 1 (9 feb 2015)
#
# --------------------------------------

# --------------------------------------
# Robust clustered SE function
# http://people.su.se/~ma/mcluster.R
# R (www.r-project.org) codes for computing multi-way clustered-standard errors
# Mahmood Arai, Jan 21, 2008. 
# See: Thompson (2006), Cameron, Gelbach and Miller (2006) and Petersen (2006).
#
# slightly modified to have it return the vcovCL object
# rather than the updated fit (since need the VC matrix)
# --------------------------------------
cl   <- function(dat,fm, cluster){
	# dat: data used to fit the model
	# fm : model fit (object)
	# cluster : vector of cluster IDs
	require(sandwich, quietly = TRUE)
	require(lmtest, quietly = TRUE)
	M <- length(unique(cluster))
	N <- length(cluster)
	K <- fm$rank
	dfc <- (M/(M-1))*((N-1)/(N-K))
	uj  <- apply(estfun(fm),2, function(x) tapply(x, cluster, sum));
	vcovCL <- dfc*sandwich(fm, meat=crossprod(uj)/N)
	return(vcovCL)
}

# --------------------------------------
# Convenience Function to run log-binomial models
# and obtain robust SEs (clusterd on hhid)
# for the beaches data
# --------------------------------------

lbreg <- function(formula,dat,beach,vcv=FALSE) {
	# log-binomial regression formula
	# dataset used to fit the model
	# beach name (character)	
	fit <- glm(formula,family=binomial(link="log"),data=dat[dat$beach==beach,])
	vcovCL <- cl(dat[dat$beach==beach,],fm=fit,cluster=dat$hhid[dat$beach==beach])
	rfit <- coeftest(fit, vcovCL) 
	if(vcv==FALSE) {
		return(rfit)
	} else {
		list(fit=rfit,vcovCL=vcovCL)
	}
}

# --------------------------------------
# Convenience Function to run modified Poisson models
# and obtain robust SEs (clusterd on hhid)
# for the beaches data
# --------------------------------------

mpreg <- function(formula,dat,beach,vcv=FALSE) {
	# modified Poisson regression formula
	# dataset used to fit the model
	# beach name (character)	
	fit <- glm(formula,family=poisson(link="log"),data=dat[dat$beach==beach,])
	vcovCL <- cl(dat[dat$beach==beach,],fm=fit,cluster=dat$hhid[dat$beach==beach])
	rfit <- coeftest(fit, vcovCL) 
	if(vcv==FALSE) {
		return(rfit)
	} else {
		list(fit=rfit,vcovCL=vcovCL)
	}
}



# --------------------------------------
# convenience functions to obtain
# point estimates and SEs from
# model fits for plotting
# --------------------------------------

# swim exposure regs w/o interactions
# (need to sum 2 coefficients: anycontact + higher level of contact)
swestci <- function(fo,vcv) {
	# fo : fit object returned from lbreg or mpreg
	# vcv: variance covariance matrix from lbreg or mpreg
	nr <- nrow(vcv)
	lc <- c(0,1,1,rep(0,nr-3)) # linear combination of betas for the estimate
	est <- exp(t(lc)%*%fo[,1])
	se <- sqrt( t(lc)%*%vcv%*%lc )
	lb <- exp(log(est)-1.96*se)
	ub <- exp(log(est)+1.96*se)
	list(est=est,lb=lb,ub=ub)
}


# swim exposure regs with binary interactions
iswestci <- function(fo,vcv) {
	# fo : fit object returned from lbreg or mpreg
	# vcv: variance covariance matrix from lbreg or mpreg
	# Important: code assumes that the interaction terms are 
	#            at the end of the list of covariates
	nr <- nrow(vcv)
	lc0 <- lc1 <- rep(0,nr)
	lc0[c(2,3)] <- 1
	lc1[c(2,3,nr-1,nr)] <- 1
	
	est0 <- exp(t(lc0)%*%fo[,1])
	se0  <- sqrt( t(lc0)%*%vcv%*%lc0 )
	lb0 <- exp(log(est0)-1.96*se0)
	ub0 <- exp(log(est0)+1.96*se0)
	
	est1 <- exp(t(lc1)%*%fo[,1])
	se1  <- sqrt( t(lc1)%*%vcv%*%lc1 )
	lb1 <- exp(log(est1)-1.96*se1)
	ub1 <- exp(log(est1)+1.96*se1)
	
	list(est0=est0,lb0=lb0,ub0=ub0,est1=est1,lb1=lb1,ub1=ub1)
}

# swim exposure regs with 3-level interactions
i3swestci <- function(fo,vcv) {
	# fo : fit object returned from lbreg or mpreg
	# vcv: variance covariance matrix from lbreg or mpreg
	# Important: code assumes that the 4 interaction terms are 
	#            at the end of the list of covariates
	nr <- nrow(vcv)
	lc0 <- lc1 <- lc2 <- rep(0,nr)
	lc0[c(2,3)] <- 1
	lc1[c(2,3,nr-3,nr-1)] <- 1
	lc2[c(2,3,nr-2,nr)] <- 1
	
	est0 <- exp(t(lc0)%*%fo[,1])
	se0  <- sqrt( t(lc0)%*%vcv%*%lc0 )
	lb0 <- exp(log(est0)-1.96*se0)
	ub0 <- exp(log(est0)+1.96*se0)
	
	est1 <- exp(t(lc1)%*%fo[,1])
	se1  <- sqrt( t(lc1)%*%vcv%*%lc1 )
	lb1 <- exp(log(est1)-1.96*se1)
	ub1 <- exp(log(est1)+1.96*se1)
	
	est2 <- exp(t(lc2)%*%fo[,1])
	se2  <- sqrt( t(lc2)%*%vcv%*%lc2 )
	lb2 <- exp(log(est2)-1.96*se2)
	ub2 <- exp(log(est2)+1.96*se2)
	
	list(est0=est0,lb0=lb0,ub0=ub0,est1=est1,lb1=lb1,ub1=ub1,est2=est2,lb2=lb2,ub2=ub2)
}



# swim exposure test of interaction
switest <- function(fo,vcv) {
	# fo : fit object returned from lbreg or mpreg
	# vcv: variance covariance matrix from lbreg or mpreg
	# this is needed because there are 2 separate exposures
	# in these regressions: anycontact and then some higher level (e.g., body immersion)
	# Important: code assumes that the 2 interaction terms are 
	#            at the end of the list of covariates
	# Calculate Z statistic and P-values for the interaction terms
	nr <- nrow(vcv)
	lc <- rep(0,nr)
	lc[c(nr-1,nr)] <- 1
	est <- lc%*%fo[,1]
	se  <- sqrt( t(lc)%*%vcv%*%lc )
	Z <- est/se
	P <- 2*pnorm(-abs(Z))
	
	cat("\nTest of interaction:\n","Z=",sprintf("%1.3f",Z),", P=",sprintf("%1.3f",P),"\n\n")
	list(Z=Z,P=P)
}


# CI calculation for indicator regs w/o interactions
estci <- function(fo) {
	# fo : fit object returned from lbreg or mpreg
	est <- exp(fo[2,1])
	lb <- exp(fo[2,1]-1.96*fo[2,2])
	ub <- exp(fo[2,1]+1.96*fo[2,2])
	list(est=est,lb=lb,ub=ub)
}

# CI calculation for indicator regs w/ interactions (Avalon, Doheny)
iestci <- function(fo,vcv) {
	# fo : fit object returned from lbreg or mpreg
	# vcv: variance covariance matrix from lbreg or mpreg
	est0 <- exp(fo[2,1])
	lb0  <- exp(fo[2,1]-1.96*fo[2,2])
	ub0  <- exp(fo[2,1]+1.96*fo[2,2])
	
	est1 <- exp(fo[2,1]+fo[4,1])
	se1  <- sqrt(vcv[2,2]+vcv[4,4]+2*vcv[2,4])
	lb1  <- exp(fo[2,1]+fo[4,1] - 1.96*se1)
	ub1  <- exp(fo[2,1]+fo[4,1] + 1.96*se1)
	
	res <- matrix(c(est0,lb0,ub0,est1,lb1,ub1),nrow=2,ncol=3,byrow=T)
	rownames(res) <- c("effmod=0","effmod=1")
	colnames(res) <- c("est","lb","ub")
	
	return(res)
}



# --------------------------------------
# Swimming exposure regressions:
# Make matrices of CIRs and 95% CIs
# from the regression model output
# note: calls swestci and iswestci from
# above. cirtab calculates
# CIRs and their CIs using linear combinations
# of the regression coefficients + VC matrix
# --------------------------------------

cirtab <- function(b,h,s,I=FALSE) {
	# b : body immersion results object from mpreg
	# h : head immersion results object from mpreg
	# s : swallowed water results object from mpreg
	# I : logical: interaction between swim cond + env (Avalon, Doheny)
	if(I==FALSE) {
		M <- rbind(
		unlist(swestci(b$fit,b$vcovCL)),  
		unlist(swestci(h$fit,h$vcovCL)),
		unlist(swestci(s$fit,s$vcovCL))	
		)
	} else {
		M <- rbind(
		unlist(iswestci(b$fit,b$vcovCL)),  
		unlist(iswestci(h$fit,h$vcovCL)),
		unlist(iswestci(s$fit,s$vcovCL))	
		)
	}
	rownames(M) <- c("Body Immersion","Head Immersion","Swallowed Water")
	return(M)
}



