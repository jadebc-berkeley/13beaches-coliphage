



# --------------------------------------
# 0-aim2-base-functions.R
# ben arnold (benarnold@berkeley.edu)
#
# description:
# base functions for 13 beaches Aim 2
# (attributable risk calculations)
#
# --------------------------------------



#-------------------------------------------
# Attributable Risk associated with body
# immersion swim exposure (non-swimmers as
# the counterfactual)
#-------------------------------------------

ARswimex <- function(fmla,data) {
	# a function to estimate a marginally adjusted attributable risk for two different parameters
	# a population attributable risk (PAR) comparing the empirical distribution vs. none exposed: P(Y|A,W) - P(Y|A=0,W)
	# and a population attributable fraction (PAF), which is the PAR / P(Y|A,W)
	#
	# in this function, the exposure (A) is the variable "bodycontact" an indicator of body immersion swimming
	# note, that we control for any water contact (because there are some people who may have swallowed water, but not
	# actually gone swimming.
	#
	# See Muller, C. J. & MacLehose, R. F. Estimating predicted probabilities from logistic regression: 
	#     different methods correspond to different target populations. Int J Epidemiol, 2014, 43, 962-970
	#     for a tutorial on marginal standardization for binary outcomes
	# See Greenland, Sander, and Karsten Drescher. Maximum likelihood estimation of the attributable fraction 
	#     from logistic models." Biometrics (1993): 865-872.
	#     for the relevant methods underlying the approach.
	# arguments:
	# fmla : a model formula for glm
	# data : data.frame that includes all of the variables in the formula, and the population used for estimation
	#         note: data must include "bodycontact" and all variables specified in the formula (fmla)
	
	# --------------------------------------
	# Fit log-linear model used to estimate
	# counterfactuals.
	# --------------------------------------
	regfit <- glm(fmla,family=poisson(link="log"),data=data)
	
	# --------------------------------------
	# create counterfactual datasets
	# and predicted P(Y|A,W)
	# --------------------------------------
	# Empirical Probability (adjusted for W)
	pY <- predict(regfit,type="response")
	
	# Counterfactual: everybody is a non-swimmer with no water contact
	cf0 <- data
	cf0$anycontact  <- "No"
	cf0$bodycontact <- "No"
	pY0 <- predict(regfit,newdata=cf0,type="response")
	
	# Counterfactual: all people are body immersion swimmers
	# not currently used to estimate parameters, but stored... for example to estimate RD
	cf1 <- data
	cf1$anycontact  <- "Yes"
	cf1$bodycontact <- "Yes"
	pY1 <- predict(regfit,newdata=cf1,type="response")
	
	# --------------------------------------
	# calculate marginally standardized means
	# for the empirical and counterfactual 
	# distributions of exposure
	# --------------------------------------
	mu  <- mean(pY)
	mu0 <- mean(pY0)
	mu1 <- mean(pY1)
	
	# --------------------------------------
	# estmimate the parameters of interest
	# PAR: P(Y|A,  W) - P(Y|A=0,W)
	# PAF: 100% * [P(Y|A,  W) - P(Y|A=0,W)] / P(Y|A,  W)
	# --------------------------------------
	thetaPAR <- mean(pY-pY0)
	thetaPAF <- 100*( thetaPAR/mean(pY) )

	# --------------------------------------
	# return marginal means and 
	# parameters of interest as a list
	# --------------------------------------
	list(thetaPAR=thetaPAR,thetaPAF=thetaPAF,mu=mu,mu0=mu0,mu1=mu1)
}


#-------------------------------------------
# Attributable Risk among swimmers
# associated with exposure
# to water with Entero EPA 1600 >35 CFU/100ml
#-------------------------------------------

AR35noswim <- function(fmla,data) {
	# a function to estimate a marginally adjusted attributable risk for two different parameters
	# a population attributable risk (PAR) comparing the empirical distribution vs. none exposed: P(Y|A,W) - P(Y|A=0,W)
	# and a population attributable fraction (PAF), which is the PAR / P(Y|A,W)
	#
	# in this function, the exposure (A) is the variable "swim35" a 3-level factor of non-swimmers, swimmers exposed Enterococcus EPA 1600 <= 35 CFU/100ml, and swimmers exposed to >35 CFU
	#
	# See Muller, C. J. & MacLehose, R. F. Estimating predicted probabilities from logistic regression: 
	#     different methods correspond to different target populations. Int J Epidemiol, 2014, 43, 962-970
	#     for a tutorial on marginal standardization for binary outcomes
	# See Greenland, Sander, and Karsten Drescher. Maximum likelihood estimation of the attributable fraction 
	#     from logistic models." Biometrics (1993): 865-872.
	#     for the relevant methods underlying the approach.
	# arguments:
	# fmla : a model formula for glm
	# data : data.frame that includes all of the variables in the formula, and the population used for estimation
	#        note: data must include "swim35"
	
	# --------------------------------------
	# Fit log-linear model used to estimate
	# counterfactuals.
	# --------------------------------------
	regfit <- glm(fmla,family=poisson(link="log"),data=data)
	
	# --------------------------------------
	# create counterfactual datasets
	# and predicted P(Y|A,W)
	# --------------------------------------
	# Empirical Probability (adjusted for W)
	pY <- predict(regfit,type="response")
	
	# Counterfactual: swimmers exposed above 35 CFU are prevented from swimming
	cf0 <- data
	cf0$swim35[cf0$swim35=="Above35cfu"] <- "Non-swimmers"
	pY0 <- predict(regfit,newdata=cf0,type="response")

	# Counterfactual: swimmers exposed above 35 CFU swim <35 CFU
	# not currently used to estimate parameters, but stored... for example compare w/ estimates from AR35cfu (below)
	cf1 <- data
	cf1$swim35[cf1$swim35=="Above35cfu"] <- "Below35cfu"
	pY1 <- predict(regfit,newdata=cf1,type="response")
	
	# --------------------------------------
	# calculate marginally standardized means
	# for the empirical and counterfactual 
	# distributions of exposure
	# --------------------------------------
	mu  <- mean(pY)
	mu0 <- mean(pY0)
	mu1 <- mean(pY1)
	
	# --------------------------------------
	# estmimate the parameters of interest
	# PAR: P(Y|A,  W) - P(Y|A=0,W)
	# PAF: 100% * [P(Y|A,  W) - P(Y|A=0,W)] / P(Y|A,  W)
	# --------------------------------------
	thetaPAR <- mean(pY-pY0)
	thetaPAF <- 100*( thetaPAR/mean(pY) )

	# --------------------------------------
	# return marginal means and 
	# parameters of interest as a list
	# --------------------------------------
	list(thetaPAR=thetaPAR,thetaPAF=thetaPAF,mu=mu,mu0=mu0,mu1=mu1)
}

#-------------------------------------------
# Attributable Risk among swimmers
# associated with exposure
# to water with Entero EPA 1600 >35 CFU/100ml
#-------------------------------------------

AR35cfu <- function(fmla,data) {
	# a function to estimate a marginally adjusted attributable risk for two different parameters
	# a population attributable risk (PAR) comparing the empirical distribution vs. none exposed: P(Y|A,W) - P(Y|A=0,W)
	# and a population attributable fraction (PAF), which is the PAR / P(Y|A,W)
	#
	# in this function, the exposure (A) is the variable "entero35" an indicator of Enterococcus EPA 1600 > 35 CFU/100ml 
	#
	# See Muller, C. J. & MacLehose, R. F. Estimating predicted probabilities from logistic regression: 
	#     different methods correspond to different target populations. Int J Epidemiol, 2014, 43, 962-970
	#     for a tutorial on marginal standardization for binary outcomes
	# See Greenland, Sander, and Karsten Drescher. Maximum likelihood estimation of the attributable fraction 
	#     from logistic models." Biometrics (1993): 865-872.
	#     for the relevant methods underlying the approach.
	# arguments:
	# fmla : a model formula for glm
	# data : data.frame that includes all of the variables in the formula, and the population used for estimation
	#         note: data must include "entero35"
	
	# --------------------------------------
	# Fit log-linear model used to estimate
	# counterfactuals.
	# --------------------------------------
	regfit <- glm(fmla,family=poisson(link="log"),data=data)
	
	# --------------------------------------
	# create counterfactual datasets
	# and predicted P(Y|A,W)
	# --------------------------------------
	# Empirical Probability (adjusted for W)
	pY <- predict(regfit,type="response")
	
	# Counterfactual: all swimmers exposed to Entero 1600 <= 35 CFU/100ml
	cf0 <- data
	cf0$entero35 <- 0
	pY0 <- predict(regfit,newdata=cf0,type="response")
	
	# Counterfactual: all swimmers exposed to Entero 1600 > 35 CFU/100ml
	# not currently used to estimate parameters, but stored... for example to estimate RD
	cf1 <- data
	cf1$entero35 <- 1
	pY1 <- predict(regfit,newdata=cf1,type="response")
	
	# --------------------------------------
	# calculate marginally standardized means
	# for the empirical and counterfactual 
	# distributions of exposure
	# --------------------------------------
	mu  <- mean(pY)
	mu0 <- mean(pY0)
	mu1 <- mean(pY1)
	
	# --------------------------------------
	# estmimate the parameters of interest
	# PAR: P(Y|A,  W) - P(Y|A=0,W)
	# PAF: 100% * [P(Y|A,  W) - P(Y|A=0,W)] / P(Y|A,  W)
	# --------------------------------------
	thetaPAR <- mean(pY-pY0)
	thetaPAF <- 100*( thetaPAR/mean(pY) )

	# --------------------------------------
	# return marginal means and 
	# parameters of interest as a list
	# --------------------------------------
	list(thetaPAR=thetaPAR,thetaPAF=thetaPAF,mu=mu,mu0=mu0,mu1=mu1)
}



#-------------------------------------------
# Stratified bootstrap resampling function
# Samples n_i households from stratum i,
# where stratum in this study is beach cohort
# 
# This function is tailored to store multiple
# statistics generated by the ARswimex and AR35cfu 
# functions above.  
#-------------------------------------------

bootAR <- function(fn,fmla,data,ID,strata,iter,dots=TRUE) {
	# fn    : an Attributable risk function to call, defined above (ARswimex, AR35cfu)
	# fmla  : a model formula for glm, which is passed to fn
	# data  : data.frame that includes all of the variables in the formula, and the population used for estimation
	#         which is passed to fn.
	# ID    : ID for resampling (e.g., household ID required b/c of repeated, potentially correlated obs w/in HH)
	# strata: stratifying variable (e.g., beach ID)
	# iter  : number of bootstrap iterations	
	# dots  : logical. print dots for iterations
	if(length(ID)!=nrow(data) | length(strata)!=nrow(data)) {
		print("data, ID, and strata args must be the same length")
		break
	}
	
	# add ID variables to data frame for convenience (+ removes empty beach strata for some subgroup analyses)
	data$ID <- ID
	data$strata <- factor(strata)
	
	if(dots) cat("\n\n\nBootstrap Iterations (",iter,") \n----|--- 1 ---|--- 2 ---|--- 3 ---|--- 4 ---| --- 5 \n",sep="") 
	start.time <- Sys.time()
	
	# create an empty matrix to store bootstrap results
	# corresponding to the PAR and PAF parameters of interest + 3 marginal means
	stats <- matrix(NA,ncol=5,nrow=iter)
	colnames(stats) <- c("thetaPAR","thetaPAF","PY","PY0","PY1")
	
	# Create the bootstrap samples of size n_i, within each stratum i
	IDlist <- as.list(tapply(data$ID,data$strata,function(x) unique(x)))
	bs <- matrix(NA,nrow=length(unique(data$ID)),ncol=iter)
	k <- 0
	for(i in 1:dim(IDlist)) {
		Nids <- length(IDlist[[i]])
		j <- k+1
		k <- j-1+Nids
		bs[j:k,] <- sample(IDlist[[i]],Nids*iter,replace=TRUE)
	}
	
	# calculate statistics for each bootstrap sample
	# call the appropriate AR function (specified in the fn argument), depending on the exposure of interest
	for (bb in 1:iter) {
		bd <- merge(data.frame(ID=bs[,bb]),data,by="ID",all.x=TRUE)
		stats[bb,] <- unlist( tryCatch( do.call(fn,args=list(fmla=fmla,data=bd)), error=function(e) rep(NA,5) )  )
		if(dots) cat(".",sep="")
    	if(dots && bb %% 50 == 0) cat(bb,"\n")
	}
	if(dots) cat("\n Bootstrap Run Time:",round(difftime(Sys.time(),start.time,units="mins"),3)," Minutes \n\n")
	
	# if there were any bootstrap samples where the model failed to converge 
	# for any reason, caught by tryCatch(), then report that.
	# also triage any estimates that are clearly from a failed model that sort-of converged
	# (due to model instability from sparse data in a bootstrap sample)
	stats[(stats[,1]>10 | stats[,1]< -10 | stats[,3]>10 | stats[,3]< -10), ] <- NA
	Nna <- sum(ifelse(is.na(stats[,1])==TRUE,1,0))
	if(Nna>0) {
		cat("\n   Warning: ",Nna," bootstrap replicates failed to converge to sensible estimates \n")
		cat("\n   Bootstrap estimates are based on the remaining ",Nboot-Nna," replicates\n\n")
	}
	
	# calculate the mean and percentile 95% confidence intervals for the statistics
	bootest  <- apply(stats,2,function(x) mean(x,na.rm=TRUE))
	boot95lb <- apply(stats,2,function(x) quantile(x,prob=0.025,na.rm=TRUE))
	boot95ub <- apply(stats,2,function(x) quantile(x,prob=0.975,na.rm=TRUE))
	
	
	# print the results
	cat("\n Bootstrap Estimates:\n\n Means    :", paste(sprintf("%1.4f",bootest)),"\n Lower 95%:", paste(sprintf("%1.4f",boot95lb)),"\n Upper 95%:",paste(sprintf("%1.4f",boot95ub)),"\n")
	
	list(bootest=bootest,boot95lb=boot95lb,boot95ub=boot95ub,stats=stats)
}

