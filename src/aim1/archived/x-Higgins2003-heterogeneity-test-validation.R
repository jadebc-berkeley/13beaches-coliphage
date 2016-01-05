

#--------------------------------------
# x-Higgins2003-heterogeneity-test-validation.R
#
# Use data from Fig 1 of
# Measuring inconsistency in meta-analyses
# Higgins et al. BMJ 2003;327:557–60
#
# to ensure that we are calculating 
# Cochrane's Q and Higgins' I^2 correctly
#
#--------------------------------------




#--------------------------------------
# input the data
#--------------------------------------
d <- matrix(c(
	16, 141, 41, 152,
	 1,  53,  8,  52, 
	 8, 136, 28, 139,
	 9,  59,  9,  51,
	32,  95, 59,  97,
	15, 107, 20,  99, 
	 2, 113, 27, 132,
	 3, 317,  5, 159),
	 nrow=8,ncol=4,byrow=T
)
colnames(d) <- c("Tn","TN","Cn","CN")
d <- data.frame(d)

# calculate study-specific ORs and 95% CIs
d$pT <- d$Tn/d$TN
d$pC <- d$Cn/d$CN
d$or <- (d$pT/(1-d$pT))/(d$pC/(1-d$pC))

d$orse <- sqrt(1/(d$Tn) + 1/(d$TN-d$Tn) + 1/(d$Cn) + 1/(d$CN-d$Cn) )
d$orlb <- exp(log(d$or) - 1.96*d$orse)
d$orub <- exp(log(d$or) + 1.96*d$orse)

# load the metafor package and confirm OR calcs are correct
library(metafor)
fit <- escalc("OR",ai=d$Tn,ci=d$Cn,n1i=d$TN,n2i=d$CN)
summary(fit)
cbind(exp(fit$yi),d$or)
cbind(summary(fit)[,3],d$orse)

# calculate Q and I^2
# should be Q(df=7) = 12.2, P=0.09, I2=44%
Qfit <- rma(yi=log(d$or),sei=d$orse)
Qfit


# test the home-made function
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

ORpooled <- exp(Qfit$b)

Qfit2 <- I2calc(subset(d,select=c("or","orlb","orub")),ORpooled)


