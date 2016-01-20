##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# Plots of the RR across the range of concentration

# Results pooled across beach
##########################################


rm(list=ls())

load("~/Documents/CRG/coliphage/results/rawoutput/regress-10day-continuous-body-joint.Rdata")


# --------------------------------------
# formatting function for CIR results
# --------------------------------------
CIRformat <- function(x) {
  # x : vector of length 3 with CIR, CIRlb, Cub
  paste(sprintf("%1.2f",x[1])," (",sprintf("%1.2f",x[2]),", ",sprintf("%1.2f",x[3]),")",sep="")
} 


# --------------------------------------
# get CIs in one object for graphing dose
# response figures
# --------------------------------------
# function to get CIR Estimates and CIs from simple stratified models
getCIR <- function(x) {
  # x : log-linear model object returned from coeftest (class=coeftest)
  # NOTE: assumes exposure of interest is the first covariate and there are no interactions
  est <- exp(x[2,1]+x[3,1]+x[23,1])
  se  <- x[2,2]  
  lb <- exp(log(est)-1.96*se)
  ub <- exp(log(est)+1.96*se)
  res <- c(est,lb,ub)
  return(res)
}

# print output 
CIRformat(getCIR(overall.fit10.fmc1601))
CIRformat(getCIR(overall.fit10.fmc1602))
CIRformat(getCIR(overall.fit10.fpc1601))
CIRformat(getCIR(overall.fit10.fpc1602))

CIRformat(getCIR(overall.fit10.fmc1602.low))
CIRformat(getCIR(overall.fit10.fmc1602.high))
CIRformat(getCIR(overall.fit10.fpc1601.low))
CIRformat(getCIR(overall.fit10.fpc1601.high))
CIRformat(getCIR(overall.fit10.fpc1602.low))
CIRformat(getCIR(overall.fit10.fpc1602.high))

