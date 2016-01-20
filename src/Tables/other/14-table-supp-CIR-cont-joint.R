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
getCIR <- function(x,coliphage) {
  # x : log-linear model object returned from coeftest (class=coeftest)
  # NOTE: assumes exposure of interest is the first covariate and there are no interactions
  row.coli=grep(coliphage,rownames(x))[1]
  row.entero=grep("entero",rownames(x))[1]
  row.joint=grep(":",rownames(x))
  est <- exp(x[row.coli,1]+x[row.entero,1]+x[row.joint,1])
  se  <- x[2,2]  
  lb <- exp(log(est)-1.96*se)
  ub <- exp(log(est)+1.96*se)
  res <- c(est,lb,ub)
  return(res)
}

# print output 
all=rbind(CIRformat(getCIR(overall.fit10.fmc1601,"fmc1601")),
  CIRformat(getCIR(overall.fit10.fmc1602,"fmc1602")),
  CIRformat(getCIR(overall.fit10.fpc1601,"fpc1601")),
  CIRformat(getCIR(overall.fit10.fpc1602,"fpc1602")))

low=rbind(NA,CIRformat(getCIR(overall.fit10.fmc1602.low,"fmc1602")),
  CIRformat(getCIR(overall.fit10.fpc1601.low,"fpc1601")),
  CIRformat(getCIR(overall.fit10.fpc1602.low,"fpc1602")))

high=rbind(NA,CIRformat(getCIR(overall.fit10.fmc1602.high,"fmc1602")),
  CIRformat(getCIR(overall.fit10.fpc1601.high,"fpc1601")),
  CIRformat(getCIR(overall.fit10.fpc1602.high,"fpc1602")))

tab=cbind(all, low, high)
label=c("Somatic coliphage 1601","Somatic coliphage 1602",
        "Male-specific coliphage 1601","Male-specific coliphage 1602")
tab.out=cbind(label,tab)

save(tab.out,file="~/Documents/CRG/coliphage/Results/Tables/CIR-10-cont-joint.RData")



