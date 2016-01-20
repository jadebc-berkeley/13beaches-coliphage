##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# Comparison of CIRs for continuous 
# coliphage indicators for main analysis
# and negative controls

##########################################


rm(list=ls())
load("~/Documents/CRG/coliphage/results/rawoutput/regress-10day-continuous-body.Rdata")
load("~/Documents/CRG/coliphage/results/rawoutput/regress-10day-continuous-body-negcontrol.Rdata")
load("~/Documents/CRG/coliphage/results/rawoutput/regress-10day-continuous-body-entero.Rdata")
load("~/Documents/CRG/coliphage/results/rawoutput/regress-10day-continuous-body-entero-negcontrol.Rdata")

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
  est <- exp(x[2,1])
  se  <- x[2,2]  
  lb <- exp(log(est)-1.96*se)
  ub <- exp(log(est)+1.96*se)
  res <- c(est,lb,ub)
  return(res)
}

# --------------------------------------
# compile estimates
# --------------------------------------

# coliphage main 
all=list(overall.fit10.fmc1601, overall.fit10.fmc1602, 
  overall.fit10.fpc1601,overall.fit10.fpc1602)

low=list(overall.fit10.fmc1602.low,overall.fit10.fpc1601.low,overall.fit10.fpc1602.low)

high=list(overall.fit10.fmc1602.high,overall.fit10.fpc1601.high,overall.fit10.fpc1602.high)

# entero main
all.entero=list(overall.fit10.entero.fmc1601,overall.fit10.entero.fmc1602,
  overall.fit10.entero.fpc1601,overall.fit10.entero.fpc1602)

high.entero=list(overall.fit10.entero.high.fmc1602,overall.fit10.entero.high.fpc1601,
  overall.fit10.entero.high.fpc1602)

low.entero=list(overall.fit10.entero.low.fmc1602,overall.fit10.entero.low.fpc1601,
  overall.fit10.entero.low.fpc1602)

# coliphage negative control
nc.all=list(nc.overall.fit10.fmc1601, nc.overall.fit10.fmc1602, 
  nc.overall.fit10.fpc1601,nc.overall.fit10.fpc1602)

nc.low=list(nc.overall.fit10.fmc1602.low,nc.overall.fit10.fpc1601.low,nc.overall.fit10.fpc1602.low)

nc.high=list(nc.overall.fit10.fmc1602.high,nc.overall.fit10.fpc1601.high,nc.overall.fit10.fpc1602.high)

# entero negative control
nc.all.entero=list(nc.overall.fit10.entero.fmc1601,nc.overall.fit10.entero.fmc1602,
  nc.overall.fit10.entero.fpc1601,nc.overall.fit10.entero.fpc1602)

nc.high.entero=list(nc.overall.fit10.entero.high.fmc1602,nc.overall.fit10.entero.high.fpc1601,
  nc.overall.fit10.entero.high.fpc1602)

nc.low.entero=list(nc.overall.fit10.entero.low.fmc1602,nc.overall.fit10.entero.low.fpc1601,
  nc.overall.fit10.entero.low.fpc1602)

# --------------------------------------
# format estimates
# --------------------------------------
est.all=unlist(lapply(all,function(x) CIRformat(getCIR(x))))
est.low=c(NA,unlist(lapply(low,function(x) CIRformat(getCIR(x)))))
est.high=c(NA,unlist(lapply(high,function(x) CIRformat(getCIR(x)))))

est.all.entero=unlist(lapply(all.entero,function(x) CIRformat(getCIR(x))))
est.low.entero=c(NA,unlist(lapply(low.entero,function(x) CIRformat(getCIR(x)))))
est.high.entero=c(NA,unlist(lapply(high.entero,function(x) CIRformat(getCIR(x)))))

est.nc.all=unlist(lapply(nc.all,function(x) CIRformat(getCIR(x))))
est.nc.low=c(NA,unlist(lapply(nc.low,function(x) CIRformat(getCIR(x)))))
est.nc.high=c(NA,unlist(lapply(nc.high,function(x) CIRformat(getCIR(x)))))

est.nc.all.entero=unlist(lapply(nc.all.entero,function(x) CIRformat(getCIR(x))))
est.nc.low.entero=c(NA,unlist(lapply(nc.low.entero,function(x) CIRformat(getCIR(x)))))
est.nc.high.entero=c(NA,unlist(lapply(nc.high.entero,function(x) CIRformat(getCIR(x)))))


# --------------------------------------
# make table
# --------------------------------------
col.tab=cbind(est.all,est.low,est.high)
col.tab.nc=cbind(est.nc.all,est.nc.low,est.nc.high)

ent.tab=cbind(est.all.entero,est.low.entero,est.high.entero)
ent.tab.nc=cbind(est.nc.all.entero,est.nc.low.entero,est.nc.high.entero)

colnames(col.tab)=c("all","low","high")
colnames(col.tab.nc)=c("all","low","high")
colnames(ent.tab)=c("all","low","high")
colnames(ent.tab.nc)=c("all","low","high")

cont.tab=rbind(col.tab[1,],ent.tab[1,],col.tab[2,],ent.tab[2,],
               col.tab[3,],ent.tab[3,],col.tab[4,],ent.tab[4,])

cont.tab.nc=rbind(col.tab.nc[1,],ent.tab.nc[1,],col.tab.nc[2,],ent.tab.nc[2,],
               col.tab.nc[3,],ent.tab.nc[3,],col.tab.nc[4,],ent.tab.nc[4,])

cont.tab.out=rbind(cont.tab,cont.tab.nc)

label=rep(c("Coliphage","Enterococcus"),8)
analysis=c(rep("Main",8),rep("Neg control",8))

cont.tab.out=cbind(label,analysis,cont.tab.out)

save(cont.tab.out,file="~/Documents/CRG/coliphage/results/tables/table-cont-negcontrol.RData")


