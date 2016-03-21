##########################################
# Coliphage analysis - 6 beaches

# Table of continuous analysis results
##########################################


rm(list=ls())

load("~/Documents/CRG/coliphage/results/rawoutput/regress-10day-body.Rdata")
load("~/Documents/CRG/coliphage/results/rawoutput/regress-10day-continuous-body-entero-n.Rdata")
load("~/Documents/CRG/coliphage/results/rawoutput/regress-10day-continuous-body-n.Rdata")
load("~/Documents/CRG/coliphage/results/rawoutput/regress-10day-continuous-body.Rdata")
load("~/Documents/CRG/coliphage/results/rawoutput/regress-10day-continuous-CIR-body.Rdata")
load("~/Documents/CRG/coliphage/results/rawoutput/regress-10day-continuous-body-entero.Rdata")
load("~/Documents/CRG/coliphage/results/rawoutput/regress-10day-continuous-CIR-body-entero.Rdata")

# --------------------------------------
# combine n's
# --------------------------------------

coli.n=c(all.n10.fmc1601.low,all.n10.fmc1601.high,
         all.n10.fmc1602.low,all.n10.fmc1602.high,
         all.n10.fpc1601.low,all.n10.fpc1601.high,
         all.n10.fpc1602.low,all.n10.fpc1602.high)

ent.n=c(all.n10.ent.fmc1601.low,all.n10.ent.fmc1601.high,
        all.n10.ent.fmc1602.low, all.n10.ent.fmc1602.high, 
        all.n10.ent.fpc1601.low, all.n10.ent.fpc1601.high,
        all.n10.ent.fpc1602.low, all.n10.ent.fpc1602.high)

# check that coli and ent n's are the same
coli.n==ent.n

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

# entero
ent.fmc1601.low=CIRformat(getCIR(overall.fit10.entero.low.fmc1601))
ent.fmc1601.high=CIRformat(getCIR(overall.fit10.entero.high.fmc1601))
ent.fmc1602.low=CIRformat(getCIR(overall.fit10.entero.low.fmc1602))
ent.fmc1602.high=CIRformat(getCIR(overall.fit10.entero.high.fmc1602))
ent.fpc1601.low=CIRformat(getCIR(overall.fit10.entero.low.fpc1601))
ent.fpc1601.high=CIRformat(getCIR(overall.fit10.entero.high.fpc1601))
ent.fpc1602.low=CIRformat(getCIR(overall.fit10.entero.low.fpc1602))
ent.fpc1602.high=CIRformat(getCIR(overall.fit10.entero.high.fpc1602))


# coliphage
fmc1601.low=CIRformat(getCIR(overall.fit10.fmc1601.low))
fmc1601.high=CIRformat(getCIR(overall.fit10.fmc1601.high))
fmc1602.low=CIRformat(getCIR(overall.fit10.fmc1602.low))
fmc1602.high=CIRformat(getCIR(overall.fit10.fmc1602.high))
fpc1601.low=CIRformat(getCIR(overall.fit10.fpc1601.low))
fpc1601.high=CIRformat(getCIR(overall.fit10.fpc1601.high))
fpc1602.low=CIRformat(getCIR(overall.fit10.fpc1602.low))
fpc1602.high=CIRformat(getCIR(overall.fit10.fpc1602.high))

# make table rows
fmc1601.row.low=c(fmc1601.low,ent.fmc1601.low)
fmc1601.row.high=c(fmc1601.high,ent.fmc1601.high)
fmc1602.row.low=c(fmc1602.low,ent.fmc1602.low)
fmc1602.row.high=c(fmc1602.high,ent.fmc1602.high)
fpc1601.row.low=c(fpc1601.low,ent.fpc1601.low)
fpc1601.row.high=c(fpc1601.high,ent.fpc1601.high)
fpc1602.row.low=c(fpc1602.low,ent.fpc1602.low)
fpc1602.row.high=c(fpc1602.high,ent.fpc1602.high)

# make table
table=data.frame(rbind(fmc1601.row.low,fmc1602.row.high,
                       fmc1602.row.low,fmc1602.row.high,
                       fpc1601.row.low,fpc1601.row.high,
                       fpc1602.row.low,fpc1602.row.high))

lab1=c("EPA Method 1601","EPA Method 1601","EPA Method 1602","EPA Method 1602",
       "EPA Method 1601","EPA Method 1601", "EPA Method 1602","EPA Method 1602")
lab2=c("Low-risk conditions","High-risk conditions",
       "Low-risk conditions","High-risk conditions",
       "Low-risk conditions","High-risk conditions",
       "Low-risk conditions","High-risk conditions")
nbeach=c(0,0,3,2,4,4,2,2)
table=cbind(lab1,lab2,nbeach,coli.n,table)

save(table,file="~/Documents/CRG/coliphage/Results/Tables/CIR-10-continuous.RData")
write.csv(table,file="~/Documents/CRG/coliphage/Results/Tables/CIR-10-continuous.csv",na="",row.names=FALSE)
