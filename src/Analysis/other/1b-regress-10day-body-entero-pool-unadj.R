##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file conducts maximum likelihood regression
# to estimate prevalence ratios
# for enterococcus

# Results pooled across beaches and assay

# 10 day gi illness
##########################################

rm(list=ls())
library(foreign)

# --------------------------------------
# load the and pre-preprocess the 
# analysis dataset
# (refer to the base functions script
# for details on the pre-processing)
# --------------------------------------
beaches13=read.csv("~/Documents/CRG/coliphage/13beaches-data/final/13beaches-analysis.csv")

# load base functions
source("~/Documents/CRG/coliphage/13beaches-coliphage/src/Analysis/0-base-functions.R")

data=preprocess.6beaches(beaches13)

# restrict to 6 beaches with coliphage data
beach.list=c("Avalon","Doheny","Malibu","Mission Bay",
             "Fairhope","Goddard")

all=data[data$beach %in% beach.list,]

# drop individuals with no water quality information
all=subset(all,nowq==0)
# subset to non-missing exposure categories
# to make the robust CI calcs work
all=subset(all,all$bodycontact=="Yes")

# --------------------------------------
# Calculate the actual Ns for each cell
# and store them for plotting and tables
# --------------------------------------
regN <- function(outcome,exposurecat) {
  sum(table(outcome,exposurecat))
}

# n's pooled by assay and beach ---------------------------------------
n10.entero35.fmc = regN(all$gici10[!is.na(all$fmc.pres)],
  all$entero35[!is.na(all$fmc.pres)])

n10.entero35.fpc = regN(all$gici10[!is.na(all$fpc.pres)],
  all$entero35[!is.na(all$fpc.pres)])

# pooled n's by risk level---------------------------------------
data=all[!is.na(all$fmc.pres),]
data.high=subset(data,data$risk=="High")
n10.entero35.fmc.high = regN(data.high$gici10,data.high$entero35)
data.low=subset(data,data$risk=="Low")
n10.entero35.fmc.low = regN(data.low$gici10,data.low$entero35)

data=all[!is.na(all$fpc.pres),]
data.high=subset(data,data$risk=="High")
n10.entero35.fpc.high = regN(data.high$gici10,data.high$entero35)
data.low=subset(data,data$risk=="Low")
n10.entero35.fpc.low = regN(data.low$gici10,data.low$entero35)


# --------------------------------------
# Estimates pooled across beach
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness
# --------------------------------------

# f- coliphage --------------------------------
all.fit10.entero.fmc <- glm(gici10~entero35,family=poisson(link="log"),data=all[!is.na(all$entero35) & 
    !is.na(all$fmc.pres),])

all.VC10.entero.fmc <- cl(all[!is.na(all$entero35) & !is.na(all$fmc.pres)],
  fm=all.fit10.entero.fmc, cluster=
      all$hhid[!is.na(all$entero35)  & !is.na(all$fmc.pres)])
overall.fit10.entero.fmc <- coeftest(all.fit10.entero.fmc, all.VC10.entero.fmc)
summary(all.fit10.entero.fmc)
overall.fit10.entero.fmc
aic.entero.fmc=AIC(all.fit10.entero.fmc)

# f+ coliphage --------------------------------
all.fit10.entero.fpc <- glm(gici10~entero35,family=poisson(link="log"),data=all[!is.na(all$entero35) & 
    !is.na(all$fpc.pres),])

all.VC10.entero.fpc <- cl(all[!is.na(all$entero35) & !is.na(all$fpc.pres)],
  fm=all.fit10.entero.fpc, cluster=
      all$hhid[!is.na(all$entero35)  & !is.na(all$fpc.pres)])
overall.fit10.entero.fpc <- coeftest(all.fit10.entero.fpc, all.VC10.entero.fpc)
summary(all.fit10.entero.fpc)
overall.fit10.entero.fpc
aic.entero.fpc=AIC(all.fit10.entero.fpc)

# --------------------------------------
# Estimates pooled across beach and stratified by conditions
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness

# all beaches ----------------
# F- coliphage #####################
# high risk conditions --------------------------------
data=all[!is.na(all$entero35) & !is.na(all$fmc.pres),]
data.high=subset(data,data$risk=="High")
all.fit10.entero.high.fmc <- glm(gici10~entero35,family=poisson(link="log"),data=data.high)

all.VC10.entero.high.fmc <- cl(data.high,fm=all.fit10.entero.high.fmc, cluster=data.high$hhid)
overall.fit10.entero.high.fmc <- coeftest(all.fit10.entero.high.fmc, all.VC10.entero.high.fmc)
summary(all.fit10.entero.high.fmc)
overall.fit10.entero.high.fmc
aic.entero.high.fmc=AIC(all.fit10.entero.high.fmc)

# low risk conditions --------------------------------
data.low=subset(data,data$risk=="Low")
all.fit10.entero.low.fmc <- glm(gici10~entero35,family=poisson(link="log"),data=data.low)

all.VC10.entero.low.fmc <- cl(data.low,fm=all.fit10.entero.low.fmc, cluster=data.low$hhid)
overall.fit10.entero.low.fmc <- coeftest(all.fit10.entero.low.fmc, all.VC10.entero.low.fmc)
summary(all.fit10.entero.low.fmc)
overall.fit10.entero.low.fmc
aic.entero.low.fmc=AIC(all.fit10.entero.low.fmc)

# F+ coliphage #####################
# high risk conditions --------------------------------
data=all[!is.na(all$entero35) & !is.na(all$fpc.pres),]
data.high=subset(data,data$risk=="High")
all.fit10.entero.high.fpc <- glm(gici10~entero35,family=poisson(link="log"),data=data.high)

all.VC10.entero.high.fpc <- cl(data.high,fm=all.fit10.entero.high.fpc, cluster=data.high$hhid)
overall.fit10.entero.high.fpc <- coeftest(all.fit10.entero.high.fpc, all.VC10.entero.high.fpc)
summary(all.fit10.entero.high.fpc)
overall.fit10.entero.high.fpc
aic.entero.high.fpc=AIC(all.fit10.entero.high.fpc)

# low risk conditions --------------------------------
data.low=subset(data,data$risk=="Low")
all.fit10.entero.low.fpc <- glm(gici10~entero35,family=poisson(link="log"),data=data.low)

all.VC10.entero.low.fpc <- cl(data.low,fm=all.fit10.entero.low.fpc, cluster=data.low$hhid)
overall.fit10.entero.low.fpc <- coeftest(all.fit10.entero.low.fpc, all.VC10.entero.low.fpc)
summary(all.fit10.entero.low.fpc)
overall.fit10.entero.low.fpc
aic.entero.low.fpc=AIC(all.fit10.entero.low.fpc)


# --------------------------------------
# save the results
# exclude glm objects and data frames
# (they are really large)
# --------------------------------------
save(
  
  n10.entero35.fmc,n10.entero35.fpc,

  n10.entero35.fmc.high,n10.entero35.fmc.low,
  n10.entero35.fpc.high,n10.entero35.fpc.low,
  
  overall.fit10.entero.fmc,overall.fit10.entero.fpc,
  
  overall.fit10.entero.high.fmc,overall.fit10.entero.high.fpc,
  overall.fit10.entero.low.fmc,overall.fit10.entero.low.fpc,
  
  all.VC10.entero.fmc,all.VC10.entero.fpc,

  all.VC10.entero.high.fmc, all.VC10.entero.low.fmc,
  all.VC10.entero.high.fpc, all.VC10.entero.low.fpc,

  aic.entero.fmc,aic.entero.fpc,
  
  aic.entero.low.fmc,aic.entero.low.fpc,
  aic.entero.high.fmc,aic.entero.high.fpc,
  
  file="~/Documents/CRG/coliphage/results/rawoutput/regress-10day-body-entero-pool-unadj.Rdata"
)


