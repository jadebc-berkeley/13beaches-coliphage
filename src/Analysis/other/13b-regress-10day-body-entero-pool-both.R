##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file conducts maximum likelihood regression
# to estimate prevalence ratios
# for enterococcus

# Results pooled across beaches, assays, and type
# of coliphage

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

# create indicator for pooled presence absence

all$pres=NA
all$pres[all$pres==1 | all$fmc.pres==1]=1
all$pres[all$pres==0 & all$fpc.pres==0]=0

# --------------------------------------
# Estimates pooled across beach
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness
# --------------------------------------

# f- coliphage --------------------------------
all.fit10.entero <- glm(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$entero35) & 
    !is.na(all$pres),])

all.VC10.entero <- cl(all[!is.na(all$entero35) & !is.na(all$pres)],
  fm=all.fit10.entero, cluster=
      all$hhid[!is.na(all$entero35)  & !is.na(all$pres)])
overall.fit10.entero <- coeftest(all.fit10.entero, all.VC10.entero)
summary(all.fit10.entero)
overall.fit10.entero
aic.entero=AIC(all.fit10.entero)


# --------------------------------------
# Estimates pooled across beach and stratified by conditions
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness

# all beaches ----------------

# high risk conditions --------------------------------
data=all[!is.na(all$entero35) & !is.na(all$pres),]
data.high=subset(data,data$risk=="High")
all.fit10.entero.high <- glm(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC10.entero.high <- cl(data.high,fm=all.fit10.entero.high, cluster=data.high$hhid)
overall.fit10.entero.high <- coeftest(all.fit10.entero.high, all.VC10.entero.high)
summary(all.fit10.entero.high)
overall.fit10.entero.high
aic.entero.high=AIC(all.fit10.entero.high)

# low risk conditions --------------------------------
data.low=subset(data,data$risk=="Low")
all.fit10.entero.low <- glm(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC10.entero.low <- cl(data.low,fm=all.fit10.entero.low, cluster=data.low$hhid)
overall.fit10.entero.low <- coeftest(all.fit10.entero.low, all.VC10.entero.low)
summary(all.fit10.entero.low)
overall.fit10.entero.low
aic.entero.low=AIC(all.fit10.entero.low)


# --------------------------------------
# save the results
# exclude glm objects and data frames
# (they are really large)
# --------------------------------------
save(
  
  overall.fit10.entero, overall.fit10.entero.high,overall.fit10.entero.low,
  
  all.VC10.entero,all.VC10.entero.high, all.VC10.entero.low,

  aic.entero,  aic.entero.low, aic.entero.high,
  
  file="~/Documents/CRG/coliphage/results/rawoutput/regress-10day-body-entero-pool-both.Rdata"
)


