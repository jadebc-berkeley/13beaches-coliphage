##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file conducts maximum likelihood regression
# to estimate prevalence ratios

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
all$pres[all$fmc.pres==1 | all$fpc.pres==1]=1
all$pres[all$fmc.pres==0 & all$fpc.pres==0]=0

# --------------------------------------
# Estimates pooled across beach and assay method
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness
# --------------------------------------

all.fit10 <- glm(gici10~pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
                           rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$pres),])

all.VC10 <- cl(all[!is.na(all$pres)],fm=all.fit10,
                       cluster=all$hhid[!is.na(all$pres)])
overall.fit10 <- coeftest(all.fit10, all.VC10)
summary(all.fit10)
overall.fit10
aic=AIC(all.fit10)


# --------------------------------------
# Estimates pooled across beach and stratified by conditions
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness

# all beaches ----------------

# high risk conditions
data=all[!is.na(all$pres),]
data.high=subset(data,data$risk=="High")
all.fit10.high <- glm(gici10~pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC10.high <- cl(data.high,fm=all.fit10.high, cluster=data.high$hhid)
overall.fit10.high <- coeftest(all.fit10.high, all.VC10.high)
summary(all.fit10.high)
overall.fit10.high
aic.high=AIC(all.fit10.high)

# low risk conditions
data.low=subset(data,data$risk=="Low")
all.fit10.low <- glm(gici10~pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC10.low <- cl(data.low,fm=all.fit10.low, cluster=data.low$hhid)
overall.fit10.low <- coeftest(all.fit10.low, all.VC10.low)
summary(all.fit10.low)
overall.fit10.low
aic.low=AIC(all.fit10.low)



# --------------------------------------
# save the results
# exclude glm objects and data frames
# (they are really large)
# --------------------------------------
save(

  all.VC10, overall.fit10,

  all.VC10.high,all.VC10.low,
  overall.fit10.high,overall.fit10.low,
  
  aic,  aic.high, aic.low,
  
  file="~/Documents/CRG/coliphage/results/rawoutput/regress-10day-body-pool-both.Rdata"
)

