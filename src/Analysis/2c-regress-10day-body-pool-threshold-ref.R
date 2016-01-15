##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file conducts maximum likelihood regression
# to estimate prevalence ratios

# Results pooled across beaches and assay

# 10 day gi illness
##########################################

rm(list=ls())
library(foreign)

setwd("~/Dropbox/Coliphage/")

# --------------------------------------
# load the and pre-preprocess the 
# analysis dataset
# (refer to the base functions script
# for details on the pre-processing)
# --------------------------------------
beaches13=read.csv("~/Dropbox/13beaches-fork-coliphage/data/final/13beaches-analysis.csv")

# load base functions
source("Programs/Analysis/0-base-functions.R")

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

# n's pooled across assay and beach ---------------------------------------
all.n10.fmc.25 = regN(all$gici10[!is.na(all$fmc_25ref)],
                       all$fmc_25ref[!is.na(all$fmc_25ref)])
all.n10.fmc.50 = regN(all$gici10[!is.na(all$fmc_50ref)],
                       all$fmc_50ref[!is.na(all$fmc_50ref)])
all.n10.fmc.75 = regN(all$gici10[!is.na(all$fmc_75ref)],
                       all$fmc_75ref[!is.na(all$fmc_75ref)])
all.n10.fpc.25 = regN(all$gici10[!is.na(all$fpc_25ref)],
                       all$fpc_25ref[!is.na(all$fpc_25ref)])
all.n10.fpc.50 = regN(all$gici10[!is.na(all$fpc_50ref)],
                       all$fpc_50ref[!is.na(all$fpc_50ref)])
all.n10.fpc.75 = regN(all$gici10[!is.na(all$fpc_75ref)],
                       all$fpc_75ref[!is.na(all$fpc_75ref)])

# pooled n's by risk level---------------------------------------
data=all[!is.na(all$fmc.pres),]
data.high=subset(data,data$risk=="High")
all.n10.fmc.high.25 = regN(data.high$gici10,data.high$fmc_25ref)
all.n10.fmc.high.50 = regN(data.high$gici10,data.high$fmc_50ref)
all.n10.fmc.high.75 = regN(data.high$gici10,data.high$fmc_75ref)
data.low=subset(data,data$risk=="Low")
all.n10.fmc.low.25 = regN(data.low$gici10,data.low$fmc_25ref)
all.n10.fmc.low.50 = regN(data.low$gici10,data.low$fmc_50ref)
all.n10.fmc.low.75 = regN(data.low$gici10,data.low$fmc_75ref)

data=all[!is.na(all$fpc.pres),]
data.high=subset(data,data$risk=="High")
all.n10.fpc.high.25 = regN(data.high$gici10,data.high$fpc_25ref)
all.n10.fpc.high.50 = regN(data.high$gici10,data.high$fpc_50ref)
all.n10.fpc.high.75 = regN(data.high$gici10,data.high$fpc_75ref)
data.low=subset(data,data$risk=="Low")
all.n10.fpc.low.25 = regN(data.low$gici10,data.low$fpc_25ref)
all.n10.fpc.low.50 = regN(data.low$gici10,data.low$fpc_50ref)
all.n10.fpc.low.75 = regN(data.low$gici10,data.low$fpc_75ref)



# --------------------------------------
# Estimates pooled across beach and assay method
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness

# f- coliphage ----------------

# >25% of samples with detectable coliphage
all.fit10.fmc.25 <- glm(gici10~fmc_25ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
                           rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc_25ref),])

all.VC10.fmc.25 <- cl(all[!is.na(all$fmc_25ref)],fm=all.fit10.fmc.25,
                       cluster=all$hhid[!is.na(all$fmc_25ref)])
overall.fit10.fmc.25 <- coeftest(all.fit10.fmc.25, all.VC10.fmc.25)
summary(all.fit10.fmc.25)
overall.fit10.fmc.25
aic.fmc.25=AIC(all.fit10.fmc.25)

# >50% of samples with detectable coliphage
all.fit10.fmc.50 <- glm(gici10~fmc_50ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc_50ref),])

all.VC10.fmc.50 <- cl(all[!is.na(all$fmc_50ref)],fm=all.fit10.fmc.50,
     cluster=all$hhid[!is.na(all$fmc_50ref)])
overall.fit10.fmc.50 <- coeftest(all.fit10.fmc.50, all.VC10.fmc.50)
summary(all.fit10.fmc.50)
overall.fit10.fmc.50
aic.fmc.50=AIC(all.fit10.fmc.50)

# >75% of samples with detectable coliphage
all.fit10.fmc.75 <- glm(gici10~fmc_75ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc_75ref),])

all.VC10.fmc.75 <- cl(all[!is.na(all$fmc_75ref)],fm=all.fit10.fmc.75,
    cluster=all$hhid[!is.na(all$fmc_75ref)])
overall.fit10.fmc.75 <- coeftest(all.fit10.fmc.75, all.VC10.fmc.75)
summary(all.fit10.fmc.75)
overall.fit10.fmc.75
aic.fmc.75=AIC(all.fit10.fmc.75)

# f+ coliphage ----------------
# >25% of samples with detectable coliphage
all.fit10.fpc.25 <- glm(gici10~fpc_25ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
                           rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc_25ref),])

all.VC10.fpc.25 <- cl(all[!is.na(all$fpc_25ref)],fm=all.fit10.fpc.25,
                       cluster=all$hhid[!is.na(all$fpc_25ref)])
overall.fit10.fpc.25 <- coeftest(all.fit10.fpc.25, all.VC10.fpc.25)
summary(all.fit10.fpc.25)
overall.fit10.fpc.25
aic.fpc.25=AIC(all.fit10.fpc.25)

# >50% of samples with detectable coliphage
all.fit10.fpc.50 <- glm(gici10~fpc_50ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc_50ref),])

all.VC10.fpc.50 <- cl(all[!is.na(all$fpc_50ref)],fm=all.fit10.fpc.50,
     cluster=all$hhid[!is.na(all$fpc_50ref)])
overall.fit10.fpc.50 <- coeftest(all.fit10.fpc.50, all.VC10.fpc.50)
summary(all.fit10.fpc.50)
overall.fit10.fpc.50
aic.fpc.50=AIC(all.fit10.fpc.50)

# >75% of samples with detectable coliphage
all.fit10.fpc.75 <- glm(gici10~fpc_75ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc_75ref),])

all.VC10.fpc.75 <- cl(all[!is.na(all$fpc_75ref)],fm=all.fit10.fpc.75,
    cluster=all$hhid[!is.na(all$fpc_75ref)])
overall.fit10.fpc.75 <- coeftest(all.fit10.fpc.75, all.VC10.fpc.75)
summary(all.fit10.fpc.75)
overall.fit10.fpc.75
aic.fpc.75=AIC(all.fit10.fpc.75)

# --------------------------------------
# Estimates pooled across beach and stratified by conditions
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness

# all beaches ----------------
# f- coliphage --------
# high risk conditions

# >25% of samples with detectable coliphage
data=all[!is.na(all$fmc_25ref),]
data.high=subset(data,data$risk=="High")
all.fit10.fmc.high.25 <- glm(gici10~fmc_25ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC10.fmc.high.25 <- cl(data.high,fm=all.fit10.fmc.high.25, cluster=data.high$hhid)
overall.fit10.fmc.high.25 <- coeftest(all.fit10.fmc.high.25, all.VC10.fmc.high.25)
summary(all.fit10.fmc.high.25)
overall.fit10.fmc.high.25
aic.fmc.high.25=AIC(all.fit10.fmc.high.25)

# >50% of samples with detectable coliphage
data=all[!is.na(all$fmc_50ref),]
data.high=subset(data,data$risk=="High")
all.fit10.fmc.high.50 <- glm(gici10~fmc_50ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC10.fmc.high.50 <- cl(data.high,fm=all.fit10.fmc.high.50, cluster=data.high$hhid)
overall.fit10.fmc.high.50 <- coeftest(all.fit10.fmc.high.50, all.VC10.fmc.high.50)
summary(all.fit10.fmc.high.50)
overall.fit10.fmc.high.50
aic.fmc.high.50=AIC(all.fit10.fmc.high.50)

# >75% of samples with detectable coliphage
data=all[!is.na(all$fmc_75ref),]
data.high=subset(data,data$risk=="High")
all.fit10.fmc.high.75 <- glm(gici10~fmc_75ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC10.fmc.high.75 <- cl(data.high,fm=all.fit10.fmc.high.75, cluster=data.high$hhid)
overall.fit10.fmc.high.75 <- coeftest(all.fit10.fmc.high.75, all.VC10.fmc.high.75)
summary(all.fit10.fmc.high.75)
overall.fit10.fmc.high.75
aic.fmc.high.75=AIC(all.fit10.fmc.high.75)

# low risk conditions

# >25% of samples with detectable coliphage
data=all[!is.na(all$fmc_25ref),]
data.low=subset(data,data$risk=="Low")
all.fit10.fmc.low.25 <- glm(gici10~fmc_25ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC10.fmc.low.25 <- cl(data.low,fm=all.fit10.fmc.low.25, cluster=data.low$hhid)
overall.fit10.fmc.low.25 <- coeftest(all.fit10.fmc.low.25, all.VC10.fmc.low.25)
summary(all.fit10.fmc.low.25)
overall.fit10.fmc.low.25
aic.fmc.low.25=AIC(all.fit10.fmc.low.25)

# >50% of samples with detectable coliphage
data=all[!is.na(all$fmc_50ref),]
data.low=subset(data,data$risk=="Low")
all.fit10.fmc.low.50 <- glm(gici10~fmc_50ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC10.fmc.low.50 <- cl(data.low,fm=all.fit10.fmc.low.50, cluster=data.low$hhid)
overall.fit10.fmc.low.50 <- coeftest(all.fit10.fmc.low.50, all.VC10.fmc.low.50)
summary(all.fit10.fmc.low.50)
overall.fit10.fmc.low.50
aic.fmc.low.50=AIC(all.fit10.fmc.low.50)

# >75% of samples with detectable coliphage
data=all[!is.na(all$fmc_75ref),]
data.low=subset(data,data$risk=="Low")
all.fit10.fmc.low.75 <- glm(gici10~fmc_75ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC10.fmc.low.75 <- cl(data.low,fm=all.fit10.fmc.low.75, cluster=data.low$hhid)
overall.fit10.fmc.low.75 <- coeftest(all.fit10.fmc.low.75, all.VC10.fmc.low.75)
summary(all.fit10.fmc.low.75)
overall.fit10.fmc.low.75
aic.fmc.low.75=AIC(all.fit10.fmc.low.75)

# f+ coliphage  --------
# >25% of samples with detectable coliphage
data=all[!is.na(all$fpc_25ref),]
data.high=subset(data,data$risk=="High")
all.fit10.fpc.high.25 <- glm(gici10~fpc_25ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC10.fpc.high.25 <- cl(data.high,fm=all.fit10.fpc.high.25, cluster=data.high$hhid)
overall.fit10.fpc.high.25 <- coeftest(all.fit10.fpc.high.25, all.VC10.fpc.high.25)
summary(all.fit10.fpc.high.25)
overall.fit10.fpc.high.25
aic.fpc.high.25=AIC(all.fit10.fpc.high.25)

# >50% of samples with detectable coliphage
data=all[!is.na(all$fpc_50ref),]
data.high=subset(data,data$risk=="High")
all.fit10.fpc.high.50 <- glm(gici10~fpc_50ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC10.fpc.high.50 <- cl(data.high,fm=all.fit10.fpc.high.50, cluster=data.high$hhid)
overall.fit10.fpc.high.50 <- coeftest(all.fit10.fpc.high.50, all.VC10.fpc.high.50)
summary(all.fit10.fpc.high.50)
overall.fit10.fpc.high.50
aic.fpc.high.50=AIC(all.fit10.fpc.high.50)

# >75% of samples with detectable coliphage
data=all[!is.na(all$fpc_75ref),]
data.high=subset(data,data$risk=="High")
all.fit10.fpc.high.75 <- glm(gici10~fpc_75ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC10.fpc.high.75 <- cl(data.high,fm=all.fit10.fpc.high.75, cluster=data.high$hhid)
overall.fit10.fpc.high.75 <- coeftest(all.fit10.fpc.high.75, all.VC10.fpc.high.75)
summary(all.fit10.fpc.high.75)
overall.fit10.fpc.high.75
aic.fpc.high.75=AIC(all.fit10.fpc.high.75)

# low risk conditions

# >25% of samples with detectable coliphage
data=all[!is.na(all$fpc_25ref),]
data.low=subset(data,data$risk=="Low")
all.fit10.fpc.low.25 <- glm(gici10~fpc_25ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC10.fpc.low.25 <- cl(data.low,fm=all.fit10.fpc.low.25, cluster=data.low$hhid)
overall.fit10.fpc.low.25 <- coeftest(all.fit10.fpc.low.25, all.VC10.fpc.low.25)
summary(all.fit10.fpc.low.25)
overall.fit10.fpc.low.25
aic.fpc.low.25=AIC(all.fit10.fpc.low.25)

# >50% of samples with detectable coliphage
data=all[!is.na(all$fpc_50ref),]
data.low=subset(data,data$risk=="Low")
all.fit10.fpc.low.50 <- glm(gici10~fpc_50ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC10.fpc.low.50 <- cl(data.low,fm=all.fit10.fpc.low.50, cluster=data.low$hhid)
overall.fit10.fpc.low.50 <- coeftest(all.fit10.fpc.low.50, all.VC10.fpc.low.50)
summary(all.fit10.fpc.low.50)
overall.fit10.fpc.low.50
aic.fpc.low.50=AIC(all.fit10.fpc.low.50)

# >75% of samples with detectable coliphage
data=all[!is.na(all$fpc_75ref),]
data.low=subset(data,data$risk=="Low")
all.fit10.fpc.low.75 <- glm(gici10~fpc_75ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC10.fpc.low.75 <- cl(data.low,fm=all.fit10.fpc.low.75, cluster=data.low$hhid)
overall.fit10.fpc.low.75 <- coeftest(all.fit10.fpc.low.75, all.VC10.fpc.low.75)
summary(all.fit10.fpc.low.75)
overall.fit10.fpc.low.75
aic.fpc.low.75=AIC(all.fit10.fpc.low.75)




# --------------------------------------
# save the results
# exclude glm objects and data frames
# (they are really large)
# --------------------------------------
save(

#   all.n10.fmc.25,all.n10.fmc.50,all.n10.fmc.75,
#   all.n10.fpc.25,all.n10.fpc.50,all.n10.fpc.75,
#   all.n10.fmc.high.25,all.n10.fmc.high.50,all.n10.fmc.high.75,
#   all.n10.fmc.low.25,all.n10.fmc.low.50,all.n10.fmc.low.75,
#   all.n10.fpc.high.25,all.n10.fpc.high.50,all.n10.fpc.high.75,
#   all.n10.fpc.low.25,all.n10.fpc.low.50,all.n10.fpc.low.75,

  all.VC10.fmc.25, all.VC10.fmc.50, all.VC10.fmc.75, 
  all.VC10.fpc.25, all.VC10.fpc.50, all.VC10.fpc.75,
  overall.fit10.fmc.25,overall.fit10.fmc.50,overall.fit10.fmc.75,
  overall.fit10.fpc.25,overall.fit10.fpc.50,overall.fit10.fpc.75,

  all.VC10.fmc.high.25,all.VC10.fmc.high.50,all.VC10.fmc.high.75,
  all.VC10.fpc.high.25,all.VC10.fpc.high.50,all.VC10.fpc.high.75,
  overall.fit10.fmc.high.25,overall.fit10.fmc.high.50,overall.fit10.fmc.high.75,
  overall.fit10.fpc.high.25,overall.fit10.fpc.high.50,overall.fit10.fpc.high.75,
  
  all.VC10.fmc.low.25,all.VC10.fmc.low.50,all.VC10.fmc.low.75,
  all.VC10.fpc.low.25,all.VC10.fpc.low.50,all.VC10.fpc.low.75,
  overall.fit10.fmc.low.25,overall.fit10.fmc.low.50,overall.fit10.fmc.low.75,
  overall.fit10.fpc.low.25,overall.fit10.fpc.low.50,overall.fit10.fpc.low.75,
  
  aic.fmc.25,aic.fmc.50,aic.fmc.75,
  aic.fpc.25,aic.fpc.50,aic.fpc.75,
  aic.fmc.high.25, aic.fmc.high.50, aic.fmc.high.75, 
  aic.fpc.high.25,aic.fpc.high.50,aic.fpc.high.75,
  aic.fmc.low.25,aic.fmc.low.50,aic.fmc.low.75, 
  aic.fpc.low.25,aic.fpc.low.50,aic.fpc.low.75,
  
  file="~/dropbox/coliphage/results/rawoutput/regress-10day-body-pool-threshold-ref.Rdata"
)

