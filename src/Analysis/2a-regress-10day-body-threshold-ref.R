##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file conducts maximum likelihood regression
# to estimate prevalence ratios

# Results pooled across beaches

# New definition of high risk

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

# 
# # n's pooled across beach ---------------------------------------
# all.n10.fmc1601 = regN(all$gici10[!is.na(all$fmc1601.pres)],
#                        all$fmc1601.pres[!is.na(all$fmc1601.pres)])
# all.n10.fmc1602 = regN(all$gici10[!is.na(all$fmc1602.pres)],
#                        all$fmc1602.pres[!is.na(all$fmc1602.pres)])
# all.n10.fpc1601 = regN(all$gici10[!is.na(all$fpc1601.pres)],
#                        all$fpc1601.pres[!is.na(all$fpc1601.pres)])
# all.n10.fpc1602 = regN(all$gici10[!is.na(all$fpc1602.pres)],
#                        all$fpc1602.pres[!is.na(all$fpc1602.pres)])
# 
# # pooled n's by risk level---------------------------------------
# data=all[!is.na(all$fmc1602.pres),]
# data.high=subset(data,data$risk=="High")
# all.n10.fmc1602.high = regN(data.high$gici10,data.high$fmc1602.pres)
# data.low=subset(data,data$risk=="Low")
# all.n10.fmc1602.low = regN(data.low$gici10,data.low$fmc1602.pres)
# 
# data=all[!is.na(all$fpc1601.pres),]
# data.high=subset(data,data$risk=="High")
# all.n10.fpc1601.high = regN(data.high$gici10,data.high$fpc1601.pres)
# data.low=subset(data,data$risk=="Low")
# all.n10.fpc1601.low = regN(data.low$gici10,data.low$fpc1601.pres)
# 
# data=all[!is.na(all$fpc1602.pres),]
# data.high=subset(data,data$risk=="High")
# all.n10.fpc1602.high = regN(data.high$gici10,data.high$fpc1602.pres)
# data.low=subset(data,data$risk=="Low")
# all.n10.fpc1602.low = regN(data.low$gici10,data.low$fpc1602.pres)
# 


# --------------------------------------
# Estimates pooled across beach
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness

# all beaches ----------------
# fmc 1601
all.fit10.fmc1601.25 <- glm(gici10~fmc1601_25ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc1601_25ref),])

all.VC10.fmc1601.25 <- cl(all[!is.na(all$fmc1601_25ref)],fm=all.fit10.fmc1601.25,
  cluster=all$hhid[!is.na(all$fmc1601_25ref)])
overall.fit10.fmc1601.25 <- coeftest(all.fit10.fmc1601.25, all.VC10.fmc1601.25)
summary(all.fit10.fmc1601.25)
overall.fit10.fmc1601.25
aic.fmc1601.25=AIC(all.fit10.fmc1601.25)

all.fit10.fmc1601.50 <- glm(gici10~fmc1601_50ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc1601_50ref),])

all.VC10.fmc1601.50 <- cl(all[!is.na(all$fmc1601_50ref)],fm=all.fit10.fmc1601.50,
  cluster=all$hhid[!is.na(all$fmc1601_50ref)])
overall.fit10.fmc1601.50 <- coeftest(all.fit10.fmc1601.50, all.VC10.fmc1601.50)
summary(all.fit10.fmc1601.50)
overall.fit10.fmc1601.50
aic.fmc1601.50=AIC(all.fit10.fmc1601.50)

all.fit10.fmc1601.75 <- glm(gici10~fmc1601_75ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc1601_75ref),])

all.VC10.fmc1601.75 <- cl(all[!is.na(all$fmc1601_75ref)],fm=all.fit10.fmc1601.75,
  cluster=all$hhid[!is.na(all$fmc1601_75ref)])
overall.fit10.fmc1601.75 <- coeftest(all.fit10.fmc1601.75, all.VC10.fmc1601.75)
summary(all.fit10.fmc1601.75)
overall.fit10.fmc1601.75
aic.fmc1601.75=AIC(all.fit10.fmc1601.75)


# fmc 1602
all.fit10.fmc1602.25 <- glm(gici10~fmc1602_25ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc1602_25ref),])

all.VC10.fmc1602.25 <- cl(all[!is.na(all$fmc1602_25ref)],fm=all.fit10.fmc1602.25,
  cluster=all$hhid[!is.na(all$fmc1602_25ref)])
overall.fit10.fmc1602.25 <- coeftest(all.fit10.fmc1602.25, all.VC10.fmc1602.25)
summary(all.fit10.fmc1602.25)
overall.fit10.fmc1602.25
aic.fmc1602.25=AIC(all.fit10.fmc1602.25)

all.fit10.fmc1602.50 <- glm(gici10~fmc1602_50ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc1602_50ref),])

all.VC10.fmc1602.50 <- cl(all[!is.na(all$fmc1602_50ref)],fm=all.fit10.fmc1602.50,
  cluster=all$hhid[!is.na(all$fmc1602_50ref)])
overall.fit10.fmc1602.50 <- coeftest(all.fit10.fmc1602.50, all.VC10.fmc1602.50)
summary(all.fit10.fmc1602.50)
overall.fit10.fmc1602.50
aic.fmc1602.50=AIC(all.fit10.fmc1602.50)

all.fit10.fmc1602.75 <- glm(gici10~fmc1602_75ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc1602_75ref),])

all.VC10.fmc1602.75 <- cl(all[!is.na(all$fmc1602_75ref)],fm=all.fit10.fmc1602.75,
  cluster=all$hhid[!is.na(all$fmc1602_75ref)])
overall.fit10.fmc1602.75 <- coeftest(all.fit10.fmc1602.75, all.VC10.fmc1602.75)
summary(all.fit10.fmc1602.75)
overall.fit10.fmc1602.75
aic.fmc1602.75=AIC(all.fit10.fmc1602.75)

# fpc 1601
all.fit10.fpc1601.25 <- glm(gici10~fpc1601_25ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1601_25ref),])

all.VC10.fpc1601.25 <- cl(all[!is.na(all$fpc1601_25ref)],fm=all.fit10.fpc1601.25,
  cluster=all$hhid[!is.na(all$fpc1601_25ref)])
overall.fit10.fpc1601.25 <- coeftest(all.fit10.fpc1601.25, all.VC10.fpc1601.25)
summary(all.fit10.fpc1601.25)
overall.fit10.fpc1601.25
aic.fpc1601.25=AIC(all.fit10.fpc1601.25)

all.fit10.fpc1601.50 <- glm(gici10~fpc1601_50ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1601_50ref),])

all.VC10.fpc1601.50 <- cl(all[!is.na(all$fpc1601_50ref)],fm=all.fit10.fpc1601.50,
  cluster=all$hhid[!is.na(all$fpc1601_50ref)])
overall.fit10.fpc1601.50 <- coeftest(all.fit10.fpc1601.50, all.VC10.fpc1601.50)
summary(all.fit10.fpc1601.50)
overall.fit10.fpc1601.50
aic.fpc1601.50=AIC(all.fit10.fpc1601.50)

all.fit10.fpc1601.75 <- glm(gici10~fpc1601_75ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1601_75ref),])

all.VC10.fpc1601.75 <- cl(all[!is.na(all$fpc1601_75ref)],fm=all.fit10.fpc1601.75,
  cluster=all$hhid[!is.na(all$fpc1601_75ref)])
overall.fit10.fpc1601.75 <- coeftest(all.fit10.fpc1601.75, all.VC10.fpc1601.75)
summary(all.fit10.fpc1601.75)
overall.fit10.fpc1601.75
aic.fpc1601.75=AIC(all.fit10.fpc1601.75)

# fpc 1602
all.fit10.fpc1602.25 <- glm(gici10~fpc1602_25ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1602_25ref),])

all.VC10.fpc1602.25 <- cl(all[!is.na(all$fpc1602_25ref)],fm=all.fit10.fpc1602.25,
  cluster=all$hhid[!is.na(all$fpc1602_25ref)])
overall.fit10.fpc1602.25 <- coeftest(all.fit10.fpc1602.25, all.VC10.fpc1602.25)
summary(all.fit10.fpc1602.25)
overall.fit10.fpc1602.25
aic.fpc1602.25=AIC(all.fit10.fpc1602.25)

all.fit10.fpc1602.50 <- glm(gici10~fpc1602_50ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1602_50ref),])

all.VC10.fpc1602.50 <- cl(all[!is.na(all$fpc1602_50ref)],fm=all.fit10.fpc1602.50,
  cluster=all$hhid[!is.na(all$fpc1602_50ref)])
overall.fit10.fpc1602.50 <- coeftest(all.fit10.fpc1602.50, all.VC10.fpc1602.50)
summary(all.fit10.fpc1602.50)
overall.fit10.fpc1602.50
aic.fpc1602.50=AIC(all.fit10.fpc1602.50)

all.fit10.fpc1602.75 <- glm(gici10~fpc1602_75ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1602_75ref),])

all.VC10.fpc1602.75 <- cl(all[!is.na(all$fpc1602_75ref)],fm=all.fit10.fpc1602.75,
  cluster=all$hhid[!is.na(all$fpc1602_75ref)])
overall.fit10.fpc1602.75 <- coeftest(all.fit10.fpc1602.75, all.VC10.fpc1602.75)
summary(all.fit10.fpc1602.75)
overall.fit10.fpc1602.75
aic.fpc1602.75=AIC(all.fit10.fpc1602.75)
# --------------------------------------
# Estimates pooled across beach and stratified by conditions
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness

# all beaches ----------------

# fmc 1602 --------
# high risk conditions
data.high=subset(all,all$risk=="High")
all.fit10.fmc1602.high.25 <- glm(gici10~fmc1602_25ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high[!is.na(data.high$fmc1602_25ref),])

all.VC10.fmc1602.high.25 <- cl(data.high[!is.na(data.high$fmc1602_25ref),],
        fm=all.fit10.fmc1602.high.25, cluster=data.high$hhid[!is.na(data.high$fmc1602_25ref)])
overall.fit10.fmc1602.high.25 <- coeftest(all.fit10.fmc1602.high.25, all.VC10.fmc1602.high.25)
summary(all.fit10.fmc1602.high.25)
overall.fit10.fmc1602.high.25
aic.fmc1602.high.25=AIC(all.fit10.fmc1602.high.25)

all.fit10.fmc1602.high.50 <- glm(gici10~fmc1602_50ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high[!is.na(data.high$fmc1602_50ref),])

all.VC10.fmc1602.high.50 <- cl(data.high[!is.na(data.high$fmc1602_50ref),],fm=all.fit10.fmc1602.high.50, 
                               cluster=data.high$hhid[!is.na(data.high$fmc1602_50ref)])
overall.fit10.fmc1602.high.50 <- coeftest(all.fit10.fmc1602.high.50, all.VC10.fmc1602.high.50)
summary(all.fit10.fmc1602.high.50)
overall.fit10.fmc1602.high.50
aic.fmc1602.high.50=AIC(all.fit10.fmc1602.high.50)

all.fit10.fmc1602.high.75 <- glm(gici10~fmc1602_75ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high[!is.na(data.high$fmc1602_75ref),])

all.VC10.fmc1602.high.75 <- cl(data.high[!is.na(data.high$fmc1602_50ref),],fm=all.fit10.fmc1602.high.75, 
                      cluster=data.high$hhid[!is.na(data.high$fmc1602_75ref)])
overall.fit10.fmc1602.high.75 <- coeftest(all.fit10.fmc1602.high.75, all.VC10.fmc1602.high.75)
summary(all.fit10.fmc1602.high.75)
overall.fit10.fmc1602.high.75
aic.fmc1602.high.75=AIC(all.fit10.fmc1602.high.75)

# low risk conditions
data.low=subset(all,all$risk=="Low")
all.fit10.fmc1602.low.25 <- glm(gici10~fmc1602_25ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low[!is.na(data.low$fmc1602_25ref),])

all.VC10.fmc1602.low.25 <- cl(data.low[!is.na(data.low$fmc1602_25ref),],
     fm=all.fit10.fmc1602.low.25, cluster=data.low$hhid[!is.na(data.low$fmc1602_25ref)])
overall.fit10.fmc1602.low.25 <- coeftest(all.fit10.fmc1602.low.25, all.VC10.fmc1602.low.25)
summary(all.fit10.fmc1602.low.25)
overall.fit10.fmc1602.low.25
aic.fmc1602.low.25=AIC(all.fit10.fmc1602.low.25)

all.fit10.fmc1602.low.50 <- glm(gici10~fmc1602_50ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low[!is.na(data.low$fmc1602_50ref),])

all.VC10.fmc1602.low.50 <- cl(data.low[!is.na(data.low$fmc1602_50ref),],
            fm=all.fit10.fmc1602.low.50, cluster=data.low$hhid[!is.na(data.low$fmc1602_50ref)])
overall.fit10.fmc1602.low.50 <- coeftest(all.fit10.fmc1602.low.50, all.VC10.fmc1602.low.50)
summary(all.fit10.fmc1602.low.50)
overall.fit10.fmc1602.low.50
aic.fmc1602.low.50=AIC(all.fit10.fmc1602.low.50)

all.fit10.fmc1602.low.75 <- glm(gici10~fmc1602_75ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low[!is.na(data.low$fmc1602_75ref),])

all.VC10.fmc1602.low.75 <- cl(data.low[!is.na(data.low$fmc1602_75ref),],
            fm=all.fit10.fmc1602.low.75, cluster=data.low$hhid[!is.na(data.low$fmc1602_75ref)])
overall.fit10.fmc1602.low.75 <- coeftest(all.fit10.fmc1602.low.75, all.VC10.fmc1602.low.75)
summary(all.fit10.fmc1602.low.75)
overall.fit10.fmc1602.low.75
aic.fmc1602.low.75=AIC(all.fit10.fmc1602.low.75)

# fpc 1601 --------
# high risk conditions
data.high=subset(all,all$risk=="High")
all.fit10.fpc1601.high.25 <- glm(gici10~fpc1601_25ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high[!is.na(data.high$fpc1601_25ref),])

all.VC10.fpc1601.high.25 <- cl(data.high[!is.na(data.high$fpc1601_25ref),],
      fm=all.fit10.fpc1601.high.25, cluster=data.high$hhid[!is.na(data.high$fpc1601_25ref)])
overall.fit10.fpc1601.high.25 <- coeftest(all.fit10.fpc1601.high.25, all.VC10.fpc1601.high.25)
summary(all.fit10.fpc1601.high.25)
overall.fit10.fpc1601.high.25
aic.fpc1601.high.25=AIC(all.fit10.fpc1601.high.25)

all.fit10.fpc1601.high.50 <- glm(gici10~fpc1601_50ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high[!is.na(data.high$fpc1601_50ref),])

all.VC10.fpc1601.high.50 <- cl(data.high[!is.na(data.high$fpc1601_50ref),],
      fm=all.fit10.fpc1601.high.50, cluster=data.high$hhid[!is.na(data.high$fpc1601_50ref)])
overall.fit10.fpc1601.high.50 <- coeftest(all.fit10.fpc1601.high.50, all.VC10.fpc1601.high.50)
summary(all.fit10.fpc1601.high.50)
overall.fit10.fpc1601.high.50
aic.fpc1601.high.50=AIC(all.fit10.fpc1601.high.50)

all.fit10.fpc1601.high.75 <- glm(gici10~fpc1601_75ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high[!is.na(data.high$fpc1601_75ref),])

all.VC10.fpc1601.high.75 <- cl(data.high[!is.na(data.high$fpc1601_75ref),],
    fm=all.fit10.fpc1601.high.75, cluster=data.high$hhid[!is.na(data.high$fpc1601_75ref)])
overall.fit10.fpc1601.high.75 <- coeftest(all.fit10.fpc1601.high.75, all.VC10.fpc1601.high.75)
summary(all.fit10.fpc1601.high.75)
overall.fit10.fpc1601.high.75
aic.fpc1601.high.75=AIC(all.fit10.fpc1601.high.75)

# low risk conditions
data.low=subset(all,all$risk=="Low")
all.fit10.fpc1601.low.25 <- glm(gici10~fpc1601_25ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low[!is.na(data.low$fpc1601_25ref),])

all.VC10.fpc1601.low.25 <- cl(data.low[!is.na(data.low$fpc1601_25ref),],
    fm=all.fit10.fpc1601.low.25, cluster=data.low$hhid[!is.na(data.low$fpc1601_25ref)])
overall.fit10.fpc1601.low.25 <- coeftest(all.fit10.fpc1601.low.25, all.VC10.fpc1601.low.25)
summary(all.fit10.fpc1601.low.25)
overall.fit10.fpc1601.low.25
aic.fpc1601.low.25=AIC(all.fit10.fpc1601.low.25)

all.fit10.fpc1601.low.50 <- glm(gici10~fpc1601_50ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low[!is.na(data.low$fpc1601_50ref),])

all.VC10.fpc1601.low.50 <- cl(data.low[!is.na(data.low$fpc1601_50ref),],
    fm=all.fit10.fpc1601.low.50, cluster=data.low$hhid[!is.na(data.low$fpc1601_50ref)])
overall.fit10.fpc1601.low.50 <- coeftest(all.fit10.fpc1601.low.50, all.VC10.fpc1601.low.50)
summary(all.fit10.fpc1601.low.50)
overall.fit10.fpc1601.low.50
aic.fpc1601.low.50=AIC(all.fit10.fpc1601.low.50)

all.fit10.fpc1601.low.75 <- glm(gici10~fpc1601_75ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low[!is.na(data.low$fpc1601_75ref),])

all.VC10.fpc1601.low.75 <- cl(data.low[!is.na(data.low$fpc1601_75ref),],
    fm=all.fit10.fpc1601.low.75, cluster=data.low$hhid[!is.na(data.low$fpc1601_75ref)])
overall.fit10.fpc1601.low.75 <- coeftest(all.fit10.fpc1601.low.75, all.VC10.fpc1601.low.75)
summary(all.fit10.fpc1601.low.75)
overall.fit10.fpc1601.low.75
aic.fpc1601.low.75=AIC(all.fit10.fpc1601.low.75)


# fpc 1602 --------
# high risk conditions
data.high=subset(all,all$risk=="High")
all.fit10.fpc1602.high.25 <- glm(gici10~fpc1602_25ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high[!is.na(data.high$fpc1602_25ref),])

all.VC10.fpc1602.high.25 <- cl(data.high[!is.na(data.high$fpc1602_25ref),],
    fm=all.fit10.fpc1602.high.25, cluster=data.high$hhid[!is.na(data.high$fpc1602_25ref)])
overall.fit10.fpc1602.high.25 <- coeftest(all.fit10.fpc1602.high.25, all.VC10.fpc1602.high.25)
summary(all.fit10.fpc1602.high.25)
overall.fit10.fpc1602.high.25
aic.fpc1602.high.25=AIC(all.fit10.fpc1602.high.25)

all.fit10.fpc1602.high.50 <- glm(gici10~fpc1602_50ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high[!is.na(data.high$fpc1602_50ref),])

all.VC10.fpc1602.high.50 <- cl(data.high[!is.na(data.high$fpc1602_50ref),],
    fm=all.fit10.fpc1602.high.50, cluster=data.high$hhid[!is.na(data.high$fpc1602_50ref)])
overall.fit10.fpc1602.high.50 <- coeftest(all.fit10.fpc1602.high.50, all.VC10.fpc1602.high.50)
summary(all.fit10.fpc1602.high.50)
overall.fit10.fpc1602.high.50
aic.fpc1602.high.50=AIC(all.fit10.fpc1602.high.50)


all.fit10.fpc1602.high.75 <- glm(gici10~fpc1602_75ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high[!is.na(data.high$fpc1602_75ref),])

all.VC10.fpc1602.high.75 <- cl(data.high[!is.na(data.high$fpc1602_75ref),],
    fm=all.fit10.fpc1602.high.75, cluster=data.high$hhid[!is.na(data.high$fpc1602_75ref)])
overall.fit10.fpc1602.high.75 <- coeftest(all.fit10.fpc1602.high.75, all.VC10.fpc1602.high.75)
summary(all.fit10.fpc1602.high.75)
overall.fit10.fpc1602.high.75
aic.fpc1602.high.75=AIC(all.fit10.fpc1602.high.75)



# low risk conditions
data.low=subset(all,all$risk=="Low")
all.fit10.fpc1602.low.25 <- glm(gici10~fpc1602_25ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low[!is.na(data.low$fpc1602_25ref),])

all.VC10.fpc1602.low.25 <- cl(data.low[!is.na(data.low$fpc1602_25ref),],
   fm=all.fit10.fpc1602.low.25, cluster=data.low$hhid[!is.na(data.low$fpc1602_25ref)])
overall.fit10.fpc1602.low.25 <- coeftest(all.fit10.fpc1602.low.25, all.VC10.fpc1602.low.25)
summary(all.fit10.fpc1602.low.25)
overall.fit10.fpc1602.low.25
aic.fpc1602.low.25=AIC(all.fit10.fpc1602.low.25)

all.fit10.fpc1602.low.50 <- glm(gici10~fpc1602_50ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low[!is.na(data.low$fpc1602_50ref),])

all.VC10.fpc1602.low.50 <- cl(data.low[!is.na(data.low$fpc1602_50ref),],
    fm=all.fit10.fpc1602.low.50, cluster=data.low$hhid[!is.na(data.low$fpc1602_50ref)])
overall.fit10.fpc1602.low.50 <- coeftest(all.fit10.fpc1602.low.50, all.VC10.fpc1602.low.50)
summary(all.fit10.fpc1602.low.50)
overall.fit10.fpc1602.low.50
aic.fpc1602.low.50=AIC(all.fit10.fpc1602.low.50)

all.fit10.fpc1602.low.75 <- glm(gici10~fpc1602_75ref+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low[!is.na(data.low$fpc1602_75ref),])

all.VC10.fpc1602.low.75 <- cl(data.low[!is.na(data.low$fpc1602_75ref),],
    fm=all.fit10.fpc1602.low.75, cluster=data.low$hhid[!is.na(data.low$fpc1602_75ref)])
overall.fit10.fpc1602.low.75 <- coeftest(all.fit10.fpc1602.low.75, all.VC10.fpc1602.low.75)
summary(all.fit10.fpc1602.low.75)
overall.fit10.fpc1602.low.75
aic.fpc1602.low.75=AIC(all.fit10.fpc1602.low.75)

# --------------------------------------
# save the results
# exclude glm objects and data frames
# (they are really large)
# --------------------------------------
save(

  #all.n10.fmc1601,all.n10.fmc1602,all.n10.fpc1601,all.n10.fpc1602,
  
  #all.n10.fmc1602.high,all.n10.fmc1602.low,all.n10.fpc1601.high,
  #all.n10.fpc1601.low,all.n10.fpc1602.high,all.n10.fpc1602.low,
  
  all.VC10.fmc1601.25,all.VC10.fmc1601.50,all.VC10.fmc1601.75,
  all.VC10.fmc1602.25,all.VC10.fmc1602.50,all.VC10.fmc1602.75,
  all.VC10.fpc1601.25,all.VC10.fpc1601.50,all.VC10.fpc1601.75,
  all.VC10.fpc1602.25,all.VC10.fpc1602.50,all.VC10.fpc1602.75,

  overall.fit10.fmc1601.25,overall.fit10.fmc1601.50,overall.fit10.fmc1601.75,
  overall.fit10.fmc1602.25,overall.fit10.fmc1602.50,overall.fit10.fmc1602.75,
  overall.fit10.fpc1601.25,overall.fit10.fpc1601.50,overall.fit10.fpc1601.75,
  overall.fit10.fpc1602.25,overall.fit10.fpc1602.50,overall.fit10.fpc1602.75,

  all.VC10.fmc1602.high.25,all.VC10.fmc1602.high.50,all.VC10.fmc1602.high.75,
  all.VC10.fpc1601.high.25,all.VC10.fpc1601.high.50,all.VC10.fpc1601.high.75,
  all.VC10.fpc1602.high.25,all.VC10.fpc1602.high.50,all.VC10.fpc1602.high.75,

  overall.fit10.fmc1602.high.25,overall.fit10.fmc1602.high.50,overall.fit10.fmc1602.high.75,
  overall.fit10.fpc1601.high.25,overall.fit10.fpc1601.high.50,overall.fit10.fpc1601.high.75,
  overall.fit10.fpc1602.high.25,overall.fit10.fpc1602.high.50,overall.fit10.fpc1602.high.75,
  
  all.VC10.fmc1602.low.25,all.VC10.fmc1602.low.50,all.VC10.fmc1602.low.75,
  all.VC10.fpc1601.low.25,all.VC10.fpc1601.low.50,all.VC10.fpc1601.low.75,
  all.VC10.fpc1602.low.25,all.VC10.fpc1602.low.50,all.VC10.fpc1602.low.75,

  overall.fit10.fmc1602.low.25,overall.fit10.fmc1602.low.50,overall.fit10.fmc1602.low.75,
  overall.fit10.fpc1601.low.25,overall.fit10.fpc1601.low.50,overall.fit10.fpc1601.low.75,
  overall.fit10.fpc1602.low.25,overall.fit10.fpc1602.low.50,overall.fit10.fpc1602.low.75,
    
  #aic.fmc1601,aic.fmc1602,aic.fpc1601,aic.fpc1602,
  #aic.fmc1602.high,aic.fmc1602.low,
  #aic.fpc1601.high,aic.fpc1601.low,aic.fpc1602.high,aic.fpc1602.low,

  file="~/dropbox/coliphage/results/rawoutput/regress-10day-body-threshold-ref.Rdata"
)

