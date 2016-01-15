##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file conducts maximum likelihood regression
# to estimate prevalence ratios

# negative control analysis

# Results pooled across beaches

# using new exposure definition

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
beaches13=read.csv("~/Dropbox/13beaches/data/final/13beaches-analysis.csv")

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
all=subset(all,all$anycontact=="No")

# --------------------------------------
# Estimates pooled across beach
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness

# all beaches ----------------
# fmc 1601
nc.fit10.fmc1601 <- glm(gici10~fmc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc1601.pres) ,])

nc.VC10.fmc1601 <- cl(all[!is.na(all$fmc1601.pres) ],fm=nc.fit10.fmc1601,
  cluster=all$hhid[!is.na(all$fmc1601.pres) ])
nc.overall.fit10.fmc1601 <- coeftest(nc.fit10.fmc1601, nc.VC10.fmc1601)
summary(nc.fit10.fmc1601)
nc.overall.fit10.fmc1601
aic.nc.fmc1601=AIC(nc.fit10.fmc1601)

# fmc 1602
nc.fit10.fmc1602 <- glm(gici10~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc1602.pres),])

nc.VC10.fmc1602 <- cl(all[!is.na(all$fmc1602.pres),],fm=nc.fit10.fmc1602,
  cluster=all$hhid[!is.na(all$fmc1602.pres)])
nc.overall.fit10.fmc1602 <- coeftest(nc.fit10.fmc1602, nc.VC10.fmc1602)
summary(nc.fit10.fmc1602)
nc.overall.fit10.fmc1602
aic.nc.fmc1602=AIC(nc.fit10.fmc1602)

# fpc 1601
nc.fit10.fpc1601 <- glm(gici10~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1601.pres) ,])

nc.VC10.fpc1601 <- cl(all[!is.na(all$fpc1601.pres) ,],fm=nc.fit10.fpc1601,
    cluster=all$hhid[!is.na(all$fpc1601.pres) ])
nc.overall.fit10.fpc1601 <- coeftest(nc.fit10.fpc1601, nc.VC10.fpc1601)
summary(nc.fit10.fpc1601)
nc.overall.fit10.fpc1601
aic.nc.fpc1601=AIC(nc.fit10.fpc1601)

# fpc 1602
nc.fit10.fpc1602 <- glm(gici10~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1602.pres) ,])

nc.VC10.fpc1602 <- cl(all[!is.na(all$fpc1602.pres) ,],fm=nc.fit10.fpc1602,
    cluster=all$hhid[!is.na(all$fpc1602.pres) ])
nc.overall.fit10.fpc1602 <- coeftest(nc.fit10.fpc1602, nc.VC10.fpc1602)
summary(nc.fit10.fpc1602)
nc.overall.fit10.fpc1602
aic.nc.fpc1602=AIC(nc.fit10.fpc1602)


# --------------------------------------
# Estimates pooled across beach and stratified by conditions
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness

# all beaches ----------------
# fmc 1602 --------
# high risk conditions
data=all[!is.na(all$fmc1602.pres),]
data.high=subset(data,data$risk=="High")
nc.fit10.fmc1602.high <- glm(gici10~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high)

nc.VC10.fmc1602.high <- cl(data.high,fm=nc.fit10.fmc1602.high, cluster=data.high$hhid)
nc.overall.fit10.fmc1602.high <- coeftest(nc.fit10.fmc1602.high, nc.VC10.fmc1602.high)
summary(nc.fit10.fmc1602.high)
nc.overall.fit10.fmc1602.high
aic.nc.fmc1602.high=AIC(nc.fit10.fmc1602.high)

# low risk conditions
data.low=subset(data,data$risk=="Low")
nc.fit10.fmc1602.low <- glm(gici10~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low)

nc.VC10.fmc1602.low <- cl(data.low,fm=nc.fit10.fmc1602.low, cluster=data.low$hhid)
nc.overall.fit10.fmc1602.low <- coeftest(nc.fit10.fmc1602.low, nc.VC10.fmc1602.low)
summary(nc.fit10.fmc1602.low)
nc.overall.fit10.fmc1602.low
aic.nc.fmc1602.low=AIC(nc.fit10.fmc1602.low)

# fpc 1601 --------
# high risk conditions
data=all[!is.na(all$fpc1601.pres) ,]
data.high=subset(data,data$risk=="High")
nc.fit10.fpc1601.high <- glm(gici10~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high)

nc.VC10.fpc1601.high <- cl(data.high,fm=nc.fit10.fpc1601.high, cluster=data.high$hhid)
nc.overall.fit10.fpc1601.high <- coeftest(nc.fit10.fpc1601.high, nc.VC10.fpc1601.high)
summary(nc.fit10.fpc1601.high)
nc.overall.fit10.fpc1601.high
aic.nc.fpc1601.high=AIC(nc.fit10.fpc1601.high)

# low risk conditions
data.low=subset(data,data$risk=="Low")
nc.fit10.fpc1601.low <- glm(gici10~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low)

nc.VC10.fpc1601.low <- cl(data.low,fm=nc.fit10.fpc1601.low, cluster=data.low$hhid)
nc.overall.fit10.fpc1601.low <- coeftest(nc.fit10.fpc1601.low, nc.VC10.fpc1601.low)
summary(nc.fit10.fpc1601.low)
nc.overall.fit10.fpc1601.low
aic.nc.fpc1601.low=AIC(nc.fit10.fpc1601.low)


# fpc 1602 --------
# high risk conditions
data=all[!is.na(all$fpc1602.pres) ,]
data.high=subset(data,data$risk=="High")
nc.fit10.fpc1602.high <- glm(gici10~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high)

nc.VC10.fpc1602.high <- cl(data.high,fm=nc.fit10.fpc1602.high, cluster=data.high$hhid)
nc.overall.fit10.fpc1602.high <- coeftest(nc.fit10.fpc1602.high, nc.VC10.fpc1602.high)
summary(nc.fit10.fpc1602.high)
nc.overall.fit10.fpc1602.high
aic.nc.fpc1602.high=AIC(nc.fit10.fpc1602.high)


# low risk conditions
data.low=subset(data,data$risk=="Low")
nc.fit10.fpc1602.low <- glm(gici10~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low)

nc.VC10.fpc1602.low <- cl(data.low,fm=nc.fit10.fpc1602.low, cluster=data.low$hhid)
nc.overall.fit10.fpc1602.low <- coeftest(nc.fit10.fpc1602.low, nc.VC10.fpc1602.low)
summary(nc.fit10.fpc1602.low)
nc.overall.fit10.fpc1602.low
aic.nc.fpc1602.low=AIC(nc.fit10.fpc1602.low)



# --------------------------------------
# save the results
# exclude glm objects and data frames
# (they are really large)
# --------------------------------------
save(

  nc.VC10.fmc1601,nc.VC10.fmc1602,nc.VC10.fpc1601,nc.VC10.fpc1602,
  nc.overall.fit10.fmc1601,nc.overall.fit10.fmc1602,nc.overall.fit10.fpc1601,
  nc.overall.fit10.fpc1602,

  nc.VC10.fmc1602.high,nc.VC10.fpc1601.high,nc.VC10.fpc1602.high,
  nc.overall.fit10.fmc1602.high,nc.overall.fit10.fpc1601.high,
  nc.overall.fit10.fpc1602.high,
  
  nc.VC10.fmc1602.low,nc.VC10.fpc1601.low,nc.VC10.fpc1602.low,
  nc.overall.fit10.fmc1602.low,nc.overall.fit10.fpc1601.low,
  nc.overall.fit10.fpc1602.low,
  
  aic.nc.fmc1601,aic.nc.fmc1602,aic.nc.fpc1601,aic.nc.fpc1602,
  aic.nc.fmc1602.high,aic.nc.fmc1602.low,
  aic.nc.fpc1601.high,aic.nc.fpc1601.low,aic.nc.fpc1602.high,aic.nc.fpc1602.low,

  file="~/dropbox/coliphage/results/rawoutput/regress-10day-body-negcontrol.Rdata"
)

