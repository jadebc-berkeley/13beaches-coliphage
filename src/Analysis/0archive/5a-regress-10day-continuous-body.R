##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# Estimate the RR across the range of concentration

# Results pooled across beaches

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
all=subset(all,nowq==0)
all=subset(all,all$bodycontact=="Yes")


#-------------------------------------------
# estimating CIR for continuous exposure
#-------------------------------------------
# fmc 1601 ---------------------------------
all.fit10.fmc1601 <- glm(gici10~fmc1601+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc1601),])

all.VC10.fmc1601 <- cl(all[!is.na(all$fmc1601),],fm=all.fit10.fmc1601,
    cluster=all$hhid[!is.na(all$fmc1601)])
overall.fit10.fmc1601 <- coeftest(all.fit10.fmc1601, all.VC10.fmc1601)

# fmc 1602 ---------------------------------
all.fit10.fmc1602 <- glm(gici10~fmc1602+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc1602),])

all.VC10.fmc1602 <- cl(all[!is.na(all$fmc1602),],fm=all.fit10.fmc1602,
    cluster=all$hhid[!is.na(all$fmc1602)])
overall.fit10.fmc1602 <- coeftest(all.fit10.fmc1602, all.VC10.fmc1602)

# high risk conditions
data=all[!is.na(all$fmc1602),]
data.high=subset(data,data$risk=="High")
all.fit10.fmc1602.high <- glm(gici10~fmc1602+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC10.fmc1602.high <- cl(data.high,fm=all.fit10.fmc1602.high, cluster=data.high$hhid)
overall.fit10.fmc1602.high <- coeftest(all.fit10.fmc1602.high, all.VC10.fmc1602.high)

# low risk conditions
data.low=subset(data,data$risk=="Low")
all.fit10.fmc1602.low <- glm(gici10~fmc1602+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC10.fmc1602.low <- cl(data.low,fm=all.fit10.fmc1602.low, cluster=data.low$hhid)
overall.fit10.fmc1602.low <- coeftest(all.fit10.fmc1602.low, all.VC10.fmc1602.low)

# fpc 1601 ---------------------------------
all.fit10.fpc1601 <- glm(gici10~fpc1601+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1601),])

all.VC10.fpc1601 <- cl(all[!is.na(all$fpc1601),],fm=all.fit10.fpc1601,
    cluster=all$hhid[!is.na(all$fpc1601)])
overall.fit10.fpc1601 <- coeftest(all.fit10.fpc1601, all.VC10.fpc1601)

# high risk conditions
data=all[!is.na(all$fpc1601),]
data.high=subset(data,data$risk=="High")
all.fit10.fpc1601.high <- glm(gici10~fpc1601+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC10.fpc1601.high <- cl(data.high,fm=all.fit10.fpc1601.high, cluster=data.high$hhid)
overall.fit10.fpc1601.high <- coeftest(all.fit10.fpc1601.high, all.VC10.fpc1601.high)

# low risk conditions
data.low=subset(data,data$risk=="Low")
all.fit10.fpc1601.low <- glm(gici10~fpc1601+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC10.fpc1601.low <- cl(data.low,fm=all.fit10.fpc1601.low, cluster=data.low$hhid)
overall.fit10.fpc1601.low <- coeftest(all.fit10.fpc1601.low, all.VC10.fpc1601.low)

# fpc 1602 ---------------------------------
all.fit10.fpc1602 <- glm(gici10~fpc1602+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1602),])

all.VC10.fpc1602 <- cl(all[!is.na(all$fpc1602),],fm=all.fit10.fpc1602,
    cluster=all$hhid[!is.na(all$fpc1602)])
overall.fit10.fpc1602 <- coeftest(all.fit10.fpc1602, all.VC10.fpc1602)

# high risk conditions
data=all[!is.na(all$fpc1602),]
data.high=subset(data,data$risk=="High")
all.fit10.fpc1602.high <- glm(gici10~fpc1602+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC10.fpc1602.high <- cl(data.high,fm=all.fit10.fpc1602.high, cluster=data.high$hhid)
overall.fit10.fpc1602.high <- coeftest(all.fit10.fpc1602.high, all.VC10.fpc1602.high)

# low risk conditions
data.low=subset(data,data$risk=="Low")
all.fit10.fpc1602.low <- glm(gici10~fpc1602+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC10.fpc1602.low <- cl(data.low,fm=all.fit10.fpc1602.low, cluster=data.low$hhid)
overall.fit10.fpc1602.low <- coeftest(all.fit10.fpc1602.low, all.VC10.fpc1602.low)

#-------------------------------------------
# estimating probability of illness over
# observed range of exposure
#-------------------------------------------
iter=1000

set.seed(92203789)

# fmc 1601 ---------------------------------
all.fmc1601.pY = boot.pY(fmla=gici10~fmc1601+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+rawfood+beach,dat=all[!is.na(all$fmc1601),],nameX="fmc1601",
    ID=all[!is.na(all$fmc1601),"hhid"],iter)

# fmc 1602 ---------------------------------
all.fmc1602.pY = boot.pY(fmla=gici10~fmc1602+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+rawfood+beach,dat=all[!is.na(all$fmc1602),],nameX="fmc1602",
    ID=all[!is.na(all$fmc1602),"hhid"],iter)
data=all[!is.na(all$fmc1602),]
data.high=subset(data,data$risk=="High")
all.fmc1602.pY.high = boot.pY(fmla=gici10~fmc1602+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+ rawfood+beach,dat=data.high,nameX="fmc1602",ID=data.high$hhid,iter)
data.low=subset(data,data$risk=="Low")
all.fmc1602.pY.low = boot.pY(fmla=gici10~fmc1602+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+ rawfood+beach,dat=data.low,nameX="fmc1602",ID=data.low$hhid,iter)

# fpc 1601 ---------------------------------
all.fpc1601.pY = boot.pY(fmla=gici10~fpc1601+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+rawfood+beach,dat=all[!is.na(all$fpc1601),],nameX="fpc1601",
    ID=all[!is.na(all$fpc1601),"hhid"],iter)
data=all[!is.na(all$fpc1601),]
data.high=subset(data,data$risk=="High")
all.fpc1601.pY.high = boot.pY(fmla=gici10~fpc1601+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+ rawfood+beach,dat=data.high,nameX="fpc1601",ID=data.high$hhid,iter)
data.low=subset(data,data$risk=="Low")
all.fpc1601.pY.low = boot.pY(fmla=gici10~fpc1601+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+ rawfood+beach,dat=data.low,nameX="fpc1601",ID=data.low$hhid,iter)

# fpc 1602 ---------------------------------
all.fpc1602.pY = boot.pY(fmla=gici10~fpc1602+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+rawfood+beach,dat=all[!is.na(all$fpc1602),],nameX="fpc1602",
    ID=all[!is.na(all$fpc1602),"hhid"],iter)
data=all[!is.na(all$fpc1602),]
data.high=subset(data,data$risk=="High")
all.fpc1602.pY.high = boot.pY(fmla=gici10~fpc1602+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+ rawfood+beach,dat=data.high,nameX="fpc1602",ID=data.high$hhid,iter)
data.low=subset(data,data$risk=="Low")
all.fpc1602.pY.low = boot.pY(fmla=gici10~fpc1602+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+ rawfood+beach,dat=data.low,nameX="fpc1602",ID=data.low$hhid,iter)


save(
  overall.fit10.fmc1601, overall.fit10.fmc1602, overall.fit10.fmc1602.high,
  overall.fit10.fmc1602.low,overall.fit10.fpc1601,overall.fit10.fpc1601.high,
  overall.fit10.fpc1601.low,overall.fit10.fpc1602,overall.fit10.fpc1602.high,
  overall.fit10.fpc1602.low,

  all.fmc1601.pY,all.fmc1602.pY,all.fmc1602.pY.high,all.fmc1602.pY.low,
  all.fpc1601.pY,all.fpc1601.pY.high,all.fpc1601.pY.low,all.fpc1602.pY,
  all.fpc1602.pY.high,all.fpc1602.pY.low,
  
  file="~/dropbox/coliphage/results/rawoutput/regress-10day-continuous-body.Rdata"
  
)












