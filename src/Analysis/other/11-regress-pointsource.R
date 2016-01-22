##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# Assess potential effect modification by 
# point source for FPC 1601 (only indicator
# with sufficient beaches of different
# point source status)

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
# assess effect modification by point source
# --------------------------------------

# f- coliphage ----------------
all.fit10.fmc <- glm(gici10~fmc.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
                       rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc.pres),])

all.VC10.fmc <- cl(all[!is.na(all$fmc.pres)],fm=all.fit10.fmc,
                   cluster=all$hhid[!is.na(all$fmc.pres)])
overall.fit10.fmc <- coeftest(all.fit10.fmc, all.VC10.fmc)

# Interaction model with point v. non-point source beaches
ps.fit10.fmc <- glm(gici10~fmc.pres*pointsource+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc.pres) ,])

summary(ps.fit10.fmc)
lr.fmc=lrtest(all.fit10.fmc,ps.fit10.fmc)
lr.fmc

# f+ coliphage ----------------
all.fit10.fpc <- glm(gici10~fpc.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
                       rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc.pres),])

all.VC10.fpc <- cl(all[!is.na(all$fpc.pres)],fm=all.fit10.fpc,
                   cluster=all$hhid[!is.na(all$fpc.pres)])
overall.fit10.fpc <- coeftest(all.fit10.fpc, all.VC10.fpc)
summary(all.fit10.fpc)
overall.fit10.fpc
aic.fpc=AIC(all.fit10.fpc)

# Interaction model with point v. non-point source beaches
ps.fit10.fpc <- glm(gici10~fpc.pres*pointsource+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc.pres) ,])

summary(ps.fit10.fpc1601)
lr.fpc=lrtest(all.fit10.fpc,ps.fit10.fpc)
lr.fpc

#################################################
# stratify by conditions
#################################################

# f- coliphage --------
# high risk conditions
data=all[!is.na(all$fmc.pres),]
data.high=subset(data,data$risk=="High")
all.fit10.fmc.high <- glm(gici10~fmc.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC10.fmc.high <- cl(data.high,fm=all.fit10.fmc.high, cluster=data.high$hhid)
overall.fit10.fmc.high <- coeftest(all.fit10.fmc.high, all.VC10.fmc.high)
summary(all.fit10.fmc.high)

# Interaction model with point v. non-point source beaches
ps.fit10.fmc.high <- glm(gici10~fmc.pres*pointsource+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+rawfood+beach,family=poisson(link="log"),data=data.high)

summary(ps.fit10.fmc.high)
lr.fmc.high=lrtest(all.fit10.fmc.high,ps.fit10.fmc.high)
lr.fmc.high

# low risk conditions
data.low=subset(data,data$risk=="Low")
all.fit10.fmc.low <- glm(gici10~fmc.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC10.fmc.low <- cl(data.low,fm=all.fit10.fmc.low, cluster=data.low$hhid)
overall.fit10.fmc.low <- coeftest(all.fit10.fmc.low, all.VC10.fmc.low)
summary(all.fit10.fmc.low)

# Interaction model with point v. non-point source beaches
ps.fit10.fmc.low <- glm(gici10~fmc.pres*pointsource+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+rawfood+beach,family=poisson(link="log"),data=data.low)

summary(ps.fit10.fmc.low)
lr.fmc.low=lrtest(all.fit10.fmc.low,ps.fit10.fmc.low)
lr.fmc.low

# f+ coliphage  --------
# high risk conditions
data=all[!is.na(all$fpc.pres),]
data.high=subset(data,data$risk=="High")
all.fit10.fpc.high <- glm(gici10~fpc.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC10.fpc.high <- cl(data.high,fm=all.fit10.fpc.high, cluster=data.high$hhid)
overall.fit10.fpc.high <- coeftest(all.fit10.fpc.high, all.VC10.fpc.high)
summary(all.fit10.fpc.high)

# Interaction model with point v. non-point source beaches
ps.fit10.fpc.high <- glm(gici10~fpc.pres*pointsource+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+rawfood+beach,family=poisson(link="log"),data=data.high)

summary(ps.fit10.fpc.high)
lr.fpc.high=lrtest(all.fit10.fpc.high,ps.fit10.fpc.high)
lr.fpc.high

# low risk conditions
data.low=subset(data,data$risk=="Low")
all.fit10.fpc.low <- glm(gici10~fpc.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC10.fpc.low <- cl(data.low,fm=all.fit10.fpc.low, cluster=data.low$hhid)
overall.fit10.fpc.low <- coeftest(all.fit10.fpc.low, all.VC10.fpc.low)
summary(all.fit10.fpc.low)

# Interaction model with point v. non-point source beaches
ps.fit10.fpc.low <- glm(gici10~fpc.pres*pointsource+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+rawfood+beach,family=poisson(link="log"),data=data.low)

summary(ps.fit10.fpc.low)
lr.fpc.low=lrtest(all.fit10.fpc.low,ps.fit10.fpc.low)
lr.fpc.low



# --------------------------------------
# save the results
# exclude glm objects and data frames
# (they are really large)
# --------------------------------------
save(
  
  lr.fmc, lr.fpc, lr.fmc.low,lr.fmc.high,
  lr.fpc.low,lr.fpc.high,
  
  file="~/Documents/CRG/coliphage/results/rawoutput/regress-10day-ps-em.Rdata"
)

