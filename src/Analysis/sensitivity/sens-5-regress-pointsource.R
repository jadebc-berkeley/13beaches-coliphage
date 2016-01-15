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
all=subset(all,all$bodycontact=="Yes")


# --------------------------------------
# assess effect modification by point source
# --------------------------------------
# fpc 1601
all.fit10.fpc1601 <- glm(gici10~fpc1601.pres+pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1601.pres) ,])

all.VC10.fpc1601 <- cl(all[!is.na(all$fpc1601.pres) ,],fm=all.fit10.fpc1601,
    cluster=all$hhid[!is.na(all$fpc1601.pres) ])
overall.fit10.fpc1601 <- coeftest(all.fit10.fpc1601, all.VC10.fpc1601)
summary(all.fit10.fpc1601)
overall.fit10.fpc1601

# Interaction model with point v. non-point source beaches
ps.fit10.fpc1601 <- glm(gici10~fpc1601.pres*pointsource+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1601.pres) ,])

summary(ps.fit10.fpc1601)
lrtest(all.fit10.fpc1601,ps.fit10.fpc1601)

