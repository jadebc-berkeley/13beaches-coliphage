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


# --------------------------------------
# Estimates pooled across beach
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness

# all beaches ----------------
# fmc 1601
all.fit10.fmc1601.25 <- glm(gici10~fmc1601_25+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc1601_25),])

all.VC10.fmc1601.25 <- cl(all[!is.na(all$fmc1601_25)],fm=all.fit10.fmc1601.25,
  cluster=all$hhid[!is.na(all$fmc1601_25)])
overall.fit10.fmc1601.25 <- coeftest(all.fit10.fmc1601.25, all.VC10.fmc1601.25)
summary(all.fit10.fmc1601.25)
overall.fit10.fmc1601.25
aic.fmc1601.25=AIC(all.fit10.fmc1601.25)

all.fit10.fmc1601.50 <- glm(gici10~fmc1601_50+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc1601_50),])

all.VC10.fmc1601.50 <- cl(all[!is.na(all$fmc1601_50)],fm=all.fit10.fmc1601.50,
  cluster=all$hhid[!is.na(all$fmc1601_50)])
overall.fit10.fmc1601.50 <- coeftest(all.fit10.fmc1601.50, all.VC10.fmc1601.50)
summary(all.fit10.fmc1601.50)
overall.fit10.fmc1601.50
aic.fmc1601.50=AIC(all.fit10.fmc1601.50)

all.fit10.fmc1601.75 <- glm(gici10~fmc1601_75+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc1601_75),])

all.VC10.fmc1601.75 <- cl(all[!is.na(all$fmc1601_75)],fm=all.fit10.fmc1601.75,
  cluster=all$hhid[!is.na(all$fmc1601_75)])
overall.fit10.fmc1601.75 <- coeftest(all.fit10.fmc1601.75, all.VC10.fmc1601.75)
summary(all.fit10.fmc1601.75)
overall.fit10.fmc1601.75
aic.fmc1601.75=AIC(all.fit10.fmc1601.75)


# fmc 1602
all.fit10.fmc1602.25 <- glm(gici10~fmc1602_25+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc1602_25),])

all.VC10.fmc1602.25 <- cl(all[!is.na(all$fmc1602_25)],fm=all.fit10.fmc1602.25,
  cluster=all$hhid[!is.na(all$fmc1602_25)])
overall.fit10.fmc1602.25 <- coeftest(all.fit10.fmc1602.25, all.VC10.fmc1602.25)
summary(all.fit10.fmc1602.25)
overall.fit10.fmc1602.25
aic.fmc1602.25=AIC(all.fit10.fmc1602.25)

all.fit10.fmc1602.50 <- glm(gici10~fmc1602_50+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc1602_50),])

all.VC10.fmc1602.50 <- cl(all[!is.na(all$fmc1602_50)],fm=all.fit10.fmc1602.50,
  cluster=all$hhid[!is.na(all$fmc1602_50)])
overall.fit10.fmc1602.50 <- coeftest(all.fit10.fmc1602.50, all.VC10.fmc1602.50)
summary(all.fit10.fmc1602.50)
overall.fit10.fmc1602.50
aic.fmc1602.50=AIC(all.fit10.fmc1602.50)

all.fit10.fmc1602.75 <- glm(gici10~fmc1602_75+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc1602_75),])

all.VC10.fmc1602.75 <- cl(all[!is.na(all$fmc1602_75)],fm=all.fit10.fmc1602.75,
  cluster=all$hhid[!is.na(all$fmc1602_75)])
overall.fit10.fmc1602.75 <- coeftest(all.fit10.fmc1602.75, all.VC10.fmc1602.75)
summary(all.fit10.fmc1602.75)
overall.fit10.fmc1602.75
aic.fmc1602.75=AIC(all.fit10.fmc1602.75)

# fpc 1601
all.fit10.fpc1601.25 <- glm(gici10~fpc1601_25+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1601_25),])

all.VC10.fpc1601.25 <- cl(all[!is.na(all$fpc1601_25)],fm=all.fit10.fpc1601.25,
  cluster=all$hhid[!is.na(all$fpc1601_25)])
overall.fit10.fpc1601.25 <- coeftest(all.fit10.fpc1601.25, all.VC10.fpc1601.25)
summary(all.fit10.fpc1601.25)
overall.fit10.fpc1601.25
aic.fpc1601.25=AIC(all.fit10.fpc1601.25)

all.fit10.fpc1601.50 <- glm(gici10~fpc1601_50+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1601_50),])

all.VC10.fpc1601.50 <- cl(all[!is.na(all$fpc1601_50)],fm=all.fit10.fpc1601.50,
  cluster=all$hhid[!is.na(all$fpc1601_50)])
overall.fit10.fpc1601.50 <- coeftest(all.fit10.fpc1601.50, all.VC10.fpc1601.50)
summary(all.fit10.fpc1601.50)
overall.fit10.fpc1601.50
aic.fpc1601.50=AIC(all.fit10.fpc1601.50)

all.fit10.fpc1601.75 <- glm(gici10~fpc1601_75+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1601_75),])

all.VC10.fpc1601.75 <- cl(all[!is.na(all$fpc1601_75)],fm=all.fit10.fpc1601.75,
  cluster=all$hhid[!is.na(all$fpc1601_75)])
overall.fit10.fpc1601.75 <- coeftest(all.fit10.fpc1601.75, all.VC10.fpc1601.75)
summary(all.fit10.fpc1601.75)
overall.fit10.fpc1601.75
aic.fpc1601.75=AIC(all.fit10.fpc1601.75)

# fpc 1602
all.fit10.fpc1602.25 <- glm(gici10~fpc1602_25+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1602_25),])

all.VC10.fpc1602.25 <- cl(all[!is.na(all$fpc1602_25)],fm=all.fit10.fpc1602.25,
  cluster=all$hhid[!is.na(all$fpc1602_25)])
overall.fit10.fpc1602.25 <- coeftest(all.fit10.fpc1602.25, all.VC10.fpc1602.25)
summary(all.fit10.fpc1602.25)
overall.fit10.fpc1602.25
aic.fpc1602.25=AIC(all.fit10.fpc1602.25)

all.fit10.fpc1602.50 <- glm(gici10~fpc1602_50+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1602_50),])

all.VC10.fpc1602.50 <- cl(all[!is.na(all$fpc1602_50)],fm=all.fit10.fpc1602.50,
  cluster=all$hhid[!is.na(all$fpc1602_50)])
overall.fit10.fpc1602.50 <- coeftest(all.fit10.fpc1602.50, all.VC10.fpc1602.50)
summary(all.fit10.fpc1602.50)
overall.fit10.fpc1602.50
aic.fpc1602.50=AIC(all.fit10.fpc1602.50)

all.fit10.fpc1602.75 <- glm(gici10~fpc1602_75+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1602_75),])

all.VC10.fpc1602.75 <- cl(all[!is.na(all$fpc1602_75)],fm=all.fit10.fpc1602.75,
  cluster=all$hhid[!is.na(all$fpc1602_75)])
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
all.fit10.fmc1602.high.25 <- glm(gici10~fmc1602_25+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high[!is.na(data.high$fmc1602_25),])

all.VC10.fmc1602.high.25 <- cl(data.high[!is.na(data.high$fmc1602_25),],
        fm=all.fit10.fmc1602.high.25, cluster=data.high$hhid[!is.na(data.high$fmc1602_25)])
overall.fit10.fmc1602.high.25 <- coeftest(all.fit10.fmc1602.high.25, all.VC10.fmc1602.high.25)
summary(all.fit10.fmc1602.high.25)
overall.fit10.fmc1602.high.25
aic.fmc1602.high.25=AIC(all.fit10.fmc1602.high.25)

all.fit10.fmc1602.high.50 <- glm(gici10~fmc1602_50+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high[!is.na(data.high$fmc1602_50),])

all.VC10.fmc1602.high.50 <- cl(data.high[!is.na(data.high$fmc1602_50),],fm=all.fit10.fmc1602.high.50, 
                               cluster=data.high$hhid[!is.na(data.high$fmc1602_50)])
overall.fit10.fmc1602.high.50 <- coeftest(all.fit10.fmc1602.high.50, all.VC10.fmc1602.high.50)
summary(all.fit10.fmc1602.high.50)
overall.fit10.fmc1602.high.50
aic.fmc1602.high.50=AIC(all.fit10.fmc1602.high.50)

all.fit10.fmc1602.high.75 <- glm(gici10~fmc1602_75+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high[!is.na(data.high$fmc1602_75),])

all.VC10.fmc1602.high.75 <- cl(data.high[!is.na(data.high$fmc1602_50),],fm=all.fit10.fmc1602.high.75, 
                      cluster=data.high$hhid[!is.na(data.high$fmc1602_75)])
overall.fit10.fmc1602.high.75 <- coeftest(all.fit10.fmc1602.high.75, all.VC10.fmc1602.high.75)
summary(all.fit10.fmc1602.high.75)
overall.fit10.fmc1602.high.75
aic.fmc1602.high.75=AIC(all.fit10.fmc1602.high.75)

# low risk conditions
data.low=subset(all,all$risk=="Low")
all.fit10.fmc1602.low.25 <- glm(gici10~fmc1602_25+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low[!is.na(data.low$fmc1602_25),])

all.VC10.fmc1602.low.25 <- cl(data.low[!is.na(data.low$fmc1602_25),],
     fm=all.fit10.fmc1602.low.25, cluster=data.low$hhid[!is.na(data.low$fmc1602_25)])
overall.fit10.fmc1602.low.25 <- coeftest(all.fit10.fmc1602.low.25, all.VC10.fmc1602.low.25)
summary(all.fit10.fmc1602.low.25)
overall.fit10.fmc1602.low.25
aic.fmc1602.low.25=AIC(all.fit10.fmc1602.low.25)

all.fit10.fmc1602.low.50 <- glm(gici10~fmc1602_50+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low[!is.na(data.low$fmc1602_50),])

all.VC10.fmc1602.low.50 <- cl(data.low[!is.na(data.low$fmc1602_50),],
            fm=all.fit10.fmc1602.low.50, cluster=data.low$hhid[!is.na(data.low$fmc1602_50)])
overall.fit10.fmc1602.low.50 <- coeftest(all.fit10.fmc1602.low.50, all.VC10.fmc1602.low.50)
summary(all.fit10.fmc1602.low.50)
overall.fit10.fmc1602.low.50
aic.fmc1602.low.50=AIC(all.fit10.fmc1602.low.50)

all.fit10.fmc1602.low.75 <- glm(gici10~fmc1602_75+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low[!is.na(data.low$fmc1602_75),])

all.VC10.fmc1602.low.75 <- cl(data.low[!is.na(data.low$fmc1602_75),],
            fm=all.fit10.fmc1602.low.75, cluster=data.low$hhid[!is.na(data.low$fmc1602_75)])
overall.fit10.fmc1602.low.75 <- coeftest(all.fit10.fmc1602.low.75, all.VC10.fmc1602.low.75)
summary(all.fit10.fmc1602.low.75)
overall.fit10.fmc1602.low.75
aic.fmc1602.low.75=AIC(all.fit10.fmc1602.low.75)

# fpc 1601 --------
# high risk conditions
data.high=subset(all,all$risk=="High")
all.fit10.fpc1601.high.25 <- glm(gici10~fpc1601_25+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high[!is.na(data.high$fpc1601_25),])

all.VC10.fpc1601.high.25 <- cl(data.high[!is.na(data.high$fpc1601_25),],
      fm=all.fit10.fpc1601.high.25, cluster=data.high$hhid[!is.na(data.high$fpc1601_25)])
overall.fit10.fpc1601.high.25 <- coeftest(all.fit10.fpc1601.high.25, all.VC10.fpc1601.high.25)
summary(all.fit10.fpc1601.high.25)
overall.fit10.fpc1601.high.25
aic.fpc1601.high.25=AIC(all.fit10.fpc1601.high.25)

all.fit10.fpc1601.high.50 <- glm(gici10~fpc1601_50+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high[!is.na(data.high$fpc1601_50),])

all.VC10.fpc1601.high.50 <- cl(data.high[!is.na(data.high$fpc1601_50),],
      fm=all.fit10.fpc1601.high.50, cluster=data.high$hhid[!is.na(data.high$fpc1601_50)])
overall.fit10.fpc1601.high.50 <- coeftest(all.fit10.fpc1601.high.50, all.VC10.fpc1601.high.50)
summary(all.fit10.fpc1601.high.50)
overall.fit10.fpc1601.high.50
aic.fpc1601.high.50=AIC(all.fit10.fpc1601.high.50)

all.fit10.fpc1601.high.75 <- glm(gici10~fpc1601_75+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high[!is.na(data.high$fpc1601_75),])

all.VC10.fpc1601.high.75 <- cl(data.high[!is.na(data.high$fpc1601_75),],
    fm=all.fit10.fpc1601.high.75, cluster=data.high$hhid[!is.na(data.high$fpc1601_75)])
overall.fit10.fpc1601.high.75 <- coeftest(all.fit10.fpc1601.high.75, all.VC10.fpc1601.high.75)
summary(all.fit10.fpc1601.high.75)
overall.fit10.fpc1601.high.75
aic.fpc1601.high.75=AIC(all.fit10.fpc1601.high.75)

# low risk conditions
data.low=subset(all,all$risk=="Low")
all.fit10.fpc1601.low.25 <- glm(gici10~fpc1601_25+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low[!is.na(data.low$fpc1601_25),])

all.VC10.fpc1601.low.25 <- cl(data.low[!is.na(data.low$fpc1601_25),],
    fm=all.fit10.fpc1601.low.25, cluster=data.low$hhid[!is.na(data.low$fpc1601_25)])
overall.fit10.fpc1601.low.25 <- coeftest(all.fit10.fpc1601.low.25, all.VC10.fpc1601.low.25)
summary(all.fit10.fpc1601.low.25)
overall.fit10.fpc1601.low.25
aic.fpc1601.low.25=AIC(all.fit10.fpc1601.low.25)

all.fit10.fpc1601.low.50 <- glm(gici10~fpc1601_50+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low[!is.na(data.low$fpc1601_50),])

all.VC10.fpc1601.low.50 <- cl(data.low[!is.na(data.low$fpc1601_50),],
    fm=all.fit10.fpc1601.low.50, cluster=data.low$hhid[!is.na(data.low$fpc1601_50)])
overall.fit10.fpc1601.low.50 <- coeftest(all.fit10.fpc1601.low.50, all.VC10.fpc1601.low.50)
summary(all.fit10.fpc1601.low.50)
overall.fit10.fpc1601.low.50
aic.fpc1601.low.50=AIC(all.fit10.fpc1601.low.50)

all.fit10.fpc1601.low.75 <- glm(gici10~fpc1601_75+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low[!is.na(data.low$fpc1601_75),])

all.VC10.fpc1601.low.75 <- cl(data.low[!is.na(data.low$fpc1601_75),],
    fm=all.fit10.fpc1601.low.75, cluster=data.low$hhid[!is.na(data.low$fpc1601_75)])
overall.fit10.fpc1601.low.75 <- coeftest(all.fit10.fpc1601.low.75, all.VC10.fpc1601.low.75)
summary(all.fit10.fpc1601.low.75)
overall.fit10.fpc1601.low.75
aic.fpc1601.low.75=AIC(all.fit10.fpc1601.low.75)


# fpc 1602 --------
# high risk conditions
data.high=subset(all,all$risk=="High")
all.fit10.fpc1602.high.25 <- glm(gici10~fpc1602_25+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high[!is.na(data.high$fpc1602_25),])

all.VC10.fpc1602.high.25 <- cl(data.high[!is.na(data.high$fpc1602_25),],
    fm=all.fit10.fpc1602.high.25, cluster=data.high$hhid[!is.na(data.high$fpc1602_25)])
overall.fit10.fpc1602.high.25 <- coeftest(all.fit10.fpc1602.high.25, all.VC10.fpc1602.high.25)
summary(all.fit10.fpc1602.high.25)
overall.fit10.fpc1602.high.25
aic.fpc1602.high.25=AIC(all.fit10.fpc1602.high.25)

all.fit10.fpc1602.high.50 <- glm(gici10~fpc1602_50+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high[!is.na(data.high$fpc1602_50),])

all.VC10.fpc1602.high.50 <- cl(data.high[!is.na(data.high$fpc1602_50),],
    fm=all.fit10.fpc1602.high.50, cluster=data.high$hhid[!is.na(data.high$fpc1602_50)])
overall.fit10.fpc1602.high.50 <- coeftest(all.fit10.fpc1602.high.50, all.VC10.fpc1602.high.50)
summary(all.fit10.fpc1602.high.50)
overall.fit10.fpc1602.high.50
aic.fpc1602.high.50=AIC(all.fit10.fpc1602.high.50)


all.fit10.fpc1602.high.75 <- glm(gici10~fpc1602_75+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high[!is.na(data.high$fpc1602_75),])

all.VC10.fpc1602.high.75 <- cl(data.high[!is.na(data.high$fpc1602_75),],
    fm=all.fit10.fpc1602.high.75, cluster=data.high$hhid[!is.na(data.high$fpc1602_75)])
overall.fit10.fpc1602.high.75 <- coeftest(all.fit10.fpc1602.high.75, all.VC10.fpc1602.high.75)
summary(all.fit10.fpc1602.high.75)
overall.fit10.fpc1602.high.75
aic.fpc1602.high.75=AIC(all.fit10.fpc1602.high.75)



# low risk conditions
data.low=subset(all,all$risk=="Low")
all.fit10.fpc1602.low.25 <- glm(gici10~fpc1602_25+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low[!is.na(data.low$fpc1602_25),])

all.VC10.fpc1602.low.25 <- cl(data.low[!is.na(data.low$fpc1602_25),],
   fm=all.fit10.fpc1602.low.25, cluster=data.low$hhid[!is.na(data.low$fpc1602_25)])
overall.fit10.fpc1602.low.25 <- coeftest(all.fit10.fpc1602.low.25, all.VC10.fpc1602.low.25)
summary(all.fit10.fpc1602.low.25)
overall.fit10.fpc1602.low.25
aic.fpc1602.low.25=AIC(all.fit10.fpc1602.low.25)

all.fit10.fpc1602.low.50 <- glm(gici10~fpc1602_50+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low[!is.na(data.low$fpc1602_50),])

all.VC10.fpc1602.low.50 <- cl(data.low[!is.na(data.low$fpc1602_50),],
    fm=all.fit10.fpc1602.low.50, cluster=data.low$hhid[!is.na(data.low$fpc1602_50)])
overall.fit10.fpc1602.low.50 <- coeftest(all.fit10.fpc1602.low.50, all.VC10.fpc1602.low.50)
summary(all.fit10.fpc1602.low.50)
overall.fit10.fpc1602.low.50
aic.fpc1602.low.50=AIC(all.fit10.fpc1602.low.50)

all.fit10.fpc1602.low.75 <- glm(gici10~fpc1602_75+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low[!is.na(data.low$fpc1602_75),])

all.VC10.fpc1602.low.75 <- cl(data.low[!is.na(data.low$fpc1602_75),],
    fm=all.fit10.fpc1602.low.75, cluster=data.low$hhid[!is.na(data.low$fpc1602_75)])
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

  file="~/Documents/CRG/coliphage/results/rawoutput/regress-10day-body-threshold.Rdata"
)

