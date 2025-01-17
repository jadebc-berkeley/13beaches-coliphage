##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file conducts maximum likelihood regression
# to estimate prevalence ratios

# Results pooled across beaches

# 3 day gi illness
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

# pooled n's ---------------------------------------
all.n3.fmc1601 = regN(all$gici3[!is.na(all$fmc1601.pres)],
                       all$fmc1601.pres[!is.na(all$fmc1601.pres)])
all.n3.fmc1602 = regN(all$gici3[!is.na(all$fmc1602.pres)],
                       all$fmc1602.pres[!is.na(all$fmc1602.pres)])
all.n3.fpc1601 = regN(all$gici3[!is.na(all$fpc1601.pres) ],
                       all$fpc1601.pres[!is.na(all$fpc1601.pres) ])
all.n3.fpc1602 = regN(all$gici3[!is.na(all$fpc1602.pres) ],
                       all$fpc1602.pres[!is.na(all$fpc1602.pres) ])

# pooled n's by risk level---------------------------------------
data=all[!is.na(all$fmc1602.pres) ,]
data.high=subset(data,data$groundwater=="Above median flow" | data$berm=="Open")
all.n3.fmc1602.high = regN(data.high$gici3,data.high$fmc1602.pres)
data.low=subset(data,data$groundwater=="Below median flow" | data$berm=="Closed")
all.n3.fmc1602.low = regN(data.low$gici3,data.low$fmc1602.pres)

data=all[!is.na(all$fpc1601.pres) ,]
data.high=subset(data,data$groundwater=="Above median flow" | data$berm=="Open")
all.n3.fpc1601.high = regN(data.high$gici3,data.high$fpc1601.pres)
data.low=subset(data,data$groundwater=="Below median flow" | data$berm=="Closed")
all.n3.fpc1601.low = regN(data.low$gici3,data.low$fpc1601.pres)

data=all[!is.na(all$fpc1602.pres) ,]
data.high=subset(data,data$groundwater=="Above median flow" | data$berm=="Open")
all.n3.fpc1602.high = regN(data.high$gici3,data.high$fpc1602.pres)
data.low=subset(data,data$groundwater=="Below median flow" | data$berm=="Closed")
all.n3.fpc1602.low = regN(data.low$gici3,data.low$fpc1602.pres)


# --------------------------------------
# Estimates pooled across beach
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness

# all beaches ----------------
# fmc 1601
all.fit3.fmc1601 <- glm(gici3~fmc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc1601.pres) ,])

all.VC3.fmc1601 <- cl(all[!is.na(all$fmc1601.pres) ],fm=all.fit3.fmc1601,
  cluster=all$hhid[!is.na(all$fmc1601.pres) ])
overall.fit3.fmc1601 <- coeftest(all.fit3.fmc1601, all.VC3.fmc1601)
summary(all.fit3.fmc1601)
overall.fit3.fmc1601


# fmc 1602
all.fit3.fmc1602 <- glm(gici3~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc1602.pres),])

all.VC3.fmc1602 <- cl(all[!is.na(all$fmc1602.pres),],fm=all.fit3.fmc1602,
  cluster=all$hhid[!is.na(all$fmc1602.pres)])
overall.fit3.fmc1602 <- coeftest(all.fit3.fmc1602, all.VC3.fmc1602)
summary(all.fit3.fmc1602)
overall.fit3.fmc1602

# fpc 1601
all.fit3.fpc1601 <- glm(gici3~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1601.pres) ,])

all.VC3.fpc1601 <- cl(all[!is.na(all$fpc1601.pres) ,],fm=all.fit3.fpc1601,
    cluster=all$hhid[!is.na(all$fpc1601.pres) ])
overall.fit3.fpc1601 <- coeftest(all.fit3.fpc1601, all.VC3.fpc1601)
summary(all.fit3.fpc1601)
overall.fit3.fpc1601

# fpc 1602
all.fit3.fpc1602 <- glm(gici3~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1602.pres) ,])

all.VC3.fpc1602 <- cl(all[!is.na(all$fpc1602.pres) ,],fm=all.fit3.fpc1602,
    cluster=all$hhid[!is.na(all$fpc1602.pres) ])
overall.fit3.fpc1602 <- coeftest(all.fit3.fpc1602, all.VC3.fpc1602)
summary(all.fit3.fpc1602)
overall.fit3.fpc1602

# --------------------------------------
# Estimates pooled across beach and stratified by conditions
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness

# all beaches ----------------
# fmc 1602 --------
# high risk conditions
data=all[!is.na(all$fmc1602.pres) ,]
data.high=subset(data,data$groundwater=="Above median flow" | data$berm=="Open")
all.fit3.fmc1602.high <- glm(gici3~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC3.fmc1602.high <- cl(data.high,fm=all.fit3.fmc1602.high, cluster=data.high$hhid)
overall.fit3.fmc1602.high <- coeftest(all.fit3.fmc1602.high, all.VC3.fmc1602.high)
summary(all.fit3.fmc1602.high)
overall.fit3.fmc1602.high

# low risk conditions
data.low=subset(data,data$groundwater=="Below median flow" | data$berm=="Closed")
all.fit3.fmc1602.low <- glm(gici3~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC3.fmc1602.low <- cl(data.low,fm=all.fit3.fmc1602.low, cluster=data.low$hhid)
overall.fit3.fmc1602.low <- coeftest(all.fit3.fmc1602.low, all.VC3.fmc1602.low)
summary(all.fit3.fmc1602.low)
overall.fit3.fmc1602.low

# fpc 1601 --------
# high risk conditions
data=all[!is.na(all$fpc1601.pres) ,]
data.high=subset(data,data$groundwater=="Above median flow" | data$berm=="Open")
all.fit3.fpc1601.high <- glm(gici3~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC3.fpc1601.high <- cl(data.high,fm=all.fit3.fpc1601.high, cluster=data.high$hhid)
overall.fit3.fpc1601.high <- coeftest(all.fit3.fpc1601.high, all.VC3.fpc1601.high)
summary(all.fit3.fpc1601.high)
overall.fit3.fpc1601.high

# low risk conditions
data.low=subset(data,data$groundwater=="Below median flow" | data$berm=="Closed")
all.fit3.fpc1601.low <- glm(gici3~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC3.fpc1601.low <- cl(data.low,fm=all.fit3.fpc1601.low, cluster=data.low$hhid)
overall.fit3.fpc1601.low <- coeftest(all.fit3.fpc1601.low, all.VC3.fpc1601.low)
summary(all.fit3.fpc1601.low)
overall.fit3.fpc1601.low


# fpc 1602 --------
# high risk conditions
data=all[!is.na(all$fpc1602.pres) ,]
data.high=subset(data,data$groundwater=="Above median flow" | data$berm=="Open")
all.fit3.fpc1602.high <- glm(gici3~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC3.fpc1602.high <- cl(data.high,fm=all.fit3.fpc1602.high, cluster=data.high$hhid)
overall.fit3.fpc1602.high <- coeftest(all.fit3.fpc1602.high, all.VC3.fpc1602.high)
summary(all.fit3.fpc1602.high)
overall.fit3.fpc1602.high

# low risk conditions
data.low=subset(data,data$groundwater=="Below median flow" | data$berm=="Closed")
all.fit3.fpc1602.low <- glm(gici3~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC3.fpc1602.low <- cl(data.low,fm=all.fit3.fpc1602.low, cluster=data.low$hhid)
overall.fit3.fpc1602.low <- coeftest(all.fit3.fpc1602.low, all.VC3.fpc1602.low)
summary(all.fit3.fpc1602.low)
overall.fit3.fpc1602.low

# --------------------------------------
# save the results
# exclude glm objects and data frames
# (they are really large)
# --------------------------------------
save(

  all.n3.fmc1601,all.n3.fmc1602,all.n3.fpc1601,all.n3.fpc1602,
  
  all.n3.fmc1602.high,all.n3.fmc1602.low,all.n3.fpc1601.high,
  all.n3.fpc1601.low,all.n3.fpc1602.high,all.n3.fpc1602.low,
  
  all.VC3.fmc1601,all.VC3.fmc1602,all.VC3.fpc1601,all.VC3.fpc1602,
  overall.fit3.fmc1601,overall.fit3.fmc1602,overall.fit3.fpc1601,
  overall.fit3.fpc1602,

  all.VC3.fmc1602.high,all.VC3.fpc1601.high,all.VC3.fpc1602.high,
  overall.fit3.fmc1602.high,overall.fit3.fpc1601.high,
  overall.fit3.fpc1602.high,
  
  all.VC3.fmc1602.low,all.VC3.fpc1601.low,all.VC3.fpc1602.low,
  overall.fit3.fmc1602.low,overall.fit3.fpc1601.low,
  overall.fit3.fpc1602.low,

  
  file="~/Documents/CRG/coliphage/results/rawoutput/regress-3day-body.Rdata"
)

