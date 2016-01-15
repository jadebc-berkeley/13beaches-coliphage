##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file conducts maximum likelihood regression
# to estimate prevalence ratios
# for enterococcus

# Results pooled across beaches

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

# n's pooled by beach ---------------------------------------
n10.entero35.fmc1601 = regN(all$gici10[!is.na(all$fmc1601.pres)],
  all$entero35[!is.na(all$fmc1601.pres)])

n10.entero35.fmc1602 = regN(all$gici10[!is.na(all$fmc1602.pres)],
   all$entero35[!is.na(all$fmc1602.pres)])

n10.entero35.fpc1601 = regN(all$gici10[!is.na(all$fpc1601.pres)],
  all$entero35[!is.na(all$fpc1601.pres)])

n10.entero35.fpc1602 = regN(all$gici10[!is.na(all$fpc1602.pres)],
  all$entero35[!is.na(all$fpc1602.pres)])


# pooled n's by risk level---------------------------------------
data=all[!is.na(all$fmc1602.pres),]
data.high=subset(data,data$risk=="High")
n10.entero35.fmc1602.high = regN(data.high$gici10,data.high$entero35)
data.low=subset(data,data$risk=="Low")
n10.entero35.fmc1602.low = regN(data.low$gici10,data.low$entero35)

data=all[!is.na(all$fpc1601.pres),]
data.high=subset(data,data$risk=="High")
n10.entero35.fpc1601.high = regN(data.high$gici10,data.high$entero35)
data.low=subset(data,data$risk=="Low")
n10.entero35.fpc1601.low = regN(data.low$gici10,data.low$entero35)

data=all[!is.na(all$fpc1602.pres),]
data.high=subset(data,data$risk=="High")
n10.entero35.fpc1602.high = regN(data.high$gici10,data.high$entero35)
data.low=subset(data,data$risk=="Low")
n10.entero35.fpc1602.low = regN(data.low$gici10,data.low$entero35)


# --------------------------------------
# Estimates pooled across beach
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness
# --------------------------------------

# all beaches --------------------------------
all.fit10.entero.fmc1601 <- glm(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$entero35) & 
    !is.na(all$fmc1601.pres),])

all.VC10.entero.fmc1601 <- cl(all[!is.na(all$entero35) & !is.na(all$fmc1601.pres)],fm=all.fit10.entero.fmc1601, cluster=
      all$hhid[!is.na(all$entero35)  & !is.na(all$fmc1601.pres)])
overall.fit10.entero.fmc1601 <- coeftest(all.fit10.entero.fmc1601, all.VC10.entero.fmc1601)
summary(all.fit10.entero.fmc1601)
overall.fit10.entero.fmc1601
aic.entero.fmc1601=AIC(all.fit10.entero.fmc1601)

all.fit10.entero.fmc1602 <- glm(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$entero35) & 
    !is.na(all$fmc1602.pres),])

all.VC10.entero.fmc1602 <- cl(all[!is.na(all$entero35) & !is.na(all$fmc1602.pres)],fm=all.fit10.entero.fmc1602, cluster=
      all$hhid[!is.na(all$entero35)  & !is.na(all$fmc1602.pres)])
overall.fit10.entero.fmc1602 <- coeftest(all.fit10.entero.fmc1602, all.VC10.entero.fmc1602)
summary(all.fit10.entero.fmc1602)
overall.fit10.entero.fmc1602
aic.entero.fmc1602=AIC(all.fit10.entero.fmc1602)

all.fit10.entero.fpc1601 <- glm(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$entero35) & 
    !is.na(all$fpc1601.pres),])

all.VC10.entero.fpc1601 <- cl(all[!is.na(all$entero35) & !is.na(all$fpc1601.pres)],fm=all.fit10.entero.fpc1601, cluster=
      all$hhid[!is.na(all$entero35)  & !is.na(all$fpc1601.pres)])
overall.fit10.entero.fpc1601 <- coeftest(all.fit10.entero.fpc1601, all.VC10.entero.fpc1601)
summary(all.fit10.entero.fpc1601)
overall.fit10.entero.fpc1601
aic.entero.fpc1601=AIC(all.fit10.entero.fpc1601)

all.fit10.entero.fpc1602 <- glm(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$entero35) & 
    !is.na(all$fpc1602.pres),])

all.VC10.entero.fpc1602 <- cl(all[!is.na(all$entero35) & !is.na(all$fpc1602.pres)],fm=all.fit10.entero.fpc1602, cluster=
      all$hhid[!is.na(all$entero35)  & !is.na(all$fpc1602.pres)])
overall.fit10.entero.fpc1602 <- coeftest(all.fit10.entero.fpc1602, all.VC10.entero.fpc1602)
summary(all.fit10.entero.fpc1602)
overall.fit10.entero.fpc1602
aic.entero.fpc1602=AIC(all.fit10.entero.fpc1602)

# FMC 1602 #####################
# high risk conditions --------------------------------
data=all[!is.na(all$entero35) & !is.na(all$fmc1602.pres),]
data.high=subset(data,data$risk=="High")
all.fit10.entero.high.fmc1602 <- glm(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC10.entero.high.fmc1602 <- cl(data.high,fm=all.fit10.entero.high.fmc1602, cluster=data.high$hhid)
overall.fit10.entero.high.fmc1602 <- coeftest(all.fit10.entero.high.fmc1602, all.VC10.entero.high.fmc1602)
summary(all.fit10.entero.high.fmc1602)
overall.fit10.entero.high.fmc1602
aic.entero.high.fmc1602=AIC(all.fit10.entero.high.fmc1602)

# low risk conditions --------------------------------
data.low=subset(data,data$risk=="Low")
all.fit10.entero.low.fmc1602 <- glm(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC10.entero.low.fmc1602 <- cl(data.low,fm=all.fit10.entero.low.fmc1602, cluster=data.low$hhid)
overall.fit10.entero.low.fmc1602 <- coeftest(all.fit10.entero.low.fmc1602, all.VC10.entero.low.fmc1602)
summary(all.fit10.entero.low.fmc1602)
overall.fit10.entero.low.fmc1602
aic.entero.low.fmc1602=AIC(all.fit10.entero.low.fmc1602)

# FPC 1601 #####################
# high risk conditions --------------------------------
data=all[!is.na(all$entero35) & !is.na(all$fpc1601.pres),]
data.high=subset(data,data$risk=="High")
all.fit10.entero.high.fpc1601 <- glm(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+
                                rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC10.entero.high.fpc1601 <- cl(data.high,fm=all.fit10.entero.high.fpc1601, cluster=data.high$hhid)
overall.fit10.entero.high.fpc1601 <- coeftest(all.fit10.entero.high.fpc1601, all.VC10.entero.high.fpc1601)
summary(all.fit10.entero.high.fpc1601)
overall.fit10.entero.high.fpc1601
aic.entero.high.fpc1601=AIC(all.fit10.entero.high.fpc1601)

# low risk conditions --------------------------------
data.low=subset(data,data$risk=="Low")
all.fit10.entero.low.fpc1601 <- glm(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+
                               rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC10.entero.low.fpc1601 <- cl(data.low,fm=all.fit10.entero.low.fpc1601, cluster=data.low$hhid)
overall.fit10.entero.low.fpc1601 <- coeftest(all.fit10.entero.low.fpc1601, all.VC10.entero.low.fpc1601)
summary(all.fit10.entero.low.fpc1601)
overall.fit10.entero.low.fpc1601
aic.entero.low.fpc1601=AIC(all.fit10.entero.low.fpc1601)

# FPC 1602 #####################
# high risk conditions --------------------------------
data=all[!is.na(all$entero35) & !is.na(all$fpc1602.pres),]
data.high=subset(data,data$risk=="High")
all.fit10.entero.high.fpc1602 <- glm(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+
                                rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC10.entero.high.fpc1602 <- cl(data.high,fm=all.fit10.entero.high.fpc1602, cluster=data.high$hhid)
overall.fit10.entero.high.fpc1602 <- coeftest(all.fit10.entero.high.fpc1602, all.VC10.entero.high.fpc1602)
summary(all.fit10.entero.high.fpc1602)
overall.fit10.entero.high.fpc1602
aic.entero.high.fpc1602=AIC(all.fit10.entero.high.fpc1602)

# low risk conditions --------------------------------
data.low=subset(data,data$risk=="Low")
all.fit10.entero.low.fpc1602 <- glm(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+
                               rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC10.entero.low.fpc1602 <- cl(data.low,fm=all.fit10.entero.low.fpc1602, cluster=data.low$hhid)
overall.fit10.entero.low.fpc1602 <- coeftest(all.fit10.entero.low.fpc1602, all.VC10.entero.low.fpc1602)
summary(all.fit10.entero.low.fpc1602)
overall.fit10.entero.low.fpc1602
aic.entero.low.fpc1602=AIC(all.fit10.entero.low.fpc1602)


# --------------------------------------
# save the results
# exclude glm objects and data frames
# (they are really large)
# --------------------------------------
save(
    
  n10.entero35.fmc1601,n10.entero35.fmc1602,n10.entero35.fpc1601,
  n10.entero35.fpc1602,

  n10.entero35.fmc1602.high,n10.entero35.fmc1602.low,
  n10.entero35.fpc1601.high,n10.entero35.fpc1601.low,
  n10.entero35.fpc1602.high,n10.entero35.fpc1602.low,
    
  overall.fit10.entero.fmc1601,overall.fit10.entero.fmc1602,
  overall.fit10.entero.fpc1601,overall.fit10.entero.fpc1602,
  
  overall.fit10.entero.high.fmc1602,overall.fit10.entero.high.fpc1601,
  overall.fit10.entero.high.fpc1602,
  
  overall.fit10.entero.low.fmc1602,overall.fit10.entero.low.fpc1601,
  overall.fit10.entero.low.fpc1602,
    
  all.VC10.entero.fmc1601,all.VC10.entero.fmc1602,all.VC10.entero.fpc1601,
  all.VC10.entero.fpc1602,
  all.VC10.entero.high.fmc1602,
  all.VC10.entero.high.fpc1601,all.VC10.entero.high.fpc1602,
  all.VC10.entero.low.fmc1602,all.VC10.entero.low.fpc1601,
  all.VC10.entero.low.fpc1602,
    
  aic.entero.fmc1601,aic.entero.fmc1602,aic.entero.fpc1601,aic.entero.fpc1602,
  aic.entero.low.fmc1602,aic.entero.low.fpc1601,aic.entero.low.fpc1602,
  aic.entero.high.fmc1602,aic.entero.high.fpc1601,aic.entero.high.fpc1602,
  
  file="~/Documents/CRG/coliphage/results/rawoutput/regress-10day-body-entero.Rdata"
)


