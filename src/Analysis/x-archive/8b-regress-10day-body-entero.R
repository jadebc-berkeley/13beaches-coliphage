##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file conducts maximum likelihood regression
# to estimate prevalence ratios
# for enterococcus

# 10 day gi illness
##########################################

rm(list=ls())
library(foreign)
library(readstata13)

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
# Calculate the actual Ns for each cell
# and store them for plotting and tables
# --------------------------------------
regN <- function(outcome,exposurecat) {
  sum(table(outcome,exposurecat))
}

# pooled n's ---------------------------------------
all.n10.entero35 = regN(all$gici10,all$entero35)

data=all[!is.na(all$entero35),]
data.high=subset(data,data$groundwater=="Above median flow" | data$berm=="Open")
all.n10.entero35.high = regN(data.high$gici10,data.high$entero35)
  
data.low=subset(data,data$groundwater=="Below median flow" | data$berm=="Closed")
all.n10.entero35.low = regN(data.low$gici10,data.low$entero35)
  

# --------------------------------------
# Estimates pooled across beach
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness

# all beaches ----------------
all.fit10.entero <- glm(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$entero35),])

all.VC10.entero <- cl(all[!is.na(all$entero35)],fm=all.fit10.entero,
                       cluster=all$hhid[!is.na(all$entero35)])
overall.fit10.entero <- coeftest(all.fit10.entero, all.VC10.entero)
summary(all.fit10.entero)
overall.fit10.entero
aic.entero=AIC(all.fit10.entero)

# high risk conditions
data=all[!is.na(all$entero35),]
data.high=subset(data,data$groundwater=="Above median flow" | data$berm=="Open")
all.fit10.entero.high <- glm(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+
                                rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC10.entero.high <- cl(data.high,fm=all.fit10.entero.high, cluster=data.high$hhid)
overall.fit10.entero.high <- coeftest(all.fit10.entero.high, all.VC10.entero.high)
summary(all.fit10.entero.high)
overall.fit10.entero.high
aic.entero.high=AIC(all.fit10.entero.high)

# low risk conditions
data.low=subset(data,data$groundwater=="Below median flow" | data$berm=="Closed")
all.fit10.entero.low <- glm(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+
                               rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC10.entero.low <- cl(data.low,fm=all.fit10.entero.low, cluster=data.low$hhid)
overall.fit10.entero.low <- coeftest(all.fit10.entero.low, all.VC10.entero.low)
summary(all.fit10.entero.low)
overall.fit10.entero.low
aic.entero.low=AIC(all.fit10.entero.low)


# --------------------------------------
# save the results
# exclude glm objects and data frames
# (they are really large)
# --------------------------------------
save(
  all.n10.entero35,all.n10.entero35.high,all.n10.entero35.low,
  
  overall.fit10.entero,overall.fit10.entero.high,overall.fit10.entero.low,
  all.VC10.entero,all.VC10.entero.high,all.VC10.entero.low,
  
  aic.entero,aic.entero.low,aic.entero.high,
  
  file="~/dropbox/coliphage/results/rawoutput/regress-10day-body-entero.Rdata"
)


