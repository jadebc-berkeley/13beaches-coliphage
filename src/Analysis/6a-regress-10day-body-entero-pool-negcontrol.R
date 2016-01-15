##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file conducts maximum likelihood regression
# to estimate prevalence ratios
# for enterococcus

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
all=subset(all,all$bodycontact=="No")

# --------------------------------------
# Estimates pooled across beach
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness
# --------------------------------------

# f- coliphage --------------------------------
all.fit10.entero.fmc <- glm(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$entero35) & 
    !is.na(all$fmc.pres),])

all.VC10.entero.fmc <- cl(all[!is.na(all$entero35) & !is.na(all$fmc.pres)],
  fm=all.fit10.entero.fmc, cluster=
      all$hhid[!is.na(all$entero35)  & !is.na(all$fmc.pres)])
overall.fit10.entero.fmc <- coeftest(all.fit10.entero.fmc, all.VC10.entero.fmc)
summary(all.fit10.entero.fmc)
overall.fit10.entero.fmc

# f+ coliphage --------------------------------
all.fit10.entero.fpc <- glm(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$entero35) & 
    !is.na(all$fpc.pres),])

all.VC10.entero.fpc <- cl(all[!is.na(all$entero35) & !is.na(all$fpc.pres)],
  fm=all.fit10.entero.fpc, cluster=
      all$hhid[!is.na(all$entero35)  & !is.na(all$fpc.pres)])
overall.fit10.entero.fpc <- coeftest(all.fit10.entero.fpc, all.VC10.entero.fpc)
summary(all.fit10.entero.fpc)
overall.fit10.entero.fpc

# --------------------------------------
# Estimates pooled across beach and stratified by conditions
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness

# all beaches ----------------
# F- coliphage #####################
# high risk conditions --------------------------------
data=all[!is.na(all$entero35) & !is.na(all$fmc.pres),]
data.high=subset(data,data$risk=="High")
all.fit10.entero.high.fmc <- glm(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC10.entero.high.fmc <- cl(data.high,fm=all.fit10.entero.high.fmc, cluster=data.high$hhid)
overall.fit10.entero.high.fmc <- coeftest(all.fit10.entero.high.fmc, all.VC10.entero.high.fmc)
summary(all.fit10.entero.high.fmc)
overall.fit10.entero.high.fmc

# low risk conditions --------------------------------
data.low=subset(data,data$risk=="Low")
all.fit10.entero.low.fmc <- glm(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC10.entero.low.fmc <- cl(data.low,fm=all.fit10.entero.low.fmc, cluster=data.low$hhid)
overall.fit10.entero.low.fmc <- coeftest(all.fit10.entero.low.fmc, all.VC10.entero.low.fmc)
summary(all.fit10.entero.low.fmc)
overall.fit10.entero.low.fmc

# F+ coliphage #####################
# high risk conditions --------------------------------
data=all[!is.na(all$entero35) & !is.na(all$fpc.pres),]
data.high=subset(data,data$risk=="High")
all.fit10.entero.high.fpc <- glm(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+
                                rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC10.entero.high.fpc <- cl(data.high,fm=all.fit10.entero.high.fpc, cluster=data.high$hhid)
overall.fit10.entero.high.fpc <- coeftest(all.fit10.entero.high.fpc, all.VC10.entero.high.fpc)
summary(all.fit10.entero.high.fpc)
overall.fit10.entero.high.fpc

# low risk conditions --------------------------------
data.low=subset(data,data$risk=="Low")
all.fit10.entero.low.fpc <- glm(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+
                               rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC10.entero.low.fpc <- cl(data.low,fm=all.fit10.entero.low.fpc, cluster=data.low$hhid)
overall.fit10.entero.low.fpc <- coeftest(all.fit10.entero.low.fpc, all.VC10.entero.low.fpc)
summary(all.fit10.entero.low.fpc)
overall.fit10.entero.low.fpc


# --------------------------------------
# save the results
# exclude glm objects and data frames
# (they are really large)
# --------------------------------------
save(

  overall.fit10.entero.fmc,overall.fit10.entero.fpc,
  
  overall.fit10.entero.high.fmc,overall.fit10.entero.high.fpc,
  overall.fit10.entero.low.fmc,overall.fit10.entero.low.fpc,
  
  all.VC10.entero.fmc,all.VC10.entero.fpc,

  all.VC10.entero.high.fmc, all.VC10.entero.low.fmc,
  all.VC10.entero.high.fpc, all.VC10.entero.low.fpc,
  
  file="~/dropbox/coliphage/results/rawoutput/regress-10day-body-entero-pool-negcontrol.Rdata"
)


