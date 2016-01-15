##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file conducts maximum likelihood regression
# to estimate prevalence ratios

# Results stratified by beach

# 3 day gi illness
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
avalon=data[data$beach %in% "Avalon",]
doheny=data[data$beach %in% "Doheny",]
malibu=data[data$beach %in% "Malibu",]
mission=data[data$beach %in% "Mission Bay",]
fairhope=data[data$beach %in% "Fairhope",]
goddard=data[data$beach %in% "Goddard",]

data.list=list(all=all,avalon=avalon,doheny=doheny,mission=mission,
               malibu=malibu,goddard=goddard,fairhope=fairhope)

data.list=lapply(data.list,function(df){
  # drop individuals with no water quality information
  df=subset(df,nowq==0)
  # subset to non-missing exposure categories
  # to make the robust CI calcs work
  df=subset(df,df$bodycontact=="Yes")
})

# convert from list back to data frames
list2env(data.list ,.GlobalEnv)

# --------------------------------------
# Calculate the actual Ns for each cell
# and store them for plotting and tables
# --------------------------------------
regN <- function(outcome,exposurecat) {
  sum(table(outcome,exposurecat))
}

#avalon
av.n3.fmc1602 = regN(avalon$gici3,avalon$fmc1602.pres)
av.n3.fpc1601 = regN(avalon$gici3,avalon$fpc1601.pres)
av.n3.fpc1602 = regN(avalon$gici3,avalon$fpc1602.pres)

#doheny
do.n3.fmc1601 = regN(doheny$gici3,doheny$fmc1601.pres)
do.n3.fmc1602 = regN(doheny$gici3,doheny$fmc1602.pres)
do.n3.fpc1601 = regN(doheny$gici3,doheny$fpc1601.pres)
do.n3.fpc1602 = regN(doheny$gici3,doheny$fpc1602.pres)

# fairhope
fa.n3.fpc1601 = regN(fairhope$gici3,fairhope$fpc1601.pres)

# goddard
# fpc 1601 always present at goddard

# malibu
# fmc 1602 and fpc 1602 always present at malibu
ma.n3.fpc1601 = regN(malibu$gici3,malibu$fpc1601.pres)

# mission bay
mb.n3.fmc1601 = regN(mission$gici3,mission$fmc1601.pres)
mb.n3.fpc1601 = regN(mission$gici3,mission$fpc1601.pres)

# n if low risk conditions ---------------------------
#avalon
av.n3.fmc1602.int0 = regN(avalon$gici3[avalon$groundwater=="Below median flow"],
  avalon$fmc1602.pres[avalon$groundwater=="Below median flow"])
av.n3.fpc1601.int0 = regN(avalon$gici3[avalon$groundwater=="Below median flow"],
  avalon$fpc1601.pres[avalon$groundwater=="Below median flow"])
av.n3.fpc1602.int0 = regN(avalon$gici3[avalon$groundwater=="Below median flow"],
  avalon$fpc1602.pres[avalon$groundwater=="Below median flow"])

#doheny
do.n3.fmc1602.int0 = regN(doheny$gici3[doheny$berm=="Closed"],
  doheny$fmc1602.pres[doheny$berm=="Closed"])
do.n3.fpc1601.int0 = regN(doheny$gici3[doheny$berm=="Closed"],
  doheny$fpc1601.pres[doheny$berm=="Closed"])
do.n3.fpc1602.int0 = regN(doheny$gici3[doheny$berm=="Closed"],
  doheny$fpc1602.pres[doheny$berm=="Closed"])

# malibu
# fmc 1602 and fpc 1602 always present at malibu
ma.n3.fpc1601.int0 = regN(malibu$gici3[malibu$berm=="Closed"],malibu$fpc1601.pres[malibu$berm=="Closed"])


# n if high risk conditions ---------------------------
#avalon
av.n3.fmc1602.int1 = regN(avalon$gici3[avalon$groundwater=="Above median flow"],
  avalon$fmc1602.pres[avalon$groundwater=="Above median flow"])
av.n3.fpc1601.int1 = regN(avalon$gici3[avalon$groundwater=="Above median flow"],
  avalon$fpc1601.pres[avalon$groundwater=="Above median flow"])
av.n3.fpc1602.int1 = regN(avalon$gici3[avalon$groundwater=="Above median flow"],
  avalon$fpc1602.pres[avalon$groundwater=="Above median flow"])

#doheny
do.n3.fmc1602.int1 = regN(doheny$gici3[doheny$berm=="Open"],
  doheny$fmc1602.pres[doheny$berm=="Open"])
do.n3.fpc1601.int1 = regN(doheny$gici3[doheny$berm=="Open"],
  doheny$fpc1601.pres[doheny$berm=="Open"])
do.n3.fpc1602.int1 = regN(doheny$gici3[doheny$berm=="Open"],
  doheny$fpc1602.pres[doheny$berm=="Open"])

# malibu
# fmc 1602 and fpc 1602 always present at malibu
ma.n3.fpc1601.int1 = regN(malibu$gici3[malibu$berm=="Open"],malibu$fpc1601.pres[malibu$berm=="Open"])


#-----------------------------------------
# 10-day illness and coliphage concentration
# pooled across berm and groundwater flow conditions
# by beach
#-----------------------------------------
# avalon
# fmc1601 -- always present at avalon 
avfit3.fmc1602 = mpreg(gici3~fmc1602.pres,
    dat=avalon,vcv=T)
avfit3.fpc1601 = mpreg(gici3~fpc1601.pres,
    dat=avalon,vcv=T)
avfit3.fpc1602 = mpreg(gici3~fpc1602.pres,
    dat=avalon,vcv=T)

# doheny
dofit3.fmc1601 = mpreg(gici3~fmc1601.pres,
     dat=doheny[!is.na(doheny$fmc1601.pres),],vcv=T)
dofit3.fmc1602 = mpreg(gici3~fmc1602.pres,
    dat=doheny,vcv=T)
dofit3.fpc1601 = mpreg(gici3~fpc1601.pres,
    dat=doheny,vcv=T)
dofit3.fpc1602 = mpreg(gici3~fpc1602.pres,
    dat=doheny,vcv=T)

# fairhope
fafit3.fpc1601 = mpreg(gici3~fpc1601.pres,
    dat=fairhope[!is.na(fairhope$fpc1601.pres),],vcv=T)

# goddard
# fpc 1601 always present at goddard

# malibu
# fmc 1602 and fpc 1602 always present at malibu
mafit3.fpc1601 = mpreg(gici3~fpc1601.pres,
    dat=malibu[!is.na(malibu$fpc1601.pres),],vcv=T)

# mission bay
mbfit3.fmc1601 = mpreg(gici3~fmc1601.pres,
    dat=mission[!is.na(mission$fmc1601.pres),],vcv=T)
mbfit3.fpc1601 = mpreg(gici3~fpc1601.pres,
    dat=mission[!is.na(mission$fpc1601.pres),],vcv=T)


#-----------------------------------------
# 3-day illness and coliphage concentration
# by berm and groundwater flow conditions
# by beach

# interaction tests
#-----------------------------------------
# reduced model for LR test
avfit3.fmc1602.glm = glm(gici3~fmc1602.pres+groundwater,
    dat=avalon,family=poisson(link="log"))
# interaction model
avfit3.fmc1602.gw = glm(gici3~fmc1602.pres*groundwater,
    family=poisson(link="log"),data=avalon)
summary(avfit3.fmc1602.gw)
lrtest(avfit3.fmc1602.glm,avfit3.fmc1602.gw)

# reduced model for LR test
avfit3.fpc1601.glm = glm(gici3~fpc1601.pres+groundwater,
    dat=avalon,family=poisson(link="log"))
# interaction model
avfit3.fpc1601.gw = glm(gici3~fpc1601.pres*groundwater,
    dat=avalon,family=poisson(link="log"))
summary(avfit3.fpc1601.gw)
lrtest(avfit3.fpc1601.glm,avfit3.fpc1601.gw)

# reduced model for LR test
avfit3.fpc1602.glm = glm(gici3~fpc1602.pres+groundwater,
    dat=avalon,family=poisson(link="log"))
# interaction model
avfit3.fpc1602.gw = glm(gici3~fpc1602.pres*groundwater,
    dat=avalon,family=poisson(link="log"))
summary(avfit3.fpc1602.gw)
lrtest(avfit3.fpc1602.glm,avfit3.fpc1602.gw)

# berm always closed when fmc1601 measured so no interaction assessment

# reduced model for LR test
dofit3.fmc1602.glm = glm(gici3~fmc1602.pres+berm,
    dat=doheny,family=poisson(link="log"))
# interaction model
dofit3.fmc1602.berm = glm(gici3~fmc1602.pres*berm,
    dat=doheny,family=poisson(link="log"))
summary(dofit3.fmc1602.berm)
lrtest(dofit3.fmc1602.glm,dofit3.fmc1602.berm)

# reduced model for LR test
dofit3.fpc1601.glm = glm(gici3~fpc1601.pres+berm,
    dat=doheny,family=poisson(link="log"))
# interaction model
dofit3.fpc1601.berm = glm(gici3~fpc1601.pres*berm,
    dat=doheny,family=poisson(link="log"))
summary(dofit3.fpc1601.berm)
lrtest(dofit3.fpc1601.glm,dofit3.fpc1601.berm)

# reduced model for LR test
dofit3.fpc1602.glm = glm(gici3~fpc1602.pres+berm,
    dat=doheny,family=poisson(link="log"))
# interaction model
dofit3.fpc1602.berm = glm(gici3~fpc1602.pres*berm,
    dat=doheny,family=poisson(link="log"))
summary(dofit3.fpc1602.berm)
lrtest(dofit3.fpc1602.glm,dofit3.fpc1602.berm)

# reduced model for LR test
mafit3.fpc1601.glm = glm(gici3~fpc1601.pres,
    dat=malibu[!is.na(malibu$fpc1601.pres),],family=poisson(link="log"))
# interaction model
mafit3.fpc1601.berm = glm(gici3~fpc1601.pres*berm,
    dat=malibu[!is.na(malibu$fpc1601.pres),],family=poisson(link="log"))
summary(mafit3.fpc1601.berm)
lrtest(mafit3.fpc1601.glm,mafit3.fpc1601.berm)

#-----------------------------------------
# 3-day illness and coliphage concentration
# by berm and groundwater flow conditions
# by beach

# stratified estimates
#-----------------------------------------
avfit3.fmc1602.gw.1 = mpreg(gici3~fmc1602.pres,
    vcv=T,dat=avalon[avalon$groundwater=="Above median flow",])
avfit3.fmc1602.gw.0 = mpreg(gici3~fmc1602.pres,
    vcv=T,dat=avalon[avalon$groundwater=="Below median flow",])

avfit3.fpc1601.gw.1 = mpreg(gici3~fpc1601.pres,
    vcv=T,dat=avalon[avalon$groundwater=="Above median flow",])
avfit3.fpc1601.gw.0 = mpreg(gici3~fpc1601.pres,
    vcv=T,dat=avalon[avalon$groundwater=="Below median flow",])

avfit3.fpc1602.gw.1 = mpreg(gici3~fpc1602.pres,
    vcv=T,dat=avalon[avalon$groundwater=="Above median flow",])
avfit3.fpc1602.gw.0 = mpreg(gici3~fpc1602.pres,
    vcv=T,dat=avalon[avalon$groundwater=="Below median flow",])

dofit3.fmc1602.berm.1 = mpreg(gici3~fmc1602.pres,
    vcv=T,dat=doheny[doheny$berm=="Open",])
dofit3.fmc1602.berm.0 = mpreg(gici3~fmc1602.pres,
    vcv=T,dat=doheny[doheny$berm=="Closed",])

dofit3.fpc1601.berm.1 = mpreg(gici3~fpc1601.pres,
    vcv=T,dat=doheny[doheny$berm=="Open",])
dofit3.fpc1601.berm.0 = mpreg(gici3~fpc1601.pres,
    vcv=T,dat=doheny[doheny$berm=="Closed",])

dofit3.fpc1602.berm.1 = mpreg(gici3~fpc1602.pres,
    vcv=T,dat=doheny[doheny$berm=="Open",])
dofit3.fpc1602.berm.0 = mpreg(gici3~fpc1602.pres,
    vcv=T,dat=doheny[doheny$berm=="Closed",])

mafit3.fpc1601.berm.1 = mpreg(gici3~fpc1601.pres,
    vcv=T,dat=malibu[malibu$berm=="Open",])
mafit3.fpc1601.berm.0 = mpreg(gici3~fpc1601.pres,
    vcv=T,dat=malibu[!is.na(malibu$fpc1601.pres) & malibu$berm=="Closed",])

# --------------------------------------
# save the results
# exclude glm objects and data frames
# (they are really large)
# --------------------------------------
save(

  av.n3.fmc1602, av.n3.fpc1601, av.n3.fpc1602, do.n3.fmc1601,
  do.n3.fmc1602, do.n3.fpc1601, do.n3.fpc1602, fa.n3.fpc1601,
  ma.n3.fpc1601, mb.n3.fmc1601, mb.n3.fpc1601,

  av.n3.fmc1602.int0, av.n3.fpc1601.int0, av.n3.fpc1602.int0, 
  do.n3.fmc1602.int0, do.n3.fpc1601.int0, do.n3.fpc1602.int0, 
  ma.n3.fpc1601.int0, 

  av.n3.fmc1602.int1, av.n3.fpc1601.int1, av.n3.fpc1602.int1, 
  do.n3.fmc1602.int1, do.n3.fpc1601.int1, do.n3.fpc1602.int1,
  ma.n3.fpc1601.int1, 

  avfit3.fmc1602, avfit3.fpc1601, avfit3.fpc1602,
  dofit3.fmc1601, dofit3.fmc1602, dofit3.fpc1601, dofit3.fpc1602,
  fafit3.fpc1601, mafit3.fpc1601,mbfit3.fmc1601,mbfit3.fpc1601,
  
  avfit3.fmc1602.gw.1, avfit3.fmc1602.gw.0, avfit3.fpc1601.gw.1, avfit3.fpc1601.gw.0, 
  avfit3.fpc1602.gw.1, avfit3.fpc1602.gw.0,
  dofit3.fmc1602.berm.1, dofit3.fmc1602.berm.0, 
  dofit3.fpc1601.berm.1, dofit3.fpc1601.berm.0, dofit3.fpc1602.berm.1, dofit3.fpc1602.berm.0,
  mafit3.fpc1601.berm.1, mafit3.fpc1601.berm.0,

  avfit3.fmc1602.gw, avfit3.fpc1601.gw, avfit3.fpc1602.gw,
  dofit3.fmc1602.berm, dofit3.fpc1601.berm, dofit3.fpc1602.berm,
  mafit3.fpc1601.berm,
  
  file="~/dropbox/coliphage/results/rawoutput/regress-3day-body-beach-unadj.Rdata"
)

