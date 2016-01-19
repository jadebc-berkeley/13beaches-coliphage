##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file conducts maximum likelihood regression
# to estimate prevalence ratios

# Results stratified by beach

# 10 day gi illness, body immersion
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
  df=subset(df,df$swallwater=="Yes")
})

# convert from list back to data frames
list2env(data.list ,.GlobalEnv)


#-----------------------------------------
# 10-day illness and coliphage concentration
# pooled across berm and groundwater flow conditions
# by beach
#-----------------------------------------
# avalon
# fmc1601 -- always present at avalon 
avfit10.fmc1602swall.unadj = mpreg(gici10~fmc1602.pres,dat=avalon,vcv=T)
avfit10.fpc1601swall.unadj = mpreg(gici10~fpc1601.pres,dat=avalon,vcv=T)
avfit10.fpc1602swall.unadj = mpreg(gici10~fpc1602.pres,dat=avalon,vcv=T)

# doheny
dofit10.fmc1601swall.unadj = mpreg(gici10~fmc1601.pres, dat=doheny[!is.na(doheny$fmc1601.pres),],vcv=T)
dofit10.fmc1602swall.unadj = mpreg(gici10~fmc1602.pres,dat=doheny,vcv=T)
dofit10.fpc1601swall.unadj = mpreg(gici10~fpc1601.pres,dat=doheny,vcv=T)
dofit10.fpc1602swall.unadj = mpreg(gici10~fpc1602.pres,dat=doheny,vcv=T)

# fairhope
fafit10.fpc1601swall.unadj = mpreg(gici10~fpc1601.pres,dat=fairhope[!is.na(fairhope$fpc1601.pres),],vcv=T)

# goddard
# fpc 1601 always present at goddard

# malibu
# fmc 1602 and fpc 1602 always present at malibu
mafit10.fpc1601swall.unadj = mpreg(gici10~fpc1601.pres,dat=malibu[!is.na(malibu$fpc1601.pres),],vcv=T)

# mission bay
mbfit10.fmc1601swall.unadj = mpreg(gici10~fmc1601.pres,dat=mission[!is.na(mission$fmc1601.pres),],vcv=T)
mbfit10.fpc1601swall.unadj = mpreg(gici10~fpc1601.pres,           dat=mission[!is.na(mission$fpc1601.pres),],vcv=T)

#-----------------------------------------
# 10-day illness and coliphage concentration
# by berm and groundwater flow conditions
# by beach

# stratified estimates
#-----------------------------------------
avfit10.fmc1602.gw.1.swall.unadj = mpreg(gici10~fmc1602.pres,vcv=T,dat=avalon[avalon$groundwater=="Above median flow",])
avfit10.fmc1602.gw.0.swall.unadj = mpreg(gici10~fmc1602.pres,vcv=T,dat=avalon[avalon$groundwater=="Below median flow",])

avfit10.fpc1601.gw.1.swall.unadj = mpreg(gici10~fpc1601.pres,vcv=T,dat=avalon[avalon$groundwater=="Above median flow",])
avfit10.fpc1601.gw.0.swall.unadj = mpreg(gici10~fpc1601.pres,vcv=T,dat=avalon[avalon$groundwater=="Below median flow",])

avfit10.fpc1602.gw.1.swall.unadj = mpreg(gici10~fpc1602.pres,vcv=T,dat=avalon[avalon$groundwater=="Above median flow",])
avfit10.fpc1602.gw.0.swall.unadj = mpreg(gici10~fpc1602.pres,vcv=T,dat=avalon[avalon$groundwater=="Below median flow",])

dofit10.fmc1602.berm.1.swall.unadj = mpreg(gici10~fmc1602.pres,vcv=T,dat=doheny[doheny$berm=="Open",])
dofit10.fmc1602.berm.0.swall.unadj = mpreg(gici10~fmc1602.pres,vcv=T,dat=doheny[doheny$berm=="Closed",])

dofit10.fpc1601.berm.1.swall.unadj = mpreg(gici10~fpc1601.pres,vcv=T,dat=doheny[doheny$berm=="Open",])
dofit10.fpc1601.berm.0.swall.unadj = mpreg(gici10~fpc1601.pres,vcv=T,dat=doheny[doheny$berm=="Closed",])

dofit10.fpc1602.berm.1.swall.unadj = mpreg(gici10~fpc1602.pres,vcv=T,dat=doheny[doheny$berm=="Open",])
dofit10.fpc1602.berm.0.swall.unadj = mpreg(gici10~fpc1602.pres,vcv=T,dat=doheny[doheny$berm=="Closed",])

mafit10.fpc1601.berm.1.swall.unadj = mpreg(gici10~fpc1601.pres,vcv=T,dat=malibu[malibu$berm=="Open",])
mafit10.fpc1601.berm.0.swall.unadj = mpreg(gici10~fpc1601.pres,vcv=T,dat=malibu[!is.na(malibu$fpc1601.pres) & malibu$berm=="Closed",])

# --------------------------------------
# save the results
# exclude glm objects and data frames
# (they are really large)
# --------------------------------------
save(

  avfit10.fmc1602swall.unadj, avfit10.fpc1601swall.unadj, avfit10.fpc1602swall.unadj,
  dofit10.fmc1601swall.unadj, dofit10.fmc1602swall.unadj, dofit10.fpc1601swall.unadj, dofit10.fpc1602swall.unadj,
  fafit10.fpc1601swall.unadj, mafit10.fpc1601swall.unadj,mbfit10.fmc1601swall.unadj,mbfit10.fpc1601swall.unadj,
  
  avfit10.fmc1602.gw.1.swall.unadj, avfit10.fmc1602.gw.0.swall.unadj, avfit10.fpc1601.gw.1.swall.unadj, avfit10.fpc1601.gw.0.swall.unadj, 
  avfit10.fpc1602.gw.1.swall.unadj, avfit10.fpc1602.gw.0.swall.unadj,
  dofit10.fmc1602.berm.1.swall.unadj, dofit10.fmc1602.berm.0.swall.unadj, 
  dofit10.fpc1601.berm.1.swall.unadj, dofit10.fpc1601.berm.0.swall.unadj, dofit10.fpc1602.berm.1.swall.unadj, dofit10.fpc1602.berm.0.swall.unadj,
  mafit10.fpc1601.berm.1.swall.unadj,mafit10.fpc1601.berm.0.swall.unadj,
  
  file="~/Documents/CRG/coliphage/results/rawoutput/regress-10day-swall-unadj.Rdata"
)



