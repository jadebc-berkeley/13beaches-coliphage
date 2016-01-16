##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file conducts maximum likelihood regression
# to estimate prevalence ratios

# Results stratified by beach

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
# Creating joint indicator variable for
# regressions
# --------------------------------------
# avalon
avalon$fmc.ent=NA
avalon$fmc.ent[avalon$fmc.pres==0 & avalon$entero35==0]=1
avalon$fmc.ent[avalon$fmc.pres==1 & avalon$entero35==0]=2
avalon$fmc.ent[avalon$fmc.pres==0 & avalon$entero35==1]=3
avalon$fmc.ent[avalon$fmc.pres==1 & avalon$entero35==1]=4
avalon$fmc.ent=as.factor(avalon$fmc.ent)

avalon$fpc.ent=NA
avalon$fpc.ent[avalon$fpc.pres==0 & avalon$entero35==0]=1
avalon$fpc.ent[avalon$fpc.pres==1 & avalon$entero35==0]=2
avalon$fpc.ent[avalon$fpc.pres==0 & avalon$entero35==1]=3
avalon$fpc.ent[avalon$fpc.pres==1 & avalon$entero35==1]=4
avalon$fpc.ent=as.factor(avalon$fpc.ent)

# doheny
doheny$fmc.ent=NA
doheny$fmc.ent[doheny$fmc.pres==0 & doheny$entero35==0]=1
doheny$fmc.ent[doheny$fmc.pres==1 & doheny$entero35==0]=2
doheny$fmc.ent[doheny$fmc.pres==0 & doheny$entero35==1]=3
doheny$fmc.ent[doheny$fmc.pres==1 & doheny$entero35==1]=4
doheny$fmc.ent=as.factor(doheny$fmc.ent)

doheny$fpc.ent=NA
doheny$fpc.ent[doheny$fpc.pres==0 & doheny$entero35==0]=1
doheny$fpc.ent[doheny$fpc.pres==1 & doheny$entero35==0]=2
doheny$fpc.ent[doheny$fpc.pres==0 & doheny$entero35==1]=3
doheny$fpc.ent[doheny$fpc.pres==1 & doheny$entero35==1]=4
doheny$fpc.ent=as.factor(doheny$fpc.ent)

# malibu
malibu$fpc.ent=NA
malibu$fpc.ent[malibu$fpc.pres==0 & malibu$entero35==0]=1
malibu$fpc.ent[malibu$fpc.pres==1 & malibu$entero35==0]=2
malibu$fpc.ent[malibu$fpc.pres==0 & malibu$entero35==1]=3
malibu$fpc.ent[malibu$fpc.pres==1 & malibu$entero35==1]=4
malibu$fpc.ent=as.factor(malibu$fpc.ent)

# mission bay
mission$fmc.ent=NA
mission$fmc.ent[mission$fmc.pres==0 & mission$entero35==0]=1
mission$fmc.ent[mission$fmc.pres==1 & mission$entero35==0]=2
mission$fmc.ent[mission$fmc.pres==0 & mission$entero35==1]=3
mission$fmc.ent[mission$fmc.pres==1 & mission$entero35==1]=4
mission$fmc.ent=as.factor(mission$fmc.ent)

mission$fpc.ent=NA
mission$fpc.ent[mission$fpc.pres==0 & mission$entero35==0]=1
mission$fpc.ent[mission$fpc.pres==1 & mission$entero35==0]=2
mission$fpc.ent[mission$fpc.pres==0 & mission$entero35==1]=3
mission$fpc.ent[mission$fpc.pres==1 & mission$entero35==1]=4
mission$fpc.ent=as.factor(mission$fpc.ent)

# fairhope
fairhope$fpc.ent=NA
fairhope$fpc.ent[fairhope$fpc.pres==0 & fairhope$entero35==0]=1
fairhope$fpc.ent[fairhope$fpc.pres==1 & fairhope$entero35==0]=2
fairhope$fpc.ent[fairhope$fpc.pres==0 & fairhope$entero35==1]=3
fairhope$fpc.ent[fairhope$fpc.pres==1 & fairhope$entero35==1]=4
fairhope$fpc.ent=as.factor(fairhope$fpc.ent)

# goddard
goddard$fpc.ent=NA
goddard$fpc.ent[goddard$fpc.pres==0 & goddard$entero35==0]=1
goddard$fpc.ent[goddard$fpc.pres==1 & goddard$entero35==0]=2
goddard$fpc.ent[goddard$fpc.pres==0 & goddard$entero35==1]=3
goddard$fpc.ent[goddard$fpc.pres==1 & goddard$entero35==1]=4
goddard$fpc.ent=as.factor(goddard$fpc.ent)

# --------------------------------------
# Calculate the actual Ns for each cell
# and store them for plotting and tables
# --------------------------------------
regN <- function(outcome,exposurecat) {
  sum(table(outcome,exposurecat))
}

#avalon
av.n10.fmc = regN(avalon$gici10,avalon$fmc.ent)
av.n10.fpc = regN(avalon$gici10,avalon$fpc.ent)

#doheny
do.n10.fmc = regN(doheny$gici10,doheny$fmc.ent)
do.n10.fpc = regN(doheny$gici10,doheny$fpc.ent)

# fairhope
fa.n10.fpc = regN(fairhope$gici10,fairhope$fpc.ent)

# goddard
# fpc 1601 always present at goddard

# malibu
ma.n10.fpc = regN(malibu$gici10,malibu$fpc.ent)

# mission bay
mb.n10.fmc = regN(mission$gici10,mission$fmc.ent)
mb.n10.fpc = regN(mission$gici10,mission$fpc.ent)

# n if low risk conditions ---------------------------
#avalon
av.n10.fmc.int0 = regN(avalon$gici10[avalon$groundwater=="Below median flow"],
  avalon$fmc.ent[avalon$groundwater=="Below median flow"])
av.n10.fpc.int0 = regN(avalon$gici10[avalon$groundwater=="Below median flow"],
  avalon$fpc.ent[avalon$groundwater=="Below median flow"])

#doheny
do.n10.fmc.int0 = regN(doheny$gici10[doheny$berm=="Closed"],
  doheny$fmc.ent[doheny$berm=="Closed"])
do.n10.fpc.int0 = regN(doheny$gici10[doheny$berm=="Closed"],
  doheny$fpc.ent[doheny$berm=="Closed"])

# malibu
# fmc 1602 and fpc 1602 always present at malibu
ma.n10.fpc.int0 = regN(malibu$gici10[malibu$berm=="Closed"],
  malibu$fpc.ent[malibu$berm=="Closed"])


# n if high risk conditions ---------------------------
#avalon
av.n10.fmc.int1 = regN(avalon$gici10[avalon$groundwater=="Above median flow"],
  avalon$fmc.ent[avalon$groundwater=="Above median flow"])
av.n10.fpc.int1 = regN(avalon$gici10[avalon$groundwater=="Above median flow"],
  avalon$fpc.ent[avalon$groundwater=="Above median flow"])

#doheny
do.n10.fmc.int1 = regN(doheny$gici10[doheny$berm=="Open"],
  doheny$fmc.ent[doheny$berm=="Open"])
do.n10.fpc.int1 = regN(doheny$gici10[doheny$berm=="Open"],
  doheny$fpc.ent[doheny$berm=="Open"])

# malibu
# fmc 1602 and fpc 1602 always present at malibu
ma.n10.fpc.int1 = regN(malibu$gici10[malibu$berm=="Open"],
  malibu$fpc.ent[malibu$berm=="Open"])


#-----------------------------------------
# 10-day illness and coliphage concentration
# pooled across berm and groundwater flow conditions
# by beach
#-----------------------------------------
# avalon
# fmc1601 -- always present at avalon 
avfit10.fmc = mpreg(gici10~fmc.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,vcv=T)
avfit10.fpc = mpreg(gici10~fpc.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,vcv=T)

# doheny
dofit10.fmc = mpreg(gici10~fmc.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
     dat=doheny[!is.na(doheny$fmc.ent),],vcv=T)
dofit10.fpc = mpreg(gici10~fpc.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,vcv=T)

# fairhope
fafit10.fpc = mpreg(gici10~fpc.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=fairhope[!is.na(fairhope$fpc.ent),],vcv=T)

# goddard
# fpc 1601 always present at goddard

# malibu
# fmc 1602 and fpc 1602 always present at malibu
mafit10.fpc = mpreg(gici10~fpc.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=malibu[!is.na(malibu$fpc.ent),],vcv=T)

# mission bay
mbfit10.fmc = mpreg(gici10~fmc.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=mission[!is.na(mission$fmc.ent),],vcv=T)
mbfit10.fpc = mpreg(gici10~fpc.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=mission[!is.na(mission$fpc.ent),],vcv=T)


#-----------------------------------------
# 3-day illness and coliphage concentration
# by berm and groundwater flow conditions
# by beach

# interaction tests
#-----------------------------------------
# reduced model for LR test
avfit10.fmc.glm = glm(gici10~fmc.ent+groundwater+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,family=poisson(link="log"))
# interaction model
avfit10.fmc.gw = glm(gici10~fmc.ent*groundwater+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    family=poisson(link="log"),data=avalon)
summary(avfit10.fmc.gw)
lrtest(avfit10.fmc.glm,avfit10.fmc.gw)

# reduced model for LR test
avfit10.fpc.glm = glm(gici10~fpc.ent+groundwater+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,family=poisson(link="log"))
# interaction model
avfit10.fpc.gw = glm(gici10~fpc.ent*groundwater+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,family=poisson(link="log"))
summary(avfit10.fpc.gw)
lrtest(avfit10.fpc.glm,avfit10.fpc.gw)


# berm always closed when fmc1601 measured so no interaction assessment

# reduced model for LR test
dofit10.fmc.glm = glm(gici10~fmc.ent+berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,family=poisson(link="log"))
# interaction model
dofit10.fmc.berm = glm(gici10~fmc.ent*berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,family=poisson(link="log"))
summary(dofit10.fmc.berm)
lrtest(dofit10.fmc.glm,dofit10.fmc.berm)

# reduced model for LR test
dofit10.fpc.glm = glm(gici10~fpc.ent+berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,family=poisson(link="log"))
# interaction model
dofit10.fpc.berm = glm(gici10~fpc.ent*berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,family=poisson(link="log"))
summary(dofit10.fpc.berm)
lrtest(dofit10.fpc.glm,dofit10.fpc.berm)

# reduced model for LR test
mafit10.fpc.glm = glm(gici10~fpc.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=malibu[!is.na(malibu$fpc.ent),],family=poisson(link="log"))
# interaction model
mafit10.fpc.berm = glm(gici10~fpc.ent*berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=malibu[!is.na(malibu$fpc.ent),],family=poisson(link="log"))
summary(mafit10.fpc.berm)
lrtest(mafit10.fpc.glm,mafit10.fpc.berm)

#-----------------------------------------
# 3-day illness and coliphage concentration
# by berm and groundwater flow conditions
# by beach

# stratified estimates
#-----------------------------------------
avfit10.fmc.gw.1 = mpreg(gici10~fmc.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=avalon[avalon$groundwater=="Above median flow",])
avfit10.fmc.gw.0 = mpreg(gici10~fmc.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=avalon[avalon$groundwater=="Below median flow",])

avfit10.fpc.gw.1 = mpreg(gici10~fpc.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=avalon[avalon$groundwater=="Above median flow",])
avfit10.fpc.gw.0 = mpreg(gici10~fpc.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=avalon[avalon$groundwater=="Below median flow",])

dofit10.fmc.berm.1 = mpreg(gici10~fmc.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=doheny[doheny$berm=="Open",])
dofit10.fmc.berm.0 = mpreg(gici10~fmc.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=doheny[doheny$berm=="Closed",])

dofit10.fpc.berm.1 = mpreg(gici10~fpc.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=doheny[doheny$berm=="Open",])
dofit10.fpc.berm.0 = mpreg(gici10~fpc.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=doheny[doheny$berm=="Closed",])

mafit10.fpc.berm.1 = mpreg(gici10~fpc.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=malibu[malibu$berm=="Open",])
mafit10.fpc.berm.0 = mpreg(gici10~fpc.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=malibu[!is.na(malibu$fpc.ent) & malibu$berm=="Closed",])


# --------------------------------------
# save the results
# exclude glm objects and data frames
# (they are really large)
# --------------------------------------
save(
  
  av.n10.fmc, av.n10.fpc, 
  do.n10.fmc, do.n10.fpc, 
  fa.n10.fpc,
  ma.n10.fpc, 
  mb.n10.fmc, mb.n10.fpc,

  av.n10.fmc.int0, av.n10.fpc.int0, 
  do.n10.fmc.int0, do.n10.fpc.int0, 
  ma.n10.fpc.int0, 

  av.n10.fmc.int1, av.n10.fpc.int1, 
  do.n10.fmc.int1, do.n10.fpc.int1, 
  ma.n10.fpc.int1, 

  avfit10.fmc, avfit10.fpc, 
  dofit10.fmc, dofit10.fpc, 
  fafit10.fpc, 
  mafit10.fpc,mbfit10.fmc,mbfit10.fpc,
  
  avfit10.fmc.gw.1, avfit10.fmc.gw.0, avfit10.fpc.gw.1, avfit10.fpc.gw.0, 
  dofit10.fmc.berm.1, dofit10.fmc.berm.0, 
  dofit10.fpc.berm.1, dofit10.fpc.berm.0, 
  mafit10.fpc.berm.1, mafit10.fpc.berm.0,

  avfit10.fmc.gw, avfit10.fpc.gw, 
  dofit10.fmc.berm, dofit10.fpc.berm, 
  mafit10.fpc.berm,

  file="~/Documents/CRG/coliphage/results/rawoutput/regress-10day-body-joint-beach.Rdata"
)

