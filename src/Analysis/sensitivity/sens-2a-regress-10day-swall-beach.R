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
  df=subset(df,df$swallwater=="Yes")
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

# pooled n's ---------------------------------------
all.n10.fmc1601swall = regN(all$gici10[!is.na(all$fmc1601.pres)],
                       all$fmc1601.pres[!is.na(all$fmc1601.pres)])
all.n10.fmc1602swall = regN(all$gici10[!is.na(all$fmc1602.pres)],
                       all$fmc1602.pres[!is.na(all$fmc1602.pres)])
all.n10.fpc1601swall = regN(all$gici10[!is.na(all$fpc1601.pres) ],
                       all$fpc1601.pres[!is.na(all$fpc1601.pres) ])
all.n10.fpc1602swall = regN(all$gici10[!is.na(all$fpc1602.pres) ],
                       all$fpc1602.pres[!is.na(all$fpc1602.pres) ])

# combined n's ---------------------------------------
#avalon
av.n10.fmc1602swall = regN(avalon$gici10,avalon$fmc1602.pres)
av.n10.fpc1601swall = regN(avalon$gici10,avalon$fpc1601.pres)
av.n10.fpc1602swall = regN(avalon$gici10,avalon$fpc1602.pres)

#doheny
do.n10.fmc1601swall = regN(doheny$gici10,doheny$fmc1601.pres)
do.n10.fmc1602swall = regN(doheny$gici10,doheny$fmc1602.pres)
do.n10.fpc1601swall = regN(doheny$gici10,doheny$fpc1601.pres)
do.n10.fpc1602swall = regN(doheny$gici10,doheny$fpc1602.pres)

# fairhope
fa.n10.fpc1601swall = regN(fairhope$gici10,fairhope$fpc1601.pres)

# goddard
# fpc 1601 always present at goddard

# malibu
# fmc 1602 and fpc 1602 always present at malibu
ma.n10.fpc1601swall = regN(malibu$gici10,malibu$fpc1601.pres)

# mission bay
mb.n10.fmc1601swall = regN(mission$gici10,mission$fmc1601.pres)
mb.n10.fpc1601swall = regN(mission$gici10,mission$fpc1601.pres)

# n if low risk conditions ---------------------------
#avalon
av.n10.fmc1602.int0swall = regN(avalon$gici10[avalon$groundwater=="Below median flow"],
  avalon$fmc1602.pres[avalon$groundwater=="Below median flow"])
av.n10.fpc1601.int0swall = regN(avalon$gici10[avalon$groundwater=="Below median flow"],
  avalon$fpc1601.pres[avalon$groundwater=="Below median flow"])
av.n10.fpc1602.int0swall = regN(avalon$gici10[avalon$groundwater=="Below median flow"],
  avalon$fpc1602.pres[avalon$groundwater=="Below median flow"])

#doheny
do.n10.fmc1602.int0swall = regN(doheny$gici10[doheny$berm=="Closed"],
  doheny$fmc1602.pres[doheny$berm=="Closed"])
do.n10.fpc1601.int0swall = regN(doheny$gici10[doheny$berm=="Closed"],
  doheny$fpc1601.pres[doheny$berm=="Closed"])
do.n10.fpc1602.int0swall = regN(doheny$gici10[doheny$berm=="Closed"],
  doheny$fpc1602.pres[doheny$berm=="Closed"])

# malibu
# fmc 1602 and fpc 1602 always present at malibu
ma.n10.fpc1601.int0swall = regN(malibu$gici10[malibu$berm=="Closed"],malibu$fpc1601.pres[malibu$berm=="Closed"])

# n if high risk conditions ---------------------------
#avalon
av.n10.fmc1602.int1swall = regN(avalon$gici10[avalon$groundwater=="Above median flow"],
  avalon$fmc1602.pres[avalon$groundwater=="Above median flow"])
av.n10.fpc1601.int1swall = regN(avalon$gici10[avalon$groundwater=="Above median flow"],
  avalon$fpc1601.pres[avalon$groundwater=="Above median flow"])
av.n10.fpc1602.int1swall = regN(avalon$gici10[avalon$groundwater=="Above median flow"],
  avalon$fpc1602.pres[avalon$groundwater=="Above median flow"])

#doheny
do.n10.fmc1602.int1swall = regN(doheny$gici10[doheny$berm=="Open"],
  doheny$fmc1602.pres[doheny$berm=="Open"])
do.n10.fpc1601.int1swall = regN(doheny$gici10[doheny$berm=="Open"],
  doheny$fpc1601.pres[doheny$berm=="Open"])
do.n10.fpc1602.int1swall = regN(doheny$gici10[doheny$berm=="Open"],
  doheny$fpc1602.pres[doheny$berm=="Open"])

# malibu
# fmc 1602 and fpc 1602 always present at malibu
ma.n10.fpc1601.int1swall = regN(malibu$gici10[malibu$berm=="Open"],malibu$fpc1601.pres[malibu$berm=="Open"])


#-----------------------------------------
# 10-day illness and coliphage concentration
# pooled across berm and groundwater flow conditions
# by beach
#-----------------------------------------
# avalon
# fmc1601 -- always present at avalon 
avfit10.fmc1602swall = mpreg(gici10~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,vcv=T)
avfit10.fpc1601swall = mpreg(gici10~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,vcv=T)
avfit10.fpc1602swall = mpreg(gici10~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,vcv=T)

# doheny
dofit10.fmc1601swall = mpreg(gici10~fmc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
     dat=doheny[!is.na(doheny$fmc1601.pres),],vcv=T)
dofit10.fmc1602swall = mpreg(gici10~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,vcv=T)
dofit10.fpc1601swall = mpreg(gici10~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,vcv=T)
dofit10.fpc1602swall = mpreg(gici10~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,vcv=T)

# fairhope
fafit10.fpc1601swall = mpreg(gici10~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=fairhope[!is.na(fairhope$fpc1601.pres),],vcv=T)

# goddard
# fpc 1601 always present at goddard

# malibu
# fmc 1602 and fpc 1602 always present at malibu
mafit10.fpc1601swall = mpreg(gici10~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=malibu[!is.na(malibu$fpc1601.pres),],vcv=T)

# mission bay
mbfit10.fmc1601swall = mpreg(gici10~fmc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=mission[!is.na(mission$fmc1601.pres),],vcv=T)
mbfit10.fpc1601swall = mpreg(gici10~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
               dat=mission[!is.na(mission$fpc1601.pres),],vcv=T)


#-----------------------------------------
# 10-day illness and coliphage concentration
# by berm and groundwater flow conditions
# by beach

# interaction tests
#-----------------------------------------
# reduced model for LR test
avfit10.fmc1602.glm.swall = glm(gici10~fmc1602.pres+groundwater+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,family=poisson(link="log"))
# interaction model
avfit10.fmc1602.gw.swall = glm(gici10~fmc1602.pres*groundwater+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    family=poisson(link="log"),data=avalon)
summary(avfit10.fmc1602.gw.swall)
lrtest(avfit10.fmc1602.glm.swall,avfit10.fmc1602.gw.swall)

# reduced model for LR test
avfit10.fpc1601.glm.swall = glm(gici10~fpc1601.pres+groundwater+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,family=poisson(link="log"))
# interaction model
avfit10.fpc1601.gw.swall = glm(gici10~fpc1601.pres*groundwater+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,family=poisson(link="log"))
summary(avfit10.fpc1601.gw.swall)
lrtest(avfit10.fpc1601.glm.swall,avfit10.fpc1601.gw.swall)

# reduced model for LR test
avfit10.fpc1602.glm.swall = glm(gici10~fpc1602.pres+groundwater+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,family=poisson(link="log"))
# interaction model
avfit10.fpc1602.gw.swall = glm(gici10~fpc1602.pres*groundwater+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,family=poisson(link="log"))
summary(avfit10.fpc1602.gw.swall)
lrtest(avfit10.fpc1602.glm.swall,avfit10.fpc1602.gw.swall)

# berm always closed when fmc1601 measured so no interaction assessment

# reduced model for LR test
dofit10.fmc1602.glm.swall = glm(gici10~fmc1602.pres+berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,family=poisson(link="log"))
# interaction model
dofit10.fmc1602.berm.swall = glm(gici10~fmc1602.pres*berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,family=poisson(link="log"))
summary(dofit10.fmc1602.berm.swall)
lrtest(dofit10.fmc1602.glm.swall,dofit10.fmc1602.berm.swall)

# reduced model for LR test
dofit10.fpc1601.glm.swall = glm(gici10~fpc1601.pres+berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,family=poisson(link="log"))
# interaction model
dofit10.fpc1601.berm.swall = glm(gici10~fpc1601.pres*berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,family=poisson(link="log"))
summary(dofit10.fpc1601.berm.swall)
lrtest(dofit10.fpc1601.glm.swall,dofit10.fpc1601.berm.swall)

# reduced model for LR test
dofit10.fpc1602.glm.swall = glm(gici10~fpc1602.pres+berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,family=poisson(link="log"))
# interaction model
dofit10.fpc1602.berm.swall = glm(gici10~fpc1602.pres*berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,family=poisson(link="log"))
summary(dofit10.fpc1602.berm.swall)
lrtest(dofit10.fpc1602.glm.swall,dofit10.fpc1602.berm.swall)

# reduced model for LR test
mafit10.fpc1601.glm.swall = glm(gici10~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=malibu[!is.na(malibu$fpc1601.pres),],family=poisson(link="log"))
# interaction model
mafit10.fpc1601.berm.swall = glm(gici10~fpc1601.pres*berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=malibu[!is.na(malibu$fpc1601.pres),],family=poisson(link="log"))
summary(mafit10.fpc1601.berm.swall)
lrtest(mafit10.fpc1601.glm.swall,mafit10.fpc1601.berm.swall)


#-----------------------------------------
# 10-day illness and coliphage concentration
# by berm and groundwater flow conditions
# by beach

# stratified estimates
#-----------------------------------------
avfit10.fmc1602.gw.1.swall = mpreg(gici10~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=avalon[avalon$groundwater=="Above median flow",])
avfit10.fmc1602.gw.0.swall = mpreg(gici10~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=avalon[avalon$groundwater=="Below median flow",])

avfit10.fpc1601.gw.1.swall = mpreg(gici10~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=avalon[avalon$groundwater=="Above median flow",])
avfit10.fpc1601.gw.0.swall = mpreg(gici10~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=avalon[avalon$groundwater=="Below median flow",])

avfit10.fpc1602.gw.1.swall = mpreg(gici10~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=avalon[avalon$groundwater=="Above median flow",])
avfit10.fpc1602.gw.0.swall = mpreg(gici10~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=avalon[avalon$groundwater=="Below median flow",])

dofit10.fmc1602.berm.1.swall = mpreg(gici10~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=doheny[doheny$berm=="Open",])
dofit10.fmc1602.berm.0.swall = mpreg(gici10~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=doheny[doheny$berm=="Closed",])

dofit10.fpc1601.berm.1.swall = mpreg(gici10~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=doheny[doheny$berm=="Open",])
dofit10.fpc1601.berm.0.swall = mpreg(gici10~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=doheny[doheny$berm=="Closed",])

dofit10.fpc1602.berm.1.swall = mpreg(gici10~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=doheny[doheny$berm=="Open",])
dofit10.fpc1602.berm.0.swall = mpreg(gici10~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=doheny[doheny$berm=="Closed",])

mafit10.fpc1601.berm.1.swall = mpreg(gici10~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=malibu[malibu$berm=="Open",])
mafit10.fpc1601.berm.0.swall = mpreg(gici10~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=malibu[!is.na(malibu$fpc1601.pres) & malibu$berm=="Closed",])

# --------------------------------------
# Pooled estimates
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness

# all beaches ----------------
# fmc 1601
all.fit10.fmc1601swall <- glm(gici10~fmc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc1601.pres) ,])

all.VC10.fmc1601swall <- cl(all[!is.na(all$fmc1601.pres) ],fm=all.fit10.fmc1601swall,
  cluster=all$hhid[!is.na(all$fmc1601.pres) ])
overall.fit10.fmc1601swall <- coeftest(all.fit10.fmc1601swall, all.VC10.fmc1601swall)
summary(all.fit10.fmc1601swall)
overall.fit10.fmc1601

# fmc 1602
all.fit10.fmc1602swall <- glm(gici10~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc1602.pres) ,])

all.VC10.fmc1602swall <- cl(all[!is.na(all$fmc1602.pres) ,],fm=all.fit10.fmc1602swall,
  cluster=all$hhid[!is.na(all$fmc1602.pres) ])
overall.fit10.fmc1602swall <- coeftest(all.fit10.fmc1602swall, all.VC10.fmc1602swall)
summary(all.fit10.fmc1602swall)
overall.fit10.fmc1602swall

# fpc 1601
all.fit10.fpc1601swall <- glm(gici10~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1601.pres) ,])

all.VC10.fpc1601swall <- cl(all[!is.na(all$fpc1601.pres) ,],fm=all.fit10.fpc1601swall,
    cluster=all$hhid[!is.na(all$fpc1601.pres) ])
overall.fit10.fpc1601swall <- coeftest(all.fit10.fpc1601swall, all.VC10.fpc1601swall)
summary(all.fit10.fpc1601swall)
overall.fit10.fpc1601swall

# fpc 1602
all.fit10.fpc1602swall <- glm(gici10~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1602.pres) ,])

all.VC10.fpc1602swall <- cl(all[!is.na(all$fpc1602.pres) ,],fm=all.fit10.fpc1602swall,
    cluster=all$hhid[!is.na(all$fpc1602.pres) ])
overall.fit10.fpc1602swall <- coeftest(all.fit10.fpc1602swall, all.VC10.fpc1602swall)
summary(all.fit10.fpc1602swall)
overall.fit10.fpc1602swall


# --------------------------------------
# save the results
# exclude glm objects and data frames
# (they are really large)
# --------------------------------------
save(

  all.n10.fmc1601swall,all.n10.fmc1602swall,all.n10.fpc1601swall,all.n10.fpc1602swall,

  av.n10.fmc1602swall, av.n10.fpc1601swall, av.n10.fpc1602swall, do.n10.fmc1601swall,
  do.n10.fmc1602swall, do.n10.fpc1601swall, do.n10.fpc1602swall, fa.n10.fpc1601swall,
  ma.n10.fpc1601swall, mb.n10.fmc1601swall, mb.n10.fpc1601swall,

  av.n10.fmc1602.int0swall, av.n10.fpc1601.int0swall, av.n10.fpc1602.int0swall, 
  do.n10.fmc1602.int0swall, do.n10.fpc1601.int0swall, do.n10.fpc1602.int0swall, 
  ma.n10.fpc1601.int0swall, 

  av.n10.fmc1602.int1swall, av.n10.fpc1601.int1swall, av.n10.fpc1602.int1swall, 
  do.n10.fmc1602.int1swall, do.n10.fpc1601.int1swall, do.n10.fpc1602.int1swall, 
  ma.n10.fpc1601.int1swall, 

  avfit10.fmc1602swall, avfit10.fpc1601swall, avfit10.fpc1602swall,
  dofit10.fmc1601swall, dofit10.fmc1602swall, dofit10.fpc1601swall, dofit10.fpc1602swall,
  fafit10.fpc1601swall, mafit10.fpc1601swall,mbfit10.fmc1601swall,mbfit10.fpc1601swall,
  
  avfit10.fmc1602.gw.1.swall, avfit10.fmc1602.gw.0.swall, avfit10.fpc1601.gw.1.swall, avfit10.fpc1601.gw.0.swall, 
  avfit10.fpc1602.gw.1.swall, avfit10.fpc1602.gw.0.swall,
  dofit10.fmc1602.berm.1.swall, dofit10.fmc1602.berm.0.swall, 
  dofit10.fpc1601.berm.1.swall, dofit10.fpc1601.berm.0.swall, dofit10.fpc1602.berm.1.swall, dofit10.fpc1602.berm.0.swall,
  mafit10.fpc1601.berm.1.swall,mafit10.fpc1601.berm.0.swall,
  
  avfit10.fmc1602.gw.swall, avfit10.fpc1601.gw.swall, avfit10.fpc1602.gw.swall,
  dofit10.fmc1602.berm.swall, dofit10.fpc1601.berm.swall, dofit10.fpc1602.berm.swall,
  mafit10.fpc1601.berm.swall,
  
  all.VC10.fmc1601swall,overall.fit10.fmc1601swall,
  all.VC10.fmc1602swall,overall.fit10.fmc1602swall,
  all.VC10.fpc1601swall,overall.fit10.fpc1601swall,
  all.VC10.fpc1602swall,overall.fit10.fpc1602swall,
  
  file="~/dropbox/coliphage/results/rawoutput/regress-10day-swall.Rdata"
)



