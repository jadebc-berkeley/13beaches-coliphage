##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file conducts maximum likelihood regression
# to estimate prevalence ratios

# Results stratified by beach

# 3 day gi illness, head immersion
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
  df=subset(df,df$headunder=="Yes")
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
all.n3.fmc1601head = regN(all$gici3[!is.na(all$fmc1601.pres)& all$beach!="Avalon"],
                       all$fmc1601.pres[!is.na(all$fmc1601.pres)& all$beach!="Avalon"])
all.n3.fmc1602head = regN(all$gici3[!is.na(all$fmc1602.pres)& all$beach!="Malibu"],
                       all$fmc1602.pres[!is.na(all$fmc1602.pres)& all$beach!="Malibu"])
all.n3.fpc1601head = regN(all$gici3[!is.na(all$fpc1601.pres) & all$beach!="Goddard"],
                       all$fpc1601.pres[!is.na(all$fpc1601.pres) & all$beach!="Goddard"])
all.n3.fpc1602head = regN(all$gici3[!is.na(all$fpc1602.pres) & all$beach!="Malibu"],
                       all$fpc1602.pres[!is.na(all$fpc1602.pres) & all$beach!="Malibu"])

# combined n's ---------------------------------------
#avalon
av.n3.fmc1602head = regN(avalon$gici3,avalon$fmc1602.pres)
av.n3.fpc1601head = regN(avalon$gici3,avalon$fpc1601.pres)
av.n3.fpc1602head = regN(avalon$gici3,avalon$fpc1602.pres)

#doheny
do.n3.fmc1601head = regN(doheny$gici3,doheny$fmc1601.pres)
do.n3.fmc1602head = regN(doheny$gici3,doheny$fmc1602.pres)
do.n3.fpc1601head = regN(doheny$gici3,doheny$fpc1601.pres)
do.n3.fpc1602head = regN(doheny$gici3,doheny$fpc1602.pres)

# fairhope
fa.n3.fpc1601head = regN(fairhope$gici3,fairhope$fpc1601.pres)

# goddard
# fpc 1601 always present at goddard

# malibu
# fmc 1602 and fpc 1602 always present at malibu
ma.n3.fpc1601head = regN(malibu$gici3,malibu$fpc1601.pres)

# mission bay
mb.n3.fmc1601head = regN(mission$gici3,mission$fmc1601.pres)
mb.n3.fpc1601head = regN(mission$gici3,mission$fpc1601.pres)

# n if low risk conditions ---------------------------
#avalon
av.n3.fmc1602.int0head = regN(avalon$gici3[avalon$groundwater=="Below median flow"],
  avalon$fmc1602.pres[avalon$groundwater=="Below median flow"])
av.n3.fpc1601.int0head = regN(avalon$gici3[avalon$groundwater=="Below median flow"],
  avalon$fpc1601.pres[avalon$groundwater=="Below median flow"])
av.n3.fpc1602.int0head = regN(avalon$gici3[avalon$groundwater=="Below median flow"],
  avalon$fpc1602.pres[avalon$groundwater=="Below median flow"])

#doheny
do.n3.fmc1602.int0head = regN(doheny$gici3[doheny$berm=="Closed"],
  doheny$fmc1602.pres[doheny$berm=="Closed"])
do.n3.fpc1601.int0head = regN(doheny$gici3[doheny$berm=="Closed"],
  doheny$fpc1601.pres[doheny$berm=="Closed"])
do.n3.fpc1602.int0head = regN(doheny$gici3[doheny$berm=="Closed"],
  doheny$fpc1602.pres[doheny$berm=="Closed"])

# malibu
# fmc 1602 and fpc 1602 always present at malibu
ma.n3.fpc1601.int0head = regN(malibu$gici3[malibu$berm=="Closed"],malibu$fpc1601.pres[malibu$berm=="Closed"])


# n if high risk conditions ---------------------------
#avalon
av.n3.fmc1602.int1head = regN(avalon$gici3[avalon$groundwater=="Above median flow"],
  avalon$fmc1602.pres[avalon$groundwater=="Above median flow"])
av.n3.fpc1601.int1head = regN(avalon$gici3[avalon$groundwater=="Above median flow"],
  avalon$fpc1601.pres[avalon$groundwater=="Above median flow"])
av.n3.fpc1602.int1head = regN(avalon$gici3[avalon$groundwater=="Above median flow"],
  avalon$fpc1602.pres[avalon$groundwater=="Above median flow"])

#doheny
do.n3.fmc1602.int1head = regN(doheny$gici3[doheny$berm=="Open"],
  doheny$fmc1602.pres[doheny$berm=="Open"])
do.n3.fpc1601.int1head = regN(doheny$gici3[doheny$berm=="Open"],
  doheny$fpc1601.pres[doheny$berm=="Open"])
do.n3.fpc1602.int1head = regN(doheny$gici3[doheny$berm=="Open"],
  doheny$fpc1602.pres[doheny$berm=="Open"])

# malibu
# fmc 1602 and fpc 1602 always present at malibu
ma.n3.fpc1601.int1head = regN(malibu$gici3[malibu$berm=="Open"],malibu$fpc1601.pres[malibu$berm=="Open"])

#-----------------------------------------
# 3-day illness and coliphage concentration
# pooled across berm and groundwater flow conditions
# by beach
#-----------------------------------------
# avalon
# fmc1601 -- always present at avalon 
avfit3.fmc1602head = mpreg(gici3~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,vcv=T)
avfit3.fpc1601head = mpreg(gici3~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,vcv=T)
avfit3.fpc1602head = mpreg(gici3~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,vcv=T)

# doheny
dofit3.fmc1601head = mpreg(gici3~fmc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
     dat=doheny[!is.na(doheny$fmc1601.pres),],vcv=T)
dofit3.fmc1602head = mpreg(gici3~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,vcv=T)
dofit3.fpc1601head = mpreg(gici3~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,vcv=T)
dofit3.fpc1602head = mpreg(gici3~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,vcv=T)

# fairhope
fafit3.fpc1601head = mpreg(gici3~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=fairhope[!is.na(fairhope$fpc1601.pres),],vcv=T)

# goddard
# fpc 1601 always present at goddard

# malibu
# fmc 1602 and fpc 1602 always present at malibu
mafit3.fpc1601head = mpreg(gici3~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=malibu[!is.na(malibu$fpc1601.pres),],vcv=T)

# mission bay
mbfit3.fmc1601head = mpreg(gici3~fmc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=mission[!is.na(mission$fmc1601.pres),],vcv=T)
mbfit3.fpc1601head = mpreg(gici3~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=mission[!is.na(mission$fpc1601.pres),],vcv=T)


#-----------------------------------------
# 3-day illness and coliphage concentration
# by berm and groundwater flow conditions
# by beach

# interaction tests
#-----------------------------------------
# reduced model for LR test
avfit3.fmc1602.glm.head = glm(gici3~fmc1602.pres+groundwater+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,family=poisson(link="log"))
# interaction model
avfit3.fmc1602.gw.head = glm(gici3~fmc1602.pres*groundwater+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    family=poisson(link="log"),data=avalon)
summary(avfit3.fmc1602.gw.head)
lrtest(avfit3.fmc1602.glm.head,avfit3.fmc1602.gw.head)

# reduced model for LR test
avfit3.fpc1601.glm.head = glm(gici3~fpc1601.pres+groundwater+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,family=poisson(link="log"))
# interaction model
avfit3.fpc1601.gw.head = glm(gici3~fpc1601.pres*groundwater+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,family=poisson(link="log"))
summary(avfit3.fpc1601.gw.head)
lrtest(avfit3.fpc1601.glm.head,avfit3.fpc1601.gw.head)

# reduced model for LR test
avfit3.fpc1602.glm.head = glm(gici3~fpc1602.pres+groundwater+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,family=poisson(link="log"))
# interaction model
avfit3.fpc1602.gw.head = glm(gici3~fpc1602.pres*groundwater+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,family=poisson(link="log"))
summary(avfit3.fpc1602.gw.head)
lrtest(avfit3.fpc1602.glm.head,avfit3.fpc1602.gw.head)

# berm always closed when fmc1601 measured so no interaction assessment

# reduced model for LR test
dofit3.fmc1602.glm.head = glm(gici3~fmc1602.pres+berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,family=poisson(link="log"))
# interaction model
dofit3.fmc1602.berm.head = glm(gici3~fmc1602.pres*berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,family=poisson(link="log"))
summary(dofit3.fmc1602.berm.head)
lrtest(dofit3.fmc1602.glm.head,dofit3.fmc1602.berm.head)

# reduced model for LR test
dofit3.fpc1601.glm.head = glm(gici3~fpc1601.pres+berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,family=poisson(link="log"))
# interaction model
dofit3.fpc1601.berm.head = glm(gici3~fpc1601.pres*berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,family=poisson(link="log"))
summary(dofit3.fpc1601.berm.head)
lrtest(dofit3.fpc1601.glm.head,dofit3.fpc1601.berm.head)

# reduced model for LR test
dofit3.fpc1602.glm.head = glm(gici3~fpc1602.pres+berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,family=poisson(link="log"))
# interaction model
dofit3.fpc1602.berm.head = glm(gici3~fpc1602.pres*berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,family=poisson(link="log"))
summary(dofit3.fpc1602.berm.head)
lrtest(dofit3.fpc1602.glm.head,dofit3.fpc1602.berm.head)

# reduced model for LR test
mafit3.fpc1601.glm.head = glm(gici3~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=malibu[!is.na(malibu$fpc1601.pres),],family=poisson(link="log"))
# interaction model
mafit3.fpc1601.berm.head = glm(gici3~fpc1601.pres*berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=malibu[!is.na(malibu$fpc1601.pres),],family=poisson(link="log"))
summary(mafit3.fpc1601.berm.head)
lrtest(mafit3.fpc1601.glm.head,mafit3.fpc1601.berm.head)


#-----------------------------------------
# 3-day illness and coliphage concentration
# by berm and groundwater flow conditions
# by beach

# stratified estimates
#-----------------------------------------
avfit3.fmc1602.gw.1.head = mpreg(gici3~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=avalon[avalon$groundwater=="Above median flow",])
avfit3.fmc1602.gw.0.head = mpreg(gici3~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=avalon[avalon$groundwater=="Below median flow",])

avfit3.fpc1601.gw.1.head = mpreg(gici3~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=avalon[avalon$groundwater=="Above median flow",])
avfit3.fpc1601.gw.0.head = mpreg(gici3~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=avalon[avalon$groundwater=="Below median flow",])

avfit3.fpc1602.gw.1.head = mpreg(gici3~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=avalon[avalon$groundwater=="Above median flow",])
avfit3.fpc1602.gw.0.head = mpreg(gici3~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=avalon[avalon$groundwater=="Below median flow",])

dofit3.fmc1602.berm.1.head = mpreg(gici3~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=doheny[doheny$berm=="Open",])
dofit3.fmc1602.berm.0.head = mpreg(gici3~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=doheny[doheny$berm=="Closed",])

dofit3.fpc1601.berm.1.head = mpreg(gici3~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=doheny[doheny$berm=="Open",])
dofit3.fpc1601.berm.0.head = mpreg(gici3~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=doheny[doheny$berm=="Closed",])

dofit3.fpc1602.berm.1.head = mpreg(gici3~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=doheny[doheny$berm=="Open",])
dofit3.fpc1602.berm.0.head = mpreg(gici3~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=doheny[doheny$berm=="Closed",])

mafit3.fpc1601.berm.1.head = mpreg(gici3~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=malibu[malibu$berm=="Open",])
mafit3.fpc1601.berm.0.head = mpreg(gici3~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=malibu[!is.na(malibu$fpc1601.pres) & malibu$berm=="Closed",])

# --------------------------------------
# Pooled estimates
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 3-day GI illness

# all beaches ----------------
# fmc 1601
all.fit3.fmc1601head <- glm(gici10~fmc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc1601.pres) & all$beach!="Avalon",])

all.VC3.fmc1601head <- cl(all[!is.na(all$fmc1601.pres) & all$beach!="Avalon"],fm=all.fit3.fmc1601head,
  cluster=all$hhid[!is.na(all$fmc1601.pres) & all$beach!="Avalon"])
overall.fit3.fmc1601head <- coeftest(all.fit3.fmc1601head, all.VC3.fmc1601head)
summary(all.fit3.fmc1601)
overall.fit3.fmc1601

# fmc 1602
all.fit3.fmc1602head <- glm(gici10~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc1602.pres) & all$beach!="Malibu",])

all.VC3.fmc1602head <- cl(all[!is.na(all$fmc1602.pres) & all$beach!="Malibu",],fm=all.fit3.fmc1602head,
  cluster=all$hhid[!is.na(all$fmc1602.pres) & all$beach!="Malibu"])
overall.fit3.fmc1602head <- coeftest(all.fit3.fmc1602head, all.VC3.fmc1602head)
summary(all.fit3.fmc1602)
overall.fit3.fmc1602

# fpc 1601
all.fit3.fpc1601head <- glm(gici10~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1601.pres) & all$beach!="Goddard",])

all.VC3.fpc1601head <- cl(all[!is.na(all$fpc1601.pres) & all$beach!="Goddard",],fm=all.fit3.fpc1601head,
    cluster=all$hhid[!is.na(all$fpc1601.pres) & all$beach!="Goddard"])
overall.fit3.fpc1601head <- coeftest(all.fit3.fpc1601head, all.VC3.fpc1601head)
summary(all.fit3.fpc1601)
overall.fit3.fpc1601

# fpc 1602
all.fit3.fpc1602head <- glm(gici10~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1602.pres) & all$beach!="Malibu",])

all.VC3.fpc1602head <- cl(all[!is.na(all$fpc1602.pres) & all$beach!="Malibu",],fm=all.fit3.fpc1602head,
    cluster=all$hhid[!is.na(all$fpc1602.pres) & all$beach!="Malibu"])
overall.fit3.fpc1602head <- coeftest(all.fit3.fpc1602head, all.VC3.fpc1602head)
summary(all.fit3.fpc1602)
overall.fit3.fpc1602



# --------------------------------------
# save the results
# exclude glm objects and data frames
# (they are really large)
# --------------------------------------
save(

  all.n3.fmc1601head,all.n3.fmc1602head,all.n3.fpc1601head,all.n3.fpc1602head,

  av.n3.fmc1602head, av.n3.fpc1601head, av.n3.fpc1602head, do.n3.fmc1601head,
  do.n3.fmc1602head, do.n3.fpc1601head, do.n3.fpc1602head, fa.n3.fpc1601head,
  ma.n3.fpc1601head, mb.n3.fmc1601head, mb.n3.fpc1601head,

  av.n3.fmc1602.int0head, av.n3.fpc1601.int0head, av.n3.fpc1602.int0head, 
  do.n3.fmc1602.int0head, do.n3.fpc1601.int0head, do.n3.fpc1602.int0head, 
  ma.n3.fpc1601.int0head, 

  av.n3.fmc1602.int1head, av.n3.fpc1601.int1head, av.n3.fpc1602.int1head, 
  do.n3.fmc1602.int1head, do.n3.fpc1601.int1head, do.n3.fpc1602.int1head, 
  ma.n3.fpc1601.int1head,

  avfit3.fmc1602head, avfit3.fpc1601head, avfit3.fpc1602head,
  dofit3.fmc1601head, dofit3.fmc1602head, dofit3.fpc1601head, dofit3.fpc1602head,
  fafit3.fpc1601head, mafit3.fpc1601head,mbfit3.fmc1601head,mbfit3.fpc1601head,

  avfit3.fmc1602.gw.1.head, avfit3.fmc1602.gw.0.head, avfit3.fpc1601.gw.1.head, avfit3.fpc1601.gw.0.head, 
  avfit3.fpc1602.gw.1.head, avfit3.fpc1602.gw.0.head,
  dofit3.fmc1602.berm.1.head, dofit3.fmc1602.berm.0.head, 
  dofit3.fpc1601.berm.1.head, dofit3.fpc1601.berm.0.head, dofit3.fpc1602.berm.1.head, dofit3.fpc1602.berm.0.head,
  mafit3.fpc1601.berm.1.head,mafit3.fpc1601.berm.0.head,
  
  avfit3.fmc1602.gw.head, avfit3.fpc1601.gw.head, avfit3.fpc1602.gw.head,
  dofit3.fmc1602.berm.head, dofit3.fpc1601.berm.head, dofit3.fpc1602.berm.head,
  mafit3.fpc1601.berm.head,
  
  all.VC3.fmc1601head,overall.fit3.fmc1601head,
  all.VC3.fmc1602head,overall.fit3.fmc1602head,
  all.VC3.fpc1601head,overall.fit3.fpc1601head,
  all.VC3.fpc1602head,overall.fit3.fpc1602head,
  
  file="~/dropbox/coliphage/results/rawoutput/regress-3day-head.Rdata"
)



