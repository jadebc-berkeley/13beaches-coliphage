##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file conducts maximum likelihood regression
# to estimate prevalence ratios

# Results stratified by beach

# 10 day gi illness, head immersion
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
all.n10.fmc1601head = regN(all$gici10[!is.na(all$fmc1601.pres)],
                       all$fmc1601.pres[!is.na(all$fmc1601.pres)])
all.n10.fmc1602head = regN(all$gici10[!is.na(all$fmc1602.pres)],
                       all$fmc1602.pres[!is.na(all$fmc1602.pres)])
all.n10.fpc1601head = regN(all$gici10[!is.na(all$fpc1601.pres) ],
                       all$fpc1601.pres[!is.na(all$fpc1601.pres) ])
all.n10.fpc1602head = regN(all$gici10[!is.na(all$fpc1602.pres) ],
                       all$fpc1602.pres[!is.na(all$fpc1602.pres) ])


# combined n's ---------------------------------------
#avalon
av.n10.fmc1602head = regN(avalon$gici10,avalon$fmc1602.pres)
av.n10.fpc1601head = regN(avalon$gici10,avalon$fpc1601.pres)
av.n10.fpc1602head = regN(avalon$gici10,avalon$fpc1602.pres)

#doheny
do.n10.fmc1601head = regN(doheny$gici10,doheny$fmc1601.pres)
do.n10.fmc1602head = regN(doheny$gici10,doheny$fmc1602.pres)
do.n10.fpc1601head = regN(doheny$gici10,doheny$fpc1601.pres)
do.n10.fpc1602head = regN(doheny$gici10,doheny$fpc1602.pres)

# fairhope
fa.n10.fpc1601head = regN(fairhope$gici10,fairhope$fpc1601.pres)

# goddard
# fpc 1601 always present at goddard

# malibu
# fmc 1602 and fpc 1602 always present at malibu
ma.n10.fpc1601head = regN(malibu$gici10,malibu$fpc1601.pres)

# mission bay
mb.n10.fmc1601head = regN(mission$gici10,mission$fmc1601.pres)
mb.n10.fpc1601head = regN(mission$gici10,mission$fpc1601.pres)

# n if low risk conditions ---------------------------
#avalon
av.n10.fmc1602.int0head = regN(avalon$gici10[avalon$groundwater=="Below median flow"],
  avalon$fmc1602.pres[avalon$groundwater=="Below median flow"])
av.n10.fpc1601.int0head = regN(avalon$gici10[avalon$groundwater=="Below median flow"],
  avalon$fpc1601.pres[avalon$groundwater=="Below median flow"])
av.n10.fpc1602.int0head = regN(avalon$gici10[avalon$groundwater=="Below median flow"],
  avalon$fpc1602.pres[avalon$groundwater=="Below median flow"])

#doheny
do.n10.fmc1602.int0head = regN(doheny$gici10[doheny$berm=="Closed"],
  doheny$fmc1602.pres[doheny$berm=="Closed"])
do.n10.fpc1601.int0head = regN(doheny$gici10[doheny$berm=="Closed"],
  doheny$fpc1601.pres[doheny$berm=="Closed"])
do.n10.fpc1602.int0head = regN(doheny$gici10[doheny$berm=="Closed"],
  doheny$fpc1602.pres[doheny$berm=="Closed"])

# malibu
# fmc 1602 and fpc 1602 always present at malibu
ma.n10.fpc1601.int0head = regN(malibu$gici10[malibu$berm=="Closed"],malibu$fpc1601.pres[malibu$berm=="Closed"])

# n if high risk conditions ---------------------------
#avalon
av.n10.fmc1602.int1head = regN(avalon$gici10[avalon$groundwater=="Above median flow"],
  avalon$fmc1602.pres[avalon$groundwater=="Above median flow"])
av.n10.fpc1601.int1head = regN(avalon$gici10[avalon$groundwater=="Above median flow"],
  avalon$fpc1601.pres[avalon$groundwater=="Above median flow"])
av.n10.fpc1602.int1head = regN(avalon$gici10[avalon$groundwater=="Above median flow"],
  avalon$fpc1602.pres[avalon$groundwater=="Above median flow"])

#doheny
do.n10.fmc1602.int1head = regN(doheny$gici10[doheny$berm=="Open"],
  doheny$fmc1602.pres[doheny$berm=="Open"])
do.n10.fpc1601.int1head = regN(doheny$gici10[doheny$berm=="Open"],
  doheny$fpc1601.pres[doheny$berm=="Open"])
do.n10.fpc1602.int1head = regN(doheny$gici10[doheny$berm=="Open"],
  doheny$fpc1602.pres[doheny$berm=="Open"])

# malibu
# fmc 1602 and fpc 1602 always present at malibu
ma.n10.fpc1601.int1head = regN(malibu$gici10[malibu$berm=="Open"],malibu$fpc1601.pres[malibu$berm=="Open"])

#-----------------------------------------
# 3-day illness and coliphage concentration
# pooled across berm and groundwater flow conditions
# by beach
#-----------------------------------------
# avalon
# fmc1601 -- always present at avalon 
avfit10.fmc1602head = mpreg(gici10~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,vcv=T)
avfit10.fpc1601head = mpreg(gici10~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,vcv=T)
avfit10.fpc1602head = mpreg(gici10~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,vcv=T)

# doheny
dofit10.fmc1601head = mpreg(gici10~fmc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
     dat=doheny[!is.na(doheny$fmc1601.pres),],vcv=T)
dofit10.fmc1602head = mpreg(gici10~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,vcv=T)
dofit10.fpc1601head = mpreg(gici10~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,vcv=T)
dofit10.fpc1602head = mpreg(gici10~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,vcv=T)

# fairhope
fafit10.fpc1601head = mpreg(gici10~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=fairhope[!is.na(fairhope$fpc1601.pres),],vcv=T)

# goddard
# fpc 1601 always present at goddard

# malibu
# fmc 1602 and fpc 1602 always present at malibu
mafit10.fpc1601head = mpreg(gici10~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=malibu[!is.na(malibu$fpc1601.pres),],vcv=T)

# mission bay
mbfit10.fmc1601head = mpreg(gici10~fmc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=mission[!is.na(mission$fmc1601.pres),],vcv=T)
mbfit10.fpc1601head = mpreg(gici10~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=mission[!is.na(mission$fpc1601.pres),],vcv=T)


#-----------------------------------------
# 3-day illness and coliphage concentration
# by berm and groundwater flow conditions
# by beach

# interaction tests
#-----------------------------------------
# reduced model for LR test
avfit10.fmc1602.glm.head = glm(gici10~fmc1602.pres+groundwater+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,family=poisson(link="log"))
# interaction model
avfit10.fmc1602.gw.head = glm(gici10~fmc1602.pres*groundwater+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    family=poisson(link="log"),data=avalon)
summary(avfit10.fmc1602.gw.head)
lrtest(avfit10.fmc1602.glm.head,avfit10.fmc1602.gw.head)

# reduced model for LR test
avfit10.fpc1601.glm.head = glm(gici10~fpc1601.pres+groundwater+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,family=poisson(link="log"))
# interaction model
avfit10.fpc1601.gw.head = glm(gici10~fpc1601.pres*groundwater+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,family=poisson(link="log"))
summary(avfit10.fpc1601.gw.head)
lrtest(avfit10.fpc1601.glm.head,avfit10.fpc1601.gw.head)

# reduced model for LR test
avfit10.fpc1602.glm.head = glm(gici10~fpc1602.pres+groundwater+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,family=poisson(link="log"))
# interaction model
avfit10.fpc1602.gw.head = glm(gici10~fpc1602.pres*groundwater+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,family=poisson(link="log"))
summary(avfit10.fpc1602.gw.head)
lrtest(avfit10.fpc1602.glm.head,avfit10.fpc1602.gw.head)

# berm always closed when fmc1601 measured so no interaction assessment

# reduced model for LR test
dofit10.fmc1602.glm.head = glm(gici10~fmc1602.pres+berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,family=poisson(link="log"))
# interaction model
dofit10.fmc1602.berm.head = glm(gici10~fmc1602.pres*berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,family=poisson(link="log"))
summary(dofit10.fmc1602.berm.head)
lrtest(dofit10.fmc1602.glm.head,dofit10.fmc1602.berm.head)

# reduced model for LR test
dofit10.fpc1601.glm.head = glm(gici10~fpc1601.pres+berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,family=poisson(link="log"))
# interaction model
dofit10.fpc1601.berm.head = glm(gici10~fpc1601.pres*berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,family=poisson(link="log"))
summary(dofit10.fpc1601.berm.head)
lrtest(dofit10.fpc1601.glm.head,dofit10.fpc1601.berm.head)

# reduced model for LR test
dofit10.fpc1602.glm.head = glm(gici10~fpc1602.pres+berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,family=poisson(link="log"))
# interaction model
dofit10.fpc1602.berm.head = glm(gici10~fpc1602.pres*berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,family=poisson(link="log"))
summary(dofit10.fpc1602.berm.head)
lrtest(dofit10.fpc1602.glm.head,dofit10.fpc1602.berm.head)

# reduced model for LR test
mafit10.fpc1601.glm.head = glm(gici10~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=malibu[!is.na(malibu$fpc1601.pres),],family=poisson(link="log"))
# interaction model
mafit10.fpc1601.berm.head = glm(gici10~fpc1601.pres*berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=malibu[!is.na(malibu$fpc1601.pres),],family=poisson(link="log"))
summary(mafit10.fpc1601.berm.head)
lrtest(mafit10.fpc1601.glm.head,mafit10.fpc1601.berm.head)


#-----------------------------------------
# 3-day illness and coliphage concentration
# by berm and groundwater flow conditions
# by beach

# stratified estimates
#-----------------------------------------
avfit10.fmc1602.gw.1.head = mpreg(gici10~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=avalon[avalon$groundwater=="Above median flow",])
avfit10.fmc1602.gw.0.head = mpreg(gici10~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=avalon[avalon$groundwater=="Below median flow",])

avfit10.fpc1601.gw.1.head = mpreg(gici10~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=avalon[avalon$groundwater=="Above median flow",])
avfit10.fpc1601.gw.0.head = mpreg(gici10~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=avalon[avalon$groundwater=="Below median flow",])

avfit10.fpc1602.gw.1.head = mpreg(gici10~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=avalon[avalon$groundwater=="Above median flow",])
avfit10.fpc1602.gw.0.head = mpreg(gici10~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=avalon[avalon$groundwater=="Below median flow",])

dofit10.fmc1602.berm.1.head = mpreg(gici10~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=doheny[doheny$berm=="Open",])
dofit10.fmc1602.berm.0.head = mpreg(gici10~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=doheny[doheny$berm=="Closed",])

dofit10.fpc1601.berm.1.head = mpreg(gici10~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=doheny[doheny$berm=="Open",])
dofit10.fpc1601.berm.0.head = mpreg(gici10~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=doheny[doheny$berm=="Closed",])

dofit10.fpc1602.berm.1.head = mpreg(gici10~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=doheny[doheny$berm=="Open",])
dofit10.fpc1602.berm.0.head = mpreg(gici10~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=doheny[doheny$berm=="Closed",])

mafit10.fpc1601.berm.1.head = mpreg(gici10~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=malibu[malibu$berm=="Open",])
mafit10.fpc1601.berm.0.head = mpreg(gici10~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=malibu[!is.na(malibu$fpc1601.pres) & malibu$berm=="Closed",])

# --------------------------------------
# Pooled estimates
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 3-day GI illness

# all beaches ----------------
# fmc 1601
all.fit10.fmc1601head <- glm(gici10~fmc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc1601.pres) ,])

all.VC10.fmc1601head <- cl(all[!is.na(all$fmc1601.pres) ],fm=all.fit10.fmc1601head,
  cluster=all$hhid[!is.na(all$fmc1601.pres) ])
overall.fit10.fmc1601head <- coeftest(all.fit10.fmc1601head, all.VC10.fmc1601head)
summary(all.fit10.fmc1601head)
overall.fit10.fmc1601head

# fmc 1602
all.fit10.fmc1602head <- glm(gici10~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc1602.pres) ,])

all.VC10.fmc1602head <- cl(all[!is.na(all$fmc1602.pres) ,],fm=all.fit10.fmc1602head,
  cluster=all$hhid[!is.na(all$fmc1602.pres) ])
overall.fit10.fmc1602head <- coeftest(all.fit10.fmc1602head, all.VC10.fmc1602head)
summary(all.fit10.fmc1602head)
overall.fit10.fmc1602head

# fpc 1601
all.fit10.fpc1601head <- glm(gici10~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1601.pres) ,])

all.VC10.fpc1601head <- cl(all[!is.na(all$fpc1601.pres) ,],fm=all.fit10.fpc1601head,
    cluster=all$hhid[!is.na(all$fpc1601.pres) ])
overall.fit10.fpc1601head <- coeftest(all.fit10.fpc1601head, all.VC10.fpc1601head)
summary(all.fit10.fpc1601head)
overall.fit10.fpc1601head

# fpc 1602
all.fit10.fpc1602head <- glm(gici10~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1602.pres) ,])

all.VC10.fpc1602head <- cl(all[!is.na(all$fpc1602.pres) ,],fm=all.fit10.fpc1602head,
    cluster=all$hhid[!is.na(all$fpc1602.pres) ])
overall.fit10.fpc1602head <- coeftest(all.fit10.fpc1602head, all.VC10.fpc1602head)
summary(all.fit10.fpc1602head)
all.VC10.fpc1602head

# --------------------------------------
# save the results
# exclude glm objects and data frames
# (they are really large)
# --------------------------------------
save(

  all.n10.fmc1601head,all.n10.fmc1602head,all.n10.fpc1601head,all.n10.fpc1602head,

  av.n10.fmc1602head, av.n10.fpc1601head, av.n10.fpc1602head, do.n10.fmc1601head,
  do.n10.fmc1602head, do.n10.fpc1601head, do.n10.fpc1602head, fa.n10.fpc1601head,
  ma.n10.fpc1601head, mb.n10.fmc1601head, mb.n10.fpc1601head,

  av.n10.fmc1602.int0head, av.n10.fpc1601.int0head, av.n10.fpc1602.int0head, 
  do.n10.fmc1602.int0head, do.n10.fpc1601.int0head, do.n10.fpc1602.int0head, 
  ma.n10.fpc1601.int0head,

  av.n10.fmc1602.int1head, av.n10.fpc1601.int1head, av.n10.fpc1602.int1head, 
  do.n10.fmc1602.int1head, do.n10.fpc1601.int1head, do.n10.fpc1602.int1head,
  ma.n10.fpc1601.int1head, 

  avfit10.fmc1602head, avfit10.fpc1601head, avfit10.fpc1602head,
  dofit10.fmc1601head, dofit10.fmc1602head, dofit10.fpc1601head, dofit10.fpc1602head,
  fafit10.fpc1601head, mafit10.fpc1601head,mbfit10.fmc1601head,mbfit10.fpc1601head,

  avfit10.fmc1602.gw.1.head, avfit10.fmc1602.gw.0.head, avfit10.fpc1601.gw.1.head, avfit10.fpc1601.gw.0.head, 
  avfit10.fpc1602.gw.1.head, avfit10.fpc1602.gw.0.head,
  dofit10.fmc1602.berm.1.head, dofit10.fmc1602.berm.0.head, 
  dofit10.fpc1601.berm.1.head, dofit10.fpc1601.berm.0.head, dofit10.fpc1602.berm.1.head, dofit10.fpc1602.berm.0.head,
  mafit10.fpc1601.berm.1.head,mafit10.fpc1601.berm.0.head,
  
  avfit10.fmc1602.gw.head, avfit10.fpc1601.gw.head, avfit10.fpc1602.gw.head,
  dofit10.fmc1602.berm.head, dofit10.fpc1601.berm.head, dofit10.fpc1602.berm.head,
  mafit10.fpc1601.berm.head,
  
  all.VC10.fmc1601head,overall.fit10.fmc1601head,
  all.VC10.fmc1602head,overall.fit10.fmc1602head,
  all.VC10.fpc1601head,overall.fit10.fpc1601head,
  all.VC10.fpc1602head,overall.fit10.fpc1602head,
  
  file="~/dropbox/coliphage/results/rawoutput/regress-10day-head.Rdata"
)



