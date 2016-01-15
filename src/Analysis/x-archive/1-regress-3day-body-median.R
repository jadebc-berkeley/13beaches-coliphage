##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file conducts maximum likelihood regression
# to estimate prevalence ratios

# 3 day gi illness, body immersion
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
beaches13=read.dta("~/Dropbox/13beaches/data/final/13beaches-analysis.dta")

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

# combined n's
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

#-----------------------------------------
# 3-day illness and coliphage concentration
# pooled across berm and groundwater flow conditions
# by beach
#-----------------------------------------
# avalon
# fmc1601 -- always present at avalon 
avfit3.fmc1602 = mpreg(gici3~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,vcv=T)
avfit3.fpc1601 = mpreg(gici3~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,vcv=T)
avfit3.fpc1602 = mpreg(gici3~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,vcv=T)

# doheny
dofit3.fmc1601 = mpreg(gici3~fmc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
     dat=doheny[!is.na(doheny$fmc1601.pres),],vcv=T)
dofit3.fmc1602 = mpreg(gici3~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,vcv=T)
dofit3.fpc1601 = mpreg(gici3~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,vcv=T)
dofit3.fpc1602 = mpreg(gici3~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,vcv=T)

# fairhope
fafit3.fpc1601 = mpreg(gici3~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=fairhope[!is.na(fairhope$fpc1601.pres),],vcv=T)

# goddard
# fpc 1601 always present at goddard

# malibu
# fmc 1602 and fpc 1602 always present at malibu
mafit3.fpc1601 = mpreg(gici3~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=malibu[!is.na(malibu$fpc1601.pres),],vcv=T)

# mission bay
mbfit3.fmc1601 = mpreg(gici3~fmc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=mission[!is.na(mission$fmc1601.pres),],vcv=T)
mbfit3.fpc1601 = mpreg(gici3~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
               dat=mission[!is.na(mission$fpc1601.pres),],vcv=T)


#-----------------------------------------
# 3-day illness and coliphage concentration
# by berm and groundwater flow conditions
# by beach

# interaction tests
#-----------------------------------------
# reduced model for LR test
avfit3.fmc1602.glm = glm(gici3~fmc1602.pres+groundwater+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,family=poisson(link="log"))
# interaction model
avfit3.fmc1602.gw = glm(gici3~fmc1602.pres*groundwater+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    family=poisson(link="log"),data=avalon)
summary(avfit3.fmc1602.gw)
lrtest(avfit3.fmc1602.glm,avfit3.fmc1602.gw)

# reduced model for LR test
avfit3.fpc1601.glm = glm(gici3~fpc1601.pres+groundwater+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,family=poisson(link="log"))
# interaction model
avfit3.fpc1601.gw = glm(gici3~fpc1601.pres*groundwater+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,family=poisson(link="log"))
summary(avfit3.fpc1601.gw)
lrtest(avfit3.fpc1601.glm,avfit3.fpc1601.gw)

# reduced model for LR test
avfit3.fpc1602.glm = glm(gici3~fpc1602.pres+groundwater+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,family=poisson(link="log"))
# interaction model
avfit3.fpc1602.gw = glm(gici3~fpc1602.pres*groundwater+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon,family=poisson(link="log"))
summary(avfit3.fpc1602.gw)
lrtest(avfit3.fpc1602.glm,avfit3.fpc1602.gw)

# berm always closed when fmc1601 measured so no interaction assessment

# reduced model for LR test
dofit3.fmc1602.glm = glm(gici3~fmc1602.pres+berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,family=poisson(link="log"))
# interaction model
dofit3.fmc1602.berm = glm(gici3~fmc1602.pres*berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,family=poisson(link="log"))
summary(dofit3.fmc1602.berm)
lrtest(dofit3.fmc1602.glm,dofit3.fmc1602.berm)

# reduced model for LR test
dofit3.fpc1601.glm = glm(gici3~fpc1601.pres+berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,family=poisson(link="log"))
# interaction model
dofit3.fpc1601.berm = glm(gici3~fpc1601.pres*berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,family=poisson(link="log"))
summary(dofit3.fpc1601.berm)
lrtest(dofit3.fpc1601.glm,dofit3.fpc1601.berm)

# reduced model for LR test
dofit3.fpc1602.glm = glm(gici3~fpc1602.pres+berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,family=poisson(link="log"))
# interaction model
dofit3.fpc1602.berm = glm(gici3~fpc1602.pres*berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny,family=poisson(link="log"))
summary(dofit3.fpc1602.berm)
lrtest(dofit3.fpc1602.glm,dofit3.fpc1602.berm)

#-----------------------------------------
# 3-day illness and coliphage concentration
# by berm and groundwater flow conditions
# by beach

# stratified estimates
#-----------------------------------------
avfit3.fmc1602.gw.1 = mpreg(gici3~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=avalon[avalon$groundwater=="Above median flow",])
avfit3.fmc1602.gw.0 = mpreg(gici3~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=avalon[avalon$groundwater=="Below median flow",])

avfit3.fpc1601.gw.1 = mpreg(gici3~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=avalon[avalon$groundwater=="Above median flow",])
avfit3.fpc1601.gw.0 = mpreg(gici3~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=avalon[avalon$groundwater=="Below median flow",])

avfit3.fpc1602.gw.1 = mpreg(gici3~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=avalon[avalon$groundwater=="Above median flow",])
avfit3.fpc1602.gw.0 = mpreg(gici3~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=avalon[avalon$groundwater=="Below median flow",])

dofit3.fmc1602.berm.1 = mpreg(gici3~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=doheny[doheny$berm=="Open",])
dofit3.fmc1602.berm.0 = mpreg(gici3~fmc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=doheny[doheny$berm=="Closed",])

dofit3.fpc1601.berm.1 = mpreg(gici3~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=doheny[doheny$berm=="Open",])
dofit3.fpc1601.berm.0 = mpreg(gici3~fpc1601.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=doheny[doheny$berm=="Closed",])

dofit3.fpc1602.berm.1 = mpreg(gici3~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=doheny[doheny$berm=="Open",])
dofit3.fpc1602.berm.0 = mpreg(gici3~fpc1602.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    vcv=T,dat=doheny[doheny$berm=="Closed",])


# --------------------------------------
# Pooled estimates and test of interaction 
# for point source

# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# all beaches ----------------
# fmc 1601
all.fit3.fmc1601 <- glm(gici3~fmc1601.pres+pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc1601.pres),])

all.VC3.fmc1601 <- cl(all[!is.na(all$fmc1601.pres)],fm=all.fit3.fmc1601,
	cluster=all$hhid[!is.na(all$fmc1601.pres)])
overall.fit3.fmc1601 <- coeftest(all.fit3.fmc1601, all.VC3.fmc1601)
summary(all.fit3.fmc1601)
overall.fit3.fmc1601

# Interaction model with point v. non-point source beaches
ps.fit3.fmc1601 <- glm(gici3~fmc1601.pres*pointsource+agecat+female+racewhite+gichron+anim_any+
                gicontactbase+rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc1601.pres),])

summary(ps.fit3.fmc1601)
lrtest(all.fit3.fmc1601,ps.fit3.fmc1601)


# fmc 1602 ----------------
all.fit3.fmc1602 <- glm(gici3~fmc1602.pres+pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc1602.pres),])

all.VC3.fmc1602 <- cl(all[!is.na(all$fmc1602.pres)],fm=all.fit3.fmc1602,
	cluster=all$hhid[!is.na(all$fmc1602.pres)])
overall.fit3.fmc1602 <- coeftest(all.fit3.fmc1602, all.VC3.fmc1602)
summary(all.fit3.fmc1602)
overall.fit3.fmc1602

# Interaction model with point v. non-point source beaches
ps.fit3.fmc1602 <- glm(gici3~fmc1602.pres*pointsource+agecat+female+racewhite+gichron+anim_any+
                gicontactbase+rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc1602.pres),])

summary(ps.fit3.fmc1602)
lrtest(all.fit3.fmc1602,ps.fit3.fmc1602)

# fpc 1601 ----------------
all.fit3.fpc1601 <- glm(gici3~fpc1601.pres+pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1601.pres),])

all.VC3.fpc1601 <- cl(all[!is.na(all$fpc1601.pres)],fm=all.fit3.fpc1601,
    cluster=all$hhid[!is.na(all$fpc1601.pres)])
overall.fit3.fpc1601 <- coeftest(all.fit3.fpc1601, all.VC3.fpc1601)
summary(all.fit3.fpc1601)
overall.fit3.fpc1601

# Interaction model with point v. non-point source beaches
ps.fit3.fpc1601 <- glm(gici3~fpc1601.pres*pointsource+agecat+female+racewhite+gichron+anim_any+
                gicontactbase+rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1601.pres),])

summary(ps.fit3.fpc1601)
lrtest(all.fit3.fpc1601,ps.fit3.fpc1601)

# fpc 1602 ----------------
all.fit3.fpc1602 <- glm(gici3~fpc1602.pres+pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1602.pres),])

all.VC3.fpc1602 <- cl(all[!is.na(all$fpc1602.pres)],fm=all.fit3.fpc1602,
    cluster=all$hhid[!is.na(all$fpc1602.pres)])
overall.fit3.fpc1602 <- coeftest(all.fit3.fpc1602, all.VC3.fpc1602)
summary(all.fit3.fpc1602)
overall.fit3.fpc1602

# Interaction model with point v. non-point source beaches
ps.fit3.fpc1602 <- glm(gici3~fpc1602.pres*pointsource+agecat+female+racewhite+gichron+anim_any+
                gicontactbase+rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1602.pres),])

summary(ps.fit3.fpc1602)
lrtest(all.fit3.fpc1602,ps.fit3.fpc1602)



# --------------------------------------
# save the results
# exclude glm objects and data frames
# (they are really large)
# --------------------------------------
save(

  av.n3.fmc1602, av.n3.fpc1601, av.n3.fpc1602, do.n3.fmc1601,
  do.n3.fmc1602, do.n3.fpc1601, do.n3.fpc1602, fa.n3.fpc1601,
  ma.n3.fpc1601, mb.n3.fmc1601, mb.n3.fpc1601,

  avfit3.fmc1602, avfit3.fpc1601, avfit3.fpc1602,
  dofit3.fmc1601, dofit3.fmc1602, dofit3.fpc1601, dofit3.fpc1602,
  fafit3.fpc1601, mafit3.fpc1601,mbfit3.fmc1601,mbfit3.fpc1601,

  avfit3.fmc1602.gw.1, avfit3.fmc1602.gw.0, avfit3.fpc1601.gw.1, avfit3.fpc1601.gw.0, 
  avfit3.fpc1602.gw.1, avfit3.fpc1602.gw.0,
  dofit3.fmc1602.berm.1, dofit3.fmc1602.berm.0, 
  dofit3.fpc1601.berm.1, dofit3.fpc1601.berm.0, dofit3.fpc1602.berm.1, dofit3.fpc1602.berm.0,
  
  avfit3.fmc1602.gw, avfit3.fpc1601.gw, avfit3.fpc1602.gw,
  dofit3.fmc1602.berm, dofit3.fpc1601.berm, dofit3.fpc1602.berm,
  
  all.fit3.fmc1601,all.VC3.fmc1601,overall.fit3.fmc1601,
  all.fit3.fmc1602,all.VC3.fmc1602,overall.fit3.fmc1602,
  all.fit3.fpc1601,all.VC3.fpc1601,overall.fit3.fpc1601,
  all.fit3.fpc1602,all.VC3.fpc1602,overall.fit3.fpc1602,
  
  ps.fit3.fmc1601,ps.fit3.fmc1602,ps.fit3.fpc1601,ps.fit3.fpc1602,

  file="~/dropbox/coliphage/results/rawoutput/regress-3day-body.Rdata"
)



