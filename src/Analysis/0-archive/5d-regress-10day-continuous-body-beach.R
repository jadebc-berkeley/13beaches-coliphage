  ##########################################
  # Coliphage analysis - 6 beaches
  # v1 by Jade 7/13/15
  
  # Estimate the RR across the range of concentration
  
  # Results stratified by beach
  
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
  

#-------------------------------------------
# estimating CIR for continuous exposure
#-------------------------------------------
# avalon
# fmc1601 -- always present at avalon 
avfit10.fmc1602.pY.gw.1 = mpreg(gici10~fmc1602+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon[avalon$groundwater=="Above median flow",], vcv=T)
avfit10.fmc1602.pY.gw.0 = mpreg(gici10~fmc1602+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon[avalon$groundwater=="Below median flow",], vcv=T)
  
avfit10.fpc1601.pY.gw.1 = mpreg(gici10~fpc1601+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon[avalon$groundwater=="Above median flow",], vcv=T)
avfit10.fpc1601.pY.gw.0 = mpreg(gici10~fpc1601+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon[avalon$groundwater=="Below median flow",], vcv=T)

avfit10.fpc1602.pY.gw.1 = mpreg(gici10~fpc1602+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon[avalon$groundwater=="Above median flow",], vcv=T)
avfit10.fpc1602.pY.gw.0 = mpreg(gici10~fpc1602+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon[avalon$groundwater=="Below median flow",], vcv=T)

# doheny
dofit10.fmc1601.pY = mpreg(gici10~fmc1601+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny[!is.na(doheny$fmc1601),],vcv=T)

dofit10.fmc1602.pY.berm.1 = mpreg(gici10~fmc1602+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny[doheny$berm=="Open",],vcv=T)
dofit10.fmc1602.pY.berm.0 = mpreg(gici10~fmc1602+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny[doheny$berm=="Closed",],vcv=T)


dofit10.fpc1601.pY.berm.1 = mpreg(gici10~fpc1601+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny[doheny$berm=="Open",],vcv=T)
dofit10.fpc1601.pY.berm.0 = mpreg(gici10~fpc1601+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny[doheny$berm=="Closed",],vcv=T)

dofit10.fpc1602.pY.berm.1 = mpreg(gici10~fpc1602+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny[doheny$berm=="Open",],vcv=T)
dofit10.fpc1602.pY.berm.0 = mpreg(gici10~fpc1602+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny[doheny$berm=="Closed",],vcv=T)

# fairhope
fafit10.fpc1601.pY = mpreg(gici10~fpc1601+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=fairhope[!is.na(fairhope$fpc1601),],vcv=T)

# goddard
gofit10.fpc1601.pY = mpreg(gici10~fpc1601+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=goddard,vcv=T)
  
# malibu
# fpc 1602 values all 0.1
  
mafit10.fmc1602.pY.berm.1 = mpreg(gici10~fmc1602+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=malibu[!is.na(malibu$fmc1602) & malibu$berm=="Open",],vcv=T)
mafit10.fmc1602.pY.berm.0 = mpreg(gici10~fmc1602+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=malibu[!is.na(malibu$fmc1602) & malibu$berm=="Closed",],vcv=T)
mafit10.fpc1601.pY.berm.1 = mpreg(gici10~fpc1601+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=malibu[!is.na(malibu$fpc1601) & malibu$berm=="Open",],vcv=T)
mafit10.fpc1601.pY.berm.0 = mpreg(gici10~fpc1601+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=malibu[!is.na(malibu$fpc1601) & malibu$berm=="Closed",],vcv=T)

# mission bay
mbfit10.fmc1601.pY = mpreg(gici10~fmc1601+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=mission[!is.na(mission$fmc1601),],vcv=T)
mbfit10.fpc1601.pY = mpreg(gici10~fpc1601+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=mission[!is.na(mission$fpc1601),],vcv=T)
  
#-------------------------------------------
# pooled
all.fit10.fpc1602 <- glm(gici10~fpc1602+pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1602) & all$beach!="Malibu",])

all.VC10.fpc1602 <- cl(all[!is.na(all$fpc1602) & all$beach!="Malibu",],fm=all.fit10.fpc1602,
    cluster=all$hhid[!is.na(all$fpc1602) & all$beach!="Malibu"])
overall.fit10.fpc1602 <- coeftest(all.fit10.fpc1602, all.VC10.fpc1602)

# high risk conditions
data=all[!is.na(all$fpc1602) & all$beach!="Malibu",]
data.high=subset(data,data$groundwater=="Above median flow" | data$berm=="Open")
all.fit10.fpc1602.high <- glm(gici10~fpc1602+pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC10.fpc1602.high <- cl(data.high,fm=all.fit10.fpc1602.high, cluster=data.high$hhid)
overall.fit10.fpc1602.high <- coeftest(all.fit10.fpc1602.high, all.VC10.fpc1602.high)

# low risk conditions
data.low=subset(data,data$groundwater=="Below median flow" | data$berm=="Closed")
all.fit10.fpc1602.low <- glm(gici10~fpc1602+pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC10.fpc1602.low <- cl(data.low,fm=all.fit10.fpc1602.low, cluster=data.low$hhid)
overall.fit10.fpc1602.low <- coeftest(all.fit10.fpc1602.low, all.VC10.fpc1602.low)

#-------------------------------------------
# estimating probability of illness over
# observed range of exposure
#-------------------------------------------
iter=1000

set.seed(92203789)

# avalon
# fmc1601 -- always present at avalon 
av.fmc1602.pY.gw.1 = boot.pY(fmla=gici10~fmc1602+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon[avalon$groundwater=="Above median flow",],nameX="fmc1602",
    ID=avalon$hhid[avalon$groundwater=="Above median flow"],iter)
av.fmc1602.pY.gw.0 = boot.pY(fmla=gici10~fmc1602+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon[avalon$groundwater=="Below median flow",],nameX="fmc1602",
    ID=avalon$hhid[avalon$groundwater=="Below median flow"],iter)  
  
av.fpc1601.pY.gw.1 = boot.pY(fmla=gici10~fpc1601+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon[avalon$groundwater=="Above median flow",],nameX="fpc1601",
    ID=avalon$hhid[avalon$groundwater=="Above median flow"],iter)
av.fpc1601.pY.gw.0 = boot.pY(fmla=gici10~fpc1601+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon[avalon$groundwater=="Below median flow",],nameX="fpc1601",
    ID=avalon$hhid[avalon$groundwater=="Below median flow"],iter)

av.fpc1602.pY.gw.1 = boot.pY(fmla=gici10~fpc1602+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon[avalon$groundwater=="Above median flow",],nameX="fpc1602",
    ID=avalon$hhid[avalon$groundwater=="Above median flow"],iter)
av.fpc1602.pY.gw.0 = boot.pY(fmla=gici10~fpc1602+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=avalon[avalon$groundwater=="Below median flow",],nameX="fpc1602",
    ID=avalon$hhid[avalon$groundwater=="Below median flow"],iter)

# doheny
do.fmc1601.pY = boot.pY(fmla=gici10~fmc1601+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny[!is.na(doheny$fmc1601),],nameX="fmc1601",ID=doheny[!is.na(doheny$fmc1601),"hhid"],iter)

do.fmc1602.pY.berm.1 = boot.pY(fmla=gici10~fmc1602+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny[doheny$berm=="Open",],nameX="fmc1602",ID=doheny$hhid[doheny$berm=="Open"],iter)
do.fmc1602.pY.berm.0 = boot.pY(fmla=gici10~fmc1602+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny[doheny$berm=="Closed",],nameX="fmc1602",ID=doheny$hhid[doheny$berm=="Closed"],iter)


do.fpc1601.pY.berm.1 = boot.pY(fmla=gici10~fpc1601+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny[doheny$berm=="Open",],nameX="fpc1601",ID=doheny$hhid[doheny$berm=="Open"],iter)
do.fpc1601.pY.berm.0 = boot.pY(fmla=gici10~fpc1601+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny[doheny$berm=="Closed",],nameX="fpc1601",ID=doheny$hhid[doheny$berm=="Closed"],iter)

do.fpc1602.pY.berm.1 = boot.pY(fmla=gici10~fpc1602+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny[doheny$berm=="Open",],nameX="fpc1602",ID=doheny$hhid[doheny$berm=="Open"],iter)
do.fpc1602.pY.berm.0 = boot.pY(fmla=gici10~fpc1602+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=doheny[doheny$berm=="Closed",],nameX="fpc1602",ID=doheny$hhid[doheny$berm=="Closed"],iter)

# fairhope
fa.fpc1601.pY = boot.pY(fmla=gici10~fpc1601+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=fairhope[!is.na(fairhope$fpc1601),],nameX="fpc1601",ID=fairhope[!is.na(fairhope$fpc1601),"hhid"],iter)

# goddard
go.fpc1601.pY = boot.pY(fmla=gici10~fpc1601+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=goddard,nameX="fpc1601",ID=goddard$hhid,iter)
  
  
# malibu
# fpc 1602 values all 0.1
ma.fmc1602.pY.berm.1 = boot.pY(fmla=gici10~fmc1602+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=malibu[!is.na(malibu$fmc1602) & malibu$berm=="Open",],nameX="fmc1602",ID=malibu[!is.na(malibu$fmc1602) & malibu$berm=="Open","hhid"],iter)
ma.fmc1602.pY.berm.0 = boot.pY(fmla=gici10~fmc1602+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=malibu[!is.na(malibu$fmc1602) & malibu$berm=="Closed",],nameX="fmc1602",ID=malibu[!is.na(malibu$fmc1602) & malibu$berm=="Closed","hhid"],iter)
ma.fpc1601.pY.berm.1 = boot.pY(fmla=gici10~fpc1601+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=malibu[!is.na(malibu$fpc1601) & malibu$berm=="Open",],nameX="fpc1601",ID=malibu[!is.na(malibu$fpc1601) & malibu$berm=="Open","hhid"],iter)
ma.fpc1601.pY.berm.0 = boot.pY(fmla=gici10~fpc1601+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=malibu[!is.na(malibu$fpc1601) & malibu$berm=="Closed",],nameX="fpc1601",ID=malibu[!is.na(malibu$fpc1601) & malibu$berm=="Closed","hhid"],iter)
  
# mission bay
mb.fmc1601.pY = boot.pY(fmla=gici10~fmc1601+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
    dat=mission[!is.na(mission$fmc1601),],nameX="fmc1601",ID=mission[!is.na(mission$fmc1601),"hhid"],iter)
mb.fpc1601.pY = boot.pY(fmla=gici10~fpc1601+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
               dat=mission[!is.na(mission$fpc1601),],nameX="fpc1601",ID=mission[!is.na(mission$fpc1601),"hhid"],iter)

# pooled
all.fpc1601.pY = boot.pY(fmla=gici10~fpc1601+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+rawfood+beach,dat=all[!is.na(all$fpc1601),],nameX="fpc1601",
    ID=all[!is.na(all$fpc1601),"hhid"],iter)
  
all.fpc1602.pY = boot.pY(fmla=gici10~fpc1602+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+rawfood+beach,dat=all[!is.na(all$fpc1602)& all$beach!="Malibu",],nameX="fpc1602",
    ID=all[!is.na(all$fpc1602)& all$beach!="Malibu","hhid"],iter)
data=all[!is.na(all$fpc1602) & all$beach!="Malibu",]
data.high=subset(data,data$groundwater=="Above median flow" | data$berm=="Open")
all.fpc1602.pY.high = boot.pY(fmla=gici10~fpc1602+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+ rawfood+beach,dat=data.high,nameX="fpc1602",ID=data.high$hhid,iter)
data.low=subset(data,data$groundwater=="Below median flow" | data$berm=="Closed")
all.fpc1602.pY.low = boot.pY(fmla=gici10~fpc1602+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+ rawfood+beach,dat=data.low,nameX="fmc1602",ID=data.low$hhid,iter)

save(
  
  avfit10.fmc1602.pY.gw.1, avfit10.fmc1602.pY.gw.0, avfit10.fpc1601.pY.gw.1, avfit10.fpc1601.pY.gw.0, 
  avfit10.fpc1602.pY.gw.1, avfit10.fpc1602.pY.gw.0,
  dofit10.fmc1601.pY,dofit10.fmc1602.pY.berm.1, dofit10.fmc1602.pY.berm.0, 
  dofit10.fpc1601.pY.berm.1, dofit10.fpc1601.pY.berm.0, dofit10.fpc1602.pY.berm.1, dofit10.fpc1602.pY.berm.0,

  fafit10.fpc1601.pY, gofit10.fpc1601.pY,mafit10.fmc1602.pY.berm.1, mafit10.fmc1602.pY.berm.0,
  mafit10.fpc1601.pY.berm.1,mafit10.fpc1601.pY.berm.0,
  mbfit10.fmc1601.pY,mbfit10.fpc1601.pY,
  
  av.fmc1602.pY.gw.1, av.fmc1602.pY.gw.0, av.fpc1601.pY.gw.1, av.fpc1601.pY.gw.0,
  av.fpc1602.pY.gw.1, av.fpc1602.pY.gw.0,
  do.fmc1601.pY, do.fmc1602.pY.berm.1, do.fmc1602.pY.berm.0, 
  do.fmc1602.pY.berm.1, do.fmc1602.pY.berm.0, do.fpc1601.pY.berm.1, do.fpc1601.pY.berm.0,
  do.fpc1602.pY.berm.1,do.fpc1602.pY.berm.0,
  fa.fpc1601.pY,go.fpc1601.pY, ma.fmc1602.pY.berm.1,ma.fmc1602.pY.berm.0,
  ma.fpc1601.pY.berm.1,ma.fpc1601.pY.berm.0,
  mb.fmc1601.pY,mb.fpc1601.pY,

  overall.fit10.fpc1602,overall.fit10.fpc1602.high,overall.fit10.fpc1602.low,
  
  file="~/dropbox/coliphage/results/rawoutput/regress-3day-continuous-body-beach.Rdata"
  
)












