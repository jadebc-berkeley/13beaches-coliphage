##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file conducts maximum likelihood regression
# to estimate prevalence ratios
# for enterococcus

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
beaches13=read.csv("~/Dropbox/13beaches/data/final/13beaches-analysis.csv")

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
  
# combined n's ---------------------------------------
#avalon
av.n10.entero35 = regN(avalon$gici10,avalon$entero35)
do.n10.entero35 = regN(doheny$gici10,doheny$entero35)
fa.n10.entero35 = regN(fairhope$gici10,fairhope$entero35)
go.n10.entero35 = regN(fairhope$gici10,fairhope$entero35)
ma.n10.entero35 = regN(malibu$gici10,malibu$entero35)
mb.n10.entero35 = regN(mission$gici10,mission$entero35)

# n if low risk conditions ---------------------------
av.n10.entero35.int0 = regN(avalon$gici3[avalon$groundwater=="Below median flow"],
                           avalon$entero35[avalon$groundwater=="Below median flow"])
do.n10.entero35.int0 = regN(doheny$gici3[doheny$berm=="Closed"],
                           doheny$entero35[doheny$berm=="Closed"])
ma.n10.entero35.int0 = regN(malibu$gici3[malibu$berm=="Closed"],
                           malibu$entero35[malibu$berm=="Closed"])

# n if high risk conditions ---------------------------
av.n10.entero35.int1 = regN(avalon$gici3[avalon$groundwater=="Above median flow"],
                           avalon$entero35[avalon$groundwater=="Above median flow"])
do.n10.entero35.int1 = regN(doheny$gici3[doheny$berm=="Open"],
                           doheny$entero35[doheny$berm=="Open"])
ma.n10.entero35.int1 = regN(malibu$gici3[malibu$berm=="Open"],
                           malibu$entero35[malibu$berm=="Open"])

#-----------------------------------------
# 10-day illness and enterococcus > 35
# pooled across berm and groundwater flow conditions
# by beach
#-----------------------------------------
avfit10.entero35 = mpreg(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
                        dat=avalon,vcv=T)
dofit10.entero35 = mpreg(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
                        dat=doheny,vcv=T)
fafit10.entero35 = mpreg(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
                        dat=fairhope,vcv=T)
# excluding goddard because entero35==1 in only 4 cases
mafit10.entero35 = mpreg(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
                        dat=malibu,vcv=T)
mbfit10.entero35 = mpreg(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
                        dat=mission[!is.na(mission$entero35),],vcv=T)


#-----------------------------------------
# 3-day illness and coliphage concentration
# by berm and groundwater flow conditions
# by beach

# interaction tests
#-----------------------------------------
# reduced model for LR test
avfit10.entero35.glm = glm(gici10~entero35+groundwater+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
                          dat=avalon,family=poisson(link="log"))
# interaction model
avfit10.entero35.gw = glm(gici10~entero35*groundwater+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
                         family=poisson(link="log"),data=avalon)
summary(avfit10.entero35.gw)
lrtest(avfit10.entero35.glm,avfit10.entero35.gw)

# berm always closed when fmc1601 measured so no interaction assessment

# reduced model for LR test
dofit10.entero35.glm = glm(gici10~entero35+berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
                          dat=doheny,family=poisson(link="log"))
# interaction model
dofit10.entero35.berm = glm(gici10~entero35*berm+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
                           dat=doheny,family=poisson(link="log"))
summary(dofit10.entero35.berm)
lrtest(dofit10.entero35.glm,dofit10.entero35.berm)

# no stratification for malibu because only 25 cases of entero35==1

#-----------------------------------------
# 3-day illness and coliphage concentration
# by berm and groundwater flow conditions
# by beach

# stratified estimates
#-----------------------------------------
avfit10.entero35.gw.1 = mpreg(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
                             vcv=T,dat=avalon[avalon$groundwater=="Above median flow",])
avfit10.entero35.gw.0 = mpreg(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
                             vcv=T,dat=avalon[avalon$groundwater=="Below median flow",])

dofit10.entero35.berm.1 = mpreg(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
                               vcv=T,dat=doheny[doheny$berm=="Open",])
dofit10.entero35.berm.0 = mpreg(gici10~entero35+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,
                               vcv=T,dat=doheny[doheny$berm=="Closed",])


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
  
  av.n10.entero35, do.n10.entero35, fa.n10.entero35,
  go.n10.entero35, ma.n10.entero35, mb.n10.entero35,

  av.n10.entero35.int0, do.n10.entero35.int0, ma.n10.entero35.int0, 

  av.n10.entero35.int1, do.n10.entero35.int1, ma.n10.entero35.int1, 
  
  avfit10.entero35, dofit10.entero35, fafit10.entero35, 
  mafit10.entero35,mbfit10.entero35,
  
  avfit10.entero35.gw.1, avfit10.entero35.gw.0,
  dofit10.entero35.berm.1, dofit10.entero35.berm.0, 
  
  avfit10.entero35.gw, dofit10.entero35.berm,
  
  overall.fit10.entero,overall.fit10.entero.high,overall.fit10.entero.low,
  all.VC10.entero,all.VC10.entero.high,all.VC10.entero.low,
  
  aic.entero,aic.entero.low,aic.entero.high,
  
  file="~/Documents/CRG/coliphage/results/rawoutput/regress-10day-body-entero-beach.Rdata"
)


