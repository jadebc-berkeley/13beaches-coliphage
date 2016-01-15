##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file conducts maximum likelihood regression
# to estimate prevalence ratios
# for enterococcus

# Results pooled across beaches

# 3 day gi illness
##########################################

rm(list=ls())
library(foreign)
library(readstata13)

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

# pooled n's ---------------------------------------
all.n3.entero35 = regN(all$gici3,all$entero35)

data=all[!is.na(all$entero35),]
data.high=subset(data,data$groundwater=="Above median flow" | data$berm=="Open")
all.n3.entero35.high = regN(data.high$gici3,data.high$entero35)
  
data.low=subset(data,data$groundwater=="Below median flow" | data$berm=="Closed")
all.n3.entero35.low = regN(data.low$gici3,data.low$entero35)
  

# --------------------------------------
# Estimates pooled across beach
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness

# all beaches ----------------
all.fit3.entero <- glm(gici3~entero35+pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$entero35),])

all.VC3.entero <- cl(all[!is.na(all$entero35)],fm=all.fit3.entero,
                       cluster=all$hhid[!is.na(all$entero35)])
overall.fit3.entero <- coeftest(all.fit3.entero, all.VC3.entero)
summary(all.fit3.entero)
overall.fit3.entero

# high risk conditions
data=all[!is.na(all$entero35),]
data.high=subset(data,data$groundwater=="Above median flow" | data$berm=="Open")
all.fit3.entero.high <- glm(gici3~entero35+pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+
                                rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC3.entero.high <- cl(data.high,fm=all.fit3.entero.high, cluster=data.high$hhid)
overall.fit3.entero.high <- coeftest(all.fit3.entero.high, all.VC3.entero.high)
summary(all.fit3.entero.high)
overall.fit3.entero.high

# low risk conditions
data.low=subset(data,data$groundwater=="Below median flow" | data$berm=="Closed")
all.fit3.entero.low <- glm(gici3~entero35+pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+
                               rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC3.entero.low <- cl(data.low,fm=all.fit3.entero.low, cluster=data.low$hhid)
overall.fit3.entero.low <- coeftest(all.fit3.entero.low, all.VC3.entero.low)
summary(all.fit3.entero.low)
overall.fit3.entero.low

# --------------------------------------
# save the results
# exclude glm objects and data frames
# (they are really large)
# --------------------------------------
save(
  all.n3.entero35,all.n3.entero35.high,all.n3.entero35.low,
  
  overall.fit3.entero,overall.fit3.entero.high,overall.fit3.entero.low,
  all.VC3.entero,all.VC3.entero.high,all.VC3.entero.low,
  
  file="~/dropbox/coliphage/results/rawoutput/regress-3day-body-entero.Rdata"
)


