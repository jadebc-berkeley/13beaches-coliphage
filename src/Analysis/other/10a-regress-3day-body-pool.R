##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file conducts maximum likelihood regression
# to estimate prevalence ratios

# Results pooled across beaches and assay

# 3 day gi illness
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

# drop individuals with no water quality information
all=subset(all,nowq==0)
# subset to non-missing exposure categories
# to make the robust CI calcs work
all=subset(all,all$bodycontact=="Yes")


# --------------------------------------
# Calculate the actual Ns for each cell
# and store them for plotting and tables
# --------------------------------------
regN <- function(outcome,exposurecat) {
  sum(table(outcome,exposurecat))
}

# n's pooled across assay and beach ---------------------------------------
all.n3.fmc = regN(all$gici3[!is.na(all$fmc.pres)],
                       all$fmc.pres[!is.na(all$fmc.pres)])
all.n3.fpc = regN(all$gici3[!is.na(all$fmc.pres)],
                       all$fmc.pres[!is.na(all$fmc.pres)])

# pooled n's by risk level---------------------------------------
data=all[!is.na(all$fmc.pres),]
data.high=subset(data,data$risk=="High")
all.n3.fmc.high = regN(data.high$gici3,data.high$fmc.pres)
data.low=subset(data,data$risk=="Low")
all.n3.fmc.low = regN(data.low$gici3,data.low$fmc.pres)

data=all[!is.na(all$fpc.pres),]
data.high=subset(data,data$risk=="High")
all.n3.fpc.high = regN(data.high$gici3,data.high$fpc.pres)
data.low=subset(data,data$risk=="Low")
all.n3.fpc.low = regN(data.low$gici3,data.low$fpc.pres)



# --------------------------------------
# Estimates pooled across beach and assay method
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness

# f- coliphage ----------------
all.fit3.fmc <- glm(gici3~fmc.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
                           rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc.pres),])

all.VC3.fmc <- cl(all[!is.na(all$fmc.pres)],fm=all.fit3.fmc,
                       cluster=all$hhid[!is.na(all$fmc.pres)])
overall.fit3.fmc <- coeftest(all.fit3.fmc, all.VC3.fmc)
summary(all.fit3.fmc)
overall.fit3.fmc
aic.fmc=AIC(all.fit3.fmc)

# f+ coliphage ----------------
all.fit3.fpc <- glm(gici3~fpc.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
                       rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc.pres),])

all.VC3.fpc <- cl(all[!is.na(all$fpc.pres)],fm=all.fit3.fpc,
                   cluster=all$hhid[!is.na(all$fpc.pres)])
overall.fit3.fpc <- coeftest(all.fit3.fpc, all.VC3.fpc)
summary(all.fit3.fpc)
overall.fit3.fpc
aic.fpc=AIC(all.fit3.fpc)


# --------------------------------------
# Estimates pooled across beach and stratified by conditions
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness

# all beaches ----------------
# f- coliphage --------
# high risk conditions
data=all[!is.na(all$fmc.pres),]
data.high=subset(data,data$risk=="High")
all.fit3.fmc.high <- glm(gici3~fmc.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC3.fmc.high <- cl(data.high,fm=all.fit3.fmc.high, cluster=data.high$hhid)
overall.fit3.fmc.high <- coeftest(all.fit3.fmc.high, all.VC3.fmc.high)
summary(all.fit3.fmc.high)
overall.fit3.fmc.high
aic.fmc.high=AIC(all.fit3.fmc.high)

# low risk conditions
data.low=subset(data,data$risk=="Low")
all.fit3.fmc.low <- glm(gici3~fmc.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC3.fmc.low <- cl(data.low,fm=all.fit3.fmc.low, cluster=data.low$hhid)
overall.fit3.fmc.low <- coeftest(all.fit3.fmc.low, all.VC3.fmc.low)
summary(all.fit3.fmc.low)
overall.fit3.fmc.low
aic.fmc.low=AIC(all.fit3.fmc.low)

# f+ coliphage  --------
# high risk conditions
data=all[!is.na(all$fpc.pres),]
data.high=subset(data,data$risk=="High")
all.fit3.fpc.high <- glm(gici3~fpc.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC3.fpc.high <- cl(data.high,fm=all.fit3.fpc.high, cluster=data.high$hhid)
overall.fit3.fpc.high <- coeftest(all.fit3.fpc.high, all.VC3.fpc.high)
summary(all.fit3.fpc.high)
overall.fit3.fpc.high
aic.fpc.high=AIC(all.fit3.fpc.high)

# low risk conditions
data.low=subset(data,data$risk=="Low")
all.fit3.fpc.low <- glm(gici3~fpc.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC3.fpc.low <- cl(data.low,fm=all.fit3.fpc.low, cluster=data.low$hhid)
overall.fit3.fpc.low <- coeftest(all.fit3.fpc.low, all.VC3.fpc.low)
summary(all.fit3.fpc.low)
overall.fit3.fpc.low
aic.fpc.low=AIC(all.fit3.fpc.low)




# --------------------------------------
# save the results
# exclude glm objects and data frames
# (they are really large)
# --------------------------------------
save(

  all.n3.fmc,all.n3.fpc,
  all.n3.fmc.high,all.n3.fmc.low,all.n3.fpc.high,all.n3.fpc.low,

  all.VC3.fmc, all.VC3.fpc,overall.fit3.fmc,overall.fit3.fpc,

  all.VC3.fmc.high,all.VC3.fpc.high,
  overall.fit3.fmc.high,overall.fit3.fpc.high,
  
  all.VC3.fmc.low,all.VC3.fpc.low,
  overall.fit3.fmc.low,overall.fit3.fpc.low,
  
  aic.fmc, aic.fpc,
  aic.fmc.high, aic.fpc.high,aic.fmc.low, aic.fpc.low,
  
  file="~/Documents/CRG/coliphage/results/rawoutput/regress-3day-body-pool.Rdata"
)

