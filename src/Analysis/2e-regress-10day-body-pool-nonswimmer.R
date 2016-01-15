##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file conducts maximum likelihood regression
# to estimate prevalence ratios

# Results pooled across beaches and assay

# Non-swimmers as reference group

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

# drop individuals with no water quality information
all=subset(all,nowq==0)
# subset to non-missing exposure categories
# to make the robust CI calcs work
all=subset(all,all$bodycontact=="Yes" | all$anycontact=="No")



# --------------------------------------
# Calculate the actual Ns for each cell
# and store them for plotting and tables
# --------------------------------------
regN <- function(outcome,exposurecat) {
  sum(table(outcome,exposurecat))
}

# n's pooled across assay and beach ---------------------------------------
all.n10.fmc = regN(all$gici10[!is.na(all$fmc.pres)],
                       all$fmc.pres[!is.na(all$fmc.pres)])
all.n10.fpc = regN(all$gici10[!is.na(all$fmc.pres)],
                       all$fmc.pres[!is.na(all$fmc.pres)])

# pooled n's by risk level---------------------------------------
data=all[!is.na(all$fmc.pres),]
data.high=subset(data,data$risk=="High")
all.n10.fmc.high = regN(data.high$gici10,data.high$fmc.pres)
data.low=subset(data,data$risk=="Low")
all.n10.fmc.low = regN(data.low$gici10,data.low$fmc.pres)

data=all[!is.na(all$fpc.pres),]
data.high=subset(data,data$risk=="High")
all.n10.fpc.high = regN(data.high$gici10,data.high$fpc.pres)
data.low=subset(data,data$risk=="Low")
all.n10.fpc.low = regN(data.low$gici10,data.low$fpc.pres)



# --------------------------------------
# Estimates pooled across beach and assay method
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness

# f- coliphage ----------------
all.fit10.fmc <- glm(gici10~swim.fmc+agecat+female+racewhite+gichron+anim_any+gicontactbase+
     rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$swim.fmc),])

all.VC10.fmc <- cl(all[!is.na(all$swim.fmc)],fm=all.fit10.fmc,
                       cluster=all$hhid[!is.na(all$swim.fmc)])
overall.fit10.fmc <- coeftest(all.fit10.fmc, all.VC10.fmc)
summary(all.fit10.fmc)
overall.fit10.fmc

# f+ coliphage ----------------
all.fit10.fpc <- glm(gici10~swim.fpc+agecat+female+racewhite+gichron+anim_any+gicontactbase+
                       rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$swim.fpc),])

all.VC10.fpc <- cl(all[!is.na(all$swim.fpc)],fm=all.fit10.fpc,
                   cluster=all$hhid[!is.na(all$swim.fpc)])
overall.fit10.fpc <- coeftest(all.fit10.fpc, all.VC10.fpc)
summary(all.fit10.fpc)
overall.fit10.fpc


# --------------------------------------
# Estimates pooled across beach and stratified by conditions
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness

# all beaches ----------------
# f- coliphage --------
# high risk conditions
data=all[!is.na(all$swim.fmc),]
data.high=subset(data,data$risk=="High")
all.fit10.fmc.high <- glm(gici10~swim.fmc+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC10.fmc.high <- cl(data.high,fm=all.fit10.fmc.high, cluster=data.high$hhid)
overall.fit10.fmc.high <- coeftest(all.fit10.fmc.high, all.VC10.fmc.high)
summary(all.fit10.fmc.high)
overall.fit10.fmc.high

# low risk conditions
data.low=subset(data,data$risk=="Low")
all.fit10.fmc.low <- glm(gici10~swim.fmc+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC10.fmc.low <- cl(data.low,fm=all.fit10.fmc.low, cluster=data.low$hhid)
overall.fit10.fmc.low <- coeftest(all.fit10.fmc.low, all.VC10.fmc.low)
summary(all.fit10.fmc.low)
overall.fit10.fmc.low

# f+ coliphage  --------
# high risk conditions
data=all[!is.na(all$swim.fpc),]
data.high=subset(data,data$risk=="High")
all.fit10.fpc.high <- glm(gici10~swim.fpc+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC10.fpc.high <- cl(data.high,fm=all.fit10.fpc.high, cluster=data.high$hhid)
overall.fit10.fpc.high <- coeftest(all.fit10.fpc.high, all.VC10.fpc.high)
summary(all.fit10.fpc.high)
overall.fit10.fpc.high

# low risk conditions
data.low=subset(data,data$risk=="Low")
all.fit10.fpc.low <- glm(gici10~swim.fpc+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC10.fpc.low <- cl(data.low,fm=all.fit10.fpc.low, cluster=data.low$hhid)
overall.fit10.fpc.low <- coeftest(all.fit10.fpc.low, all.VC10.fpc.low)
summary(all.fit10.fpc.low)
overall.fit10.fpc.low



# --------------------------------------
# save the results
# exclude glm objects and data frames
# (they are really large)
# --------------------------------------
save(

  all.n10.fmc,all.n10.fpc,
  all.n10.fmc.high,all.n10.fmc.low,all.n10.fpc.high,all.n10.fpc.low,

  all.VC10.fmc, all.VC10.fpc,overall.fit10.fmc,overall.fit10.fpc,

  all.VC10.fmc.high,all.VC10.fpc.high,
  overall.fit10.fmc.high,overall.fit10.fpc.high,
  
  all.VC10.fmc.low,all.VC10.fpc.low,
  overall.fit10.fmc.low,overall.fit10.fpc.low,

  file="~/Documents/CRG/coliphage/results/rawoutput/regress-10day-body-pool-nonswimmer.Rdata"
)

