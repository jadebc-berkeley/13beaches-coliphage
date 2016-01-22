##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file conducts maximum likelihood regression
# to estimate prevalence ratios

# Results pooled across beaches and assay

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
all=subset(all,all$swallwater=="Yes")


# --------------------------------------
# Calculate the actual Ns for each cell
# and store them for plotting and tables
# --------------------------------------
regN <- function(outcome,exposurecat) {
  sum(table(outcome,exposurecat))
}

# n's pooled across assay and beach ---------------------------------------
all.n10.fmc.swall = regN(all$gici10[!is.na(all$fmc.pres)],
                       all$fmc.pres[!is.na(all$fmc.pres)])
all.n10.fpc.swall = regN(all$gici10[!is.na(all$fmc.pres)],
                       all$fmc.pres[!is.na(all$fmc.pres)])

# pooled n's by risk level---------------------------------------
data=all[!is.na(all$fmc.pres),]
data.high=subset(data,data$risk=="High")
all.n10.fmc.high.swall = regN(data.high$gici10,data.high$fmc.pres)
data.low=subset(data,data$risk=="Low")
all.n10.fmc.low.swall = regN(data.low$gici10,data.low$fmc.pres)

data=all[!is.na(all$fpc.pres),]
data.high=subset(data,data$risk=="High")
all.n10.fpc.high.swall = regN(data.high$gici10,data.high$fpc.pres)
data.low=subset(data,data$risk=="Low")
all.n10.fpc.low.swall = regN(data.low$gici10,data.low$fpc.pres)



# --------------------------------------
# Estimates pooled across beach and assay method
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness

# f- coliphage ----------------
all.fit10.fmc.swall <- glm(gici10~fmc.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
                           rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc.pres),])

all.VC10.fmc.swall <- cl(all[!is.na(all$fmc.pres)],fm=all.fit10.fmc.swall,
                       cluster=all$hhid[!is.na(all$fmc.pres)])
overall.fit10.fmc.swall <- coeftest(all.fit10.fmc.swall, all.VC10.fmc.swall)
summary(all.fit10.fmc.swall)
overall.fit10.fmc.swall
aic.fmc.swall=AIC(all.fit10.fmc.swall)

# f+ coliphage ----------------
all.fit10.fpc.swall <- glm(gici10~fpc.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
                       rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc.pres),])

all.VC10.fpc.swall <- cl(all[!is.na(all$fpc.pres)],fm=all.fit10.fpc.swall,
                   cluster=all$hhid[!is.na(all$fpc.pres)])
overall.fit10.fpc.swall <- coeftest(all.fit10.fpc.swall, all.VC10.fpc.swall)
summary(all.fit10.fpc.swall)
overall.fit10.fpc.swall
aic.fpc.swall=AIC(all.fit10.fpc.swall)


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
all.fit10.fmc.high.swall <- glm(gici10~fmc.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC10.fmc.high.swall <- cl(data.high,fm=all.fit10.fmc.high.swall, cluster=data.high$hhid)
overall.fit10.fmc.high.swall <- coeftest(all.fit10.fmc.high.swall, all.VC10.fmc.high.swall)
summary(all.fit10.fmc.high.swall)
overall.fit10.fmc.high.swall
aic.fmc.high.swall=AIC(all.fit10.fmc.high.swall)

# low risk conditions
data.low=subset(data,data$risk=="Low")
all.fit10.fmc.low.swall <- glm(gici10~fmc.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC10.fmc.low.swall <- cl(data.low,fm=all.fit10.fmc.low.swall, cluster=data.low$hhid)
overall.fit10.fmc.low.swall <- coeftest(all.fit10.fmc.low.swall, all.VC10.fmc.low.swall)
summary(all.fit10.fmc.low.swall)
overall.fit10.fmc.low.swall
aic.fmc.low.swall=AIC(all.fit10.fmc.low.swall)

# f+ coliphage  --------
# high risk conditions
data=all[!is.na(all$fpc.pres),]
data.high=subset(data,data$risk=="High")
all.fit10.fpc.high.swall <- glm(gici10~fpc.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC10.fpc.high.swall <- cl(data.high,fm=all.fit10.fpc.high.swall, cluster=data.high$hhid)
overall.fit10.fpc.high.swall <- coeftest(all.fit10.fpc.high.swall, all.VC10.fpc.high.swall)
summary(all.fit10.fpc.high.swall)
overall.fit10.fpc.high.swall
aic.fpc.high.swall=AIC(all.fit10.fpc.high.swall)

# low risk conditions
data.low=subset(data,data$risk=="Low")
all.fit10.fpc.low.swall <- glm(gici10~fpc.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC10.fpc.low.swall <- cl(data.low,fm=all.fit10.fpc.low.swall, cluster=data.low$hhid)
overall.fit10.fpc.low.swall <- coeftest(all.fit10.fpc.low.swall, all.VC10.fpc.low.swall)
summary(all.fit10.fpc.low.swall)
overall.fit10.fpc.low.swall
aic.fpc.low.swall=AIC(all.fit10.fpc.low.swall)




# --------------------------------------
# save the results
# exclude glm objects and data frames
# (they are really large)
# --------------------------------------
save(

  all.n10.fmc.swall,all.n10.fpc.swall,
  all.n10.fmc.high.swall,all.n10.fmc.low.swall,all.n10.fpc.high.swall,all.n10.fpc.low.swall,

  all.VC10.fmc.swall, all.VC10.fpc.swall,overall.fit10.fmc.swall,overall.fit10.fpc.swall,

  all.VC10.fmc.high.swall,all.VC10.fpc.high.swall,
  overall.fit10.fmc.high.swall,overall.fit10.fpc.high.swall,
  
  all.VC10.fmc.low.swall,all.VC10.fpc.low.swall,
  overall.fit10.fmc.low.swall,overall.fit10.fpc.low.swall,
  
  aic.fmc.swall, aic.fpc.swall,
  aic.fmc.high.swall, aic.fpc.high.swall,aic.fmc.low.swall, aic.fpc.low.swall,
  
  file="~/Documents/CRG/coliphage/results/rawoutput/regress-10day-swall-pool.Rdata"
)

