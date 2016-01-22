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
all=subset(all,all$headunder=="Yes")


# --------------------------------------
# Calculate the actual Ns for each cell
# and store them for plotting and tables
# --------------------------------------
regN <- function(outcome,exposurecat) {
  sum(table(outcome,exposurecat))
}

# n's pooled across assay and beach ---------------------------------------
all.n10.fmc.head = regN(all$gici10[!is.na(all$fmc.pres)],
                       all$fmc.pres[!is.na(all$fmc.pres)])
all.n10.fpc.head = regN(all$gici10[!is.na(all$fmc.pres)],
                       all$fmc.pres[!is.na(all$fmc.pres)])

# pooled n's by risk level---------------------------------------
data=all[!is.na(all$fmc.pres),]
data.high=subset(data,data$risk=="High")
all.n10.fmc.high.head = regN(data.high$gici10,data.high$fmc.pres)
data.low=subset(data,data$risk=="Low")
all.n10.fmc.low.head = regN(data.low$gici10,data.low$fmc.pres)

data=all[!is.na(all$fpc.pres),]
data.high=subset(data,data$risk=="High")
all.n10.fpc.high.head = regN(data.high$gici10,data.high$fpc.pres)
data.low=subset(data,data$risk=="Low")
all.n10.fpc.low.head = regN(data.low$gici10,data.low$fpc.pres)


# --------------------------------------
# Estimates pooled across beach and assay method
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness

# f- coliphage ----------------
all.fit10.fmc.head <- glm(gici10~fmc.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
                           rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc.pres),])

all.VC10.fmc.head <- cl(all[!is.na(all$fmc.pres)],fm=all.fit10.fmc.head,
                       cluster=all$hhid[!is.na(all$fmc.pres)])
overall.fit10.fmc.head <- coeftest(all.fit10.fmc.head, all.VC10.fmc.head)
summary(all.fit10.fmc.head)
overall.fit10.fmc.head
aic.fmc.head=AIC(all.fit10.fmc.head)

# f+ coliphage ----------------
all.fit10.fpc.head <- glm(gici10~fpc.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
                       rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc.pres),])

all.VC10.fpc.head <- cl(all[!is.na(all$fpc.pres)],fm=all.fit10.fpc.head,
                   cluster=all$hhid[!is.na(all$fpc.pres)])
overall.fit10.fpc.head <- coeftest(all.fit10.fpc.head, all.VC10.fpc.head)
summary(all.fit10.fpc.head)
overall.fit10.fpc.head
aic.fpc.head=AIC(all.fit10.fpc.head)


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
all.fit10.fmc.high.head <- glm(gici10~fmc.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC10.fmc.high.head <- cl(data.high,fm=all.fit10.fmc.high.head, cluster=data.high$hhid)
overall.fit10.fmc.high.head <- coeftest(all.fit10.fmc.high.head, all.VC10.fmc.high.head)
summary(all.fit10.fmc.high.head)
overall.fit10.fmc.high.head
aic.fmc.high.head=AIC(all.fit10.fmc.high.head)

# low risk conditions
data.low=subset(data,data$risk=="Low")
all.fit10.fmc.low.head <- glm(gici10~fmc.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC10.fmc.low.head <- cl(data.low,fm=all.fit10.fmc.low.head, cluster=data.low$hhid)
overall.fit10.fmc.low.head <- coeftest(all.fit10.fmc.low.head, all.VC10.fmc.low.head)
summary(all.fit10.fmc.low.head)
overall.fit10.fmc.low.head
aic.fmc.low.head=AIC(all.fit10.fmc.low.head)

# f+ coliphage  --------
# high risk conditions
data=all[!is.na(all$fpc.pres),]
data.high=subset(data,data$risk=="High")
all.fit10.fpc.high.head <- glm(gici10~fpc.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC10.fpc.high.head <- cl(data.high,fm=all.fit10.fpc.high.head, cluster=data.high$hhid)
overall.fit10.fpc.high.head <- coeftest(all.fit10.fpc.high.head, all.VC10.fpc.high.head)
summary(all.fit10.fpc.high.head)
overall.fit10.fpc.high.head
aic.fpc.high.head=AIC(all.fit10.fpc.high.head)

# low risk conditions
data.low=subset(data,data$risk=="Low")
all.fit10.fpc.low.head <- glm(gici10~fpc.pres+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC10.fpc.low.head <- cl(data.low,fm=all.fit10.fpc.low.head, cluster=data.low$hhid)
overall.fit10.fpc.low.head <- coeftest(all.fit10.fpc.low.head, all.VC10.fpc.low.head)
summary(all.fit10.fpc.low.head)
overall.fit10.fpc.low.head
aic.fpc.low.head=AIC(all.fit10.fpc.low.head)




# --------------------------------------
# save the results
# exclude glm objects and data frames
# (they are really large)
# --------------------------------------
save(

  all.n10.fmc.head, all.n10.fpc.head,
  all.n10.fmc.high.head, all.n10.fmc.low.head, all.n10.fpc.high.head, all.n10.fpc.low.head,

  all.VC10.fmc.head, all.VC10.fpc.head, overall.fit10.fmc.head, overall.fit10.fpc.head,

  all.VC10.fmc.high.head, all.VC10.fpc.high.head,
  overall.fit10.fmc.high.head, overall.fit10.fpc.high.head,
  
  all.VC10.fmc.low.head, all.VC10.fpc.low.head,
  overall.fit10.fmc.low.head, overall.fit10.fpc.low.head,
  
  aic.fmc.head, aic.fpc.head,
  aic.fmc.high.head, aic.fpc.high.head, aic.fmc.low.head, aic.fpc.low.head,
  
  file="~/Documents/CRG/coliphage/results/rawoutput/regress-10day-head-pool.Rdata"
)

