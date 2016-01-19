##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file conducts maximum likelihood regression
# to estimate prevalence ratios

# Results pooled across beaches, assays, and type
# of coliphage

# 10 day gi illness
##########################################

rm(list=ls())

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
all=subset(all,nowq==0)
all=subset(all,all$bodycontact=="Yes")

# subset to observations with no missing enterococcus information
all=subset(all,!is.na(all$entero35))

# create indicator for pooled presence absence

all$pres=NA
all$pres[all$fmc.pres==1 | all$fpc.pres==1]=1
all$pres[all$fmc.pres==0 & all$fpc.pres==0]=0

# --------------------------------------
# Creating joint indicator variable for
# regressions
# --------------------------------------
all$ent=NA
all$ent[all$pres==0 & all$entero35==0]=1
all$ent[all$pres==1 & all$entero35==0]=2
all$ent[all$pres==0 & all$entero35==1]=3
all$ent[all$pres==1 & all$entero35==1]=4
all$ent=as.factor(all$ent)


# --------------------------------------
# All conditions

# Estimates pooled across beach
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness

# all beaches ----------------

# -------------------------------------

all.fit10 <- glm(gici10~ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$ent),])

all.VC10 <- cl(all[!is.na(all$ent)],fm=all.fit10,
    cluster=all$hhid[!is.na(all$ent)])
overall.fit10.int <- coeftest(all.fit10, all.VC10)
summary(all.fit10)
overall.fit10.int
aic.int=AIC(all.fit10)


# --------------------------------------
# Estimates pooled across beach and stratified by conditions

# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness

# all beaches ----------------

# -------------------------------------
# f- coliphage
# -------------------------------------
data=all[!is.na(all$ent),]

# high risk conditions ----------------
data.high=subset(data,data$risk=="High")

all.fit10.high <- glm(gici10~ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+
          rawfood+beach,family=poisson(link="log"),data=data.high)
all.VC10.high <- cl(data.high,fm=all.fit10.high, cluster=data.high$hhid)
overall.fit10.high.int <- coeftest(all.fit10.high, all.VC10.high)
summary(all.fit10.high)
overall.fit10.high.int
aic.high.int=AIC(all.fit10.high)


# low risk conditions ----------------
data.low=subset(data,data$risk=="Low")

all.fit10.low <- glm(gici10~ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+
          rawfood+beach,family=poisson(link="log"),data=data.low)
all.VC10.low <- cl(data.low,fm=all.fit10.low, cluster=data.low$hhid)
overall.fit10.low.int <- coeftest(all.fit10.low, all.VC10.low)
summary(all.fit10.low)
overall.fit10.low.int
aic.low.int=AIC(all.fit10.low)




# --------------------------------------
# save the results
# exclude glm objects and data frames
# (they are really large)
# --------------------------------------
save(
  
  overall.fit10.int,

  overall.fit10.low.int,overall.fit10.high.int,

  aic.int, 
  aic.low.int, aic.high.int,

  file="~/Documents/CRG/coliphage/results/rawoutput/regress-10day-body-joint-pool-both.Rdata"
)


