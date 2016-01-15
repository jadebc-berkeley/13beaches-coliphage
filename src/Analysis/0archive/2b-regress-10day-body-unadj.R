##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file conducts maximum likelihood regression
# to estimate prevalence ratios

# Results pooled across beaches

# unadjusted analyses

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

# pooled n's ---------------------------------------
all.n10.fmc1601 = regN(all$gici10[!is.na(all$fmc1601.pres)],
                       all$fmc1601.pres[!is.na(all$fmc1601.pres)])
all.n10.fmc1602 = regN(all$gici10[!is.na(all$fmc1602.pres)],
                       all$fmc1602.pres[!is.na(all$fmc1602.pres)])
all.n10.fpc1601 = regN(all$gici10[!is.na(all$fpc1601.pres)],
                       all$fpc1601.pres[!is.na(all$fpc1601.pres)])
all.n10.fpc1602 = regN(all$gici10[!is.na(all$fpc1602.pres)],
                       all$fpc1602.pres[!is.na(all$fpc1602.pres)])

# pooled n's by risk level---------------------------------------
data=all[!is.na(all$fmc1602.pres),]
data.high=subset(data,data$groundwater=="Above median flow" | data$berm=="Open")
all.n10.fmc1602.high = regN(data.high$gici10,data.high$fmc1602.pres)
data.low=subset(data,data$groundwater=="Below median flow" | data$berm=="Closed")
all.n10.fmc1602.low = regN(data.low$gici10,data.low$fmc1602.pres)

data=all[!is.na(all$fpc1601.pres),]
data.high=subset(data,data$groundwater=="Above median flow" | data$berm=="Open")
all.n10.fpc1601.high = regN(data.high$gici10,data.high$fpc1601.pres)
data.low=subset(data,data$groundwater=="Below median flow" | data$berm=="Closed")
all.n10.fpc1601.low = regN(data.low$gici10,data.low$fpc1601.pres)

data=all[!is.na(all$fpc1602.pres),]
data.high=subset(data,data$groundwater=="Above median flow" | data$berm=="Open")
all.n10.fpc1602.high = regN(data.high$gici10,data.high$fpc1602.pres)
data.low=subset(data,data$groundwater=="Below median flow" | data$berm=="Closed")
all.n10.fpc1602.low = regN(data.low$gici10,data.low$fpc1602.pres)


# --------------------------------------
# Estimates pooled across beach
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness

# all beaches ----------------
# fmc 1601
all.fit10.fmc1601 <- glm(gici10~fmc1601.pres,family=poisson(link="log"),data=all[!is.na(all$fmc1601.pres),])

all.VC10.fmc1601 <- cl(all[!is.na(all$fmc1601.pres)],fm=all.fit10.fmc1601,
  cluster=all$hhid[!is.na(all$fmc1601.pres)])
overall.fit10.fmc1601 <- coeftest(all.fit10.fmc1601, all.VC10.fmc1601)
summary(all.fit10.fmc1601)
overall.fit10.fmc1601
aic.fmc1601=AIC(all.fit10.fmc1601)

# fmc 1602
all.fit10.fmc1602 <- glm(gici10~fmc1602.pres,family=poisson(link="log"),data=all[!is.na(all$fmc1602.pres),])

all.VC10.fmc1602 <- cl(all[!is.na(all$fmc1602.pres),],fm=all.fit10.fmc1602,
  cluster=all$hhid[!is.na(all$fmc1602.pres)])
overall.fit10.fmc1602 <- coeftest(all.fit10.fmc1602, all.VC10.fmc1602)
summary(all.fit10.fmc1602)
overall.fit10.fmc1602
aic.fmc1602=AIC(all.fit10.fmc1602)

# fpc 1601
all.fit10.fpc1601 <- glm(gici10~fpc1601.pres,family=poisson(link="log"),data=all[!is.na(all$fpc1601.pres),])

all.VC10.fpc1601 <- cl(all[!is.na(all$fpc1601.pres),],fm=all.fit10.fpc1601,
    cluster=all$hhid[!is.na(all$fpc1601.pres)])
overall.fit10.fpc1601 <- coeftest(all.fit10.fpc1601, all.VC10.fpc1601)
summary(all.fit10.fpc1601)
overall.fit10.fpc1601
aic.fpc1601=AIC(all.fit10.fpc1601)

# fpc 1602
all.fit10.fpc1602 <- glm(gici10~fpc1602.pres,family=poisson(link="log"),data=all[!is.na(all$fpc1602.pres),])

all.VC10.fpc1602 <- cl(all[!is.na(all$fpc1602.pres),],fm=all.fit10.fpc1602,
    cluster=all$hhid[!is.na(all$fpc1602.pres)])
overall.fit10.fpc1602 <- coeftest(all.fit10.fpc1602, all.VC10.fpc1602)
summary(all.fit10.fpc1602)
overall.fit10.fpc1602
aic.fpc1602=AIC(all.fit10.fpc1602)


# --------------------------------------
# Estimates pooled across beach and stratified by conditions
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness

# all beaches ----------------
# fmc 1602 --------
# high risk conditions
data=all[!is.na(all$fmc1602.pres),]
data.high=subset(data,data$groundwater=="Above median flow" | data$berm=="Open")
all.fit10.fmc1602.high <- glm(gici10~fmc1602.pres,family=poisson(link="log"),data=data.high)

all.VC10.fmc1602.high <- cl(data.high,fm=all.fit10.fmc1602.high, cluster=data.high$hhid)
overall.fit10.fmc1602.high <- coeftest(all.fit10.fmc1602.high, all.VC10.fmc1602.high)
summary(all.fit10.fmc1602.high)
overall.fit10.fmc1602.high
aic.fmc1602.high=AIC(all.fit10.fmc1602.high)

# low risk conditions
data.low=subset(data,data$groundwater=="Below median flow" | data$berm=="Closed")
all.fit10.fmc1602.low <- glm(gici10~fmc1602.pres,family=poisson(link="log"),data=data.low)

all.VC10.fmc1602.low <- cl(data.low,fm=all.fit10.fmc1602.low, cluster=data.low$hhid)
overall.fit10.fmc1602.low <- coeftest(all.fit10.fmc1602.low, all.VC10.fmc1602.low)
summary(all.fit10.fmc1602.low)
overall.fit10.fmc1602.low
aic.fmc1602.low=AIC(all.fit10.fmc1602.low)

# fpc 1601 --------
# high risk conditions
data=all[!is.na(all$fpc1601.pres),]
data.high=subset(data,data$groundwater=="Above median flow" | data$berm=="Open")
all.fit10.fpc1601.high <- glm(gici10~fpc1601.pres,family=poisson(link="log"),data=data.high)

all.VC10.fpc1601.high <- cl(data.high,fm=all.fit10.fpc1601.high, cluster=data.high$hhid)
overall.fit10.fpc1601.high <- coeftest(all.fit10.fpc1601.high, all.VC10.fpc1601.high)
summary(all.fit10.fpc1601.high)
overall.fit10.fpc1601.high
aic.fpc1601.high=AIC(all.fit10.fpc1601.high)

# low risk conditions
data.low=subset(data,data$groundwater=="Below median flow" | data$berm=="Closed")
all.fit10.fpc1601.low <- glm(gici10~fpc1601.pres,family=poisson(link="log"),data=data.low)

all.VC10.fpc1601.low <- cl(data.low,fm=all.fit10.fpc1601.low, cluster=data.low$hhid)
overall.fit10.fpc1601.low <- coeftest(all.fit10.fpc1601.low, all.VC10.fpc1601.low)
summary(all.fit10.fpc1601.low)
overall.fit10.fpc1601.low
aic.fpc1601.low=AIC(all.fit10.fpc1601.low)


# fpc 1602 --------
# high risk conditions
data=all[!is.na(all$fpc1602.pres),]
data.high=subset(data,data$groundwater=="Above median flow" | data$berm=="Open")
all.fit10.fpc1602.high <- glm(gici10~fpc1602.pres,family=poisson(link="log"),data=data.high)

all.VC10.fpc1602.high <- cl(data.high,fm=all.fit10.fpc1602.high, cluster=data.high$hhid)
overall.fit10.fpc1602.high <- coeftest(all.fit10.fpc1602.high, all.VC10.fpc1602.high)
summary(all.fit10.fpc1602.high)
overall.fit10.fpc1602.high
aic.fpc1602.high=AIC(all.fit10.fpc1602.high)


# low risk conditions
data.low=subset(data,data$groundwater=="Below median flow" | data$berm=="Closed")
all.fit10.fpc1602.low <- glm(gici10~fpc1602.pres,family=poisson(link="log"),data=data.low)

all.VC10.fpc1602.low <- cl(data.low,fm=all.fit10.fpc1602.low, cluster=data.low$hhid)
overall.fit10.fpc1602.low <- coeftest(all.fit10.fpc1602.low, all.VC10.fpc1602.low)
summary(all.fit10.fpc1602.low)
overall.fit10.fpc1602.low
aic.fpc1602.low=AIC(all.fit10.fpc1602.low)



# --------------------------------------
# save the results
# exclude glm objects and data frames
# (they are really large)
# --------------------------------------
save(

  all.n10.fmc1601,all.n10.fmc1602,all.n10.fpc1601,all.n10.fpc1602,
  
  all.n10.fmc1602.high,all.n10.fmc1602.low,all.n10.fpc1601.high,
  all.n10.fpc1601.low,all.n10.fpc1602.high,all.n10.fpc1602.low,

  all.VC10.fmc1601,all.VC10.fmc1602,all.VC10.fpc1601,all.VC10.fpc1602,
  overall.fit10.fmc1601,overall.fit10.fmc1602,overall.fit10.fpc1601,
  overall.fit10.fpc1602,

  all.VC10.fmc1602.high,all.VC10.fpc1601.high,all.VC10.fpc1602.high,
  overall.fit10.fmc1602.high,overall.fit10.fpc1601.high,
  overall.fit10.fpc1602.high,
  
  all.VC10.fmc1602.low,all.VC10.fpc1601.low,all.VC10.fpc1602.low,
  overall.fit10.fmc1602.low,overall.fit10.fpc1601.low,
  overall.fit10.fpc1602.low,
  
  aic.fmc1601,aic.fmc1602,aic.fpc1601,aic.fpc1602,
  aic.fmc1602.high,aic.fmc1602.low,
  aic.fpc1601.high,aic.fpc1601.low,aic.fpc1602.high,aic.fpc1602.low,

  file="~/dropbox/coliphage/results/rawoutput/regress-10day-body-unadj.Rdata"
)

