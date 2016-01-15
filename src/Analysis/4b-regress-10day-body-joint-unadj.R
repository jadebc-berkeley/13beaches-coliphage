##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file conducts maximum likelihood regression
# to estimate prevalence ratios

# Results pooled across beaches

# unadjusted analysis

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
source("~/Documents/CRG/coliphage/13beaches-data/src/Analysis/0-base-functions.R")

data=preprocess.6beaches(beaches13)

# restrict to 6 beaches with coliphage data
beach.list=c("Avalon","Doheny","Malibu","Mission Bay",
             "Fairhope","Goddard")

all=data[data$beach %in% beach.list,]
all=subset(all,nowq==0)
all=subset(all,all$bodycontact=="Yes")

# subset to observations with no missing enterococcus information
all=subset(all,!is.na(all$entero35))

# --------------------------------------
# Creating joint indicator variable for
# regressions
# --------------------------------------
all$fmc1601.ent=NA
all$fmc1601.ent[all$fmc1601.pres==0 & all$entero35==0]=1
all$fmc1601.ent[all$fmc1601.pres==1 & all$entero35==0]=2
all$fmc1601.ent[all$fmc1601.pres==0 & all$entero35==1]=3
all$fmc1601.ent[all$fmc1601.pres==1 & all$entero35==1]=4
all$fmc1601.ent=as.factor(all$fmc1601.ent)

all$fmc1602.ent=NA
all$fmc1602.ent[all$fmc1602.pres==0 & all$entero35==0]=1
all$fmc1602.ent[all$fmc1602.pres==1 & all$entero35==0]=2
all$fmc1602.ent[all$fmc1602.pres==0 & all$entero35==1]=3
all$fmc1602.ent[all$fmc1602.pres==1 & all$entero35==1]=4
all$fmc1602.ent=as.factor(all$fmc1602.ent)

all$fpc1601.ent=NA
all$fpc1601.ent[all$fpc1601.pres==0 & all$entero35==0]=1
all$fpc1601.ent[all$fpc1601.pres==1 & all$entero35==0]=2
all$fpc1601.ent[all$fpc1601.pres==0 & all$entero35==1]=3
all$fpc1601.ent[all$fpc1601.pres==1 & all$entero35==1]=4
all$fpc1601.ent=as.factor(all$fpc1601.ent)

all$fpc1602.ent=NA
all$fpc1602.ent[all$fpc1602.pres==0 & all$entero35==0]=1
all$fpc1602.ent[all$fpc1602.pres==1 & all$entero35==0]=2
all$fpc1602.ent[all$fpc1602.pres==0 & all$entero35==1]=3
all$fpc1602.ent[all$fpc1602.pres==1 & all$entero35==1]=4
all$fpc1602.ent=as.factor(all$fpc1602.ent)

# --------------------------------------
# Calculate the actual Ns for each cell
# and store them for plotting and tables
# --------------------------------------
regN <- function(outcome,exposurecat) {
  sum(table(outcome,exposurecat))
}

# all conditions
all.n10.fmc1601 = data.frame(table(all$fmc1601.ent))[,2]
all.n10.fmc1602 = data.frame(table(all$fmc1602.ent))[,2]
all.n10.fpc1601 = data.frame(table(all$fpc1601.ent))[,2]
all.n10.fpc1602 = data.frame(table(all$fpc1602.ent))[,2]

# stratified by conditions
data=all[!is.na(all$fmc1602.pres),]
data.high=subset(data,data$risk=="High")
all.n10.fmc1602.high = data.frame(table(data.high$fmc1602.ent))[,2]
data.low=subset(data,data$risk=="Low")
all.n10.fmc1602.low = data.frame(table(data.low$fmc1602.ent))[,2]

data=all[!is.na(all$fpc1601.pres),]
data.high=subset(data,data$risk=="High")
all.n10.fpc1601.high = data.frame(table(data.high$fpc1601.ent))[,2]
data.low=subset(data,data$risk=="Low")
all.n10.fpc1601.low = data.frame(table(data.low$fpc1601.ent))[,2]

data=all[!is.na(all$fpc1602.pres),]
data.high=subset(data,data$risk=="High")
all.n10.fpc1602.high = data.frame(table(data.high$fpc1602.ent))[,2]
data.low=subset(data,data$risk=="Low")
all.n10.fpc1602.low = data.frame(table(data.low$fpc1602.ent))[,2]

# --------------------------------------
# All conditions

# Estimates pooled across beach
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness

# all beaches ----------------

# -------------------------------------
# fmc 1601
# -------------------------------------
all.fit10.fmc1601 <- glm(gici10~fmc1601.ent,family=poisson(link="log"),data=all[!is.na(all$fmc1601.pres),])

all.VC10.fmc1601 <- cl(all[!is.na(all$fmc1601.pres)],fm=all.fit10.fmc1601,
                       cluster=all$hhid[!is.na(all$fmc1601.pres)])
overall.fit10.fmc1601.int <- coeftest(all.fit10.fmc1601, all.VC10.fmc1601)
summary(all.fit10.fmc1601)
overall.fit10.fmc1601.int
aic.fmc1601.int=AIC(all.fit10.fmc1601)


# -------------------------------------
# fmc 1602
# -------------------------------------
all.fit10.fmc1602 <- glm(gici10~fmc1602.ent,family=poisson(link="log"),data=all[!is.na(all$fmc1602.pres),])

all.VC10.fmc1602 <- cl(all[!is.na(all$fmc1602.pres),],fm=all.fit10.fmc1602,
                       cluster=all$hhid[!is.na(all$fmc1602.pres)])
overall.fit10.fmc1602.int <- coeftest(all.fit10.fmc1602, all.VC10.fmc1602)
summary(all.fit10.fmc1602)
overall.fit10.fmc1602.int
aic.fmc1602.int=AIC(all.fit10.fmc1602)


# -------------------------------------
# fpc 1601
# -------------------------------------
all.fit10.fpc1601 <- glm(gici10~fpc1601.ent,family=poisson(link="log"),data=all[!is.na(all$fpc1601.pres),])

all.VC10.fpc1601 <- cl(all[!is.na(all$fpc1601.pres),],fm=all.fit10.fpc1601,
                       cluster=all$hhid[!is.na(all$fpc1601.pres)])
overall.fit10.fpc1601.int <- coeftest(all.fit10.fpc1601, all.VC10.fpc1601)
summary(all.fit10.fpc1601)
overall.fit10.fpc1601.int
aic.fpc1601.int=AIC(all.fit10.fpc1601)

# -------------------------------------
# fpc 1602
# -------------------------------------
all.fit10.fpc1602 <- glm(gici10~fpc1602.ent,family=poisson(link="log"),data=all[!is.na(all$fpc1602.pres),])

all.VC10.fpc1602 <- cl(all[!is.na(all$fpc1602.pres),],fm=all.fit10.fpc1602,
                       cluster=all$hhid[!is.na(all$fpc1602.pres)])
overall.fit10.fpc1602.int <- coeftest(all.fit10.fpc1602, all.VC10.fpc1602)
summary(all.fit10.fpc1602)
overall.fit10.fpc1602.int
aic.fpc1602.int=AIC(all.fit10.fpc1602)

# -------------------------------------
# fmc 1601
# -------------------------------------

# --------------------------------------
# Estimates pooled across beach and stratified by conditions

# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness

# all beaches ----------------

# -------------------------------------
# fmc 1602
# -------------------------------------
data=all[!is.na(all$fmc1602.pres),]

# high risk conditions ----------------
data.high=subset(data,data$risk=="High")

all.fit10.fmc1602.high <- glm(gici10~fmc1602.ent,family=poisson(link="log"),data=data.high)
all.VC10.fmc1602.high <- cl(data.high,fm=all.fit10.fmc1602.high, cluster=data.high$hhid)
overall.fit10.fmc1602.high.int <- coeftest(all.fit10.fmc1602.high, all.VC10.fmc1602.high)
summary(all.fit10.fmc1602.high)
overall.fit10.fmc1602.high.int
aic.fmc1602.high.int=AIC(all.fit10.fmc1602.high)


# low risk conditions ----------------
data.low=subset(data,data$risk=="Low")

all.fit10.fmc1602.low <- glm(gici10~fmc1602.ent,family=poisson(link="log"),data=data.low)
all.VC10.fmc1602.low <- cl(data.low,fm=all.fit10.fmc1602.low, cluster=data.low$hhid)
overall.fit10.fmc1602.low.int <- coeftest(all.fit10.fmc1602.low, all.VC10.fmc1602.low)
summary(all.fit10.fmc1602.low)
overall.fit10.fmc1602.low.int
aic.fmc1602.low.int=AIC(all.fit10.fmc1602.low)


# -------------------------------------
# fpc 1601
# -------------------------------------
data=all[!is.na(all$fpc1601.pres),]

# high risk conditions ----------------
data.high=subset(data,data$risk=="High")

all.fit10.fpc1601.high <- glm(gici10~fpc1601.ent,family=poisson(link="log"),data=data.high)
all.VC10.fpc1601.high <- cl(data.high,fm=all.fit10.fpc1601.high, cluster=data.high$hhid)
overall.fit10.fpc1601.high.int <- coeftest(all.fit10.fpc1601.high, all.VC10.fpc1601.high)
summary(all.fit10.fpc1601.high)
overall.fit10.fpc1601.high.int
aic.fpc1601.high.int=AIC(all.fit10.fpc1601.high)


# low risk conditions ----------------
data.low=subset(data,data$risk=="Low")

all.fit10.fpc1601.low <- glm(gici10~fpc1601.ent,family=poisson(link="log"),data=data.low)
all.VC10.fpc1601.low <- cl(data.low,fm=all.fit10.fpc1601.low, 
  cluster=data.low$hhid)
overall.fit10.fpc1601.low.int <- coeftest(all.fit10.fpc1601.low, all.VC10.fpc1601.low)
summary(all.fit10.fpc1601.low)
overall.fit10.fpc1601.low.int
aic.fpc1601.low.int=AIC(all.fit10.fpc1601.low)


# -------------------------------------
# fpc 1602
# -------------------------------------
data=all[!is.na(all$fpc1602.pres),]

# high risk conditions ----------------
data.high=subset(data,data$risk=="High")

all.fit10.fpc1602.high <- glm(gici10~fpc1602.ent,family=poisson(link="log"),data=data.high)
all.VC10.fpc1602.high <- cl(data.high,fm=all.fit10.fpc1602.high, cluster=data.high$hhid)
overall.fit10.fpc1602.high.int <- coeftest(all.fit10.fpc1602.high, all.VC10.fpc1602.high)
summary(all.fit10.fpc1602.high)
overall.fit10.fpc1602.high.int
aic.fpc1602.high.int=AIC(all.fit10.fpc1602.high)


# low risk conditions ----------------
data.low=subset(data,data$risk=="Low")

all.fit10.fpc1602.low <- glm(gici10~fpc1602.ent,family=poisson(link="log"),data=data.low)
all.VC10.fpc1602.low <- cl(data.low,fm=all.fit10.fpc1602.low, cluster=data.low$hhid)
overall.fit10.fpc1602.low.int <- coeftest(all.fit10.fpc1602.low, all.VC10.fpc1602.low)
summary(all.fit10.fpc1602.low)
overall.fit10.fpc1602.low.int
aic.fpc1602.low.int=AIC(all.fit10.fpc1602.low)


# --------------------------------------
# save the results
# exclude glm objects and data frames
# (they are really large)
# --------------------------------------
save(

  all.n10.fmc1601,all.n10.fmc1602,all.n10.fpc1601,all.n10.fpc1602,

  all.n10.fmc1602.low, all.n10.fmc1602.high,  all.n10.fpc1601.low,
  all.n10.fpc1601.high, all.n10.fpc1602.low,  all.n10.fpc1602.high,
  
  overall.fit10.fmc1601.int,overall.fit10.fmc1602.int,
  overall.fit10.fpc1601.int,overall.fit10.fpc1602.int,

  overall.fit10.fmc1602.low.int,overall.fit10.fmc1602.high.int,
  overall.fit10.fpc1601.low.int,overall.fit10.fpc1601.high.int,
  overall.fit10.fpc1602.low.int,overall.fit10.fpc1602.high.int,

  aic.fmc1601.int, aic.fmc1602.int, aic.fpc1601.int,aic.fpc1602.int,

  aic.fmc1602.low.int, aic.fmc1602.high.int,
  aic.fpc1601.low.int, aic.fpc1601.high.int,
  aic.fpc1602.low.int, aic.fpc1602.high.int, 
  aic.fpc1602.low.int, aic.fpc1602.high.int,

  file="~/Documents/CRG/coliphage/results/rawoutput/regress-10day-body-joint-unadj.Rdata"
)


