##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file conducts maximum likelihood regression
# to estimate prevalence ratios

# Results pooled across beaches

# 10 day gi illness

# negative control analysis
##########################################

rm(list=ls())

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
all=subset(all,nowq==0)
all=subset(all,all$anycontact=="No")

# subset to observations with no missing enterococcus information
all=subset(all,!is.na(all$entero35))

# --------------------------------------
# Creating joint indicator variable for
# regressions
# --------------------------------------
all$fmc1601.ent=NA
all$fmc1601.ent[all$fmc1601.pres==0]=1
all$fmc1601.ent[all$fmc1601.pres==1]=2
all$fmc1601.ent[all$fmc1601.pres==0 & all$entero35==1]=3
all$fmc1601.ent[all$fmc1601.pres==1 & all$entero35==1]=4
all$fmc1601.ent=as.factor(all$fmc1601.ent)

all$fmc1602.ent=NA
all$fmc1602.ent[all$fmc1602.pres==0]=1
all$fmc1602.ent[all$fmc1602.pres==1]=2
all$fmc1602.ent[all$fmc1602.pres==0 & all$entero35==1]=3
all$fmc1602.ent[all$fmc1602.pres==1 & all$entero35==1]=4
all$fmc1602.ent=as.factor(all$fmc1602.ent)

all$fpc1601.ent=NA
all$fpc1601.ent[all$fpc1601.pres==0]=1
all$fpc1601.ent[all$fpc1601.pres==1]=2
all$fpc1601.ent[all$fpc1601.pres==0 & all$entero35==1]=3
all$fpc1601.ent[all$fpc1601.pres==1 & all$entero35==1]=4
all$fpc1601.ent=as.factor(all$fpc1601.ent)

all$fpc1602.ent=NA
all$fpc1602.ent[all$fpc1602.pres==0]=1
all$fpc1602.ent[all$fpc1602.pres==1]=2
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
data=all[!is.na(all$fmc1602.pres) ,]
data.high=subset(data,data$groundwater=="Above median flow" | data$berm=="Open")
all.n10.fmc1602.high = data.frame(table(data.high$fmc1602.ent))[,2]
data.low=subset(data,data$groundwater=="Below median flow" | data$berm=="Closed")
all.n10.fmc1602.low = data.frame(table(data.low$fmc1602.ent))[,2]

data=all[!is.na(all$fpc1601.pres) ,]
data.high=subset(data,data$groundwater=="Above median flow" | data$berm=="Open")
all.n10.fpc1601.high = data.frame(table(data.high$fpc1601.ent))[,2]
data.low=subset(data,data$groundwater=="Below median flow" | data$berm=="Closed")
all.n10.fpc1601.low = data.frame(table(data.low$fpc1601.ent))[,2]

data=all[!is.na(all$fpc1602.pres) ,]
data.high=subset(data,data$groundwater=="Above median flow" | data$berm=="Open")
all.n10.fpc1602.high = data.frame(table(data.high$fpc1602.ent))[,2]
data.low=subset(data,data$groundwater=="Below median flow" | data$berm=="Closed")
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
nc.fit10.fmc1601 <- glm(gici10~fmc1601.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+
                           rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc1601.pres) ,])

nc.VC10.fmc1601 <- cl(all[!is.na(all$fmc1601.pres) ],fm=nc.fit10.fmc1601,
                       cluster=all$hhid[!is.na(all$fmc1601.pres) ])
nc.overall.fit10.fmc1601.int <- coeftest(nc.fit10.fmc1601, nc.VC10.fmc1601)
summary(nc.fit10.fmc1601)
nc.overall.fit10.fmc1601.int
nc.aic.fmc1601.int=AIC(nc.fit10.fmc1601)


# -------------------------------------
# fmc 1602
# -------------------------------------
nc.fit10.fmc1602 <- glm(gici10~fmc1602.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+
        rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fmc1602.pres),])

nc.VC10.fmc1602 <- cl(all[!is.na(all$fmc1602.pres),],fm=nc.fit10.fmc1602,
                       cluster=all$hhid[!is.na(all$fmc1602.pres)])
nc.overall.fit10.fmc1602.int <- coeftest(nc.fit10.fmc1602, nc.VC10.fmc1602)
summary(nc.fit10.fmc1602)
nc.overall.fit10.fmc1602.int
nc.aic.fmc1602.int=AIC(nc.fit10.fmc1602)


# -------------------------------------
# fpc 1601
# -------------------------------------
nc.fit10.fpc1601 <- glm(gici10~fpc1601.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+
                           rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1601.pres) ,])

nc.VC10.fpc1601 <- cl(all[!is.na(all$fpc1601.pres) ,],fm=nc.fit10.fpc1601,
                       cluster=all$hhid[!is.na(all$fpc1601.pres) ])
nc.overall.fit10.fpc1601.int <- coeftest(nc.fit10.fpc1601, nc.VC10.fpc1601)
summary(nc.fit10.fpc1601)
nc.overall.fit10.fpc1601.int
nc.aic.fpc1601.int=AIC(nc.fit10.fpc1601)

# -------------------------------------
# fpc 1602
# -------------------------------------
nc.fit10.fpc1602 <- glm(gici10~fpc1602.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+
                           rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$fpc1602.pres) & 
                           all$beach!="Malibu",])

nc.VC10.fpc1602 <- cl(all[!is.na(all$fpc1602.pres) ,],fm=nc.fit10.fpc1602,
                       cluster=all$hhid[!is.na(all$fpc1602.pres) ])
nc.overall.fit10.fpc1602.int <- coeftest(nc.fit10.fpc1602, nc.VC10.fpc1602)
summary(nc.fit10.fpc1602)
nc.overall.fit10.fpc1602.int
nc.aic.fpc1602.int=AIC(nc.fit10.fpc1602)

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
data=all[!is.na(all$fmc1602.pres) ,]

# high risk conditions ----------------
data.high=subset(data,data$risk=="High")

nc.fit10.fmc1602.high <- glm(gici10~fmc1602.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+
          rawfood+beach,family=poisson(link="log"),data=data.high)
nc.VC10.fmc1602.high <- cl(data.high,fm=nc.fit10.fmc1602.high, cluster=data.high$hhid)
nc.overall.fit10.fmc1602.high.int <- coeftest(nc.fit10.fmc1602.high, nc.VC10.fmc1602.high)
summary(nc.fit10.fmc1602.high)
nc.overall.fit10.fmc1602.high.int
nc.aic.fmc1602.high.int=AIC(nc.fit10.fmc1602.high)


# low risk conditions ----------------
data.low=subset(data,data$risk=="Low")

nc.fit10.fmc1602.low <- glm(gici10~fmc1602.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+
          rawfood+beach,family=poisson(link="log"),data=data.low)
nc.VC10.fmc1602.low <- cl(data.low,fm=nc.fit10.fmc1602.low, cluster=data.low$hhid)
nc.overall.fit10.fmc1602.low.int <- coeftest(nc.fit10.fmc1602.low, nc.VC10.fmc1602.low)
summary(nc.fit10.fmc1602.low)
nc.overall.fit10.fmc1602.low.int
nc.aic.fmc1602.low.int=AIC(nc.fit10.fmc1602.low)


# -------------------------------------
# fpc 1601
# -------------------------------------
data=all[!is.na(all$fpc1601.pres) ,]

# high risk conditions ----------------
data.high=subset(data,data$risk=="High")

nc.fit10.fpc1601.high <- glm(gici10~fpc1601.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+
          rawfood+beach,family=poisson(link="log"),data=data.high)
nc.VC10.fpc1601.high <- cl(data.high,fm=nc.fit10.fpc1601.high, cluster=data.high$hhid)
nc.overall.fit10.fpc1601.high.int <- coeftest(nc.fit10.fpc1601.high, nc.VC10.fpc1601.high)
summary(nc.fit10.fpc1601.high)
nc.overall.fit10.fpc1601.high.int
nc.aic.fpc1601.high.int=AIC(nc.fit10.fpc1601.high)


# low risk conditions ----------------
data.low=subset(data,data$risk=="Low")

nc.fit10.fpc1601.low <- glm(gici10~fpc1601.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+
          rawfood+beach,family=poisson(link="log"),data=data.low)
nc.VC10.fpc1601.low <- cl(data.low,fm=nc.fit10.fpc1601.low, 
  cluster=data.low$hhid)
nc.overall.fit10.fpc1601.low.int <- coeftest(nc.fit10.fpc1601.low, nc.VC10.fpc1601.low)
summary(nc.fit10.fpc1601.low)
nc.overall.fit10.fpc1601.low.int
nc.aic.fpc1601.low.int=AIC(nc.fit10.fpc1601.low)


# -------------------------------------
# fpc 1602
# -------------------------------------
data=all[!is.na(all$fpc1602.pres) ,]

# high risk conditions ----------------
data.high=subset(data,data$risk=="High")

nc.fit10.fpc1602.high <- glm(gici10~fpc1602.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+
          rawfood+beach,family=poisson(link="log"),data=data.high)
nc.VC10.fpc1602.high <- cl(data.high,fm=nc.fit10.fpc1602.high, cluster=data.high$hhid)
nc.overall.fit10.fpc1602.high.int <- coeftest(nc.fit10.fpc1602.high, nc.VC10.fpc1602.high)
summary(nc.fit10.fpc1602.high)
nc.overall.fit10.fpc1602.high.int
nc.aic.fpc1602.high.int=AIC(nc.fit10.fpc1602.high)


# low risk conditions ----------------
data.low=subset(data,data$risk=="Low")

nc.fit10.fpc1602.low <- glm(gici10~fpc1602.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+
          rawfood+beach,family=poisson(link="log"),data=data.low)
nc.VC10.fpc1602.low <- cl(data.low,fm=nc.fit10.fpc1602.low, cluster=data.low$hhid)
nc.overall.fit10.fpc1602.low.int <- coeftest(nc.fit10.fpc1602.low, nc.VC10.fpc1602.low)
summary(nc.fit10.fpc1602.low)
nc.overall.fit10.fpc1602.low.int
nc.aic.fpc1602.low.int=AIC(nc.fit10.fpc1602.low)


# --------------------------------------
# save the results
# exclude glm objects and data frames
# (they are really large)
# --------------------------------------
save(

  all.n10.fmc1601,all.n10.fmc1602,all.n10.fpc1601,all.n10.fpc1602,

  all.n10.fmc1602.low, all.n10.fmc1602.high,  all.n10.fpc1601.low,
  all.n10.fpc1601.high, all.n10.fpc1602.low,  all.n10.fpc1602.high,
  
  nc.overall.fit10.fmc1601.int,nc.overall.fit10.fmc1602.int,
  nc.overall.fit10.fpc1601.int,nc.overall.fit10.fpc1602.int,

  nc.overall.fit10.fmc1602.low.int,nc.overall.fit10.fmc1602.high.int,
  nc.overall.fit10.fpc1601.low.int,nc.overall.fit10.fpc1601.high.int,
  nc.overall.fit10.fpc1602.low.int,nc.overall.fit10.fpc1602.high.int,

  nc.aic.fmc1601.int, nc.aic.fmc1602.int, nc.aic.fpc1601.int,nc.aic.fpc1602.int,

  nc.aic.fmc1602.low.int, nc.aic.fmc1602.high.int,
  nc.aic.fpc1601.low.int, nc.aic.fpc1601.high.int,
  nc.aic.fpc1602.low.int, nc.aic.fpc1602.high.int, 
  nc.aic.fpc1602.low.int, nc.aic.fpc1602.high.int,

  file="~/dropbox/coliphage/results/rawoutput/regress-10day-body-joint-negcontrol.Rdata"
)


