##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file conducts maximum likelihood regression
# to estimate prevalence ratios
# for enterococcus concentration

# Results pooled across beaches

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
all=subset(all,all$bodycontact=="Yes")

# --------------------------------------
# Calculate the actual Ns for each cell
# and store them for plotting and tables
# --------------------------------------
regN <- function(outcome,exposurecat) {
  sum(table(outcome,exposurecat))
}


# n's pooled across beach ---------------------------------------
all.n10.ent.fmc1601 = regN(all$gici10[!is.na(all$entero) & !is.na(all$fmc1601.pres)],
                       all$fmc1601.pres[!is.na(all$entero) & !is.na(all$fmc1601.pres)])
all.n10.ent.fmc1602 = regN(all$gici10[!is.na(all$entero) & !is.na(all$fmc1602.pres)],
                       all$fmc1602.pres[!is.na(all$entero) & !is.na(all$fmc1602.pres)])
all.n10.ent.fpc1601 = regN(all$gici10[!is.na(all$entero) & !is.na(all$fpc1601.pres)],
                       all$fpc1601.pres[!is.na(all$entero) & !is.na(all$fpc1601.pres)])
all.n10.ent.fpc1602 = regN(all$gici10[!is.na(all$entero) & !is.na(all$fpc1602.pres)],
                       all$fpc1602.pres[!is.na(all$entero) & !is.na(all$fpc1602.pres)])

# pooled n's by risk level---------------------------------------
data=all[!is.na(all$entero) & !is.na(all$fmc1601.pres),]
data.high=subset(data,data$risk=="High")
all.n10.ent.fmc1601.high = regN(data.high$gici10,data.high$fmc1601.pres)
data.low=subset(data,data$risk=="Low")
all.n10.ent.fmc1601.low = regN(data.low$gici10,data.low$fmc1601.pres)

data=all[!is.na(all$entero) & !is.na(all$fmc1602.pres),]
data.high=subset(data,data$risk=="High")
all.n10.ent.fmc1602.high = regN(data.high$gici10,data.high$fmc1602.pres)
data.low=subset(data,data$risk=="Low")
all.n10.ent.fmc1602.low = regN(data.low$gici10,data.low$fmc1602.pres)

data=all[!is.na(all$entero) & !is.na(all$fpc1601.pres),]
data.high=subset(data,data$risk=="High")
all.n10.ent.fpc1601.high = regN(data.high$gici10,data.high$fpc1601.pres)
data.low=subset(data,data$risk=="Low")
all.n10.ent.fpc1601.low = regN(data.low$gici10,data.low$fpc1601.pres)

data=all[!is.na(all$entero) & !is.na(all$fpc1602.pres),]
data.high=subset(data,data$risk=="High")
all.n10.ent.fpc1602.high = regN(data.high$gici10,data.high$fpc1602.pres)
data.low=subset(data,data$risk=="Low")
all.n10.ent.fpc1602.low = regN(data.low$gici10,data.low$fpc1602.pres)

# --------------------------------------
# Estimates pooled across beach
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness
# --------------------------------------

# all beaches --------------------------------
all.fit10.entero.fmc1601 <- glm(gici10~entero+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$entero) & 
    !is.na(all$fmc1601) ,])

all.VC10.entero.fmc1601 <- cl(all[!is.na(all$entero) & !is.na(all$fmc1601) & 
    all$beach!="Avalon"],fm=all.fit10.entero.fmc1601, cluster=
      all$hhid[!is.na(all$entero)  & !is.na(all$fmc1601) ])
overall.fit10.entero.fmc1601 <- coeftest(all.fit10.entero.fmc1601, all.VC10.entero.fmc1601)
summary(all.fit10.entero.fmc1601)
overall.fit10.entero.fmc1601
aic.entero.fmc1601=AIC(all.fit10.entero.fmc1601)

all.fit10.entero.fmc1602 <- glm(gici10~entero+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$entero) & 
    !is.na(all$fmc1602),])

all.VC10.entero.fmc1602 <- cl(all[!is.na(all$entero) & !is.na(all$fmc1602)],fm=all.fit10.entero.fmc1602, cluster=
      all$hhid[!is.na(all$entero)  & !is.na(all$fmc1602)])
overall.fit10.entero.fmc1602 <- coeftest(all.fit10.entero.fmc1602, all.VC10.entero.fmc1602)
summary(all.fit10.entero.fmc1602)
overall.fit10.entero.fmc1602
aic.entero.fmc1602=AIC(all.fit10.entero.fmc1602)

all.fit10.entero.fpc1601 <- glm(gici10~entero+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$entero) & 
    !is.na(all$fpc1601) ,])

all.VC10.entero.fpc1601 <- cl(all[!is.na(all$entero) & !is.na(all$fpc1601) & 
    all$beach!="Goddard"],fm=all.fit10.entero.fpc1601, cluster=
      all$hhid[!is.na(all$entero)  & !is.na(all$fpc1601) ])
overall.fit10.entero.fpc1601 <- coeftest(all.fit10.entero.fpc1601, all.VC10.entero.fpc1601)
summary(all.fit10.entero.fpc1601)
overall.fit10.entero.fpc1601
aic.entero.fpc1601=AIC(all.fit10.entero.fpc1601)

all.fit10.entero.fpc1602 <- glm(gici10~entero+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$entero) & 
    !is.na(all$fpc1602) ,])

all.VC10.entero.fpc1602 <- cl(all[!is.na(all$entero) & !is.na(all$fpc1602) & 
    all$beach!="Malibu"],fm=all.fit10.entero.fpc1602, cluster=
      all$hhid[!is.na(all$entero)  & !is.na(all$fpc1602) ])
overall.fit10.entero.fpc1602 <- coeftest(all.fit10.entero.fpc1602, all.VC10.entero.fpc1602)
summary(all.fit10.entero.fpc1602)
overall.fit10.entero.fpc1602
aic.entero.fpc1602=AIC(all.fit10.entero.fpc1602)

# FMC 1601 #####################
# high risk conditions --------------------------------
data=all[!is.na(all$entero) & !is.na(all$fmc1601) ,]
data.high=subset(data,data$risk=="High")
all.fit10.entero.high.fmc1601 <- glm(gici10~entero+agecat+female+racewhite+gichron+anim_any+gicontactbase+
                                       rawfood,family=poisson(link="log"),data=data.high)

all.VC10.entero.high.fmc1601 <- cl(data.high,fm=all.fit10.entero.high.fmc1601, cluster=data.high$hhid)
overall.fit10.entero.high.fmc1601 <- coeftest(all.fit10.entero.high.fmc1601, all.VC10.entero.high.fmc1601)
summary(all.fit10.entero.high.fmc1601)
overall.fit10.entero.high.fmc1601
aic.entero.high.fmc1601=AIC(all.fit10.entero.high.fmc1601)

# low risk conditions --------------------------------
data.low=subset(data,data$risk=="Low")
all.fit10.entero.low.fmc1601 <- glm(gici10~entero+agecat+female+racewhite+gichron+anim_any+gicontactbase+
                                      rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC10.entero.low.fmc1601 <- cl(data.low,fm=all.fit10.entero.low.fmc1601, cluster=data.low$hhid)
overall.fit10.entero.low.fmc1601 <- coeftest(all.fit10.entero.low.fmc1601, all.VC10.entero.low.fmc1601)
summary(all.fit10.entero.low.fmc1601)
overall.fit10.entero.low.fmc1601
aic.entero.low.fmc1601=AIC(all.fit10.entero.low.fmc1601)


# FMC 1602 #####################
# high risk conditions --------------------------------
data=all[!is.na(all$entero) & !is.na(all$fmc1602) ,]
data.high=subset(data,data$risk=="High")
all.fit10.entero.high.fmc1602 <- glm(gici10~entero+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC10.entero.high.fmc1602 <- cl(data.high,fm=all.fit10.entero.high.fmc1602, cluster=data.high$hhid)
overall.fit10.entero.high.fmc1602 <- coeftest(all.fit10.entero.high.fmc1602, all.VC10.entero.high.fmc1602)
summary(all.fit10.entero.high.fmc1602)
overall.fit10.entero.high.fmc1602
aic.entero.high.fmc1602=AIC(all.fit10.entero.high.fmc1602)

# low risk conditions --------------------------------
data.low=subset(data,data$risk=="Low")
all.fit10.entero.low.fmc1602 <- glm(gici10~entero+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC10.entero.low.fmc1602 <- cl(data.low,fm=all.fit10.entero.low.fmc1602, cluster=data.low$hhid)
overall.fit10.entero.low.fmc1602 <- coeftest(all.fit10.entero.low.fmc1602, all.VC10.entero.low.fmc1602)
summary(all.fit10.entero.low.fmc1602)
overall.fit10.entero.low.fmc1602
aic.entero.low.fmc1602=AIC(all.fit10.entero.low.fmc1602)

# FPC 1601 #####################
# high risk conditions --------------------------------
data=all[!is.na(all$entero) & !is.na(all$fpc1601) ,]
data.high=subset(data,data$risk=="High")
all.fit10.entero.high.fpc1601 <- glm(gici10~entero+agecat+female+racewhite+gichron+anim_any+gicontactbase+
                                rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC10.entero.high.fpc1601 <- cl(data.high,fm=all.fit10.entero.high.fpc1601, cluster=data.high$hhid)
overall.fit10.entero.high.fpc1601 <- coeftest(all.fit10.entero.high.fpc1601, all.VC10.entero.high.fpc1601)
summary(all.fit10.entero.high.fpc1601)
overall.fit10.entero.high.fpc1601
aic.entero.high.fpc1601=AIC(all.fit10.entero.high.fpc1601)

# low risk conditions --------------------------------
data.low=subset(data,data$risk=="Low")
all.fit10.entero.low.fpc1601 <- glm(gici10~entero+agecat+female+racewhite+gichron+anim_any+gicontactbase+
                               rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC10.entero.low.fpc1601 <- cl(data.low,fm=all.fit10.entero.low.fpc1601, cluster=data.low$hhid)
overall.fit10.entero.low.fpc1601 <- coeftest(all.fit10.entero.low.fpc1601, all.VC10.entero.low.fpc1601)
summary(all.fit10.entero.low.fpc1601)
overall.fit10.entero.low.fpc1601
aic.entero.low.fpc1601=AIC(all.fit10.entero.low.fpc1601)

# FPC 1602 #####################
# high risk conditions --------------------------------
data=all[!is.na(all$entero) & !is.na(all$fpc1602) ,]
data.high=subset(data,data$risk=="High")
all.fit10.entero.high.fpc1602 <- glm(gici10~entero+agecat+female+racewhite+gichron+anim_any+gicontactbase+
                                rawfood+beach,family=poisson(link="log"),data=data.high)

all.VC10.entero.high.fpc1602 <- cl(data.high,fm=all.fit10.entero.high.fpc1602, cluster=data.high$hhid)
overall.fit10.entero.high.fpc1602 <- coeftest(all.fit10.entero.high.fpc1602, all.VC10.entero.high.fpc1602)
summary(all.fit10.entero.high.fpc1602)
overall.fit10.entero.high.fpc1602
aic.entero.high.fpc1602=AIC(all.fit10.entero.high.fpc1602)

# low risk conditions --------------------------------
data.low=subset(data,data$risk=="Low")
all.fit10.entero.low.fpc1602 <- glm(gici10~entero+agecat+female+racewhite+gichron+anim_any+gicontactbase+
                               rawfood+beach,family=poisson(link="log"),data=data.low)

all.VC10.entero.low.fpc1602 <- cl(data.low,fm=all.fit10.entero.low.fpc1602, cluster=data.low$hhid)
overall.fit10.entero.low.fpc1602 <- coeftest(all.fit10.entero.low.fpc1602, all.VC10.entero.low.fpc1602)
summary(all.fit10.entero.low.fpc1602)
overall.fit10.entero.low.fpc1602
aic.entero.low.fpc1602=AIC(all.fit10.entero.low.fpc1602)

#-------------------------------------------
# estimating probability of illness over
# observed range of exposure
#-------------------------------------------

iter=1000

set.seed(92203789)

# fmc 1601 ---------------------------------
all.fmc1601.pY = boot.pY(fmla=gici10~entero+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+rawfood+beach,dat=all[!is.na(all$fmc1601),],nameX="entero",
    ID=all[!is.na(all$fmc1601),"hhid"],iter)
data=all[!is.na(all$fmc1601),]
data.high=subset(data,data$risk=="High")
all.fmc1601.pY.high = boot.pY(fmla=gici10~entero+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+ rawfood,dat=data.high,nameX="entero",ID=data.high$hhid,iter)
data.low=subset(data,data$risk=="Low")
all.fmc1601.pY.low = boot.pY(fmla=gici10~entero+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+ rawfood+beach,dat=data.low,nameX="entero",ID=data.low$hhid,iter)

# fmc 1602 ---------------------------------
all.fmc1602.pY = boot.pY(fmla=gici10~entero+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+rawfood+beach,dat=all[!is.na(all$fmc1602),],nameX="entero",
    ID=all[!is.na(all$fmc1602),"hhid"],iter)
data=all[!is.na(all$fmc1602),]
data.high=subset(data,data$risk=="High")
all.fmc1602.pY.high = boot.pY(fmla=gici10~entero+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+ rawfood+beach,dat=data.high,nameX="entero",ID=data.high$hhid,iter)
data.low=subset(data,data$risk=="Low")
all.fmc1602.pY.low = boot.pY(fmla=gici10~entero+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+ rawfood+beach,dat=data.low,nameX="entero",ID=data.low$hhid,iter)

# fpc 1601 ---------------------------------
all.fpc1601.pY = boot.pY(fmla=gici10~entero+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+rawfood+beach,dat=all[!is.na(all$fpc1601),],nameX="entero",
    ID=all[!is.na(all$fpc1601),"hhid"],iter)
data=all[!is.na(all$fpc1601),]
data.high=subset(data,data$risk=="High")
all.fpc1601.pY.high = boot.pY(fmla=gici10~entero+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+ rawfood+beach,dat=data.high,nameX="entero",ID=data.high$hhid,iter)
data.low=subset(data,data$risk=="Low")
all.fpc1601.pY.low = boot.pY(fmla=gici10~entero+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+ rawfood+beach,dat=data.low,nameX="entero",ID=data.low$hhid,iter)

# fpc 1602 ---------------------------------
all.fpc1602.pY = boot.pY(fmla=gici10~entero+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+rawfood+beach,dat=all[!is.na(all$fpc1602),],nameX="entero",
    ID=all[!is.na(all$fpc1602),"hhid"],iter)
data=all[!is.na(all$fpc1602),]
data.high=subset(data,data$risk=="High")
all.fpc1602.pY.high = boot.pY(fmla=gici10~entero+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+ rawfood+beach,dat=data.high,nameX="entero",ID=data.high$hhid,iter)
data.low=subset(data,data$risk=="Low")
all.fpc1602.pY.low = boot.pY(fmla=gici10~entero+agecat+female+racewhite+gichron+anim_any+
    gicontactbase+ rawfood+beach,dat=data.low,nameX="entero",ID=data.low$hhid,iter)


# --------------------------------------
# save the results
# exclude glm objects and data frames
# (they are really large)
# --------------------------------------

save(
all.n10.ent.fmc1601,all.n10.ent.fmc1602,all.n10.ent.fpc1601,all.n10.ent.fpc1602,
all.n10.ent.fmc1601.high, all.n10.ent.fmc1601.low,
all.n10.ent.fmc1602.high, all.n10.ent.fmc1602.low,
all.n10.ent.fpc1601.high, all.n10.ent.fpc1601.low,
all.n10.ent.fpc1602.high, all.n10.ent.fpc1602.low,

file="~/Documents/CRG/coliphage/results/rawoutput/regress-10day-continuous-body-entero-n.Rdata")

save(

  overall.fit10.entero.fmc1601,overall.fit10.entero.fmc1602,
  overall.fit10.entero.fpc1601,overall.fit10.entero.fpc1602,
  
  overall.fit10.entero.high.fmc1601,overall.fit10.entero.high.fmc1602,
  overall.fit10.entero.high.fpc1601,overall.fit10.entero.high.fpc1602,
  
  overall.fit10.entero.low.fmc1601,overall.fit10.entero.low.fmc1602,
  overall.fit10.entero.low.fpc1601,overall.fit10.entero.low.fpc1602,
  
  all.VC10.entero.fmc1601,all.VC10.entero.fmc1602,all.VC10.entero.fpc1601,
  all.VC10.entero.fpc1602,
  all.VC10.entero.high.fmc1602,
  all.VC10.entero.high.fpc1601,all.VC10.entero.high.fpc1602,
  all.VC10.entero.low.fmc1602,all.VC10.entero.low.fpc1601,
  all.VC10.entero.low.fpc1602,

    file="~/Documents/CRG/coliphage/results/rawoutput/regress-10day-continuous-CIR-body-entero.Rdata"
)

save(
  all.fmc1601.pY,all.fmc1601.pY.high,all.fmc1601.pY.low,
  all.fmc1602.pY,all.fmc1602.pY.high,all.fmc1602.pY.low,
  all.fpc1601.pY,all.fpc1601.pY.high,all.fpc1601.pY.low,all.fpc1602.pY,
  all.fpc1602.pY.high,all.fpc1602.pY.low,

  file="~/Documents/CRG/coliphage/results/rawoutput/regress-10day-continuous-body-entero.Rdata"
)



