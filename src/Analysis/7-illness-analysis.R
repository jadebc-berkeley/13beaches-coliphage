##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 11/3/15

# Illness figures for main text
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

# drop individuals with no water quality information
all=subset(all,nowq==0)

# --------------------------------------
# % swimmers in waters with coliphage
# --------------------------------------
sum(table(all$fmc.pres[all$bodycontact=="Yes"]))
sum(table(all$fpc.pres[all$bodycontact=="Yes"]))

prop.table(table(all$fmc.pres[all$bodycontact=="Yes"]))
prop.table(table(all$fpc.pres[all$bodycontact=="Yes"]))

# --------------------------------------
# illness rates
# --------------------------------------
# illness rates among all enrollees
table(all$gici10)
prop.table(table(all$gici10))
prop.table(table(all$gici10,all$beach),2)

# illness rates among swimmers
table(all$gici10[all$bodycontact=="Yes"])
prop.table(table(all$gici10[all$bodycontact=="Yes"]))
prop.table(table(all$gici10[all$bodycontact=="Yes"],all$beach[all$bodycontact=="Yes"]),2)

# illness rates among non-swimmers
table(all$gici10[all$anycontact=="No"])
prop.table(table(all$gici10[all$anycontact=="No"]))
prop.table(table(all$gici10[all$anycontact=="No"],all$beach[all$anycontact=="No"]),2)

# --------------------------------------
# regression for swimming only
# --------------------------------------
swim.reg=glm(gici10~bodycontact+agecat+female+racewhite+gichron+anim_any+gicontactbase+
  rawfood+beach,family=poisson(link="log"),data=all)
# CIR for swimming
exp(swim.reg$coef[["bodycontactYes"]])
exp(swim.reg$coef[["bodycontactYes"]]-
      qnorm(.975)*summary(swim.reg)$coefficients[2,2])
exp(swim.reg$coef[["bodycontactYes"]]+
      qnorm(.975)*summary(swim.reg)$coefficients[2,2])

swimrisk.reg=glm(gici10~bodycontact*risk+agecat+female+racewhite+gichron+anim_any+gicontactbase+
  rawfood+beach,family=poisson(link="log"),data=all)
# CIR for swimming under high risk
exp(swimrisk.reg$coef[["bodycontactYes"]]+
    swimrisk.reg$coef[["bodycontactYes:riskHigh"]])

# CIR for swimming under low risk
exp(swimrisk.reg$coef[["bodycontactYes"]])

lrtest(swim.reg,swimrisk.reg)

save(swim.reg,file="~/Documents/CRG/coliphage/results/rawoutput/regress-10day-swim.Rdata")
