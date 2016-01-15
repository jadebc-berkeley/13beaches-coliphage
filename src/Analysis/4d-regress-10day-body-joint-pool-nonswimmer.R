##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file conducts maximum likelihood regression
# to estimate prevalence ratios

# Results pooled across beaches

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
# subset to non-missing exposure categories
# to make the robust CI calcs work
all=subset(all,all$bodycontact=="Yes" | all$anycontact=="No")


# subset to observations with no missing enterococcus information
all=subset(all,!is.na(all$entero35))


# --------------------------------------
# Calculate the actual Ns for each cell
# and store them for plotting and tables
# --------------------------------------
regN <- function(outcome,exposurecat) {
  sum(table(outcome,exposurecat))
}

all.n10.fmc.joint=regN(all$gici10[!is.na(all$bodyfmc.ent)],
     all$bodyfmc.ent[!is.na(all$bodyfmc.ent)])
all.n10.fpc.joint=regN(all$gici10[!is.na(all$bodyfmc.ent)],
     all$bodyfmc.ent[!is.na(all$bodyfmc.ent)])

data=all[!is.na(all$bodyfmc.ent),]
data.high=subset(data,data$risk=="High")
all.n10.fmc.high.joint = regN(data.high$gici10,data.high$bodyfmc.ent)
data.low=subset(data,data$risk=="Low")
all.n10.fmc.low.joint = regN(data.low$gici10,data.low$bodyfmc.ent)

data=all[!is.na(all$bodyfpc.ent),]
data.high=subset(data,data$risk=="High")
all.n10.fpc.high.joint = regN(data.high$gici10,data.high$bodyfpc.ent)
data.low=subset(data,data$risk=="Low")
all.n10.fpc.low.joint = regN(data.low$gici10,data.low$bodyfpc.ent)

# --------------------------------------
# All conditions

# Estimates pooled across beach
# (can't use the mpreg fn because we 
# need the actual glm returned object 
# for the LR tests)

# 10-day GI illness

# all beaches ----------------

# -------------------------------------
# f- coliphage
# -------------------------------------
all.fit10.fmc <- glm(gici10~swim.fmc.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$swim.fmc.ent),])

all.VC10.fmc <- cl(all[!is.na(all$swim.fmc.ent)],fm=all.fit10.fmc,
    cluster=all$hhid[!is.na(all$swim.fmc.ent)])
overall.fit10.fmc.int <- coeftest(all.fit10.fmc, all.VC10.fmc)
summary(all.fit10.fmc)
overall.fit10.fmc.int

# -------------------------------------
# f+ coliphage
# -------------------------------------
all.fit10.fpc <- glm(gici10~swim.fpc.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+
    rawfood+beach,family=poisson(link="log"),data=all[!is.na(all$swim.fpc.ent),])

all.VC10.fpc <- cl(all[!is.na(all$swim.fpc.ent),],fm=all.fit10.fpc,
    cluster=all$hhid[!is.na(all$swim.fpc.ent)])
overall.fit10.fpc.int <- coeftest(all.fit10.fpc, all.VC10.fpc)
summary(all.fit10.fpc)
overall.fit10.fpc.int

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
data=all[!is.na(all$swim.fmc.ent),]

# high risk conditions ----------------
data.high=subset(data,data$risk=="High")

all.fit10.fmc.high <- glm(gici10~swim.fmc.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+
          rawfood+beach,family=poisson(link="log"),data=data.high)
all.VC10.fmc.high <- cl(data.high,fm=all.fit10.fmc.high, cluster=data.high$hhid)
overall.fit10.fmc.high.int <- coeftest(all.fit10.fmc.high, all.VC10.fmc.high)
summary(all.fit10.fmc.high)
overall.fit10.fmc.high.int

# low risk conditions ----------------
data.low=subset(data,data$risk=="Low")

all.fit10.fmc.low <- glm(gici10~swim.fmc.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+
          rawfood+beach,family=poisson(link="log"),data=data.low)
all.VC10.fmc.low <- cl(data.low,fm=all.fit10.fmc.low, cluster=data.low$hhid)
overall.fit10.fmc.low.int <- coeftest(all.fit10.fmc.low, all.VC10.fmc.low)
summary(all.fit10.fmc.low)
overall.fit10.fmc.low.int

# -------------------------------------
# f+ coliphage
# -------------------------------------
data=all[!is.na(all$swim.fpc.ent),]

# high risk conditions ----------------
data.high=subset(data,data$risk=="High")

all.fit10.fpc.high <- glm(gici10~swim.fpc.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+
          rawfood+beach,family=poisson(link="log"),data=data.high)
all.VC10.fpc.high <- cl(data.high,fm=all.fit10.fpc.high, cluster=data.high$hhid)
overall.fit10.fpc.high.int <- coeftest(all.fit10.fpc.high, all.VC10.fpc.high)
summary(all.fit10.fpc.high)
overall.fit10.fpc.high.int


# low risk conditions ----------------
data.low=subset(data,data$risk=="Low")

all.fit10.fpc.low <- glm(gici10~swim.fpc.ent+agecat+female+racewhite+gichron+anim_any+gicontactbase+
          rawfood+beach,family=poisson(link="log"),data=data.low)
all.VC10.fpc.low <- cl(data.low,fm=all.fit10.fpc.low, 
  cluster=data.low$hhid)
overall.fit10.fpc.low.int <- coeftest(all.fit10.fpc.low, all.VC10.fpc.low)
summary(all.fit10.fpc.low)
overall.fit10.fpc.low.int


# --------------------------------------
# save the results
# exclude glm objects and data frames
# (they are really large)
# --------------------------------------
save(

  all.n10.fmc.joint,all.n10.fpc.joint,

  all.n10.fmc.low.joint, all.n10.fmc.high.joint,
  all.n10.fpc.low.joint, all.n10.fpc.high.joint, 
  
  overall.fit10.fmc.int, overall.fit10.fpc.int,

  overall.fit10.fmc.low.int,overall.fit10.fmc.high.int,
  overall.fit10.fpc.low.int,overall.fit10.fpc.high.int,

  file="~/Documents/CRG/coliphage/13beaches-coliphage/results/rawoutput/regress-10day-body-joint-pool-nonswimmer.Rdata"
)


