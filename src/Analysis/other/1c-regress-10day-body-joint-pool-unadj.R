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
all=subset(all,all$bodycontact=="Yes")

# subset to observations with no missing enterococcus information
all=subset(all,!is.na(all$entero35))

# --------------------------------------
# Creating joint indicator variable for
# regressions
# --------------------------------------
all$fmc.ent=NA
all$fmc.ent[all$fmc.pres==0 & all$entero35==0]=1
all$fmc.ent[all$fmc.pres==1 & all$entero35==0]=2
all$fmc.ent[all$fmc.pres==0 & all$entero35==1]=3
all$fmc.ent[all$fmc.pres==1 & all$entero35==1]=4
all$fmc.ent=as.factor(all$fmc.ent)

all$fpc.ent=NA
all$fpc.ent[all$fpc.pres==0 & all$entero35==0]=1
all$fpc.ent[all$fpc.pres==1 & all$entero35==0]=2
all$fpc.ent[all$fpc.pres==0 & all$entero35==1]=3
all$fpc.ent[all$fpc.pres==1 & all$entero35==1]=4
all$fpc.ent=as.factor(all$fpc.ent)

# --------------------------------------
# Calculate the actual Ns for each cell
# and store them for plotting and tables
# --------------------------------------
regN <- function(outcome,exposurecat) {
  sum(table(outcome,exposurecat))
}

# all conditions
all.n10.fmc.joint=regN(all$gici10[!is.na(all$fmc.ent)],
     all$fmc.ent[!is.na(all$fmc.ent)])
all.n10.fpc.joint=regN(all$gici10[!is.na(all$fmc.ent)],
     all$fmc.ent[!is.na(all$fmc.ent)])

# stratified by risk conditions
data=all[!is.na(all$fmc.ent),]
data.high=subset(data,data$risk=="High")
all.n10.fmc.high.joint = regN(data.high$gici10,data.high$fmc.ent)
data.low=subset(data,data$risk=="Low")
all.n10.fmc.low.joint = regN(data.low$gici10,data.low$fmc.ent)

data=all[!is.na(all$fpc.ent),]
data.high=subset(data,data$risk=="High")
all.n10.fpc.high.joint = regN(data.high$gici10,data.high$fpc.ent)
data.low=subset(data,data$risk=="Low")
all.n10.fpc.low.joint = regN(data.low$gici10,data.low$fpc.ent)

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
all.fit10.fmc <- glm(gici10~fmc.ent,family=poisson(link="log"),data=all[!is.na(all$fmc.pres),])

all.VC10.fmc <- cl(all[!is.na(all$fmc.pres)],fm=all.fit10.fmc,
    cluster=all$hhid[!is.na(all$fmc.pres)])
overall.fit10.fmc.int <- coeftest(all.fit10.fmc, all.VC10.fmc)
summary(all.fit10.fmc)
overall.fit10.fmc.int
aic.fmc.int=AIC(all.fit10.fmc)


# -------------------------------------
# f+ coliphage
# -------------------------------------
all.fit10.fpc <- glm(gici10~fpc.ent,family=poisson(link="log"),data=all[!is.na(all$fpc.pres),])

all.VC10.fpc <- cl(all[!is.na(all$fpc.pres),],fm=all.fit10.fpc,
    cluster=all$hhid[!is.na(all$fpc.pres)])
overall.fit10.fpc.int <- coeftest(all.fit10.fpc, all.VC10.fpc)
summary(all.fit10.fpc)
overall.fit10.fpc.int
aic.fpc.int=AIC(all.fit10.fpc)

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
data=all[!is.na(all$fmc.pres),]

# high risk conditions ----------------
data.high=subset(data,data$risk=="High")

all.fit10.fmc.high <- glm(gici10~fmc.ent,family=poisson(link="log"),data=data.high)
all.VC10.fmc.high <- cl(data.high,fm=all.fit10.fmc.high, cluster=data.high$hhid)
overall.fit10.fmc.high.int <- coeftest(all.fit10.fmc.high, all.VC10.fmc.high)
summary(all.fit10.fmc.high)
overall.fit10.fmc.high.int
aic.fmc.high.int=AIC(all.fit10.fmc.high)


# low risk conditions ----------------
data.low=subset(data,data$risk=="Low")

all.fit10.fmc.low <- glm(gici10~fmc.ent,family=poisson(link="log"),data=data.low)
all.VC10.fmc.low <- cl(data.low,fm=all.fit10.fmc.low, cluster=data.low$hhid)
overall.fit10.fmc.low.int <- coeftest(all.fit10.fmc.low, all.VC10.fmc.low)
summary(all.fit10.fmc.low)
overall.fit10.fmc.low.int
aic.fmc.low.int=AIC(all.fit10.fmc.low)


# -------------------------------------
# f+ coliphage
# -------------------------------------
data=all[!is.na(all$fpc.pres),]

# high risk conditions ----------------
data.high=subset(data,data$risk=="High")

all.fit10.fpc.high <- glm(gici10~fpc.ent,family=poisson(link="log"),data=data.high)
all.VC10.fpc.high <- cl(data.high,fm=all.fit10.fpc.high, cluster=data.high$hhid)
overall.fit10.fpc.high.int <- coeftest(all.fit10.fpc.high, all.VC10.fpc.high)
summary(all.fit10.fpc.high)
overall.fit10.fpc.high.int
aic.fpc.high.int=AIC(all.fit10.fpc.high)

##### HERE

# low risk conditions ----------------
data.low=subset(data,data$risk=="Low")

all.fit10.fpc.low <- glm(gici10~fpc.ent,family=poisson(link="log"),data=data.low)
all.VC10.fpc.low <- cl(data.low,fm=all.fit10.fpc.low, 
  cluster=data.low$hhid)
overall.fit10.fpc.low.int <- coeftest(all.fit10.fpc.low, all.VC10.fpc.low)
summary(all.fit10.fpc.low)
overall.fit10.fpc.low.int
aic.fpc.low.int=AIC(all.fit10.fpc.low)



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

  aic.fmc.int, aic.fpc.int,
  aic.fmc.low.int, aic.fmc.high.int,
  aic.fpc.low.int, aic.fpc.high.int, 

  file="~/Documents/CRG/coliphage/results/rawoutput/regress-10day-body-joint-pool-unadj.Rdata"
)


