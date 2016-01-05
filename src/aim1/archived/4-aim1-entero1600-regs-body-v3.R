



# --------------------------------------
# 4-aim1-entero1600-regs-body.R
# ben arnold (benarnold@berkeley.edu)
#
# description:
# estimate the association between water
# quality indicator concentrations and
# the risk of Diarrhea among swimmers
# for the 13 beaches study
#
# Analyses are conducted for EPA 1600
# Among Swimmers with Body Immmersion
#
# version 3 (9 mar 2015)
# updated to estimate analysis for quartiles of Entero
#
# version 2 (21 feb 2015)
#
# version 1 (28 feb 2014)
#
# --------------------------------------

# --------------------------------------
# preamble
# --------------------------------------

rm(list=ls())
library(sandwich)
library(lmtest)

# source the base functions
source("~/dropbox/13beaches/src/aim1/0-aim1-base-functions.R")

# --------------------------------------
# load the and pre-preprocess the 
# analysis dataset
# (refer to the base functions script
# for details on the pre-processing)
# --------------------------------------

ad <- preprocess.13beaches("~/dropbox/13beaches/data/final/13beaches-analysis.csv")

# --------------------------------------
# esimate risk of GI illness associated
# with exposure to EPA 1600 quintiles
# Body Immersion
# All ages
# --------------------------------------

# subset to non-missing exposure categories
# to make the robust CI calcs work
ah <- subset(ad,ad$bodycontact=="Yes" & is.na(ad$qentero1600)==FALSE)


# tests of interaction by environmental conditions for Avalon and Doheny
av.h.noint <-glm(diarrheaci3~ qentero1600+groundwater +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,family=poisson(link="log"),data=ah[ah$beach=="Avalon",])
av.h <-glm(diarrheaci3~ qentero1600*groundwater +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,family=poisson(link="log"),data=ah[ah$beach=="Avalon",])
	lrtest(av.h.noint,av.h)

dh.h.noint <-glm(diarrheaci3~qentero1600+berm +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,family=poisson(link="log"),data=ah[ah$beach=="Doheny",])
dh.h <-glm(diarrheaci3~ qentero1600*berm +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,family=poisson(link="log"),data=ah[ah$beach=="Doheny",])
	lrtest(dh.h.noint,dh.h)


# --------------------------------------
# Freshwater beaches

# Huntington
hufit <- mpreg(diarrheaci3~qentero1600 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Huntington",vcv=TRUE)

# Silver
sifit <- mpreg(diarrheaci3~qentero1600 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Silver",vcv=TRUE)

# Washington Park
wpfit <- mpreg(diarrheaci3~qentero1600 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Washington Park",vcv=TRUE)

# West
wefit <- mpreg(diarrheaci3~qentero1600 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="West",vcv=TRUE)

# --------------------------------------
# Marine beaches

# Avalon
avfit <- mpreg(diarrheaci3~qentero1600*groundwater +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Avalon",vcv=TRUE)

# Boqueron
bofit <- mpreg(diarrheaci3~qentero1600 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Boqueron",vcv=TRUE)

# Doheny
dhfit <- mpreg(diarrheaci3~qentero1600*berm +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Doheny",vcv=TRUE)

# Edgewater
edfit <- mpreg(diarrheaci3~qentero1600 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Edgewater",vcv=TRUE)

# Fairhope
fafit <- mpreg(diarrheaci3~qentero1600 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Fairhope",vcv=TRUE)

# Goddard
gdfit <- mpreg(diarrheaci3~qentero1600 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Goddard",vcv=TRUE)

# Malibu
mafit <- mpreg(diarrheaci3~qentero1600 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Malibu",vcv=TRUE)

# Mission Bay
mbfit <- mpreg(diarrheaci3~qentero1600 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Mission Bay",vcv=TRUE)

# Surfside
sufit <- mpreg(diarrheaci3~qentero1600 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Surfside",vcv=TRUE)


# --------------------------------------
# Pooled estimates

# all beaches
all.fit <- glm(diarrheaci3~qentero1600 +marine+pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ah)
	all.VC <- cl(ah,fm=all.fit,cluster=ah$hhid)
	all.body <- coeftest(all.fit, all.VC) 

# Interaction model with fresh v. marine beaches
mf.fit <- glm(diarrheaci3~qentero1600*marine +pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ah)
	# commented out robust SE calcs b/c not used (no interaction)
	# mf.VC <- cl(ah,fm=mf.fit,cluster=ah$hhid)
	# mf.head <- coeftest(mf.fit, mf.VC) 
	lrtest(all.fit,mf.fit)

# Interaction model with point v. non-point source beaches
ps.fit <- glm(diarrheaci3~qentero1600*pointsource +marine+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ah)
	# commented out robust SE calcs b/c not used (no interaction)
	ps.VC <- cl(ah,fm=ps.fit,cluster=ah$hhid)
	ps.head <- coeftest(ps.fit, ps.VC) 
	lrtest(all.fit,ps.fit)

# --------------------------------------
# Age-stratified estimates and LR tests of
# interaction

#subset dataset to those without missing ages to facilitate programming below
table(ah$agestrat)
aha <- subset(ah,agestrat!="")

# reduced models for LR tests of indicator x age interactions
noage.fit <- glm(diarrheaci3~qentero1600 +agestrat+marine+pointsource+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=aha)

# Pooled estimate (Age 0-4, 5-10, >10), Head Immersion
agestrat.fit <- glm(diarrheaci3~ qentero1600*agestrat +marine+pointsource+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=aha)
	agestrat.VC <- cl(aha,fm=agestrat.fit,cluster=aha$hhid)
	agestrat.body <- coeftest(agestrat.fit, agestrat.VC) 
	lrtest(noage.fit,agestrat.fit)

# 3-way interaction model of entero x age x pointsource
ageps.fit <- glm(diarrheaci3~ (qentero1600*agestrat*pointsource) +marine+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=aha)
	ageps.VC <- cl(aha,fm=ageps.fit,cluster=aha$hhid)
	ageps.head <- coeftest(ageps.fit,ageps.VC)
	lrtest(noage.fit,ageps.fit)
	lrtest(agestrat.fit,ageps.fit)
	

# --------------------------------------
# Estimate adjusted CIRs
# From the overall pooled model and
# From the Entero X Age interaction model
# --------------------------------------

# function to get Estimates and CIs from a linear combination of regression coefficients
lccalc <- function(lc,x,vcv) {
	# lc : linear combination of coefficients
	# x : log-linear model object returned from coeftest (class=coeftest)
	# vcv : variance-covariance matrix of coefficients for robust SEs
	est <- exp(t(lc)%*%x[,1])
	se  <- sqrt( t(lc)%*%vcv%*%lc )
	lb <- exp(log(est)-1.96*se)
	ub <- exp(log(est)+1.96*se)
	res <- c(est,lb,ub)
	return(res)
}


# --------------------------------------
# Overall adjusted CIRs
# lc1 <- c(1,rep(0,nrow(all.body)-1))
lc2 <- c(0,1,rep(0,nrow(all.body)-2))
lc3 <- c(0,0,1,rep(0,nrow(all.body)-3))
lc4 <- c(0,0,0,1,rep(0,nrow(all.body)-4))

cir.all <- cbind(
	lccalc(lc2,all.body,all.VC),
	lccalc(lc3,all.body,all.VC),
	lccalc(lc4,all.body,all.VC))


# --------------------------------------
# specify the linear combinations of 
# regression coefficients for each estimate
# Reference = Q1 Entero within each age group
 
lc <- rep(0,nrow(agestrat.body))
vs <- rownames(agestrat.body)

# Ages >10, Q1 - Q4
lc10n1 <- lc10n2 <- lc10n3 <- lc10n4 <- lc
# lc10n1[vs=="(Intercept)"] <- 1
lc10n2[vs=="qentero16002"] <- 1
lc10n3[vs=="qentero16003"] <- 1
lc10n4[vs=="qentero16004"] <- 1
lc.age11plus <- list(lc10n2,lc10n3,lc10n4)

# Ages 5 - 10, Q1 - Q4
lc5n1 <- lc5n2 <- lc5n3 <- lc5n4 <- lc
# lc5n1[vs=="agestrat(4, 10]"]<- 1
lc5n2[vs=="qentero16002"|vs=="qentero16002:agestrat(4, 10]"] <- 1
lc5n3[vs=="qentero16003"|vs=="qentero16003:agestrat(4, 10]"] <- 1
lc5n4[vs=="qentero16004"|vs=="qentero16004:agestrat(4, 10]"] <- 1
lc.age5to10 <- list(lc5n2,lc5n3,lc5n4)

# Ages 0 - 4,  Q1 - Q4 
lc0n1 <- lc0n2 <- lc0n3 <- lc0n4 <- lc
# lc0n1[vs=="agestrat(0, 4]"]<- 1
lc0n2[vs=="qentero16002"|vs=="qentero16002:agestrat(0, 4]"] <- 1
lc0n3[vs=="qentero16003"|vs=="qentero16003:agestrat(0, 4]"] <- 1
lc0n4[vs=="qentero16004"|vs=="qentero16004:agestrat(0, 4]"] <- 1
lc.age0to4 <- list(lc0n2,lc0n3,lc0n4)


# calculate all of the estimates
cir.age11plus <- sapply(lc.age11plus,lccalc,x=agestrat.body,vcv=agestrat.VC)
cir.age5to10  <- sapply(lc.age5to10,lccalc,x=agestrat.body,vcv=agestrat.VC)
cir.age0to4   <- sapply(lc.age0to4,lccalc,x=agestrat.body,vcv=agestrat.VC)

colnames(cir.all) <- colnames(cir.age11plus) <- colnames(cir.age5to10) <- colnames(cir.age0to4) <- paste("Entero1600-Q",2:4,sep="")

rownames(cir.all) <- rownames(cir.age11plus) <- rownames(cir.age5to10) <- rownames(cir.age0to4) <- c("CIR","CIRlb","CIRub")



# --------------------------------------
# Estimate adjusted cumulative incidence 
# of Diarrhea
# From the Entero x Age interaction model
# --------------------------------------

# --------------------------------------
# specify the linear combinations of 
# regression coefficients for each estimate
# Reference = Entero Q1 within each Age and point/non-point condition 
lc <- rep(0,nrow(agestrat.body))
vs <- rownames(agestrat.body)

# Ages >10, Q1 - Q4
lc10n1 <- lc10n2 <- lc10n3 <- lc10n4 <- lc
lc10n1[vs=="(Intercept)"] <- 1
lc10n2[vs=="(Intercept)"|vs=="qentero16002"] <- 1
lc10n3[vs=="(Intercept)"|vs=="qentero16003"] <- 1
lc10n4[vs=="(Intercept)"|vs=="qentero16004"] <- 1
lc.age11plusP <- list(lc10n1,lc10n2,lc10n3,lc10n4)

# Ages 5 - 10,  Q1 - Q4
lc5n1 <- lc5n2 <- lc5n3 <- lc5n4 <- lc
lc5n1[vs=="(Intercept)"|vs=="agestrat(4, 10]"]<- 1
lc5n2[vs=="(Intercept)"|vs=="agestrat(4, 10]"|vs=="qentero16002"|vs=="qentero16002:agestrat(4, 10]"] <- 1
lc5n3[vs=="(Intercept)"|vs=="agestrat(4, 10]"|vs=="qentero16003"|vs=="qentero16003:agestrat(4, 10]"] <- 1
lc5n4[vs=="(Intercept)"|vs=="agestrat(4, 10]"|vs=="qentero16004"|vs=="qentero16004:agestrat(4, 10]"] <- 1
lc.age5to10P <- list(lc5n1,lc5n2,lc5n3,lc5n4)

# Ages 0 - 4, Q1 - Q4 
lc0n1 <- lc0n2 <- lc0n3 <- lc0n4 <- lc
lc0n1[vs=="(Intercept)"|vs=="agestrat(0, 4]"]<- 1
lc0n2[vs=="(Intercept)"|vs=="agestrat(0, 4]"|vs=="qentero16002"|vs=="qentero16002:agestrat(0, 4]"] <- 1
lc0n3[vs=="(Intercept)"|vs=="agestrat(0, 4]"|vs=="qentero16003"|vs=="qentero16003:agestrat(0, 4]"] <- 1
lc0n4[vs=="(Intercept)"|vs=="agestrat(0, 4]"|vs=="qentero16004"|vs=="qentero16004:agestrat(0, 4]"] <- 1
lc.age0to4P <- list(lc0n1,lc0n2,lc0n3,lc0n4)


# calculate all of the estimates
p.age11plus <- sapply(lc.age11plusP, lccalc,x=agestrat.body,vcv=agestrat.VC)
p.age5to10  <- sapply(lc.age5to10P, lccalc,x=agestrat.body,vcv=agestrat.VC)
p.age0to4   <- sapply(lc.age0to4P, lccalc,x=agestrat.body,vcv=agestrat.VC)

colnames(p.age11plus) <- colnames(p.age5to10) <- colnames(p.age0to4) <- paste("Entero1600-Q",1:4,sep="")

rownames(p.age11plus) <- rownames(p.age5to10) <- rownames(p.age0to4) <- c("P(GI)","P(GI)lb","P(GI)ub")



# --------------------------------------
# Calculate unadjusted cumulative 
# incidence per 1000 overall and by age
# category
# --------------------------------------
ci.all <- tapply(ah$diarrheaci3,ah$qentero1600,mean)*1000
ci.age0to4 <- tapply(ah$diarrheaci3[ah$agestrat=="(0, 4]"],ah$qentero1600[ah$agestrat=="(0, 4]"],mean)*1000
ci.age5to10 <- tapply(ah$diarrheaci3[ah$agestrat=="(4, 10]"],ah$qentero1600[ah$agestrat=="(4, 10]"],mean)*1000
ci.age11plus <- tapply(ah$diarrheaci3[ah$agestrat==">10"],ah$qentero1600[ah$agestrat==">10"],mean)*1000


# --------------------------------------
# compare unadjusted and adjusted cumulative estimates
# --------------------------------------
round(rbind(p.age11plus[1,]*1000,ci.age11plus))

round(rbind(p.age5to10[1,]*1000,ci.age5to10))

round(rbind(p.age0to4[1,]*1000,ci.age0to4))



# --------------------------------------
# save the objects
# --------------------------------------
save.image("~/dropbox/13beaches/aim1-results/rawoutput/aim1-entero1600-regs-body.Rdata")









