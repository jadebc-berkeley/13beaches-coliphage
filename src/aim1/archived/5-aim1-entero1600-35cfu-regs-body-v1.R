



# --------------------------------------
# 5-aim1-entero1600-35cfu-regs-body.R
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
# The exposure categories are <=35 and >35 CFU/100ml
#
# version 1 (18 mar 2015)
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
# esimate risk of Diarrhea associated
# with exposure to EPA 1600 above/below 35 CFU/100ml
# Body Immersion
# All ages
# --------------------------------------

# subset to non-missing exposure categories
# to make the robust CI calcs work
ah <- subset(ad,ad$bodycontact=="Yes")
	dim(ad)
	dim(ah)
ah <- subset(ah,is.na(ah$entero35)==FALSE)
	dim(ah)

# tests of interaction by environmental conditions for Avalon and Doheny
av.h.noint <-glm(diarrheaci3~ entero35+groundwater +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,family=poisson(link="log"),data=ah[ah$beach=="Avalon",])
av.h <-glm(diarrheaci3~ entero35*groundwater +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,family=poisson(link="log"),data=ah[ah$beach=="Avalon",])
	lrtest(av.h.noint,av.h)

dh.h.noint <-glm(diarrheaci3~ entero35+berm +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,family=poisson(link="log"),data=ah[ah$beach=="Doheny",])
dh.h <-glm(diarrheaci3~ entero35*berm +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,family=poisson(link="log"),data=ah[ah$beach=="Doheny",])
	lrtest(dh.h.noint,dh.h)


# --------------------------------------
# Freshwater beaches

# Huntington
hufit <- mpreg.beach(diarrheaci3~entero35 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Huntington",vcv=TRUE)

# Silver
sifit <- mpreg.beach(diarrheaci3~entero35 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Silver",vcv=TRUE)

# Washington Park
wpfit <- mpreg.beach(diarrheaci3~entero35 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Washington Park",vcv=TRUE)

# West
wefit <- mpreg.beach(diarrheaci3~entero35 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="West",vcv=TRUE)

# --------------------------------------
# Marine beaches

# Avalon
avfit <- mpreg.beach(diarrheaci3~entero35*groundwater +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Avalon",vcv=TRUE)

# Boqueron
bofit <- mpreg.beach(diarrheaci3~entero35 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Boqueron",vcv=TRUE)

# Doheny
dhfit <- mpreg.beach(diarrheaci3~entero35*berm +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Doheny",vcv=TRUE)

# Edgewater
edfit <- mpreg.beach(diarrheaci3~entero35 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Edgewater",vcv=TRUE)

# Fairhope
fafit <- mpreg.beach(diarrheaci3~entero35 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Fairhope",vcv=TRUE)

# Goddard
gdfit <- mpreg.beach(diarrheaci3~entero35 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Goddard",vcv=TRUE)

# Malibu
mafit <- mpreg.beach(diarrheaci3~entero35 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Malibu",vcv=TRUE)

# Mission Bay
mbfit <- mpreg.beach(diarrheaci3~entero35 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Mission Bay",vcv=TRUE)

# Surfside
sufit <- mpreg.beach(diarrheaci3~entero35 +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Surfside",vcv=TRUE)


# --------------------------------------
# Pooled estimates

# all beaches
all.fit <- glm(diarrheaci3~entero35 +marine+pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ah)
	all.VC <- cl(ah,fm=all.fit,cluster=ah$hhid)
	all.body <- coeftest(all.fit, all.VC) 

# Interaction model with fresh v. marine beaches
mf.fit <- glm(diarrheaci3~entero35*marine +pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ah)
	# commented out robust SE calcs b/c not used (no interaction)
	# mf.VC <- cl(ah,fm=mf.fit,cluster=ah$hhid)
	# mf.body <- coeftest(mf.fit, mf.VC) 
	lrtest(all.fit,mf.fit)

# Interaction model with point v. non-point source beaches
ps.fit <- glm(diarrheaci3~entero35*pointsource +marine+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ah)
	# commented out robust SE calcs b/c not used (no interaction)
	ps.VC <- cl(ah,fm=ps.fit,cluster=ah$hhid)
	ps.body <- coeftest(ps.fit, ps.VC) 
	lrtest(all.fit,ps.fit)

# --------------------------------------
# Age-stratified estimates and LR tests of
# interaction

# reduced models for LR tests of indicator x age interactions
noage.fit <- glm(diarrheaci3~entero35 +agestrat+marine+pointsource+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ah)

# Pooled estimate (Age 0-4, 5-10, >10), Body Immersion
agestrat.fit <- glm(diarrheaci3~ entero35*agestrat +marine+pointsource+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ah)
	agestrat.VC <- cl(ah,fm=agestrat.fit,cluster=ah$hhid)
	agestrat.body <- coeftest(agestrat.fit, agestrat.VC) 
	lrtest(noage.fit,agestrat.fit)

# 3-way interaction model of entero x age x pointsource
ageps.fit <- glm(diarrheaci3~ (entero35*agestrat*pointsource) +marine+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,family=poisson(link="log"),data=ah)
	ageps.VC <- cl(ah,fm=ageps.fit,cluster=ah$hhid)
	ageps.body <- coeftest(ageps.fit,ageps.VC)
	lrtest(noage.fit,ageps.fit)
	lrtest(agestrat.fit,ageps.fit)
	

# --------------------------------------
# Stratified Models 
# based on tests of interaction (above)
# stratify the results by non-point and
# point source conditions
# --------------------------------------

# Pooled Non-point source and point-source
nps <- mpreg(diarrheaci3~entero35 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ah[ah$pointsource=="No",],vcv=T)
ps <- mpreg(diarrheaci3~entero35 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ah[ah$pointsource=="Yes",],vcv=T)

# Non-point source estimates by age group
nps0to4 <- mpreg(diarrheaci3~entero35 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ah[ah$pointsource=="No" & ah$agestrat=="(0, 4]",],vcv=T)

nps5to10 <- mpreg(diarrheaci3~entero35 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ah[ah$pointsource=="No" & ah$agestrat=="(4, 10]",],vcv=T)

nps11plus <- mpreg(diarrheaci3~entero35 +female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ah[ah$pointsource=="No" & ah$agestrat==">10",],vcv=T)

# Point-source estimates by age group
ps0to4 <- mpreg(diarrheaci3~entero35 +marine+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ah[ah$pointsource=="Yes" & ah$agestrat=="(0, 4]",],vcv=T)

ps5to10 <- mpreg(diarrheaci3~entero35 +marine+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ah[ah$pointsource=="Yes" & ah$agestrat=="(4, 10]",],vcv=T)

ps11plus <- mpreg(diarrheaci3~entero35 +marine+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach,dat=ah[ah$pointsource=="Yes" & ah$agestrat==">10",],vcv=T)


# --------------------------------------
# Estimate adjusted CIRs
# From the overall pooled model and
# From stratified models
# --------------------------------------

# function to get Estimates and CIs from a linear combination of regression coefficients
# more useful if you need to do lincoms -- here we actually don't with the fully stratified models
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


cir.all <- rbind(
	lccalc(c(0,1,rep(0,nrow(nps$fit)-2)),nps$fit,nps$vcovCL),
	lccalc(c(0,1,rep(0,nrow(ps$fit)-2)),ps$fit,ps$vcovCL),
	lccalc(c(0,1,rep(0,nrow(all.body)-2)),all.body,all.VC)
	)



# --------------------------------------
# Linear combinations for CIRs

# overall, all ages
lc.all <- c(0,1,rep(0,nrow(all.body)-2))
# overall, stratified by age
vs <- rownames(agestrat.body)
lc.11plus <- lc.5to10 <- lc.0to4 <- rep(0,nrow(agestrat.body))
lc.11plus[vs=="entero35"] <- 1
lc.5to10[vs=="entero35"|vs=="entero35:agestrat(4, 10]"] <- 1
lc.0to4[vs=="entero35"|vs=="entero35:agestrat(0, 4]"] <- 1
# non-point source and point source, all ages
vs <- rownames(ps.body)
lc.nps <- lc.ps <- rep(0,nrow(ps.body))
lc.nps[vs=="entero35"] <- 1
lc.ps[vs=="entero35"|vs=="entero35:pointsourceYes"] <- 1
# non-point source and point source, stratified by age
vs <- rownames(ageps.body)
lc.nps.11plus <- lc.nps.5to10 <- lc.nps.0to4 <- lc.ps.11plus <- lc.ps.5to10 <- lc.ps.0to4 <- rep(0,nrow(ageps.body))
lc.nps.11plus[vs=="entero35"]   <- 1
lc.nps.5to10[vs=="entero35"|vs=="entero35:agestrat(4, 10]"] <- 1
lc.nps.0to4[vs=="entero35"|vs=="entero35:agestrat(0, 4]"]  <- 1
lc.ps.11plus[vs=="entero35"|vs=="entero35:pointsourceYes"]   <- 1
lc.ps.5to10[vs=="entero35"|vs=="entero35:pointsourceYes"|vs=="entero35:agestrat(4, 10]"|vs=="entero35:agestrat(4, 10]:pointsourceYes"] <- 1
lc.ps.0to4[vs=="entero35"|vs=="entero35:pointsourceYes"|vs=="entero35:agestrat(0, 4]"|vs=="entero35:agestrat(0, 4]:pointsourceYes"]  <- 1

# bind all of the CIRs together, organized by
# non-point source vs. point source and then age, then overall
cir.all <- rbind(
	lccalc(lc.nps,ps.body,ps.VC),
	lccalc(lc.ps,ps.body,ps.VC),
	lccalc(lc.all,all.body,all.VC)
	)

cir.age0to4 <- rbind(
	lccalc(lc.nps.0to4,ageps.body,ageps.VC),
	lccalc(lc.ps.0to4,ageps.body,ageps.VC),
	lccalc(lc.0to4,agestrat.body,agestrat.VC)
	)
	
cir.age5to10 <- rbind(
	lccalc(lc.nps.5to10,ageps.body,ageps.VC),
	lccalc(lc.ps.5to10,ageps.body,ageps.VC),
	lccalc(lc.5to10,agestrat.body,agestrat.VC)
	)
	
cir.age11plus <- rbind(
	lccalc(lc.nps.11plus,ageps.body,ageps.VC),
	lccalc(lc.ps.11plus,ageps.body,ageps.VC),
	lccalc(lc.11plus,agestrat.body,agestrat.VC)
	)
	
colnames(cir.all) <- colnames(cir.age11plus) <- colnames(cir.age5to10) <- colnames(cir.age0to4) <- c("CIR","CIRlb","CIRub")
rownames(cir.all) <- rownames(cir.age11plus) <- rownames(cir.age5to10) <- rownames(cir.age0to4) <-  c("Non-point source","Point source","Overall")

#### LEFT OFF HERE ####
# need to calculate stratified estimates for avalon + doheny by point/nonpoint
# for forest plot summaries by point/non-point conditions

# --------------------------------------
# Estimate adjusted cumulative incidence 
# of Diarrhea
# From the Entero x Point source x Age interaction model
# --------------------------------------

# --------------------------------------
# specify the linear combinations of 
# regression coefficients for each estimate
# Reference = Entero <=35 CFU/100ml within each Age and point/non-point condition

# overall, all ages
lc.all.ci1 <- c(1,0,rep(0,nrow(all.body)-2))
lc.all.ci2 <- c(1,1,rep(0,nrow(all.body)-2))

# overall, stratified by age
vs <- rownames(agestrat.body)
lc.11plus.ci1 <- lc.11plus.ci2 <- lc.5to10.ci1 <- lc.5to10.ci2 <- lc.0to4.ci1 <- lc.0to4.ci2 <- rep(0,nrow(agestrat.body))
lc.11plus.ci1[vs=="(Intercept)"] <- 1
lc.11plus.ci2[vs=="(Intercept)"|vs=="entero35"] <- 1
lc.5to10.ci1[vs=="(Intercept)"|vs=="agestrat(4, 10]"] <- 1
lc.5to10.ci2[vs=="(Intercept)"|vs=="agestrat(4, 10]"|vs=="entero35"|vs=="entero35:agestrat(4, 10]"] <- 1
lc.0to4.ci1[vs=="(Intercept)"|vs=="agestrat(0, 4]"] <- 1
lc.0to4.ci2[vs=="(Intercept)"|vs=="agestrat(0, 4]"|vs=="entero35"|vs=="entero35:agestrat(0, 4]"] <- 1

# non-point source and point source, all ages
vs <- rownames(ps.body)
lc.nps.ci1 <- lc.nps.ci2 <- lc.ps.ci1 <- lc.ps.ci2 <- rep(0,nrow(ps.body))
lc.nps.ci1[vs=="(Intercept)"] <- 1
lc.nps.ci2[vs=="(Intercept)"|vs=="entero35"] <- 1
lc.ps.ci1[vs=="(Intercept)"|vs=="pointsourceYes"] <- 1
lc.ps.ci2[vs=="(Intercept)"|vs=="pointsourceYes"|vs=="entero35"|vs=="entero35:pointsourceYes"] <- 1

# age stratified

LEFT OFF HERE
NEED TO UPDATE NAMES OF THE LINEAR COMBINATION OBJECTS, THEN FINISH THE RESULTS COMPILATION BELOW

vs <- rownames(ageps.body)
lc10n1 <- lc10n2 <- lc5n1 <- lc5n2 <- lc0n1 <- lc0n2 <- rep(0,nrow(ageps.body))

# non-point source conditions
lc10n1[vs=="(Intercept)"] <- 1
lc10n2[vs=="(Intercept)"|vs=="entero35"] <- 1 
lc5n1[vs=="(Intercept)"|vs=="agestrat(4, 10]"] <- 1 
lc5n2[vs=="(Intercept)"|vs=="agestrat(4, 10]"|vs=="entero35"|vs=="entero35:agestrat(4, 10]"] <- 1
lc0n1[vs=="(Intercept)"|vs=="agestrat(0, 4]"] <- 1 
lc0n2[vs=="(Intercept)"|vs=="agestrat(0, 4]"|vs=="entero35"|vs=="entero35:agestrat(0, 4]"] <- 1

# point source conditions
lc10p1 <- lc10p2 <- lc5p1 <- lc5p2 <- lc0p1 <- lc0p2 <- rep(0,nrow(ageps.body))
lc10p1[vs=="(Intercept)"|vs=="pointsourceYes"] <- 1
lc10p2[vs=="(Intercept)"|vs=="pointsourceYes"|vs=="entero35"|vs=="entero35:pointsourceYes"] <- 1 
lc5p1[vs=="(Intercept)"|vs=="pointsourceYes"|vs=="agestrat(4, 10]"|vs=="agestrat(4, 10]:pointsourceYes"] <- 1 
lc5p2[vs=="(Intercept)"|vs=="pointsourceYes"|vs=="agestrat(4, 10]"|vs=="agestrat(4, 10]:pointsourceYes"|vs=="entero35"|vs=="entero35:pointsourceYes"|vs=="entero35:agestrat(4, 10]:pointsourceYes"] <- 1 
lc0p1[vs=="(Intercept)"|vs=="pointsourceYes"|vs=="agestrat(0, 4]"|vs=="agestrat(0, 4]:pointsourceYes"] <- 1
lc0p2[vs=="(Intercept)"|vs=="pointsourceYes"|vs=="agestrat(0, 4]"|vs=="agestrat(0, 4]:pointsourceYes"|vs=="entero35"|vs=="entero35:pointsourceYes"|vs=="entero35:agestrat(0, 4]:pointsourceYes"] <- 1



# bind all of the adjusted CIs together, organized by
# non-point source vs. point source and then age, then overall
ci.all <- rbind(
	lccalc(lcn1,ps.body,ps.VC),
	lccalc(lc.ps,ps.body,ps.VC),
	lccalc(lc.all,all.body,all.VC)
	)

ci.age0to4 <- rbind(
	lccalc(lc.nps.0to4,ageps.body,ageps.VC),
	lccalc(lc.ps.0to4,ageps.body,ageps.VC),
	lccalc(lc.0to4,agestrat.body,agestrat.VC)
	)
	
ci.age5to10 <- rbind(
	lccalc(lc.nps.5to10,ageps.body,ageps.VC),
	lccalc(lc.ps.5to10,ageps.body,ageps.VC),
	lccalc(lc.5to10,agestrat.body,agestrat.VC)
	)
	
ci.age11plus <- rbind(
	lccalc(lc.nps.11plus,ageps.body,ageps.VC),
	lccalc(lc.ps.11plus,ageps.body,ageps.VC),
	lccalc(lc.11plus,agestrat.body,agestrat.VC)
	)
	
colnames(cir.all) <- colnames(cir.age11plus) <- colnames(cir.age5to10) <- colnames(cir.age0to4) <- c("CIR","CIRlb","CIRub")
rownames(cir.all) <- rownames(cir.age11plus) <- rownames(cir.age5to10) <- rownames(cir.age0to4) <-  c("Non-point source","Point source","Overall")


# Ages >10 
lc.age11plusP <- list(lc10n1,lc10n2,lc10p1,lc10p2)

# Ages 5 - 10
lc5n1 <- lc5n2  <- lc
lc5n1[vs=="(Intercept)"|vs=="agestrat(4, 10]"]<- 1
lc5n2[vs=="(Intercept)"|vs=="agestrat(4, 10]"|vs=="entero35"|vs=="entero35:agestrat(4, 10]"] <- 1
lc.age5to10P <- list(lc5n1,lc5n2)

# Ages 0 - 4 
lc0n1 <- lc0n2  <- lc
lc0n1[vs=="(Intercept)"|vs=="agestrat(0, 4]"]<- 1
lc0n2[vs=="(Intercept)"|vs=="agestrat(0, 4]"|vs=="entero35"|vs=="entero35:agestrat(0, 4]"] <- 1
lc.age0to4P <- list(lc0n1,lc0n2)


# calculate all of the estimates
p.age11plus <- sapply(lc.age11plusP, lccalc,x=agestrat.body,vcv=agestrat.VC)
p.age5to10  <- sapply(lc.age5to10P, lccalc,x=agestrat.body,vcv=agestrat.VC)
p.age0to4   <- sapply(lc.age0to4P, lccalc,x=agestrat.body,vcv=agestrat.VC)

colnames(p.age11plus) <- colnames(p.age5to10) <- colnames(p.age0to4) <- c("Entero1600<=35cfu","Entero1600>35cfu")

rownames(p.age11plus) <- rownames(p.age5to10) <- rownames(p.age0to4) <- c("P(GI)","P(GI)lb","P(GI)ub")



# --------------------------------------
# Calculate unadjusted cumulative 
# incidence per 1000 overall and by age
# category
# --------------------------------------
ci.all <- tapply(ah$diarrheaci3,ah$entero35,mean)*1000
ci.age0to4 <- tapply(ah$diarrheaci3[ah$agestrat=="(0, 4]"],ah$entero35[ah$agestrat=="(0, 4]"],mean)*1000
ci.age5to10 <- tapply(ah$diarrheaci3[ah$agestrat=="(4, 10]"],ah$entero35[ah$agestrat=="(4, 10]"],mean)*1000
ci.age11plus <- tapply(ah$diarrheaci3[ah$agestrat==">10"],ah$entero35[ah$agestrat==">10"],mean)*1000


# --------------------------------------
# compare unadjusted and adjusted cumulative estimates
# (barely any difference due to adjustment)
# --------------------------------------
round(rbind(p.age11plus[1,]*1000,ci.age11plus))

round(rbind(p.age5to10[1,]*1000,ci.age5to10))

round(rbind(p.age0to4[1,]*1000,ci.age0to4))



# --------------------------------------
# save the objects
# --------------------------------------
save.image("~/dropbox/13beaches/aim1-results/rawoutput/aim1-entero1600-35cfu-regs-body.Rdata")









