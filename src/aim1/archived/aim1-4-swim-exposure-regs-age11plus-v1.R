



# --------------------------------------
# aim1-4-swim-exposure-regs-age11plus.R
# ben arnold (benarnold@berkeley.edu)
#
# description:
# estimate the association between 
# water exposure and the risk of GI illness
# for the 13 beaches study
#
# Ages 11 years and older
#
# version 1 (10 feb 2015)
#
# --------------------------------------

# --------------------------------------
# preamble
# --------------------------------------

rm(list=ls())
library(sandwich)
library(lmtest)


# source the base functions
source("~/dropbox/13beaches/src/aim1/aim1-0-base-functions.R")


# --------------------------------------
# load the analysis dataset
# --------------------------------------

d <- read.csv("~/dropbox/13beaches/data/final/13beaches-analysis.csv")

# convert ID variables from factors to strings
d$hhid <- as.character(d$hhid)
d$indid <- as.character(d$indid)

# --------------------------------------
# subset to observations for analysis
# --------------------------------------

# restrict the dataset to ages 11 years and older
table(d$agestrat)
d <- subset(d,age>10)

# drop individuals at Silver beach Venfest
table(d$venfest)
d <- subset(d,venfest==0|is.na(venfest)==TRUE)

# drop individuals with baseline GI illness
table(d$gibase)
d <- subset(d,gibase=="No")
dim(d)


# --------------------------------------
# subset dataset to variables of interest
# to speed up computations / processing
# --------------------------------------

# included: agef, female, racewhite, gichron, anim_any, gicontactbase, rawfood
# excluded: allergies (not relevant), frequency of beach visits (not measured at all beaches), digging in the sand (not measured at all beaches)

ad <- subset(d,select=c("beach","pointsource","marine","hhid","indid","groundwater","berm","anycontact","bodycontact","headunder","swallwater","gici3","gici10","age","agecat","agestrat","female","racewhite","gichron","anim_any","gicontactbase","rawfood"))

# for exposure variables, replace missing factor level with missing values and recode the factors
ad$anycontact[ad$anycontact==""] <- NA
	ad$anycontact <- factor(ad$anycontact)
ad$bodycontact[ad$bodycontact==""] <- NA
	ad$bodycontact <- factor(ad$bodycontact)
ad$headunder[ad$headunder==""] <- NA
	ad$headunder <- factor(ad$headunder)
ad$swallwater[ad$swallwater==""] <- NA
	ad$swallwater <- factor(ad$swallwater)

# for some covariates, move missing category to the last factor category
levels(ad$racewhite) <- c("Missing","No","Yes")
	ad$racewhite <- factor(ad$racewhite,levels=c("No","Yes","Missing"))
levels(ad$gichron) <- c("Missing","No","Yes")
	ad$gichron <- factor(ad$gichron,levels=c("No","Yes","Missing"))
ad$gicontactbase <- factor(ad$gicontactbase,levels=c("No","Yes","Missing"))

# reorder groundwater factor for convenience w/ reg estimates
ad$groundwater <- factor(ad$groundwater,levels=c("Below median flow","Above median flow",""))

# create a factor variable for fixed effects that includes separate indicators for avalon/doheny
# effect modification conditions
ad$beachfi <- factor(ad$beach,levels=c(levels(ad$beach),"Avalon, Low groundwater","Avalon, High groundwater","Doheny, Berm open","Doheny, Berm closed"))
ad$beachfi[ad$beach=="Avalon" & ad$groundwater=="Above median flow"] <- "Avalon, High groundwater"
ad$beachfi[ad$beach=="Avalon" & ad$groundwater=="Below median flow"] <- "Avalon, Low groundwater"
ad$beachfi[ad$beach=="Doheny" & ad$berm=="Open"]   <- "Doheny, Berm open"
ad$beachfi[ad$beach=="Doheny" & ad$berm=="Closed"] <- "Doheny, Berm closed"
ad$beachfi <- relevel(ad$beachfi,"Boqueron")


# --------------------------------------
# esimate risk of GI illness in 3 days
# associated with ocean exposure
# adjusted estimates only
# --------------------------------------

# adjustment covariates:
# +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood

# for beach-specific estimates, cannot adjust for gichron (sparse data)

# --------------------------------------
# 3 day follow-up
# Body Immersion
# --------------------------------------

# subset to non-missing exposure categories
# to make the robust CI calcs work
ab <- subset(ad,is.na(ad$anycontact)==FALSE)
ab <- subset(ab,is.na(ab$bodycontact)==FALSE)


# Freshwater beaches

# Huntington
hu.3body <- mpreg(gici3~anycontact+bodycontact +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ab,beach="Huntington",vcv=TRUE)

# Silver
si.3body <- mpreg(gici3~anycontact+bodycontact +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ab,beach="Silver",vcv=TRUE)

# Washington Park
wp.3body <- mpreg(gici3~anycontact+bodycontact +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ab,beach="Washington Park",vcv=TRUE)

# West
we.3body <- mpreg(gici3~anycontact+bodycontact +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ab,beach="West",vcv=TRUE)

# Marine beaches

# Avalon
av.3body <- mpreg(gici3~anycontact+bodycontact +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ab,beach="Avalon",vcv=TRUE)

# Boqueron
bo.3body <- mpreg(gici3~anycontact+bodycontact +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ab,beach="Boqueron",vcv=TRUE)

# Doheny
dh.3body <- mpreg(gici3~anycontact+bodycontact +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ab,beach="Doheny",vcv=TRUE)

# Edgewater
ed.3body <- mpreg(gici3~anycontact+bodycontact +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ab,beach="Edgewater",vcv=TRUE)

# Fairhope
fa.3body <- mpreg(gici3~anycontact+bodycontact +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ab,beach="Fairhope",vcv=TRUE)

# Goddard
gd.3body <- mpreg(gici3~anycontact+bodycontact +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ab,beach="Goddard",vcv=TRUE)

# Malibu
ma.3body <- mpreg(gici3~anycontact+bodycontact +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ab,beach="Malibu",vcv=TRUE)

# Mission Bay
mb.3body <- mpreg(gici3~anycontact+bodycontact +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ab,beach="Mission Bay",vcv=TRUE)

# Surfside
su.3body <- mpreg(gici3~anycontact+bodycontact +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ab,beach="Surfside",vcv=TRUE)

### Pooled estimate, Body Immersion
all.3body.fit <- glm(gici3~anycontact+bodycontact +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beachfi,family=poisson(link="log"),data=ab)
all.3body.VC <- cl(ab,fm=all.3body.fit,cluster=ab$hhid)
all.3body <- coeftest(all.3body.fit, all.3body.VC) 

### Pooled estimate (point vs. non-point source conditions), Body Immersion
ps.3body.fit <- glm(gici3~(anycontact+bodycontact)*pointsource +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beachfi,family=poisson(link="log"),data=ab)
ps.3body.VC <- cl(ab,fm=ps.3body.fit,cluster=ab$hhid)
ps.3body <- coeftest(ps.3body.fit, ps.3body.VC) 

### Pooled estimate (fresh vs. marine), Body Immersion
fm.3body.fit <- glm(gici3~(anycontact+bodycontact)*marine +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beachfi,family=poisson(link="log"),data=ab)
fm.3body.VC <- cl(ab,fm=fm.3body.fit,cluster=ab$hhid)
fm.3body <- coeftest(fm.3body.fit, fm.3body.VC) 


# --------------------------------------
# 3 day follow-up
# Head Immersion
# --------------------------------------

# subset to non-missing exposure categories
# to make the robust CI calcs work
ah <- subset(ad,is.na(ad$anycontact)==FALSE)
ah <- subset(ah,is.na(ah$headunder)==FALSE)


# Freshwater beaches

# Huntington
hu.3head <- mpreg(gici3~anycontact+headunder +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Huntington",vcv=TRUE)

# Silver
si.3head <- mpreg(gici3~anycontact+headunder +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Silver",vcv=TRUE)

# Washington Park
wp.3head <- mpreg(gici3~anycontact+headunder +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Washington Park",vcv=TRUE)

# West
we.3head <- mpreg(gici3~anycontact+headunder +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="West",vcv=TRUE)

# Marine beaches

# Avalon
av.3head <- mpreg(gici3~anycontact+headunder +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Avalon",vcv=TRUE)

# Boqueron
bo.3head <- mpreg(gici3~anycontact+headunder +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Boqueron",vcv=TRUE)

# Doheny
dh.3head <- mpreg(gici3~anycontact+headunder +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Doheny",vcv=TRUE)

# Edgewater
ed.3head <- mpreg(gici3~anycontact+headunder +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Edgewater",vcv=TRUE)

# Fairhope
fa.3head <- mpreg(gici3~anycontact+headunder +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Fairhope",vcv=TRUE)

# Goddard
gd.3head <- mpreg(gici3~anycontact+headunder +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Goddard",vcv=TRUE)

# Malibu
ma.3head <- mpreg(gici3~anycontact+headunder +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Malibu",vcv=TRUE)

# Mission Bay
mb.3head <- mpreg(gici3~anycontact+headunder +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Mission Bay",vcv=TRUE)

# Surfside
su.3head <- mpreg(gici3~anycontact+headunder +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=ah,beach="Surfside",vcv=TRUE)

### Pooled estimate, Head Immersion
all.3head.fit <- glm(gici3~anycontact+headunder +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beachfi,family=poisson(link="log"),data=ah)
all.3head.VC <- cl(ah,fm=all.3head.fit,cluster=ah$hhid)
all.3head <- coeftest(all.3head.fit, all.3head.VC) 

### Pooled estimate (point vs. non-point source conditions), Head Immersion
ps.3head.fit <- glm(gici3~(anycontact+headunder)*pointsource +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beachfi,family=poisson(link="log"),data=ah)
ps.3head.VC <- cl(ah,fm=ps.3head.fit,cluster=ah$hhid)
ps.3head <- coeftest(ps.3head.fit, ps.3head.VC) 

### Pooled estimate (fresh vs. marine), Head Immersion
fm.3head.fit <- glm(gici3~(anycontact+headunder)*marine +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beachfi,family=poisson(link="log"),data=ah)
fm.3head.VC <- cl(ah,fm=fm.3head.fit,cluster=ah$hhid)
fm.3head <- coeftest(fm.3head.fit, fm.3head.VC) 

# --------------------------------------
# 3 day follow-up
# Swallowed Water
# --------------------------------------

# subset to non-missing exposure categories
# to make the robust CI calcs work
as <- subset(ad,is.na(ad$anycontact)==FALSE)
as <- subset(as,is.na(as$swallwater)==FALSE)


# Freshwater beaches

# Huntington
hu.3swal <- mpreg(gici3~anycontact+swallwater +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=as,beach="Huntington",vcv=TRUE)

# Silver
si.3swal <- mpreg(gici3~anycontact+swallwater +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=as,beach="Silver",vcv=TRUE)

# Washington Park
wp.3swal <- mpreg(gici3~anycontact+swallwater +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=as,beach="Washington Park",vcv=TRUE)

# West
we.3swal <- mpreg(gici3~anycontact+swallwater +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=as,beach="West",vcv=TRUE)

# Marine beaches

# Avalon
av.3swal <- mpreg(gici3~anycontact+swallwater +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=as,beach="Avalon",vcv=TRUE)

# Boqueron
bo.3swal <- mpreg(gici3~anycontact+swallwater +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=as,beach="Boqueron",vcv=TRUE)

# Doheny
dh.3swal <- mpreg(gici3~anycontact+swallwater +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=as,beach="Doheny",vcv=TRUE)

# Edgewater
ed.3swal <- mpreg(gici3~anycontact+swallwater +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=as,beach="Edgewater",vcv=TRUE)

# Fairhope
fa.3swal <- mpreg(gici3~anycontact+swallwater +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=as,beach="Fairhope",vcv=TRUE)

# Goddard
gd.3swal <- mpreg(gici3~anycontact+swallwater +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=as,beach="Goddard",vcv=TRUE)

# Malibu
ma.3swal <- mpreg(gici3~anycontact+swallwater +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=as,beach="Malibu",vcv=TRUE)

# Mission Bay
mb.3swal <- mpreg(gici3~anycontact+swallwater +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=as,beach="Mission Bay",vcv=TRUE)

# Surfside
su.3swal <- mpreg(gici3~anycontact+swallwater +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood,dat=as,beach="Surfside",vcv=TRUE)


### Pooled estimate, Swallowed Water
all.3swal.fit <- glm(gici3~anycontact+swallwater +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beachfi,family=poisson(link="log"),data=as)
all.3swal.VC <- cl(as,fm=all.3swal.fit,cluster=as$hhid)
all.3swal <- coeftest(all.3swal.fit, all.3swal.VC) 

### Pooled estimate (point vs. non-point source conditions), Swallowed Water
ps.3swal.fit <- glm(gici3~(anycontact+swallwater)*pointsource +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beachfi,family=poisson(link="log"),data=as)
ps.3swal.VC <- cl(as,fm=ps.3swal.fit,cluster=as$hhid)
ps.3swal <- coeftest(ps.3swal.fit, ps.3swal.VC) 

### Pooled estimate (fresh vs. marine), Swallowed Water
fm.3swal.fit <- glm(gici3~(anycontact+swallwater)*marine +agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beachfi,family=poisson(link="log"),data=as)
fm.3swal.VC <- cl(as,fm=fm.3swal.fit,cluster=as$hhid)
fm.3swal <- coeftest(fm.3swal.fit, fm.3swal.VC) 


# --------------------------------------
# Calculate CIRs and 95% CIs
# from the regression model output
# note: calls swestci and iswestci from
# the Aim 1 base functions. These calculate
# CIRs and their CIs using linear combinations
# of the regression coefficients + VC matrix
# --------------------------------------

cirtab <- function(b,h,s,I=FALSE) {
	# b : body immersion results object from mpreg
	# h : head immersion results object from mpreg
	# s : swallowed water results object from mpreg
	# I : logical: interaction between swim cond + env (Avalon, Doheny)
	if(I==FALSE) {
		M <- rbind(
		unlist(swestci(b$fit,b$vcovCL)),  
		unlist(swestci(h$fit,h$vcovCL)),
		unlist(swestci(s$fit,s$vcovCL))	
		)
	} else {
		M <- rbind(
		unlist(iswestci(b$fit,b$vcovCL)),  
		unlist(iswestci(h$fit,h$vcovCL)),
		unlist(iswestci(s$fit,s$vcovCL))	
		)
	}
	rownames(M) <- c("Body Immersion","Head Immersion","Swallowed Water")
	return(M)
}

hu.cir <- cirtab(hu.3body,hu.3head,hu.3swal)
si.cir <- cirtab(si.3body,si.3head,si.3swal)
wp.cir <- cirtab(wp.3body,wp.3head,wp.3swal)
we.cir <- cirtab(we.3body,we.3head,we.3swal)
av.cir <- cirtab(av.3body,av.3head,av.3swal)
bo.cir <- cirtab(bo.3body,bo.3head,bo.3swal)
dh.cir <- cirtab(dh.3body,dh.3head,dh.3swal)
ed.cir <- cirtab(ed.3body,ed.3head,ed.3swal)
fa.cir <- cirtab(fa.3body,fa.3head,fa.3swal)
gd.cir <- cirtab(gd.3body,gd.3head,gd.3swal)
ma.cir <- cirtab(ma.3body,ma.3head,ma.3swal)
mb.cir <- cirtab(mb.3body,mb.3head,mb.3swal)
su.cir <- cirtab(su.3body,su.3head,su.3swal)

all.cir <- cirtab(
			b=list(fit=all.3body,vcovCL=all.3body.VC),
			h=list(fit=all.3head,vcovCL=all.3head.VC),
			s=list(fit=all.3swal,vcovCL=all.3swal.VC)
			)

ps.cir <- cirtab(
			b=list(fit=ps.3body,vcovCL=ps.3body.VC),
			h=list(fit=ps.3head,vcovCL=ps.3head.VC),
			s=list(fit=ps.3swal,vcovCL=ps.3swal.VC), I=TRUE
			)
	ps.3body.test <- switest(ps.3body,ps.3body.VC)
	ps.3head.test <- switest(ps.3head,ps.3head.VC)
	ps.3swal.test <- switest(ps.3swal,ps.3swal.VC)
			
fm.cir <- cirtab(
			b=list(fit=fm.3body,vcovCL=fm.3body.VC),
			h=list(fit=fm.3head,vcovCL=fm.3head.VC),
			s=list(fit=fm.3swal,vcovCL=fm.3swal.VC), I=TRUE
			)
	fm.3body.test <- switest(fm.3body,fm.3body.VC)
	fm.3head.test <- switest(fm.3head,fm.3head.VC)
	fm.3swal.test <- switest(fm.3swal,fm.3swal.VC)
	

# --------------------------------------
# save the objects
# --------------------------------------
save.image("~/dropbox/13beaches/aim1-results/rawoutput/aim1-swim-exposure-regs-age11plus.Rdata")



