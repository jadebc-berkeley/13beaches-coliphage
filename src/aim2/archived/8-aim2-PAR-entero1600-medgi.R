



# --------------------------------------
# 8-aim2-PAR-entero1600-medgi.R
# ben arnold (benarnold@berkeley.edu)
#
# description:
# estimate the Diarrhea risk attributable
# to swimming in water that exceeds
# 35 CFU/100ml based on the Entero EPA 1600
# indicator.  Estimate the effect for 
# the overall population, E(Y|R,W) - E(Y|R=0,W),
# which is a PIM parameter, and also
# just among those exposed, E(Y|R=1,W) - E(Y|R=0,W),
# which is analogous to the aim 1 analysis
#
# Analyses are conducted for EPA 1600
# Among Swimmers with Body Immmersion
#
# The exposure categories are <=35 and >35 CFU/100ml
#
# version 1 (9 apr 2015)
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
source("~/dropbox/13beaches/src/aim2/0-aim2-base-functions.R")

# --------------------------------------
# load the and pre-preprocess the 
# analysis dataset
# (refer to the base functions script
# for details on the pre-processing)
# --------------------------------------

ad <- preprocess.13beaches("~/dropbox/13beaches/data/final/13beaches-analysis.csv")

# drop individuals with no water quality information
table(ad$nowq)
ad <- subset(ad,nowq==0)
dim(ad)

# subset to non-missing exposure categories
table(ad$bodycontact)
ad <- subset(ad,ad$bodycontact=="Yes")
	dim(ad)
ad <- subset(ad,is.na(ad$entero35)==FALSE)
	dim(ad)


# --------------------------------------
# Calculate the actual Ns for each cell
# and store them for plotting and tables
# --------------------------------------

# --------------------------------------
# Summarize 
# the number of swimmers and number
# of cumulative incident episodes by
# Entero 1600 category
# --------------------------------------
calcNs <- function(outcome,exposurecat) {
  n <- tapply(outcome,exposurecat,sum)
  N <- tapply(outcome,exposurecat,function(x) length(x))
  cbind(n,N)
}

# overall
N.all <- calcNs(ad$medgi,ad$entero35)
# age stratified
N.age0to4   <- calcNs(ad$medgi[ad$agestrat=="(0, 4]"],ad$entero35[ad$agestrat=="(0, 4]"] )
N.age5to10  <- calcNs(ad$medgi[ad$agestrat=="(4, 10]"],ad$entero35[ad$agestrat=="(4, 10]"] )
N.age11plus <- calcNs(ad$medgi[ad$agestrat==">10"],ad$entero35[ad$agestrat==">10"] )
# point source vs. non-point source stratified
N.nps <- calcNs(ad$medgi[ad$pointsource=="No"], ad$entero35[ad$pointsource=="No"])
N.ps  <- calcNs(ad$medgi[ad$pointsource=="Yes"],ad$entero35[ad$pointsource=="Yes"])

# --------------------------------------
# Estimate and bootstrap the parameters
# of interest
# see the aim2 base functions for
# details of the bootAR35cfu function
# --------------------------------------
set.seed(335324)
fmla <- formula(medgi~entero35+marine+pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)
AR35cfu.medgi <- bootAR(fn="AR35cfu",fmla=fmla,data=ad,ID=ad$hhid,strata=ad$beach,iter=1000)

# save results down
rm(ad)
save.image(file="~/dropbox/13beaches/aim2-results/rawoutput/aim2-PARentero1600-medgi.RData")







