



# --------------------------------------
# 4-aim2-PAR-swimex-medgi.R
# ben arnold (benarnold@berkeley.edu)
#
# description:
# estimate the risk attributable
# to swimming in water using swimmers and
# non-swimmers for the number of medical 
# visits (including phone, in-person and 
# emergency department)
# due to GI illness
#
# Estimate the effect for 
# the overall population, E(Y|R,W) - E(Y|R=0,W),
# which is the population attributable risk (PAR).
# also divide the PAR by the baseline incidence to
# estimate the population attributable fraction (PAF)
#
#
# The exposure categories are Body Immersion
# vs. no water contact.
#
# version 1 (23 apr 2015)
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



# --------------------------------------
# Calculate the actual Ns for each cell
# and store them for plotting and tables
# --------------------------------------

# --------------------------------------
# Summarize 
# the number of people and number
# of cumulative incident episodes by
# Swim exposure category
# --------------------------------------
N.swimex <- function(data) {
	# data: data frame that includes the variables referenced below
	#       (designed to be subset to a specific swim exposure category:
	#        non-swimmers, body immersion, head immersion, swallowed water)
	calcNs <- function(x) {
  		n <- sum(x)
  		N <- length(x)
  		cbind(n,N)
		}

	N.all       <- calcNs(data$medgi)
	N.age0to4   <- calcNs(data$medgi[data$agestrat=="(0, 4]"])
	N.age5to10  <- calcNs(data$medgi[data$agestrat=="(4, 10]"])
	N.age11plus <- calcNs(data$medgi[data$agestrat==">10"])
	Ns <- rbind(N.all,N.age0to4,N.age5to10,N.age11plus)
	rownames(Ns) <- c("All Ages","Age 0 to 4","Age 5 to 10","Age >10")
	return(Ns)
}

# Total population, Non-swimmers and body immersion swimmers
N.total  <- N.swimex(ad)
N.noswim <- N.swimex(subset(ad,anycontact=="No"))
N.body   <- N.swimex(subset(ad,bodycontact=="Yes"))

# --------------------------------------
# Estimate and bootstrap the parameters
# of interest
# see the aim2 base functions for
# details of the bootAR35cfu function
# --------------------------------------
set.seed(23426)
fmla <- formula(medgi~anycontact+bodycontact+marine+pointsource+agecat+female+racewhite+gichron+anim_any+gicontactbase+rawfood+beach)


PARswimex.medgi <- bootAR(fn="ARswimex",fmla=fmla,data=ad,ID=ad$hhid,strata=ad$beach,iter=1000)

# save results down
rm(ad)
save.image(file="~/dropbox/13beaches/aim2-results/rawoutput/aim2-PARswimex-medgi.RData")







