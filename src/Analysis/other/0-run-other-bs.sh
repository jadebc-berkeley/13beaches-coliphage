#! /bin/sh
##########################################

R CMD BATCH 1d-regress-10day-continuous-body-unadj.R
R CMD BATCH 4d-regress-10day-continuous-body-beach.R
R CMD BATCH 6a-regress-10day-continuous-body-risk.R
R CMD BATCH 6b-regress-10day-continuous-body-risk-nonswimmer.R