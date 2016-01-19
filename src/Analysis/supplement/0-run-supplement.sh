#! /bin/sh
##########################################
R CMD BATCH 1a-regress-10day-body.R
R CMD BATCH 1b-regress-10day-body-entero.R
R CMD BATCH 1c-regress-10day-body-joint.R
R CMD BATCH 2a-regress-10day-body-pool-negcontrol.R
R CMD BATCH 2b-regress-10day-body-entero-pool-negcontrol.R
R CMD BATCH 2c-regress-10day-body-joint-pool-negcontrol.R

