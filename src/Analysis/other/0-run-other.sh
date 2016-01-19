#! /bin/sh
##########################################
R CMD BATCH 1a-regress-10day-body-unadj.R
R CMD BATCH 1b-regress-10day-body-entero-unadj.R
R CMD BATCH 1c-regress-10day-body-joint-unadj.R
R CMD BATCH 2-regress-10day-head-pool.R
R CMD BATCH 3-regress-10day-swall-pool.R
R CMD BATCH 4a-regress-10day-body-beach.R
R CMD BATCH 4b-regress-10day-body-entero-beach.R
R CMD BATCH 4c-regress-10day-body-joint-beach.R
R CMD BATCH 5a-regress-10day-body-pool-threshold.R
R CMD BATCH 5b-regress-10day-body-threshold.R
R CMD BATCH 7a-regress-10day-body-negcontrol.R
R CMD BATCH 7b-regress-10day-body-entero-negcontrol.R
R CMD BATCH 7c-regress-10day-body-joint-negcontrol.R
R CMD BATCH 8a-regress-10day-continuous-body-negcontrol.R
R CMD BATCH 8b-regress-10day-continuous-body-entero-negcontrol.R
R CMD BATCH 9-regress-10day-continuous-body-entero.R
R CMD BATCH 10a-regress-3day-body.R
R CMD BATCH 10b-regress-3day-body-entero.R
R CMD BATCH 11-regress-pointsource.R
R CMD BATCH 12a-regress-10day-body-pool-diarrhea.R
R CMD BATCH 12b-regress-10day-body-entero-pool-diarrhea.R
R CMD BATCH 12c-regress-10day-body-joint-pool-diarrhea.R
R CMD BATCH 13a-regress-10day-body-pool-both.R
R CMD BATCH 13b-regress-10day-body-entero-pool-both.R
R CMD BATCH 13c-regress-10day-body-joint-pool-both.R
R CMD BATCH 14-regress-10day-continuous-body-joint.R