#! /bin/sh
##########################################
R CMD BATCH 1-append-wq.R
R CMD BATCH 2a-regress-10day-body.R
R CMD BATCH 2b-regress-10day-body-unadj.R
R CMD BATCH 2c-regress-10day-body-threshold.R
R CMD BATCH 2d-regress-10day-body-beach.R
R CMD BATCH 2e-regress-10day-body-pool.R
R CMD BATCH 2f-regress-10day-body-pool-threshold.R
R CMD BATCH 2g-regress-10day-body-pool-nonswimmer.R
R CMD BATCH 3a-regress-10day-body-entero.R
R CMD BATCH 3b-regress-10day-body-entero-unadj.R
R CMD BATCH 3c-regress-10day-body-entero-pool.R
R CMD BATCH 3d-regress-10day-body-entero-beach.R
R CMD BATCH 3e-3e-regress-10day-body-entero-pool-nonswimmer.R
R CMD BATCH 4a-regress-10day-body-joint.R
R CMD BATCH 4b-regress-10day-body-joint-unadj.R
R CMD BATCH 4c-regress-10day-body-joint-pool.R
R CMD BATCH 4d-regress-10day-body-joint-pool-nonswimmer.R
R CMD BATCH 5c-regress-10day-continuous-body-entero.R 
R CMD BATCH 6a-regress-10day-body-negcontrol.R
R CMD BATCH 6b-regress-10day-body-entero-negcontrol.R
R CMD BATCH 6c-regress-10day-body-joint-negcontrol.R
R CMD BATCH 6d-regress-10day-continuous-body-negcontrol.R
R CMD BATCH 6e-regress-10day-continuous-body-entero-negcontrol.R

cd sensitivity

R CMD BATCH sens-1a-regress-10day-head-beach.R
R CMD BATCH sens-1b-regress-10day-head-beach-unadj.R
R CMD BATCH sens-2a-regress-10day-swall-beach.R
R CMD BATCH sens-2b-regress-10day-swall-beach-unadj.R
R CMD BATCH sens-3a-regress-3day-head-beach.R
R CMD BATCH sens-3b-regress-3day-head-beach-unadj.R
R CMD BATCH sens-4a-regress-3day-swall-beach.R
R CMD BATCH sens-4b-regress-3day-swall-beach-unadj.R
R CMD BATCH sens-5-regress-pointsource.R

cd ../supplement

R CMD BATCH sup-1a-regress-3day-body.R
R CMD BATCH sup-1b-regress-3day-body-beach.R
R CMD BATCH sup-1c-regress-3day-body-beach-unadj.R
R CMD BATCH sup-2a-regress-3day-body-entero.R

cd ..