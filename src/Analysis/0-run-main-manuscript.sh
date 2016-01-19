#! /bin/sh
##########################################
R CMD BATCH 1-append-wq.R
R CMD BATCH 2a-regress-10day-body-pool.R
R CMD BATCH 2b-regress-10day-body-pool-nonswimmer.R
R CMD BATCH 3a-regress-10day-body-entero-pool.R
R CMD BATCH 3b-regress-10day-body-entero-pool-nonswimmer.R
R CMD BATCH 4a-regress-10day-body-joint-pool.R
R CMD BATCH 4b-regress-10day-body-joint-pool-nonswimmer.R
R CMD BATCH 6-wq-analysis.R
R CMD BATCH 7-illness-analysis.R
R CMD BATCH 8-bargraph-analysis.R
