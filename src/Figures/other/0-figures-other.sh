
#!/bin/bash
##########################################


R CMD BATCH 1-figure-dose-response-unadj.R 
R CMD BATCH 5a-figure-forest-pool-threshold.R 
R CMD BATCH 5b-figure-forest-threshold.R
R CMD BATCH 6-figure-dose-response-risk.R
R CMD BATCH 13-figure-forest-pool-both.R
R CMD BATCH 16-figure-scatter.R

