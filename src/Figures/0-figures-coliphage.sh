
#!/bin/bash
##########################################


R CMD BATCH 1a-figure-forest-pool.R
R CMD BATCH 1b-figure-forest.R
R CMD BATCH 2-figure-forest-negcontrol.R
R CMD BATCH 3a-figure-dose-response.R
R CMD BATCH 3b-figure-dose-response-unadj.R
R CMD BATCH 4-figure-scatter.R

