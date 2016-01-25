
#!/bin/bash
##########################################


R CMD BATCH 1a-figure-forest-pool-nonswimmer.R
R CMD BATCH 1b-figure-forest.R
R CMD BATCH 2a-figure-forest-pool-negcontrol.R
R CMD BATCH 2b-figure-forest-negcontrol.R

