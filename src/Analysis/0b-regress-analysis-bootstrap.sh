
#!/bin/bash
##########################################
R CMD BATCH 5a-regress-10day-continuous-body.R
R CMD BATCH 5b-regress-10day-continuous-body-unadj.R
R CMD BATCH 5c-regress-10day-continuous-body-entero.R
R CMD BATCH 5d-regress-10day-continuous-body-beach.R
