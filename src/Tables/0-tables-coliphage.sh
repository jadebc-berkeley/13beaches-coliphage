
#!/bin/bash
##########################################


R CMD BATCH table-concen.R
R CMD BATCH table-continuous-negcontrol.R
R CMD BATCH table-supp-AIC-pool.R
R CMD BATCH table-supp-CIR-10-pool.R
R CMD BATCH table-supp-CIR-10.R
R CMD BATCH table-supp-CIR-beach.R
R CMD BATCH table-supp-desc.R

cd unadjusted
R CMD BATCH table-CIR-10-joint-unadj.R
R CMD BATCH table-CIR-10-single-unadj.R
R CMD BATCH table-supp-CIR-10-beach-unadj.R
R CMD BATCH table-supp-CIR-10-single-unadj.R

cd ../sensitivity
R CMD BATCH table-CIR-3-exp.R
R CMD BATCH table-CIR-10-exp.R

cd ..