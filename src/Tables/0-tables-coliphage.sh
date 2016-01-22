
#!/bin/bash
##########################################


R CMD BATCH 1-table-concen.R
R CMD BATCH 2-table-concen-risk.R
R CMD BATCH 3-table-supp-desc.R
R CMD BATCH 4-table-supp-CIR-10-pool.R

cd other
R CMD BATCH table-CIR-3-exp.R
R CMD BATCH table-CIR-10-exp.R
R CMD BATCH table-CIR-10-joint-unadj.R
R CMD BATCH table-CIR-10-single-unadj.R
R CMD BATCH table-supp-CIR-10-beach-unadj.R
R CMD BATCH table-supp-CIR-10-single-unadj.R
R CMD BATCH table-continuous-negcontrol.R
R CMD BATCH table-supp-AIC-pool.R
R CMD BATCH table-supp-CIR-10.R
R CMD BATCH table-supp-CIR-beach.R

cd ..