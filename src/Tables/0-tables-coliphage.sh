
#!/bin/bash
##########################################


R CMD BATCH 1-table-concen.R
R CMD BATCH 2-table-concen-risk.R
R CMD BATCH 3-table-supp-desc.R
R CMD BATCH 4-table-supp-CIR-10-pool.R

cd other
R CMD BATCH 1-table-supp-CIR-10-pool-unadj.R
R CMD BATCH 2-3-table-CIR-10-exp.R
R CMD BATCH 4-table-supp-CIR-beach.R
R CMD BATCH 5-table-threshold.R
R CMD BATCH 8-table-continuous-negcontrol.R
R CMD BATCH 10a-table-supp-CIR-3-pool.R
R CMD BATCH 10b-table-supp-CIR-3.R
R CMD BATCH 11-table-point-source-em.R
R CMD BATCH 12-table-supp-CIR-10-pool-diarrhea.R
R CMD BATCH 13-table-supp-CIR-10-pool-both.R
R CMD BATCH 15-table-supp-AIC-pool.R 

cd ../supplement

R CMD BATCH table-supp-CIR-10.R 

cd ..

