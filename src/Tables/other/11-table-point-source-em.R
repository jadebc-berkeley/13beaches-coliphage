##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file makes a table with regression output
# For a table with the combined and interaction
# results for body immersion

# Single & joint indicators
##########################################


rm(list=ls())

load("~/Documents/CRG/coliphage/results/rawoutput/regress-10day-ps-em.Rdata")

pval=c(lr.fmc[5][2,],lr.fpc[5][2,],
       lr.fmc.low[5][2,],lr.fpc.low[5][2,],
       lr.fmc.high[5][2,],lr.fpc.high[5][2,])
pvalf=sprintf("%0.3f",pval)

lab=c("Somatic coliphage - all conditions","Male-specific coliphage - all conditions",
      "Somatic coliphage - low risk", "Male-specific coliphage - low risk",
      "Somatic coliphage - high risk", "Male-specific coliphage - high risk")

tab=data.frame(cbind(lab,pvalf))

save(tab,file="~/Documents/CRG/coliphage/Results/Tables/Pointsource-em.RData")




