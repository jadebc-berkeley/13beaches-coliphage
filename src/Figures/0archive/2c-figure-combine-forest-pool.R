##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# Combine forest plots for main analysis

# Results pooled by assay and beach
##########################################
library(ggplot2)
library(gridExtra)



#-----------------------------------------
# Main analysis
#-----------------------------------------
rm(list=ls())
load("~/dropbox/coliphage/results/rawoutput/gg_forest_pool.RData")
load("~/dropbox/coliphage/results/rawoutput/gg_forest_joint_pool.RData")

# Get the ggplot grobs
gA <- ggplotGrob(forest)  
gB <- ggplotGrob(forest.joint)

# Arrange the four charts
grid.arrange(gA,gB,nrow=2)

# Combine the plots   
g = rbind(gA, gB, size = "last")

# draw it
grid.newpage()

pdf("~/dropbox/coliphage/results/figures/forestplots_pool.pdf",height=7,width=11)
grid.draw(g)
dev.off()
