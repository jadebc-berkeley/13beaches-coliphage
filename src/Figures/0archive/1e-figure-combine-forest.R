##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# Combine forest plots for main analysis
# and negative control analysis
##########################################

library(ggplot2)
library(gridExtra)


#-----------------------------------------
# Main analysis
#-----------------------------------------
rm(list=ls())
load("~/dropbox/coliphage/results/rawoutput/gg_forest.RData")
load("~/dropbox/coliphage/results/rawoutput/gg_forest_joint.RData")

# Get the ggplot grobs
gA <- ggplotGrob(forest)  
gB <- ggplotGrob(forest.joint)

# Arrange the four charts
grid.arrange(gA,gB,nrow=2)

# Combine the plots   
g = rbind(gA, gB, size = "last")

# draw it
grid.newpage()

pdf("~/dropbox/coliphage/results/figures/forestplots.pdf",height=8,width=16)
grid.draw(g)
dev.off()

#-----------------------------------------
# Negative control analysis
#-----------------------------------------
rm(list=ls())
load("~/dropbox/coliphage/results/rawoutput/gg_forest_nc.RData")
load("~/dropbox/coliphage/results/rawoutput/gg_forest_joint_nc.RData")


# Get the ggplot grobs
gA <- ggplotGrob(forest.nc)  
gB <- ggplotGrob(forest.joint.nc)

# Arrange the four charts
grid.arrange(gA,gB,nrow=2)

# Combine the plots   
g = rbind(gA, gB, size = "last")

# draw it
grid.newpage()

pdf("~/dropbox/coliphage/results/figures/forestplots_negcontrol.pdf",height=8,width=16)
grid.draw(g)
dev.off()


