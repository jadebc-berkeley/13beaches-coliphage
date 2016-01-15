##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# Threshold table

# 10 day gi illness
##########################################

rm(list=ls())
library(foreign)

setwd("~/Dropbox/Coliphage/")

# --------------------------------------
# load the and pre-preprocess the 
# analysis dataset
# (refer to the base functions script
# for details on the pre-processing)
# --------------------------------------
beaches13=read.csv("~/Dropbox/13beaches-fork-coliphage/data/final/13beaches-analysis.csv")

# load base functions
source("Programs/Analysis/0-base-functions.R")

data=preprocess.6beaches(beaches13)

# restrict to 6 beaches with coliphage data
beach.list=c("Avalon","Doheny","Malibu","Mission Bay",
             "Fairhope","Goddard")

all=data[data$beach %in% beach.list,]

# drop individuals with no water quality information
all=subset(all,nowq==0)
# subset to non-missing exposure categories
# to make the robust CI calcs work
all=subset(all,all$bodycontact=="Yes")

all$beach=as.character(all$beach)

fmctab=cbind(sprintf("%0.2f",prop.table(table(all$beach,all$fmc.pres),1)[,2]),
             sprintf("%0.2f",prop.table(table(all$beach,all$fmc_25ref),1)[,2]),
             sprintf("%0.2f",prop.table(table(all$beach,all$fmc_50ref),1)[,2]),
             sprintf("%0.2f",prop.table(table(all$beach,all$fmc_75ref),1)[,2]))

fmctab.out.ref=data.frame(cbind(names(table(all$beach)),fmctab))

fpctab=cbind(sprintf("%0.2f",prop.table(table(all$beach,all$fpc.pres),1)[,2]),
             sprintf("%0.2f",prop.table(table(all$beach,all$fpc_25ref),1)[,2]),
             sprintf("%0.2f",prop.table(table(all$beach,all$fpc_50ref),1)[,2]),
             sprintf("%0.2f",prop.table(table(all$beach,all$fpc_75ref),1)[,2]))

fpctab.out.ref=data.frame(cbind(names(table(all$beach)),fpctab))

fmctab.out.ref$X2=as.character(fmctab.out.ref$X2)
fmctab.out.ref$X3=as.character(fmctab.out.ref$X3)
fmctab.out.ref$X4=as.character(fmctab.out.ref$X4)
fmctab.out.ref$X5=as.character(fmctab.out.ref$X5)

fmctab.out.ref[c(3:4),c(2:5)]="N/A"

save(fmctab.out.ref,fpctab.out.ref,file="~/Dropbox/Coliphage/Data/Temp/table-threshold.RData")
