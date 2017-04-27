# number of samples analyzed for coliphage
# across 6 beaches

rm(list=ls())

source("~/Documents/CRG/coliphage/13beaches-coliphage/src/figures/theme_complete_bw.R")

# --------------------------------------
# load the and pre-preprocess the 
# analysis dataset
# (refer to the base functions script
# for details on the pre-processing)
# --------------------------------------
wq=read.csv("~/Documents/CRG/coliphage/13beaches-data/temp/beaches-coli-ent-wq.csv",stringsAsFactors=TRUE)

wq$coli=ifelse(!is.na(wq$fpc1601) | !is.na(wq$fpc1602) |
                 !is.na(wq$fmc1601) | !is.na(wq$fmc1602) ,1,0)

wqc=wq[wq$coli==1,]

wqc=wqc[wqc$beachcode!="Avalon-D" & wqc$beachcode!="Doheny-E",]

nrow(wqc)