##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 11/3/15

# Analysis of water quality associations for main text
##########################################


rm(list=ls())
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)
source("~/Documents/CRG/coliphage/13beaches-coliphage/src/figures/theme_complete_bw.R")

# --------------------------------------
# load the and pre-preprocess the 
# analysis dataset
# (refer to the base functions script
# for details on the pre-processing)
# --------------------------------------
wq=read.csv("~/Documents/CRG/coliphage/13beaches-data/temp/beaches-coli-ent-wq.csv",stringsAsFactors=TRUE)

wq$fmc1601[wq$fmc1601_nd=="Below detection"]=0.1
wq$fmc1602[wq$fmc1602_nd=="Below detection"]=0.1
wq$fpc1601[wq$fpc1601_nd=="Below detection"]=0.1
wq$fpc1602[wq$fpc1602_nd=="Below detection"]=0.1
wq$entero[wq$entero_nd=="Below detection"]=0.1

wq$logfmc1601=log(wq$fmc1601,base=10)
wq$logfmc1602=log(wq$fmc1602,base=10)
wq$logfpc1601=log(wq$fpc1601,base=10)
wq$logfpc1602=log(wq$fpc1602,base=10)
wq$logentero=log(wq$entero,base=10)

wq=wq[,c("beach","logentero","logfmc1601","logfmc1602",
         "logfpc1601","logfpc1602","risk",
         "fmc1601","fmc1602","fpc1601","fpc1602")]



# --------------------------------------
# summary statistics for correlations
# --------------------------------------
fmc1601=subset(wq,!is.na(wq$logfmc1601))
fmc1601=fmc1601[!is.na(fmc1601$logentero),]
cor(fmc1601$logentero,fmc1601$logfmc1601,method="spearman")
cor.test(fmc1601$logentero,fmc1601$logfmc1601,method="spearman")

fmc1602=subset(wq,!is.na(wq$logfmc1602))
fmc1602=fmc1602[!is.na(fmc1602$logentero),]
cor(fmc1602$logentero,fmc1602$logfmc1602,method="spearman")
cor.test(fmc1602$logentero,fmc1602$logfmc1602,method="spearman")

fpc1601=subset(wq,!is.na(wq$logfpc1601))
fpc1601=fpc1601[!is.na(fpc1601$logentero),]
cor(fpc1601$logentero,fpc1601$logfpc1601,method="spearman")
cor.test(fpc1601$logentero,fpc1601$logfpc1601,method="spearman")

fpc1602=subset(wq,!is.na(wq$logfpc1602))
fpc1602=fpc1602[!is.na(fpc1602$logentero),]
cor(fpc1602$logentero,fpc1602$logfpc1602,method="spearman")
cor.test(fpc1602$logentero,fpc1602$logfpc1602,method="spearman")

# --------------------------------------
# compare presence absence of entero and coli
# --------------------------------------
wq$ent.pres=ifelse(wq$logentero==-1,0,1)
wq$ent.pres[is.na(wq$logentero)]=NA
wq$ent35=ifelse(wq$logentero>log(35, base=10),1,0)
wq$ent35[is.na(wq$logentero)]=NA
wq$fmc1601.pres=ifelse(wq$logfmc1601==-1,0,1)
wq$fmc1601.pres[is.na(wq$logfmc1601)]=NA
wq$fmc1602.pres=ifelse(wq$logfmc1602==-1,0,1)
wq$fmc1602.pres[is.na(wq$logfmc1602)]=NA
wq$fpc1601.pres=ifelse(wq$logfpc1601==-1,0,1)
wq$fpc1601.pres[is.na(wq$logfpc1601)]=NA
wq$fpc1602.pres=ifelse(wq$logfpc1602==-1,0,1)
wq$fpc1602.pres[is.na(wq$logfpc1602)]=NA

prop.table(table(wq$fmc1601.pres[wq$ent.pres==1]))
prop.table(table(wq$fmc1602.pres[wq$ent.pres==1]))
prop.table(table(wq$fpc1601.pres[wq$ent.pres==1]))
prop.table(table(wq$fpc1602.pres[wq$ent.pres==1]))

prop.table(table(wq$fmc1601.pres[wq$ent35==1]))
prop.table(table(wq$fmc1602.pres[wq$ent35==1]))
prop.table(table(wq$fpc1601.pres[wq$ent35==1]))
prop.table(table(wq$fpc1602.pres[wq$ent35==1]))

prop.table(table(wq$fmc1601.pres[wq$ent35==0]))
prop.table(table(wq$fmc1602.pres[wq$ent35==0]))
prop.table(table(wq$fpc1601.pres[wq$ent35==0]))
prop.table(table(wq$fpc1602.pres[wq$ent35==0]))

wq$fmc.pres=NA
wq$fmc.pres[wq$fmc1601.pres==1 | wq$fmc1602.pres==1]=1 
wq$fmc.pres[wq$fmc1601.pres==0 & wq$fmc1602.pres==0]=0

wq$fpc.pres=NA
wq$fpc.pres[wq$fpc1601.pres==1 | wq$fpc1602.pres==1]=1 
wq$fpc.pres[wq$fpc1601.pres==0 & wq$fpc1602.pres==0]=0

prop.table(table(wq$fmc.pres[wq$ent35==1]))
prop.table(table(wq$fpc.pres[wq$ent35==1]))
prop.table(table(wq$fmc.pres[wq$ent35==0]))
prop.table(table(wq$fpc.pres[wq$ent35==0]))


