##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# Box plots of coliphage
##########################################


rm(list=ls())
library(ggplot2)
library(gridExtra)
library(grid)

setwd("~/Dropbox/Coliphage/")
source("~/dropbox/coliphage/programs/figures/theme_complete_bw.R")


# --------------------------------------
# load the and pre-preprocess the 
# analysis dataset
# (refer to the base functions script
# for details on the pre-processing)
# --------------------------------------
wq=read.csv("~/Dropbox/13beaches/data/final/13beaches-wq-samples.csv",stringsAsFactors=FALSE)

# restrict to 6 beaches with coliphage data
adm.list=c("Avalon","Doheny","Malibu")
other.list=c("Mission Bay","Fairhope","Goddard")

wq.adm=wq[wq$beach %in% adm.list,]
wq.oth=wq[wq$beach %in% other.list,]

# subset to relevant indicators
ind.list=c("14FPC051AV071","14FPC051DO071","14FPC051AV071","16FPC081AV071",
           "14FMC051AV071","16FMC081AV071","12ENT041AV071","12ENT041AV081","14FPC051DO071",
           "16FPC081DO071","14FMC051DO071","16FMC081DO071","15ENT041DO071","15ENT041DO081",
           "16FMC081MA091","14FPC051MA091","16FPC081MA091","12ENT041MA091")

wq.adm=wq.adm[wq.adm$groupindex %in% ind.list,]

all.inds=rbind(wq.adm,wq.oth)

all.inds$label=NULL
all.inds$label[all.inds$groupindex=="14FMC051AV071"]="EPA 1601 F- Coliphage"
all.inds$label[all.inds$groupindex=="14FMC051DO071"]="EPA 1601 F- Coliphage"

all.inds$label[all.inds$groupindex=="16FMC081AV071"]="EPA 1602 F- Coliphage"
all.inds$label[all.inds$groupindex=="16FMC081DO071"]="EPA 1602 F- Coliphage"
all.inds$label[all.inds$groupindex=="16FMC081MA091"]="EPA 1602 F- Coliphage"

all.inds$label[all.inds$groupindex=="14FPC051AV071"]="EPA 1601 F+ Coliphage"
all.inds$label[all.inds$groupindex=="14FPC051DO071"]="EPA 1601 F+ Coliphage"
all.inds$label[all.inds$groupindex=="14FPC051MA091"]="EPA 1601 F+ Coliphage"

all.inds$label[all.inds$groupindex=="16FPC081AV071"]="EPA 1602 F+ Coliphage"
all.inds$label[all.inds$groupindex=="16FPC081DO071"]="EPA 1602 F+ Coliphage"
all.inds$label[all.inds$groupindex=="16FPC081MA091"]="EPA 1602 F+ Coliphage"

all.inds$label[all.inds$groupindex=="12ENT041AV071"]="EPA 1600 Enterococcus 1"
all.inds$label[all.inds$groupindex=="12ENT041AV081"]="EPA 1600 Enterococcus 2"
all.inds$label[all.inds$groupindex=="15ENT041DO071"]="EPA 1600 Enterococcus 1"
all.inds$label[all.inds$groupindex=="15ENT041DO081"]="EPA 1600 Enterococcus 2"
all.inds$label[all.inds$groupindex=="12ENT041MA091"]="EPA 1600 Enterococcus 1"

all.inds$label[all.inds$analysismethod=="EPA 1601" & 
                 all.inds$indicator=="Fminus coliphage"]="EPA 1601 F- Coliphage"
all.inds$label[all.inds$analysismethod=="EPA 1601" & 
                 all.inds$indicator=="Fplus coliphage"]="EPA 1601 F+ Coliphage"
all.inds$label[all.inds$analysismethod=="EPA 1602" & 
                 all.inds$indicator=="Fminus coliphage"]="EPA 1602 F- Coliphage"
all.inds$label[all.inds$analysismethod=="EPA 1602" & 
                 all.inds$indicator=="Fplus coliphage"]="EPA 1602 F+ Coliphage"


# --------------------------------------
# Boxplots of log10 concentration
# --------------------------------------
all.inds.bx=subset(all.inds,all.inds$label!="EPA 1600 Enterococcus 1")
all.inds.bx=subset(all.inds.bx,all.inds.bx$label!="EPA 1600 Enterococcus 2")
all.inds.bx$label=as.factor(all.inds.bx$label)
all.inds.bx$label=factor(all.inds.bx$label,levels=c("EPA 1602 F+ Coliphage","EPA 1601 F+ Coliphage",
      "EPA 1602 F- Coliphage","EPA 1601 F- Coliphage"))

# Drop non-detects
#all.inds.bx$log10[all.inds.bx$log10=="-1"]=NA

pdf(file="~/Dropbox/Coliphage/Results/Figures/boxplot-coliphage.pdf",
    onefile=TRUE,width=10,height=6)
ggplot(all.inds.bx,aes(x=label,y=log10))+geom_boxplot(fill="lightgrey")+
  facet_wrap(~beach)+coord_flip()+
  ylab(bquote(Log[10] ~ "Concentration (PFU/100 ml)"))+xlab("")+
  theme_complete_bw()
dev.off()


