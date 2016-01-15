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
all.inds$label[all.inds$groupindex=="14FMC051AV071"]="F- Coliphage (EPA 1601)"
all.inds$label[all.inds$groupindex=="14FMC051DO071"]="F- Coliphage (EPA 1601)"

all.inds$label[all.inds$groupindex=="16FMC081AV071"]="F- Coliphage (EPA 1602)"
all.inds$label[all.inds$groupindex=="16FMC081DO071"]="F- Coliphage (EPA 1602)"
all.inds$label[all.inds$groupindex=="16FMC081MA091"]="F- Coliphage (EPA 1602)"

all.inds$label[all.inds$groupindex=="14FPC051AV071"]="F+ Coliphage (EPA 1601)"
all.inds$label[all.inds$groupindex=="14FPC051DO071"]="F+ Coliphage (EPA 1601)"
all.inds$label[all.inds$groupindex=="14FPC051MA091"]="F+ Coliphage (EPA 1601)"

all.inds$label[all.inds$groupindex=="16FPC081AV071"]="F+ Coliphage (EPA 1602)"
all.inds$label[all.inds$groupindex=="16FPC081DO071"]="F+ Coliphage (EPA 1602)"
all.inds$label[all.inds$groupindex=="16FPC081MA091"]="F+ Coliphage (EPA 1602)"

all.inds$label[all.inds$groupindex=="12ENT041AV071"]="EPA 1600 Enterococcus 1"
all.inds$label[all.inds$groupindex=="12ENT041AV081"]="EPA 1600 Enterococcus 2"
all.inds$label[all.inds$groupindex=="15ENT041DO071"]="EPA 1600 Enterococcus 1"
all.inds$label[all.inds$groupindex=="15ENT041DO081"]="EPA 1600 Enterococcus 2"
all.inds$label[all.inds$groupindex=="12ENT041MA091"]="EPA 1600 Enterococcus 1"

all.inds$label[all.inds$analysismethod=="EPA 1601" & 
                 all.inds$indicator=="Fminus coliphage"]="F- Coliphage (EPA 1601)"
all.inds$label[all.inds$analysismethod=="EPA 1601" & 
                 all.inds$indicator=="Fplus coliphage"]="F+ Coliphage (EPA 1601)"
all.inds$label[all.inds$analysismethod=="EPA 1602" & 
                 all.inds$indicator=="Fminus coliphage"]="F- Coliphage (EPA 1602)"
all.inds$label[all.inds$analysismethod=="EPA 1602" & 
                 all.inds$indicator=="Fplus coliphage"]="F+ Coliphage (EPA 1602)"


# --------------------------------------
# Histograms of log10 concentration
# --------------------------------------
all.inds.bx=subset(all.inds,all.inds$label!="EPA 1600 Enterococcus 1")
all.inds.bx=subset(all.inds.bx,all.inds.bx$label!="EPA 1600 Enterococcus 2")
all.inds.bx$label=as.factor(all.inds.bx$label)
all.inds.bx$label=factor(all.inds.bx$label,levels=c("F+ Coliphage (EPA 1602)","F+ Coliphage (EPA 1601)",
      "F- Coliphage (EPA 1602)","F- Coliphage (EPA 1601)"))

fmc1601.plot=ggplot(all.inds.bx[all.inds.bx$label=="F- Coliphage (EPA 1601)",],
  aes(x=log10))+scale_y_continuous(limits=c(0,250))+
  scale_x_continuous(limits=c(-1.5,4),breaks=c(-1,0,1,2,3,4))+
  geom_histogram(binwidth=0.1,color="black",fill="white")+
  facet_grid(beach~.)+theme_complete_bw()+
  xlab("Log10 Coliphage Concentration")+
  ylab("Number of samples")+ggtitle("F- Coliphage (EPA 1601)")+
  geom_vline(x=-.95,linetype="dashed")

fmc1602.plot=ggplot(all.inds.bx[all.inds.bx$label=="F- Coliphage (EPA 1602)",],
  aes(x=log10))+scale_y_continuous(limits=c(0,250))+
  scale_x_continuous(limits=c(-1.5,4),breaks=c(-1,0,1,2,3,4))+
  geom_histogram(binwidth=0.1,color="black",fill="white")+
  facet_grid(beach~.)+theme_complete_bw()+
  xlab("Log10 Coliphage Concentration")+
  ylab("Number of samples")+ggtitle("F- Coliphage (EPA 1602)")+
  geom_vline(x=-.95,linetype="dashed")

fpc1601.plot=ggplot(all.inds.bx[all.inds.bx$label=="F+ Coliphage (EPA 1601)",],
  aes(x=log10))+scale_y_continuous(limits=c(0,250))+
  scale_x_continuous(limits=c(-1.5,4),breaks=c(-1,0,1,2,3,4))+
  geom_histogram(binwidth=0.1,color="black",fill="white")+
  facet_grid(beach~.)+theme_complete_bw()+
  xlab("Log10 Coliphage Concentration")+
  ylab("Number of samples")+ggtitle("F+ Coliphage (EPA 1601)")+
  geom_vline(aes(xintercept=-.95,linetype="Imputed value for observations below detection limit"),
    show_guide=TRUE)+
  scale_linetype_manual(name="",values=c("dashed","dashed"))+theme(legend.position="bottom")

fpc1602.plot=ggplot(all.inds.bx[all.inds.bx$label=="F+ Coliphage (EPA 1602)",],
  aes(x=log10))+scale_y_continuous(limits=c(0,250))+
  scale_x_continuous(limits=c(-1.5,4),breaks=c(-1,0,1,2,3,4))+
  geom_histogram(binwidth=0.1,color="black",fill="white")+
  facet_grid(beach~.)+theme_complete_bw()+
  xlab("Log10 Coliphage Concentration")+
  ylab("Number of samples")+ggtitle("F+ Coliphage (EPA 1602)")+
  geom_vline(x=-.95,linetype="dashed")

pdf("~/Dropbox/Coliphage/Results/Figures/histograms.pdf",onefile=TRUE,
    width=8,height=12)
grid.arrange(fpc1601.plot,arrangeGrob(fmc1601.plot,fmc1602.plot,fpc1602.plot,ncol=1),
             ncol=2,heights=c(1,2))
dev.off()



