#######################################
# Graphs for coliphage paper
#######################################
rm(list=ls())
library(foreign)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

avalon=read.dta("/Users/jadebc/Dropbox/Coliphage/Data/Untouched/avalon3.dta",convert.factors=FALSE)
doheny=read.dta("/Users/jadebc/Dropbox/Coliphage/Data/Untouched/doheny3.dta",convert.factors=FALSE)
malibu=read.dta("/Users/jadebc/Dropbox/Coliphage/Data/Untouched/malibu3.dta",convert.factors=FALSE)

a.cols=data.frame(id=avalon$psid,avalon$sspda14FMC0511,avalon$sspda16FMC0811,
                  avalon$sspda14FPC0511,avalon$sspda16FPC0811)

a.cols.dl=a.cols
a.cols.dl$avalon.sspda14FMC0511[a.cols.dl$avalon.sspda14FMC0511<0]=NA
a.cols.dl$avalon.sspda14FPC0511[a.cols.dl$avalon.sspda14FPC0511<0]=NA
a.cols.dl$avalon.sspda16FPC0811[a.cols.dl$avalon.sspda16FPC0811<0]=NA
a.cols.dl$avalon.sspda16FMC0811[a.cols.dl$avalon.sspda16FMC0811<0]=NA

d.cols=data.frame(doheny$sspda14FMC0511,doheny$sspda16FMC0811,
                  doheny$sspda14FPC0511,doheny$sspda16FPC0811)
d.cols.dl=d.cols
d.cols.dl$doheny.sspda14FMC0511[d.cols.dl$doheny.sspda14FMC0511<0]=NA
d.cols.dl$doheny.sspda14FPC0511[d.cols.dl$doheny.sspda14FPC0511<0]=NA
d.cols.dl$doheny.sspda16FPC0811[d.cols.dl$doheny.sspda16FPC0811<0]=NA
d.cols.dl$doheny.sspda16FMC0811[d.cols.dl$doheny.sspda16FMC0811<0]=NA

m.cols=data.frame(malibu$sspda16FMC0811,malibu$sspda14FPC0511)

#######################################
# Boxplots of log10 concentration
#######################################
avalon.ind=read.dta("/Users/jadebc/Dropbox/Coliphage/Data/Untouched/avalon_inddata.dta")
doheny.ind=read.dta("/Users/jadebc/Dropbox/Coliphage/Data/Untouched/doheny_inddata.dta")
malibu.ind=read.dta("/Users/jadebc/Dropbox/Coliphage/Data/Untouched/malibu_inddata.dta")

a.inds=subset(avalon.ind,avalon.ind$groupindex=="14FPC051AV071" | 
                avalon.ind$groupindex=="16FPC081AV071" |
                avalon.ind$groupindex=="14FMC051AV071" |
                avalon.ind$groupindex=="16FMC081AV071"| 
                avalon.ind$groupindex=="12ENT041AV071"|
                avalon.ind$groupindex=="12ENT041AV081")

d.inds=subset(doheny.ind,doheny.ind$groupindex=="14FPC051DO071" | 
                doheny.ind$groupindex=="16FPC081DO071" |
                doheny.ind$groupindex=="14FMC051DO071" |
                doheny.ind$groupindex=="16FMC081DO071" |
                doheny.ind$groupindex=="15ENT041DO071"|
                doheny.ind$groupindex=="15ENT041DO081")

m.inds=subset(malibu.ind,malibu.ind$groupindex=="16FMC081MA091" | 
                malibu.ind$groupindex=="14FPC051MA091" |
                malibu.ind$groupindex=="16FPC081MA091"|
                malibu.ind$groupindex=="12ENT041MA091")


all.inds=rbind(a.inds,d.inds,m.inds)
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


all.inds.bx=subset(all.inds,all.inds$label!="EPA 1600 Enterococcus 1")
all.inds.bx=subset(all.inds.bx,all.inds.bx$label!="EPA 1600 Enterococcus 2")
all.inds.bx$label=as.factor(all.inds.bx$label)
all.inds.bx$label=factor(all.inds.bx$label,levels=c("EPA 1602 F+ Coliphage","EPA 1601 F+ Coliphage",
    "EPA 1602 F- Coliphage","EPA 1601 F- Coliphage"))

pdf(file="/Users/jadebc/Dropbox/Coliphage/Results/Figures/boxplot-coliphage.pdf",
    onefile=TRUE,width=10,height=6)
ggplot(all.inds.bx,aes(x=label,y=log10))+geom_boxplot(fill="lightgrey")+coord_flip()+
  facet_wrap(~beach)+
  ylab(bquote(Log[10] ~ "Concentration (PFU/100 ml)"))+xlab("")+theme_bw()+
  theme(panel.grid.major = element_line(colour = "lightgrey",size=0.25),
        panel.border = element_rect(colour="black"),
        strip.background = element_rect(fill="grey"))
dev.off()

#######################################
# Boxplots: coliphage dist if the enterococcus criteria met /failed
#######################################
#---------------------------------------
#Reorganizing data for plotting
#---------------------------------------
a.epa=data.frame(id=avalon$psid,avalon$epacrit_12ENT0411)
a.epa.long=melt(a.epa,id.vars="id")
a.cols.long=melt(a.cols,id.vars="id")
a.all=merge(a.epa.long,a.cols.long,by="id")
colnames(a.all)=c("id","variable.x","epacrit","indicator","concen")
a.all=subset(a.all,!is.na(a.all$concen))
a.all$indicator=as.factor(a.all$indicator)
levels(a.all$indicator)=c("F- coliphage (EPA 1601)","F+ coliphage (EPA 1601)",
      "F")
a.all$epacrit=as.factor(a.all$epacrit)
levels(a.all$epacrit)=c("Meets enterococcus criteria","Does not meet criteria")
a.all=a.all[,c(1,3,4,5)]

ggplot(a.all, aes(x=indicator, y=concen, fill=epacrit)) + geom_boxplot()+
  scale_fill_manual(name = "", values = c("white", "grey"), 
  labels = c("0" = "Foo", "1" = "Bar")) +
  coord_flip()+ylab("Log10 coliphage concentration per 100 ml")+
  xlab("")+theme_bw()

#######################################
# Probabiliy of illness plots - code from Ben
#######################################

#-------------------------------------------
# plotting function to make probability of
# illness plots (dose-reponse + histogram)
#-------------------------------------------
pillplot <- function(pdata,ptitle,ytics=seq(0,0.2,by=0.02)) {
  # pdata  : plot data subset to the indicator / conditions of interest
  # ptitle : Plot title indicator information
  # ytics  : customize y-axis tick marks (if necessary)
  
  orest <- paste("aOR = ",sprintf("%1.2f",pdata$or[1])," (",sprintf("%1.2f",pdata$orlb[1]),", ",sprintf("%1.2f",pdata$orub[1]),")",sep="")
  
  # X axis ticks and labels
  xtics <- c(0.1,1,10,100)
  xticlabs <- c("0.1","1","10","100")
  xminor <- log(c(sapply(xtics, function(x) seq(0, x, x/10))), 10)
  
  # P(illness) by indicator concentration
  op <- par(mar=c(2,5,8,2)+0.1)
  plot(pdata$lev,pdata$pr,type="n",
       log="x",xaxt="n",xlim=c(0.02,max(xtics)),xlab="",
       yaxt="n",ylim=range(ytics),ylab="Probability of GI Illness",
       las=1
  )
  abline(v=(xtics), col="gray90")
  abline(v=(10^xminor), col="gray90")
  abline(h=(ytics), col="gray90")
  axis(1,at=xtics,labels=xticlabs)
  axis(1,at=10^xminor,labels=FALSE)
  axis(2,at=ytics,las=1)
  
  lines(pdata$lev,pdata$pr,lwd=1.25)
  lines(pdata$lev,pdata$prlb,lty="dashed")
  lines(pdata$lev,pdata$prub,lty="dashed")
  
  mtext(paste(ptitle,orest,sep=""),side=3,line=0,adj=0)
  
  par(op)
  
}

pillhist <- function(pdata,xlab="Concentration MPN/100 ml") {
  # pdata : plot data
  # xlab  : X-axis label (for different units)
  
  op <- par(mar=c(5,5,0,2)+0.1)
  histytics <- seq(0,200,by=25)
  col <- brewer.pal(3,"Set1")[2]
  
  xtics <- c(0.1,1,10,100)
  xticlabs <- c("0.1","1","10","100")
  xminor <- log(c(sapply(xtics, function(x) seq(0, x, x/10))), 10)
  
  plot(pdata$lev,pdata$n,type="n",
       log="x",xaxt="n",xlim=c(0.02,max(xtics)),xlab=xlab,
       ylim=range(histytics),ylab="N Individuals\nExposed",
       bty="l",
       las=1
  )
  abline(v=(xtics), col="gray90")
  abline(v=(10^xminor), col="gray90")
  abline(h=(histytics), col="gray90")
  axis(1,at=xtics,labels=xticlabs)
  axis(1,at=10^xminor,labels=FALSE)
  
  lines(pdata$lev,pdata$n,type="h",lwd=2,col=col)
  
  par(op)
  
}

#-------------------------------------------
# Avalon figures
#-------------------------------------------
av <- read.dta("~/Dropbox/Coliphage/Data/Temp/avalon-coliphage-pillness.dta")

#-------------------------------------------
# F+ EPA 1601
#-------------------------------------------
fpc1601 <- subset(av,indicator=="sspda14FPC0511" & cond=="Combined" & !is.na(pr))
fpc1601gwhgh <- subset(av,indicator=="sspda14FPC0511" & cond=="Groundwater above median" & !is.na(pr))
fpc1601gwlow <- subset(av,indicator=="sspda14FPC0511" & cond=="Groundwater below median" & !is.na(pr))

pdf("/Users/jadebc/Dropbox/Coliphage/Results/Figures/avalon-Fplus-coliphage-epa1601.pdf",width=6,height=6)
lo <- layout(mat=matrix(c(1,2),nrow=2,ncol=1),heights=c(1,0.5))
pillplot(pdata=fpc1601,ptitle="Avalon beach\nEPA 1601 F+ Coliphage\n")
pillhist(pdata=fpc1601)
dev.off()

#-------------------------------------------
# F- EPA 1601
#-------------------------------------------
fmc1601 <- subset(av,indicator=="sspda14FMC0511" & cond=="Combined" & !is.na(pr))
fmc1601gwhgh <- subset(av,indicator=="sspda14FMC0511" & cond=="Groundwater above median" & !is.na(pr))
fmc1601gwlow <- subset(av,indicator=="sspda14FMC0511" & cond=="Groundwater below median" & !is.na(pr))

pdf("/Users/jadebc/Dropbox/Coliphage/Results/Figures/avalon-Fminus-coliphage-epa1601.pdf",width=6,height=6)
lo <- layout(mat=matrix(1:2,nrow=2,ncol=1),heights=c(1,0.5))
pillplot(pdata=fmc1601,ptitle="Avalon beach\nEPA 1601 F- Coliphage\n")
pillhist(pdata=fmc1601)
dev.off()

#-------------------------------------------
# F+ EPA 1602
# (data are so sparse - only present combined)
#-------------------------------------------
fpc1602 <- subset(av,indicator=="sspda16FPC0811" & cond=="Combined" & !is.na(pr))
fpc1602gwhgh <- subset(av,indicator=="sspda16FPC0811" & cond=="Groundwater above median" & !is.na(pr))
fpc1602gwlow <- subset(av,indicator=="sspda16FPC0811" & cond=="Groundwater below median" & !is.na(pr))

pdf("/Users/jadebc/Dropbox/Coliphage/Results/Figures/avalon-Fplus-coliphage-epa1602.pdf",width=6,height=6)
lo <- layout(mat=matrix(1:2,nrow=2,ncol=1),heights=c(1,0.5))
pillplot(pdata=fpc1602,ptitle="Avalon beach\nEPA 1602 F+ Coliphage\n",ytics=seq(0,0.35,by=0.05))
pillhist(pdata=fpc1602,xlab="Concentration PFU/100 ml")
dev.off()

#-------------------------------------------
# F- EPA 1602
#-------------------------------------------
fmc1602 <- subset(av,indicator=="sspda16FMC0811" & cond=="Combined" & !is.na(pr))
fmc1602gwhgh <- subset(av,indicator=="sspda16FMC0811" & cond=="Groundwater above median" & !is.na(pr))
fmc1602gwlow <- subset(av,indicator=="sspda16FMC0811" & cond=="Groundwater below median" & !is.na(pr))

pdf("/Users/jadebc/Dropbox/Coliphage/Results/Figures/avalon-Fminus-coliphage-epa1602.pdf",width=6,height=6)
lo <- layout(mat=matrix(1:2,nrow=2,ncol=1),heights=c(1,0.5))
pillplot(pdata=fmc1602,ptitle="Avalon beach\nEPA 1602 F- Coliphage\n")
pillhist(pdata=fmc1602,xlab="Concentration PFU/100 ml")
dev.off()

#-------------------------------------------
# Doheny figures
#-------------------------------------------
doh <- read.dta("/Users/jadebc/Dropbox/Coliphage/Data/Temp/doheny-coliphage-pillness.dta")

#-------------------------------------------
# F+ EPA 1601
#-------------------------------------------
fpc1601 <- subset(doh,indicator=="sspda14FPC0511" & cond=="Combined" & !is.na(pr))
fpc1601bermop <- subset(doh,indicator=="sspda14FPC0511" & cond=="Berm open" & !is.na(pr))
fpc1601bermcl <- subset(doh,indicator=="sspda14FPC0511" & cond=="Berm closed" & !is.na(pr))

pdf("/Users/jadebc/Dropbox/Coliphage/Results/Figures/doheny-Fplus-coliphage-epa1601.pdf",width=6,height=6)
lo <- layout(mat=matrix(1:2,nrow=2,ncol=1),heights=c(1,0.5))
pillplot(pdata=fpc1601,ptitle="Doheny beach\nEPA 1601 F+ Coliphage\n")
pillhist(pdata=fpc1601)
dev.off()


#-------------------------------------------
# F- EPA 1601
#-------------------------------------------
fmc1601 <- subset(doh,indicator=="sspda14FMC0511" & cond=="Combined" & !is.na(pr))
fmc1601bermop <- subset(doh,indicator=="sspda14FMC0511" & cond=="Berm open" & !is.na(pr))
fmc1601bermcl <- subset(doh,indicator=="sspda14FMC0511" & cond=="Berm closed" & !is.na(pr))

pdf("/Users/jadebc/Dropbox/Coliphage/Results/Figures/doheny-Fminus-coliphage-epa1601.pdf",width=6,height=6)
lo <- layout(mat=matrix(1:2,nrow=2,ncol=1),heights=c(1,0.5))
pillplot(pdata=fmc1601,ptitle="Doheny beach\nEPA 1601 F- Coliphage\n")
pillhist(pdata=fmc1601)
dev.off()

#-------------------------------------------
# F+ EPA 1602
#-------------------------------------------
fpc1602 <- subset(doh,indicator=="sspda16FPC0811" & cond=="Combined" & !is.na(pr))
fpc1602bermop <- subset(doh,indicator=="sspda16FPC0811" & cond=="Berm open" & !is.na(pr))
fpc1602bermcl <- subset(doh,indicator=="sspda16FPC0811" & cond=="Berm closed" & !is.na(pr))

pdf("/Users/jadebc/Dropbox/Coliphage/Results/Figures/doheny-Fplus-coliphage-epa1602.pdf",width=6,height=6)
lo <- layout(mat=matrix(1:2,nrow=2,ncol=1),heights=c(1,0.5))
pillplot(pdata=fpc1602,ptitle="Doheny beach\nEPA 1602 F+ Coliphage\n",ytics=seq(0,0.5,by=0.05))
pillhist(pdata=fpc1602,xlab="Concentration PFU/100 ml")
dev.off()

#-------------------------------------------
# F- EPA 1602
#-------------------------------------------
fmc1602 <- subset(doh,indicator=="sspda16FMC0811" & cond=="Combined" & !is.na(pr))
fmc1602bermop <- subset(doh,indicator=="sspda16FMC0811" & cond=="Berm open" & !is.na(pr))
fmc1602bermcl <- subset(doh,indicator=="sspda16FMC0811" & cond=="Berm closed" & !is.na(pr))

pdf("/Users/jadebc/Dropbox/Coliphage/Results/Figures/doheny-Fminus-coliphage-epa1602.pdf",width=6,height=6)
lo <- layout(mat=matrix(1:2,nrow=2,ncol=1),heights=c(1,0.5))
pillplot(pdata=fmc1602,ptitle="Doheny beach\nEPA 1602 F- Coliphage\n")
pillhist(pdata=fmc1602,xlab="Concentration PFU/100 ml")
dev.off()

#-------------------------------------------
# Malibu figures
#-------------------------------------------
ma <- read.dta("/Users/jadebc/Dropbox/Coliphage/Data/Temp/malibu-coliphage-pillness.dta")


#-------------------------------------------
# F+ EPA 1601
#-------------------------------------------
fpc1601 <- subset(ma,indicator=="sspda14FPC0511" & cond=="Combined" & !is.na(pr))
fpc1601bermop <- subset(ma,indicator=="sspda14FPC0511" & cond=="Berm open" & !is.na(pr))
fpc1601bermcl <- subset(ma,indicator=="sspda14FPC0511" & cond=="Berm closed" & !is.na(pr))

pdf("/Users/jadebc/Dropbox/Coliphage/Results/Figures/malibu-Fplus-coliphage-epa1601.pdf",width=6,height=6)
lo <- layout(mat=matrix(1:2,nrow=2,ncol=1),heights=c(1,0.5))
pillplot(pdata=fpc1601,ptitle="Malibu beach\nEPA 1601 F+ Coliphage\n")
pillhist(pdata=fpc1601)
dev.off()


#-------------------------------------------
# F- EPA 1601 NOT MEASURED IN MALIBU
#-------------------------------------------

#-------------------------------------------
# F+ EPA 1602 NOT MEASURED IN MALIBU
#-------------------------------------------

#-------------------------------------------
# F- EPA 1602
#-------------------------------------------
fmc1602 <- subset(ma,indicator=="sspda16FMC0811" & cond=="Combined" & !is.na(pr))
fmc1602bermop <- subset(ma,indicator=="sspda16FMC0811" & cond=="Berm open" & !is.na(pr))
fmc1602bermcl <- subset(ma,indicator=="sspda16FMC0811" & cond=="Berm closed" & !is.na(pr))

pdf("/Users/jadebc/Dropbox/Coliphage/Results/Figures/malibu-Fminus-coliphage-epa1602.pdf",width=6,height=6)
lo <- layout(mat=matrix(1:2,nrow=2,ncol=1),heights=c(1,0.5))
pillplot(pdata=fmc1602,ptitle="Malibu beach\nEPA 1602 F- Coliphage\n")
pillhist(pdata=fmc1602,xlab="Concentration PFU/100 ml")
dev.off()






