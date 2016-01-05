



# --------------------------------------
# aim1-regs-QPCR.R
# ben arnold (benarnold@berkeley.edu)
#
# description:
# estimate the association between water
# quality indicator concentrations and
# the risk of GI illness among swimmers
# for the 13 beaches study
#
# Analyses are conducted for EPA QPCR
# for now, among swimmers with head immersion
#
# version 1 (28 feb 2014)
#
# --------------------------------------

# --------------------------------------
# preamble
# --------------------------------------

rm(list=ls())
library(foreign)
library(sandwich)
library(lmtest)
library(ggplot2)
library(scales)
library(RColorBrewer)

# --------------------------------------
# Robust clustered SE function
# http://people.su.se/~ma/mcluster.R
# R (www.r-project.org) codes for computing multi-way clustered-standard errors
# Mahmood Arai, Jan 21, 2008. 
# See: Thompson (2006), Cameron, Gelbach and Miller (2006) and Petersen (2006).
#
# slightly modified to have it return the vcovCL object
# rather than the updated fit (since need the VC matrix)
# --------------------------------------
cl   <- function(dat,fm, cluster){
	# dat: data used to fit the model
	# fm : model fit (object)
	# cluster : vector of cluster IDs
	require(sandwich, quietly = TRUE)
	require(lmtest, quietly = TRUE)
	M <- length(unique(cluster))
	N <- length(cluster)
	K <- fm$rank
	dfc <- (M/(M-1))*((N-1)/(N-K))
	uj  <- apply(estfun(fm),2, function(x) tapply(x, cluster, sum));
	vcovCL <- dfc*sandwich(fm, meat=crossprod(uj)/N)
	return(vcovCL)
}

# --------------------------------------
# Convenience Function to run log-binomial models
# and obtain robust SEs (clusterd on hhid)
# for the beaches data
# --------------------------------------

lbreg <- function(formula,dat,beach,vcv=FALSE) {
	# log-binomial regression formula
	# dataset used to fit the model
	# beach name (character)	
	fit <- glm(formula,family=binomial(link="log"),data=dat[dat$beach==beach,])
	vcovCL <- cl(dat[dat$beach==beach,],fm=fit,cluster=dat$hhid[dat$beach==beach])
	rfit <- coeftest(fit, vcovCL) 
	if(vcv==FALSE) {
		return(rfit)
	} else {
		list(fit=rfit,vcovCL=vcovCL)
	}
}



# --------------------------------------
# convenience functions to obtain
# point estimates and SEs from
# model fits for plotting
# --------------------------------------
# regs w/o interactions
estci <- function(fo) {
	# fo : fit object returned from lbreg
	est <- exp(fo[2,1])
	lb <- exp(fo[2,1]-1.96*fo[2,2])
	ub <- exp(fo[2,1]+1.96*fo[2,2])
	list(est=est,lb=lb,ub=ub)
}

# regs w/ interactions (Avalon, Doheny)
iestci <- function(fo,vcv) {
	# fo : fit object returned from lbreg
	# vcv: variance covariance matrix from lbreg
	est0 <- exp(fo[2,1])
	lb0  <- exp(fo[2,1]-1.96*fo[2,2])
	ub0  <- exp(fo[2,1]+1.96*fo[2,2])
	
	est1 <- exp(fo[2,1]+fo[4,1])
	se1  <- sqrt(vcv[2,2]+vcv[4,4]+2*vcv[2,4])
	lb1  <- exp(fo[2,1]+fo[4,1] - 1.96*se1)
	ub1  <- exp(fo[2,1]+fo[4,1] + 1.96*se1)
	
	res <- matrix(c(est0,lb0,ub0,est1,lb1,ub1),nrow=2,ncol=3,byrow=T)
	rownames(res) <- c("effmod=0","effmod=1")
	colnames(res) <- c("est","lb","ub")
	
	return(res)
}



# --------------------------------------
# load the analysis dataset
# --------------------------------------

d <- read.dta("~/dropbox/13beaches/data/final/13beaches-analysis.dta")

# --------------------------------------
# subset to observations for analysis
# --------------------------------------

# drop non-swimmers
table(d$swimmer)
d <- subset(d,swimmer=="Yes")

# restrict to swimmwers with head immersion or more
d <- subset(d,headunder=="Yes"|swallwater=="Yes")

# drop individuals at Silver beach Venfest
table(d$venfest)
d <- subset(d,venfest==0|is.na(venfest)==TRUE)

# drop individuals with no water quality information
table(d$nowq)
d <- subset(d,nowq==0)

# for EPA QPCR analyses, there are 903 indivs at Doheny
# with no water quality information.  Need to drop them.
table(d$beach,is.na(d$avgdyenteropcr))
d <- subset(d,!is.na(avgdyenteropcr))

# drop individuals with baseline GI illness
table(d$gibase)
d <- subset(d,gibase=="No")
dim(d)

# create a dataset for children <=10
dch <- subset(d,d$age<=10)
dim(dch)

# --------------------------------------
# esimate risk of GI illness associated
# with a 1-log10 increase in Entero
# EPA 1600
# estimates for all ages and children<=10
# ordered by fresh and marine beaches
# --------------------------------------


# --------------------------------------
# Freshwater beaches
# all ages
# --------------------------------------

# Huntington
hufit <- lbreg(gici10~avgdyenteropcr,dat=d,beach="Huntington")

# Silver
sifit <- lbreg(gici10~avgdyenteropcr,dat=d,beach="Silver")

# Washington Park
wpfit <- lbreg(gici10~avgdyenteropcr,dat=d,beach="Washington Park")

# West
wefit <- lbreg(gici10~avgdyenteropcr,dat=d,beach="West")


# --------------------------------------
# Marine beaches
# all ages
# --------------------------------------

# Avalon
avfit <- lbreg(gici10~avgdyenteropcr*groundwater,dat=d,beach="Avalon",vcv=TRUE)

# Boqueron
bofit <- lbreg(gici10~avgdyenteropcr,dat=d,beach="Boqueron")

# Doheny
dhfit <- lbreg(gici10~avgdyenteropcr*berm,dat=d,beach="Doheny",vcv=TRUE)

# Edgewater
edfit <- lbreg(gici10~avgdyenteropcr,dat=d,beach="Edgewater")

# Fairhope
fafit <- lbreg(gici10~avgdyenteropcr,dat=d,beach="Fairhope")

# Goddard
gdfit <- lbreg(gici10~avgdyenteropcr,dat=d,beach="Goddard")

# Malibu
mafit <- lbreg(gici10~avgdyenteropcr,dat=d,beach="Malibu")

# Surfside
sufit <- lbreg(gici10~avgdyenteropcr,dat=d,beach="Surfside")


# --------------------------------------
# Freshwater beaches
# children <= 10
# --------------------------------------

# Huntington
hucfit <- lbreg(gici10~avgdyenteropcr,dat=dch,beach="Huntington")

# Silver
sicfit <- lbreg(gici10~avgdyenteropcr,dat=dch,beach="Silver")

# Washington Park
wpcfit <- lbreg(gici10~avgdyenteropcr,dat=dch,beach="Washington Park")

# West
wecfit <- lbreg(gici10~avgdyenteropcr,dat=dch,beach="West")


# --------------------------------------
# Marine beaches
# children <=10
# --------------------------------------

# Avalon
avcfit <- lbreg(gici10~avgdyenteropcr*groundwater,dat=dch,beach="Avalon",vcv=TRUE)

# Boqueron
bocfit <- lbreg(gici10~avgdyenteropcr,dat=dch,beach="Boqueron")

# Doheny
dhcfit <- lbreg(gici10~avgdyenteropcr*berm,dat=dch,beach="Doheny",vcv=TRUE)

# Edgewater
edcfit <- lbreg(gici10~avgdyenteropcr,dat=dch,beach="Edgewater")

# Fairhope
facfit <- lbreg(gici10~avgdyenteropcr,dat=dch,beach="Fairhope")

# Goddard
gdcfit <- lbreg(gici10~avgdyenteropcr,dat=dch,beach="Goddard")

# Malibu
macfit <- lbreg(gici10~avgdyenteropcr,dat=dch,beach="Malibu")

# Surfside
sucfit <- lbreg(gici10~avgdyenteropcr,dat=dch,beach="Surfside")

# --------------------------------------
# save the objects
# --------------------------------------
save.image("~/dropbox/13beaches/data/temp/aim1-regs.Rdata")


# --------------------------------------
# aggregate the results for plotting
# --------------------------------------

freshall  <- list(hufit,sifit,wpfit,wefit)
marineall <- list(bofit,edfit,fafit,gdfit,mafit,sufit)

freshch  <- list(hucfit,sicfit,wpcfit,wecfit)
marinech <- list(bocfit,edcfit,facfit,gdcfit,macfit,sucfit)

fallests <- t(sapply(freshall,estci))
fallests <- matrix(as.numeric(fallests),nrow=4,ncol=3)
mallests <- t(sapply(marineall,estci))
mallests <- matrix(as.numeric(mallests),nrow=6,ncol=3)

fchests <- t(sapply(freshch,estci))
fchests <- matrix(as.numeric(fchests),nrow=4,ncol=3)
mchests <- t(sapply(marinech,estci))
mchests <- matrix(as.numeric(mchests),nrow=6,ncol=3)

# get the pesky interaction results
avres <- iestci(avfit$fit,avfit$vcovCL)
dores <- iestci(dhfit$fit,dhfit$vcovCL)
avresch <- iestci(avcfit$fit,avcfit$vcovCL)
doresch <- iestci(dhcfit$fit,dhcfit$vcovCL)
mallests2 <- rbind(avres,dores,mallests)
mchests2 <- rbind(avresch,doresch,mchests)

# label the beaches / conditions
beachnmsf <- unique(d$beach[d$marine==0])
beachnmsm <- unique(d$beach[d$marine==1])
beachnmsm2 <- c("Avalon, low groundwater","Avalon, high groundwater","Doheny, berm closed","Doheny, berm open",beachnmsm[c(2,4:length(beachnmsm))])

rownames(fallests) <- beachnmsf
rownames(fchests) <- beachnmsf
rownames(mallests2) <- beachnmsm2
rownames(mchests2) <- beachnmsm2


# --------------------------------------
# plot results
# --------------------------------------

# datasets for plotting
# total population
allpdata <- data.frame(rbind(fallests,mallests2))
allpdata$type<-factor(c(rep("Fresh",4),rep("Marine",10)))
bl <- rev(rownames(allpdata))
allpdata$beach<-factor(rownames(allpdata),levels=bl)
rownames(allpdata) <- 1:nrow(allpdata)
allpdata <- allpdata[,c("beach","type","est","lb","ub")]

# children <=10
chpdata <- data.frame(rbind(fchests,mchests2))
chpdata$type<-factor(c(rep("Fresh",4),rep("Marine",10)))
bl <- rev(rownames(chpdata))
chpdata$beach<-factor(rownames(chpdata),levels=bl)
rownames(chpdata) <- 1:nrow(chpdata)
chpdata <- chpdata[,c("beach","type","est","lb","ub")]

# global plot parameters
ytics <- c(0.5,1,2,4)
ylims <- c(0.3,4)
cols <- brewer.pal(3,"Set1")

# plot of CIRs for total population
pdf("~/dropbox/13beaches/aim1-results/figs/cirplot-epaQPCR-all.pdf",width=6,height=8)
dp <- ggplot(data=allpdata,aes(order=type,color=type) ) 
dp  + labs(x="",y="CIR for GI Illness") + 
scale_y_continuous(breaks=ytics,labels=round(ytics,2),trans=log_trans()) +
coord_flip(ylim=ylims) +
scale_x_discrete(breaks=allpdata$beach,labels=allpdata$beach) +
geom_hline(xintercept=1, color="grey40") +
geom_errorbar( aes(x=beach,y=est,ymin=lb,ymax=ub,width=0 ) ) +
stat_identity( aes(x=beach,y=est ) ) +
scale_colour_manual(values=cols) 
dev.off()

# plot of CIRs for children <=10
ylims <- c(0.1,4)
pdf("~/dropbox/13beaches/aim1-results/figs/cirplot-epaQPCR-child.pdf",width=6,height=8)
dp <- ggplot(data=chpdata,aes(order=type,color=type) ) 
dp  + labs(x="",y="CIR for GI Illness") + 
scale_y_continuous(breaks=ytics,labels=round(ytics,2),trans=log_trans()) +
coord_flip(ylim=ylims) +
scale_x_discrete(breaks=chpdata$beach,labels=chpdata$beach) +
geom_hline(xintercept=1, color="grey40") +
geom_errorbar( aes(x=beach,y=est,ymin=lb,ymax=ub,width=0 ) ) +
stat_identity( aes(x=beach,y=est ) ) +
scale_colour_manual(values=cols) 
dev.off()


