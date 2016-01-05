
#-------------------------------------
# x-mission-bay-entero-summary.R
#
# ben arnold
#
# summarize the mission bay entero
# data using EPA 1600 and Colilert
#
#-------------------------------------


#-------------------------------------
# input files:
#   mb-wq-samples.csv
#
# output files:
#   xxx
#
#-------------------------------------


#-------------------------------------
# preamble
#-------------------------------------

rm(list=ls())
library(RColorBrewer)
library(scales)


#-------------------------------------
# load the MB water sample data
#-------------------------------------

d <- read.csv("~/dropbox/13beaches/data/final/mb-wq-samples.csv")
d$coldate <- as.Date(d$coldate,"%d%b%Y")

# drop empty data in May 24-26
d <- subset(d,coldate>as.Date("2003-05-26"))

#-------------------------------------
# impute non-detects and calculate
# log10 values
#-------------------------------------

d$entero1600cfu[d$entero1600cfu_nd=="Below detection" | d$entero1600cfu<=0] <- 0.1
d$enteroCLTmpn[d$enteroCLTmpn_nd=="Below detection" | d$enteroCLTmpn<=0] <- 0.1

d$log10entero1600 <- log10(d$entero1600cfu)
d$log10enteroCLT  <- log10(d$enteroCLTmpn)


#-------------------------------------
# summarize the number of samples of
# each assay
#-------------------------------------

d1600 <- subset(d,!is.na(d$entero1600cfu))
# ensure just one sample per beach / day
dim(d1600)
table(duplicated(d1600[,c("beachcode","coldate")]))

# tabulate number of NDs
table(d1600$entero1600cfu_nd)


dCLT <- subset(d,!is.na(d$enteroCLTmpn))
# ensure just one sample per beach / day / station / samplenumber
dim(dCLT)
table(duplicated(dCLT[,c("beachcode","coldate","stationid","samplenum")]))

# tabulate number of NDs
table(dCLT$enteroCLTmpn_nd)


#-------------------------------------
# calculate the daily mean log10
# values for each assay
#-------------------------------------

# overall geometric means
10^mean(d$log10entero1600,na.rm=T)
10^mean(d$log10enteroCLT,na.rm=T)

mu1600 <- tapply(d$log10entero1600,list(d$coldate,d$beachcode),function(x) mean(x,na.rm=T),simplify=TRUE)
muCLT <- tapply(d$log10enteroCLT,list(d$coldate,d$beachcode),function(x) mean(x,na.rm=T))

# now count non-detects after taking the averages
table(mu1600==-1)
table(muCLT==-1)

# cross-tab of values >35 / 100ml
table(ifelse(10^mu1600>35," >35 EPA 1600","<=35 EPA 1600"),ifelse(10^muCLT>35,">35 Colilert","<=35 Colilert"))


# calculate Spearman's rank correlation coeff (rho)
sp.rho <- cor.test(x=mu1600,y=muCLT,method="spearman")
rho.text <- substitute(paste("Spearman's ",rho," = ",rho.txt ),list(rho.txt=sprintf("%1.2f",sp.rho$estimate)))


pdf("~/dropbox/13beaches/aim1-results/figs/aim1-Mission-Bay-Enterococcus-1600-v-Colilert-scatter.pdf")
cols <- c(brewer.pal(9,"YlGnBu")[6], brewer.pal(9,"YlOrBr")[4]   )
ytics <- xtics <- seq(-1,4,by=1)
# combined plot
plot(mu1600,muCLT,type="n",
	xaxt="n",xlab="Log10 Enterococcus EPA 1600 CFU/100 ml",xlim=range(xtics),
	yaxt="n",ylab="",ylim=range(ytics),
	bty="n",las=1,main=""
	)
	axis(1,at=xtics,las=1)
	axis(2,at=ytics,las=1)
	mtext("Log10 Enterococcus Colilert MPN/100 ml",side=2,line=2)
	segments(x0=log10(35),y0=min(ytics),y1=max(ytics),lty=2,lwd=1.5,col="gray40")
	segments(y0=log10(35),x0=min(xtics),x1=max(xtics),lty=2,lwd=1.5,col="gray40")
	text(x=max(xtics),y=log10(35)+0.2,"35 MPN",adj=1,col="gray20")
	text(x=log10(35),y=max(ytics)-0.2,"35 CFU",pos=4,col="gray20")
	
	points(mu1600,muCLT,pch=16,cex=0.7,col=alpha(cols[1],0.5))

	mtext("Colilert vs. EPA 1600\nMission Bay Beach-level Daily Averages",side=3,line=1,cex=1.2,font=2)
	text(3,0,rho.text)
dev.off()


pdf("~/dropbox/13beaches/aim1-results/figs/aim1-Mission-Bay-Enterococcus-1600-v-Colilert-density.pdf")
plot(density(muCLT[!is.na(muCLT)]),
	main="Distribution of Mission Bay Daily Averages",
	xlim=c(-2,4), xlab="Average log10 Enterococcus per 100 ml",
	ylim=c(0,1),
	las=1,
	bty="n"
)
lines(density(mu1600[!is.na(mu1600)]),col=cols[1])
text(0,0.5,"Colilert")
text(2.5,0.5,"EPA 1600",col=cols[1])
dev.off()



