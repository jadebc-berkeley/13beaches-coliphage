##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# Plots of the RR for coliphage across the range of 
# concentration Enterococcus to assess
# continuous x continuous interaction

# Couldn't get the delta method function to 
# work inside a function, so coding this in an
# inefficient way

# NOTE these SEs are NOT adjusted for clustering

# Results pooled across beach
##########################################

library(msm)
library(car)
rm(list=ls())

load("~/Documents/CRG/coliphage/results/rawoutput/regress-10day-continuous-body-joint.Rdata")

# --------------------------------------
# All conditions
# --------------------------------------
entero.levels=seq(-1,4,length=30)

####### FMC 1601 ####### 
# extract coefficients
reg=all.fit10.fmc1601
x=reg$coef
inter.row=grep(":",names(x))
slopes <- x[2] + x[inter.row]*entero.levels
estmean<-coef(reg)
var<-vcov(reg)

# estimate standard error
SEs <- rep(NA, length(entero.levels))

for (i in 1:length(entero.levels)){
  j <- entero.levels[i]
  SEs[i] <- deltamethod (~ (x2) + (x23)*j, estmean, var)
}

upper <- slopes + 1.96*SEs
lower <- slopes - 1.96*SEs

df=data.frame(cbind(entero.levels, exp(slopes), exp(lower), exp(upper)))
colnames(df)=c("entero","cir.coli","lb","ub")

pdf("~/Documents/CRG/coliphage/results/figures/Entero_X_Coli_intrxn_fmc1601.pdf",height=4,width=3.5)
plot(df$entero, df$cir.coli, type = "l", lty = 1, ylim = c(0,3), xlab = "Enterococcus level", 
     ylab = "CIR for Somatic Coliphage 1601")
points(df$entero, df$lb, type = "l", lty = 2)
points(df$entero, df$ub, type = "l", lty = 2)
points(df$entero, rep(1, length(df$entero)),type = "l", col = "red",lty=3)
title("All conditions")
dev.off()

####### FMC 1602 ####### 
# extract coefficients
reg=all.fit10.fmc1602
x=reg$coef
inter.row=grep(":",names(x))
slopes <- x[2] + x[inter.row]*entero.levels
estmean<-coef(reg)
var<-vcov(reg)

# estimate standard error
SEs <- rep(NA, length(entero.levels))

for (i in 1:length(entero.levels)){
  j <- entero.levels[i]
  SEs[i] <- deltamethod (~ (x2) + (x22)*j, estmean, var)
}

upper <- slopes + 1.96*SEs
lower <- slopes - 1.96*SEs

df=data.frame(cbind(entero.levels, exp(slopes), exp(lower), exp(upper)))
colnames(df)=c("entero","cir.coli","lb","ub")

pdf("~/Documents/CRG/coliphage/results/figures/Entero_X_Coli_intrxn_fmc1602.pdf",height=4,width=3.5)
plot(df$entero, df$cir.coli, type = "l", lty = 1, ylim = c(0,3), xlab = "Enterococcus level", 
     ylab = "CIR for Somatic Coliphage 1602")
points(df$entero, df$lb, type = "l", lty = 2)
points(df$entero, df$ub, type = "l", lty = 2)
points(df$entero, rep(1, length(df$entero)),type = "l", col = "red",lty=3)
title("All conditions")
dev.off()

####### FPC 1601 ####### 
# extract coefficients
reg=all.fit10.fpc1601
x=reg$coef
inter.row=grep(":",names(x))
slopes <- x[2] + x[inter.row]*entero.levels
estmean<-coef(reg)
var<-vcov(reg)

# estimate standard error
SEs <- rep(NA, length(entero.levels))

for (i in 1:length(entero.levels)){
  j <- entero.levels[i]
  SEs[i] <- deltamethod (~ (x2) + (x26)*j, estmean, var)
}

upper <- slopes + 1.96*SEs
lower <- slopes - 1.96*SEs

df=data.frame(cbind(entero.levels, exp(slopes), exp(lower), exp(upper)))
colnames(df)=c("entero","cir.coli","lb","ub")

pdf("~/Documents/CRG/coliphage/results/figures/Entero_X_Coli_intrxn_fpc1601.pdf",height=4,width=3.5)
plot(df$entero, df$cir.coli, type = "l", lty = 1, ylim = c(0,3), xlab = "Enterococcus level", 
     ylab = "CIR for Male-Specific Coliphage 1601")
points(df$entero, df$lb, type = "l", lty = 2)
points(df$entero, df$ub, type = "l", lty = 2)
points(df$entero, rep(1, length(df$entero)),type = "l", col = "red",lty=3)
title("All conditions")
dev.off()

####### FPC 1602 ####### 
# extract coefficients
reg=all.fit10.fpc1602
x=reg$coef
inter.row=grep(":",names(x))
slopes <- x[2] + x[inter.row]*entero.levels
estmean<-coef(reg)
var<-vcov(reg)

# estimate standard error
SEs <- rep(NA, length(entero.levels))

for (i in 1:length(entero.levels)){
  j <- entero.levels[i]
  SEs[i] <- deltamethod (~ (x2) + (x22)*j, estmean, var)
}

upper <- slopes + 1.96*SEs
lower <- slopes - 1.96*SEs

df=data.frame(cbind(entero.levels, exp(slopes), exp(lower), exp(upper)))
colnames(df)=c("entero","cir.coli","lb","ub")

pdf("~/Documents/CRG/coliphage/results/figures/Entero_X_Coli_intrxn_fpc1602.pdf",height=4,width=3.5)
plot(df$entero, df$cir.coli, type = "l", lty = 1, ylim = c(-2,12), xlab = "Enterococcus level", 
     ylab = "CIR for Male-Specific Coliphage 1602")
points(df$entero, df$lb, type = "l", lty = 2)
points(df$entero, df$ub, type = "l", lty = 2)
points(df$entero, rep(1, length(df$entero)),type = "l", col = "red",lty=3)
title("All conditions")
dev.off()


# --------------------------------------
# Low risk conditions
# --------------------------------------
####### FMC 1602 ####### 
# extract coefficients
reg=all.fit10.fmc1602.low
x=reg$coef
inter.row=grep(":",names(x))
slopes <- x[2] + x[inter.row]*entero.levels
estmean<-coef(reg)
var<-vcov(reg)

# estimate standard error
SEs <- rep(NA, length(entero.levels))

for (i in 1:length(entero.levels)){
  j <- entero.levels[i]
  SEs[i] <- deltamethod (~ (x2) + (x22)*j, estmean, var)
}

upper <- slopes + 1.96*SEs
lower <- slopes - 1.96*SEs

df=data.frame(cbind(entero.levels, exp(slopes), exp(lower), exp(upper)))
colnames(df)=c("entero","cir.coli","lb","ub")

pdf("~/Documents/CRG/coliphage/results/figures/Entero_X_Coli_intrxn_fmc1602_low.pdf",height=4,width=3.5)
plot(df$entero, df$cir.coli, type = "l", lty = 1, ylim = c(0,3), xlab = "Enterococcus level", 
     ylab = "CIR for Somatic Coliphage 1602")
points(df$entero, df$lb, type = "l", lty = 2)
points(df$entero, df$ub, type = "l", lty = 2)
points(df$entero, rep(1, length(df$entero)),type = "l", col = "red",lty=3)
title("Low risk conditions")
dev.off()

####### FPC 1601 ####### 
# extract coefficients
reg=all.fit10.fpc1601.low
x=reg$coef
inter.row=grep(":",names(x))
slopes <- x[2] + x[inter.row]*entero.levels
estmean<-coef(reg)
var<-vcov(reg)

# estimate standard error
SEs <- rep(NA, length(entero.levels))

for (i in 1:length(entero.levels)){
  j <- entero.levels[i]
  SEs[i] <- deltamethod (~ (x2) + (x24)*j, estmean, var)
}

upper <- slopes + 1.96*SEs
lower <- slopes - 1.96*SEs

df=data.frame(cbind(entero.levels, exp(slopes), exp(lower), exp(upper)))
colnames(df)=c("entero","cir.coli","lb","ub")

pdf("~/Documents/CRG/coliphage/results/figures/Entero_X_Coli_intrxn_fpc1601_low.pdf",height=4,width=3.5)
plot(df$entero, df$cir.coli, type = "l", lty = 1, ylim = c(0,3), xlab = "Enterococcus level", 
     ylab = "CIR for Male-Specific Coliphage 1601")
points(df$entero, df$lb, type = "l", lty = 2)
points(df$entero, df$ub, type = "l", lty = 2)
points(df$entero, rep(1, length(df$entero)),type = "l", col = "red",lty=3)
title("Low risk conditions")
dev.off()


####### FPC 1602 ####### 
# extract coefficients
reg=all.fit10.fpc1602.low
x=reg$coef
inter.row=grep(":",names(x))
slopes <- x[2] + x[inter.row]*entero.levels
estmean<-coef(reg)
var<-vcov(reg)

# estimate standard error
SEs <- rep(NA, length(entero.levels))

for (i in 1:length(entero.levels)){
  j <- entero.levels[i]
  SEs[i] <- deltamethod (~ (x2) + (x22)*j, estmean, var)
}

upper <- slopes + 1.96*SEs
lower <- slopes - 1.96*SEs

df=data.frame(cbind(entero.levels, exp(slopes), exp(lower), exp(upper)))
colnames(df)=c("entero","cir.coli","lb","ub")

pdf("~/Documents/CRG/coliphage/results/figures/Entero_X_Coli_intrxn_fpc1602_low.pdf",height=4,width=3.5)
plot(df$entero, df$cir.coli, type = "l", lty = 1, ylim = c(-2,12), xlab = "Enterococcus level", 
     ylab = "CIR for Male-Specific Coliphage 1602")
points(df$entero, df$lb, type = "l", lty = 2)
points(df$entero, df$ub, type = "l", lty = 2)
points(df$entero, rep(1, length(df$entero)),type = "l", col = "red",lty=3)
title("Low risk conditions")
dev.off()


# --------------------------------------
# High risk conditions conditions
# --------------------------------------
####### FMC 1602 ####### 
# extract coefficients
reg=all.fit10.fmc1602.high
x=reg$coef
inter.row=grep(":",names(x))
slopes <- x[2] + x[inter.row]*entero.levels
estmean<-coef(reg)
var<-vcov(reg)

# estimate standard error
SEs <- rep(NA, length(entero.levels))

for (i in 1:length(entero.levels)){
  j <- entero.levels[i]
  SEs[i] <- deltamethod (~ (x2) + (x22)*j, estmean, var)
}

upper <- slopes + 1.96*SEs
higher <- slopes - 1.96*SEs

df=data.frame(cbind(entero.levels, exp(slopes), exp(higher), exp(upper)))
colnames(df)=c("entero","cir.coli","lb","ub")

pdf("~/Documents/CRG/coliphage/results/figures/Entero_X_Coli_intrxn_fmc1602_high.pdf",height=4,width=3.5)
plot(df$entero, df$cir.coli, type = "l", lty = 1, ylim = c(0,3), xlab = "Enterococcus level", 
     ylab = "CIR for Somatic Coliphage 1602")
points(df$entero, df$lb, type = "l", lty = 2)
points(df$entero, df$ub, type = "l", lty = 2)
points(df$entero, rep(1, length(df$entero)),type = "l", col = "red",lty=3)
title("High risk conditions")
dev.off()

####### FPC 1601 ####### 
# extract coefficients
reg=all.fit10.fpc1601.high
x=reg$coef
inter.row=grep(":",names(x))
slopes <- x[2] + x[inter.row]*entero.levels
estmean<-coef(reg)
var<-vcov(reg)

# estimate standard error
SEs <- rep(NA, length(entero.levels))

for (i in 1:length(entero.levels)){
  j <- entero.levels[i]
  SEs[i] <- deltamethod (~ (x2) + (x24)*j, estmean, var)
}

upper <- slopes + 1.96*SEs
higher <- slopes - 1.96*SEs

df=data.frame(cbind(entero.levels, exp(slopes), exp(higher), exp(upper)))
colnames(df)=c("entero","cir.coli","lb","ub")

pdf("~/Documents/CRG/coliphage/results/figures/Entero_X_Coli_intrxn_fpc1601_high.pdf",height=4,width=3.5)
plot(df$entero, df$cir.coli, type = "l", lty = 1, ylim = c(0,3), xlab = "Enterococcus level", 
     ylab = "CIR for Male-Specific Coliphage 1601")
points(df$entero, df$lb, type = "l", lty = 2)
points(df$entero, df$ub, type = "l", lty = 2)
points(df$entero, rep(1, length(df$entero)),type = "l", col = "red",lty=3)
title("High risk conditions")
dev.off()

####### FPC 1602 ####### 
# extract coefficients
reg=all.fit10.fpc1602.high
x=reg$coef
inter.row=grep(":",names(x))
slopes <- x[2] + x[inter.row]*entero.levels
estmean<-coef(reg)
var<-vcov(reg)

# estimate standard error
SEs <- rep(NA, length(entero.levels))

for (i in 1:length(entero.levels)){
  j <- entero.levels[i]
  SEs[i] <- deltamethod (~ (x2) + (x22)*j, estmean, var)
}

upper <- slopes + 1.96*SEs
higher <- slopes - 1.96*SEs

df=data.frame(cbind(entero.levels, exp(slopes), exp(higher), exp(upper)))
colnames(df)=c("entero","cir.coli","lb","ub")

pdf("~/Documents/CRG/coliphage/results/figures/Entero_X_Coli_intrxn_fpc1602_high.pdf",height=4,width=3.5)
plot(df$entero, df$cir.coli, type = "l", lty = 1, ylim = c(-2,12), xlab = "Enterococcus level", 
     ylab = "CIR for Male-Specific Coliphage 1602")
points(df$entero, df$lb, type = "l", lty = 2)
points(df$entero, df$ub, type = "l", lty = 2)
points(df$entero, rep(1, length(df$entero)),type = "l", col = "red",lty=3)
title("High risk conditions")
dev.off()

