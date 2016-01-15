##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 11/10/15

# Analyses for bar graph comparing adjusted risk of illness
# under different conditions
##########################################

rm(list=ls())

# --------------------------------------
# load the and pre-preprocess the 
# analysis dataset
# (refer to the base functions script
# for details on the pre-processing)
# --------------------------------------
beaches13=read.csv("~/Dropbox/13beaches/data/final/13beaches-analysis.csv")

# load base functions
source("Programs/Analysis/0-base-functions.R")

data=preprocess.6beaches(beaches13)

# restrict to 6 beaches with coliphage data
beach.list=c("Avalon","Doheny","Malibu","Mission Bay",
             "Fairhope","Goddard")

all=data[data$beach %in% beach.list,]

# drop individuals with no water quality information
all=subset(all,nowq==0)

# ------------------------------------------------
# function to make table row with exponentiated point estimate
# and indicator of p-value significant at the 0.01 level
# ------------------------------------------------
est.pY=function(form,data,nameX,xlevel){
  # fit the model
  fit <- glm(form,family=poisson(link="log"),data=data)
  
  # get the marginal average predicted probability if all are exposed
  pd <- data
  pd[,c(nameX)]<- xlevel
  pY <- mean(predict(fit,newdata=pd,type="response"))
  return(pY)
}

boot.pY <- function(fmla,data,nameX,ID,iter,xlevel,dots=TRUE) {
  # jade modified ben's code since ben's code stratified by beach
  # fmla  : a model formula for glm, which is used to fit the exposure-outcome relationship
  # data  : data.frame that includes all of the variables in the formula, and the population used for estimation
  # nameX : string name of the exposure variable (entero1600, enteroQPCR)
  # ID    : ID for resampling (e.g., household ID required b/c of repeated, potentially correlated obs w/in HH)
  # iter  : number of bootstrap iterations  
  # dots  : logical. print dots for iterations
  if(length(ID)!=nrow(data)) {
    print("data and ID args must be the same length")
    break
  }
  
  # add ID variables to data frame for convenience 
  data$ID <- ID
  
  if(dots) cat("\n\n\nBootstrap Iterations (",iter,") \n----|--- 1 ---|--- 2 ---|--- 3 ---|--- 4 ---| --- 5 \n",sep="") 
  start.time <- Sys.time()
  
  # create an empty matrix to store bootstrap results
  # corresponding to the PAR and PAF parameters of interest + 3 marginal means
  pY <- matrix(NA,ncol=1,nrow=iter)
  
  # Create the bootstrap samples of size n_i, within each stratum i
  IDlist <- unique(data$ID)
  bs <- matrix(sample(IDlist, length(IDlist)*iter, replace=TRUE), ncol=iter)
  
  # calculate P(Y|X,W) for each bootstrap sample
  for (bb in 1:iter) {   
    bd <- merge(data.frame(ID=bs[,bb]),data,by="ID",all.x=TRUE)
    pY[bb,] <- tryCatch( do.call(est.pY,args=list(form=fmla,data=bd,nameX=nameX,xlevel=xlevel)), 
                         error=function(e) rep(NA,1)) 
    if(dots) cat(".",sep="")
    if(dots && bb %% 50 == 0) cat(bb,"\n")
  }
  if(dots) cat("\n Bootstrap Run Time:",round(difftime(Sys.time(),start.time,units="mins"),3)," Minutes \n\n")
  
  # if there were any bootstrap samples where the model failed to converge 
  # for any reason, caught by tryCatch(), then report that.
  # also triage any estimates that are clearly from a failed model that sort-of converged
  # (due to model instability from sparse data in a bootstrap sample)
  pY[(pY[,1]>10 | pY[,1]< -10), ] <- NA
  Nna <- sum(ifelse(is.na(pY[,1])==TRUE,1,0))
  if(Nna>0) {
    cat("\n   Warning: ",Nna," bootstrap replicates failed to converge to sensible estimates \n")
    cat("\n   Bootstrap estimates are based on the remaining ",iter-Nna," replicates\n\n")
  }
  
  # calculate the mean and percentile 95% confidence intervals for the statistics
  bootest  <- apply(pY,2,function(x) mean(x,na.rm=TRUE))
  bootse <- apply(pY,2,function(x) sd(x,na.rm=TRUE))
  boot95lb <- apply(pY,2,function(x) quantile(x,prob=0.025,na.rm=TRUE))
  boot95ub <- apply(pY,2,function(x) quantile(x,prob=0.975,na.rm=TRUE))
  
  list(bootest=bootest,boot95lb=boot95lb,boot95ub=boot95ub,bootse=bootse,pY=pY)
}

iter = 2

#-----------------------------------------
# Non-swimmers
#-----------------------------------------
all$nobodycontact=ifelse(all$anycontact=="Yes","No","Yes")
all$nobodycontact=as.factor(all$nobodycontact)
ns.bs=boot.pY(fmla=gici10~nobodycontact+agecat+female+racewhite+gichron+anim_any+
   gicontactbase+rawfood+beach,dat=all,nameX="nobodycontact",xlevel="Yes",
   ID=all[,"hhid"],iter)

#-----------------------------------------
# Swimmers not exposed to coliphage
#-----------------------------------------
unexposed=subset(all,all$fmc.pres==0 & all$fpc.pres==0)
s.bs=boot.pY(fmla=gici10~bodycontact+agecat+female+racewhite+gichron+anim_any+
   gicontactbase+rawfood+beach,dat=unexposed,nameX="nocoli",xlevel=1,
   ID=unexposed[,"hhid"],iter)

#-----------------------------------------
# Swimmers exposed to somatic coliphage
#-----------------------------------------
swimmers=subset(all,all$bodycontact=="Yes")
fmc.bs=boot.pY(fmla=gici10~as.factor(fmc.pres)+agecat+female+racewhite+gichron+anim_any+
  gicontactbase+rawfood+beach,dat=swimmers[!is.na(swimmers$fmc.pres),],nameX="fmc.pres",xlevel=1,
  ID=swimmers[!is.na(swimmers$fmc.pres),"hhid"],iter)
fmc.bs2=boot.pY(fmla=gici10~as.factor(fmc.pres)+agecat+female+racewhite+gichron+anim_any+
                 gicontactbase+rawfood+beach,dat=all[!is.na(all$fmc.pres),],nameX="fmc.pres",xlevel=1,
               ID=all[!is.na(all$fmc.pres),"hhid"],iter)
#-----------------------------------------
# Swimmers exposed to somatic coliphage
#-----------------------------------------
fpc.bs=boot.pY(fmla=gici10~as.factor(fpc.pres)+agecat+female+racewhite+gichron+anim_any+
  gicontactbase+rawfood+beach,dat=swimmers[!is.na(swimmers$fpc.pres),],nameX="fpc.pres",xlevel=1,
  ID=swimmers[!is.na(swimmers$fpc.pres),"hhid"],iter)

#-----------------------------------------
# Swimmers exposed to somatic coliphage + 
# enterococcus > 35
#-----------------------------------------
swimmers$fmc.ent=NA
swimmers$fmc.ent[swimmers$fmc.pres==0 & swimmers$entero35==0]=1
swimmers$fmc.ent[swimmers$fmc.pres==1 & swimmers$entero35==0]=2
swimmers$fmc.ent[swimmers$fmc.pres==0 & swimmers$entero35==1]=3
swimmers$fmc.ent[swimmers$fmc.pres==1 & swimmers$entero35==1]=4
swimmers$fmc.ent=as.factor(swimmers$fmc.ent)

fmc.ent.bs=boot.pY(fmla=gici10~as.factor(fmc.ent)+agecat+female+racewhite+gichron+anim_any+
  gicontactbase+rawfood+beach,dat=swimmers[!is.na(swimmers$fmc.ent),],nameX="fmc.ent",xlevel=4,
  ID=swimmers[!is.na(swimmers$fmc.ent),"hhid"],iter)

#-----------------------------------------
# Swimmers exposed to male-specific coliphage + 
# enterococcus > 35
#-----------------------------------------
swimmers$fpc.ent=NA
swimmers$fpc.ent[swimmers$fpc.pres==0 & swimmers$entero35==0]=1
swimmers$fpc.ent[swimmers$fpc.pres==1 & swimmers$entero35==0]=2
swimmers$fpc.ent[swimmers$fpc.pres==0 & swimmers$entero35==1]=3
swimmers$fpc.ent[swimmers$fpc.pres==1 & swimmers$entero35==1]=4
swimmers$fpc.ent=as.factor(swimmers$fpc.ent)

fpc.ent.bs=boot.pY(fmla=gici10~as.factor(fpc.ent)+agecat+female+racewhite+gichron+anim_any+
  gicontactbase+rawfood+beach,dat=swimmers[!is.na(swimmers$fpc.ent),],nameX="fpc.ent",xlevel=1,
  ID=swimmers[!is.na(swimmers$fpc.ent),"hhid"],iter)

save(ns.bs, s.bs, fmc.bs, fpc.bs, fmc.ent.bs, fpc.ent.bs,
     file="~/dropbox/coliphage/results/rawoutput/bargraph-output.Rdata")



