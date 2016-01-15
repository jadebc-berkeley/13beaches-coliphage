

# --------------------------------------
# base-functions.R
# jade benjamin-chung
# modified code by benjamin arnold
#
# description:
# base functions used in analyses
#
# version 2 (sep 16, 2015) 
# added pooling across epa 1601 and 1602
# for coliphage indicators
#
# version 1 (july 13, 2015)
# --------------------------------------



# --------------------------------------
# Standard routine to load the 6 beaches
# analysis data and do any final processing
# required in R for the final water 
# quality analyses
# --------------------------------------


preprocess.6beaches <- function(dataset) {

  # --------------------------------------
  # load the analysis dataset
  # --------------------------------------
  d <- dataset
  
  # convert ID variables from factors to strings
  d$hhid <- as.character(d$hhid)
  d$indid <- as.character(d$indid)
  
  # --------------------------------------
  # subset to observations for analysis
  # --------------------------------------
  
  cat("\nSuccessfully loaded the data\n","Total sample size =",nrow(d),"\n")
  
  # drop individuals with baseline GI illness
  cat("\nDropping individuals with GI illness at enrollment\n","N =",table(d$gibase)["Yes"],"\n")
  d <- subset(d,gibase!="Yes")
  
  
  cat("\nFinal sample size =",nrow(d),"\n")
  
  
  # --------------------------------------
  # subset dataset to variables of interest
  # to speed up computations / processing
  # --------------------------------------
  
  # included: agecat, female, racewhite, gichron, anim_any, gicontactbase, rawfood
  # excluded: allergies (not relevant), frequency of beach visits (not measured at all beaches), digging in the sand (not measured at all beaches)
  cat("\nSubsetting the data to relevant variables and completing final variable pre-processing")
  ad <- subset(d,select=c("beach","pointsource","marine","hhid","indid",
    "groundwater","berm","anycontact","bodycontact","headunder","swallwater",
    "diarrheaci1","diarrheaci2","diarrheaci3","diarrheaci4","diarrheaci5",
    "diarrheaci6","diarrheaci7","diarrheaci8","diarrheaci9","diarrheaci10",
    "gici3","gici10","dailygi","workgi","medgi","medvisits","age","agecat",
    "agestrat","female","race","gichron","anim_any","gicontactbase","rawfood",
    "nowq","avgdyentero1600","qavgdyentero1600","avgdyenteropcr","qavgdyenteropcr",
    "avgdyfmc1601","avgdyfmc1602","avgdyfpc1601","avgdyfpc1602","avgdyenteroELT"))
  
  # create a race=white indicator
  ad$racewhite <- factor(ifelse(ad$race=="white","Yes",NA),
     levels=c("No","Yes","Missing"))
  ad$racewhite[ad$race=="missing"]<-"Missing"
  ad$racewhite[ad$race!="white" & ad$race!="missing"] <- "No"
  
  # for some covariates, move missing category to the last factor category
  levels(ad$gichron) <- c("Missing","No","Yes")
  ad$gichron <- factor(ad$gichron,levels=c("No","Yes","Missing"))
  ad$gicontactbase <- factor(ad$gicontactbase,levels=c("No","Yes","Missing"))
  
  # reorder groundwater factor for convenience w/ reg estimates
  ad$groundwater <- factor(ad$groundwater,levels=c("Below median flow",
    "Above median flow",""))
  
  # reorder agecat factor so they are in correct order
  ad$agecat <- factor(ad$agecat,levels=c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65-74","75+","Missing"))
  
  # re-order the age stratification factor so that older ages are reference (for convenience, since younger children have higher CIRs)
  ad$agestrat <- factor(ad$agestrat,levels=c(">10","(4, 10]","(0, 4]",""),labels=c(">10","(4, 10]","(0, 4]","Missing"))
  
  # create an indicator of >35 CFU/100ml for Entero EPA 1600
  ad$entero <- ad$avgdyentero1600
  ad$entero[ad$beach=="Mission Bay"]= NA
  ad$entero[ad$beach=="Mission Bay"] = ad$avgdyenteroELT[ad$beach=="Mission Bay"]
    
  ad$entero35 <- as.factor(ifelse(ad$entero>log10(35),1,0))
  
  # create an indicator of > 470
  ad$entero470 <- ifelse(ad$avgdyenteropcr>log10(470),1,0)
  
  # create a shorter variable name for Entero for convenience
  ad$entero1600 <- ad$avgdyentero1600
  ad$qentero1600 <- factor(ad$qavgdyentero1600)
  ad$enteroQPCR <- ad$avgdyenteropcr
  ad$qenteroQPCR <- factor(ad$qavgdyenteropcr)
  
  drops <- c("avgdyentero1600","qavgdyentero1600","avgdyenteropcr",
             "qavgdyenteropcr")
  ad=ad[,!(names(ad) %in% drops)]
  
  # create a shorter variable name for coliphage for convenience
  ad$fmc1601 <- ad$avgdyfmc1601
  ad$fmc1602 <- ad$avgdyfmc1602
  ad$fpc1601 <- ad$avgdyfpc1601
  ad$fpc1602 <- ad$avgdyfpc1602
  
  drops <- c("avgdyfmc1601","avgdyfmc1602","avgdyfpc1601",
             "avgdyfpc1602")
  ad=ad[,!(names(ad) %in% drops)]
  
  # generate indicator for presence/absence of coliphage
  ad$fmc1601.pres <- as.factor(ifelse(ad$fmc1601==-1,0,1))
  ad$fmc1602.pres <- as.factor(ifelse(ad$fmc1602==-1,0,1))
  ad$fpc1601.pres <- as.factor(ifelse(ad$fpc1601==-1,0,1))
  ad$fpc1602.pres <- as.factor(ifelse(ad$fpc1602==-1,0,1))

  # generate indicator for above/below median coliphage conc
  ad$fmc1601.med <- as.factor(ifelse(ad$fmc1601==-1,0,1))
  ad$fmc1602.med <- as.factor(ifelse(ad$fmc1602==-1,0,1))
  ad$fpc1601.med <- as.factor(ifelse(ad$fpc1601==-1,0,1))
  ad$fpc1602.med <- as.factor(ifelse(ad$fpc1602==-1,0,1))
  
  # generate pooled indicator for presence/absence of coliphage
  ad$fmc.pres <- 0
  ad$fmc.pres[ad$fmc1601.pres==1 | ad$fmc1602.pres==1] <- 1
  ad$fmc.pres[is.na(ad$fmc1601.pres) & is.na(ad$fmc1602.pres)] <- NA
  
  ad$fpc.pres <- 0
  ad$fpc.pres[ad$fpc1601.pres==1 | ad$fpc1602.pres==1] <- 1
  ad$fpc.pres[is.na(ad$fpc1601.pres) & is.na(ad$fpc1602.pres)] <- NA
  
  # new risk variable
  # high risk 
  ad$risk=NA
  ad$risk[ad$beach=="Fairhope"| ad$beach=="Goddard"|
              (ad$beach=="Doheny" & ad$berm=="Open")|
              (ad$beach=="Avalon" & ad$groundwater=="Above median flow")]=1
  ad$risk[ad$beach=="Malibu"|ad$beach=="Mission Bay"|
              (ad$beach=="Doheny" & ad$berm=="Closed")|
              (ad$beach=="Avalon" & ad$groundwater=="Below median flow")]=0
  ad$risk=factor(ad$risk)
  levels(ad$risk)=c("Low","High")
  
  # drop the complete dataset to save memory and speed up computation
  # return the analysis-ready data: ad
  rm(d)
  return(ad)
  
}


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
# Convenience wrapper function to run 
# modified Poisson models and obtain 
# robust SEs (clusterd on hhid)
# this is the work horse of all the 
# regressions run in this analysis
# --------------------------------------

mpreg <- function(formula,dat,vcv=FALSE) {
  # modified Poisson regression formula
  # dataset used to fit the model	
  fit <- glm(formula,family=poisson(link="log"),data=dat)
  vcovCL <- cl(dat,fm=fit,cluster=dat$hhid)
  rfit <- coeftest(fit, vcovCL)
  print(summary(fit))
  cat("\n\nRobust, Sandwich Standard Errors Account for Clustering:\n")
  print(rfit) 
  if(vcv==FALSE) {
    return(rfit)
  } else {
    list(fit=rfit,vcovCL=vcovCL)
  }
}



#-------------------------------------------
# Function to calculate a marginally adjusted
# dose-response curve for continuous 
# Enterococcus exposures 
#-------------------------------------------
marginal.pY <- function(form,data,nameX,pX) {
  # calculates a marginally adjusted dose-response curve for
  # an exposure (e.g., entero1600, enteroQPCR), adjusted for covariates
  #
  # arguments:
  # form: glm model formula, must include "entero1600" as a covariate
  # data: data for glm estimation, including the var specified in nameX
  # nameX : string of the name of the exposure variable: "entero1600" or "enteroQPCR"
  # pX    : levels of the exposure for predicted probabilities
  #
  # returned objects:
  # pY : vector of predicted probabilities at pX
  
  # fit the model
  fit <- glm(form,family=poisson(link="log"),data=data)
  
  # for each level of the exposure, get the marginal average predicted probability
  pY <- as.numeric(rep(NA,length(pX)))
  pd <- data
  for (i in 1:length(pY)) {
    pd[,c(nameX)]<- pX[i]
    pY[i] <- mean(predict(fit,newdata=pd,type="response"))
  }
  
  return(pY)
}


#-------------------------------------------
# Stratified bootstrap resampling function
# Samples n_i households from stratum i,
# where stratum in this study is beach cohort
# 
# This function is tailored to store 
# exposure-outcome dose response curves 
#-------------------------------------------

boot.pY <- function(fmla,data,nameX,ID,iter,dots=TRUE) {
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
  
  # get a range of the exposure
  pX <- seq(round(min(data[,c(nameX)]),2),round(max(data[,c(nameX)]),2),by=0.1)
  
  # create an empty matrix to store bootstrap results
  # corresponding to the PAR and PAF parameters of interest + 3 marginal means
  pY <- matrix(NA,ncol=length(pX),nrow=iter)
  
  # Create the bootstrap samples of size n_i, within each stratum i
  IDlist <- unique(data$ID)
  bs <- matrix(sample(IDlist, length(IDlist)*iter, replace=TRUE), ncol=iter)
  
  # calculate P(Y|X,W) for each bootstrap sample
  for (bb in 1:iter) {   
    bd <- merge(data.frame(ID=bs[,bb]),data,by="ID",all.x=TRUE)
    pY[bb,] <- tryCatch( do.call(marginal.pY,args=list(form=fmla,data=bd,nameX=nameX,pX=pX)), 
          error=function(e) rep(NA,length(pX))) 
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
  boot95lb <- apply(pY,2,function(x) quantile(x,prob=0.025,na.rm=TRUE))
  boot95ub <- apply(pY,2,function(x) quantile(x,prob=0.975,na.rm=TRUE))
  
  list(bootest=bootest,boot95lb=boot95lb,boot95ub=boot95ub,pX=pX,pY=pY)
}


# function to get Estimates and CIs from a linear combination of regression coefficients
lccalc <- function(lc,x,vcv) {
  # lc : linear combination of coefficients
  # x : log-linear model object returned from coeftest (class=coeftest)
  # vcv : variance-covariance matrix of coefficients for robust SEs
  est <- exp(t(lc)%*%x[,1])
  se  <- sqrt( t(lc)%*%vcv%*%lc )
  lb <- exp(log(est)-1.96*se)
  ub <- exp(log(est)+1.96*se)
  res <- c(est,lb,ub)
  return(res)
}




