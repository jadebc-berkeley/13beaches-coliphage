

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
    "groundwater","berm","anycontact","bodycontact","headunder","swallwater","watertime",
    "diarrheaci1","diarrheaci2","diarrheaci3","diarrheaci4","diarrheaci5",
    "diarrheaci6","diarrheaci7","diarrheaci8","diarrheaci9","diarrheaci10",
    "gici3","gici10","dailygi","workgi","medgi","medvisits","age","agecat",
    "agestrat","female","race","gichron","anim_any","gicontactbase","rawfood",
    "nowq","avgdyentero1600","qavgdyentero1600","avgdyenteropcr","qavgdyenteropcr",
    "avgdyfmc1601","avgdyfmc1602","avgdyfpc1601","avgdyfpc1602","avgdyenteroELT",
    "fmc1601perdet","fmc1602perdet","fpc1601perdet","fpc1602perdet",
    "fpcperdet","fmcperdet"))
  
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
  
  # generate indicator for >25%, >50%, >75% of samples with detectable
  # coliphage each day
  ad$fmc1601_25 <- 0
  ad$fmc1601_25[ad$fmc1601perdet>0.25]  <- 1
  ad$fmc1601_25[is.na(ad$fmc1601perdet)] <- NA
  
  ad$fmc1601_50 <- 0
  ad$fmc1601_50[ad$fmc1601perdet>0.50]  <- 1
  ad$fmc1601_50[is.na(ad$fmc1601perdet)] <- NA
  
  ad$fmc1601_75 <- 0
  ad$fmc1601_75[ad$fmc1601perdet>0.75]  <- 1
  ad$fmc1601_75[is.na(ad$fmc1601perdet)] <- NA

  ad$fmc1602_50 <- 0
  ad$fmc1602_50[ad$fmc1602perdet>0.50]  <- 1
  ad$fmc1602_50[is.na(ad$fmc1602perdet)] <- NA
  
  ad$fmc1602_75 <- 0
  ad$fmc1602_75[ad$fmc1602perdet>0.75]  <- 1
  ad$fmc1602_75[is.na(ad$fmc1602perdet)] <- NA

  ad$fpc1601_25 <- 0
  ad$fpc1601_25[ad$fpc1601perdet>0.25]  <- 1
  ad$fpc1601_25[is.na(ad$fpc1601perdet)] <- NA
  
  ad$fpc1601_50 <- 0
  ad$fpc1601_50[ad$fpc1601perdet>0.50]  <- 1
  ad$fpc1601_50[is.na(ad$fpc1601perdet)] <- NA
  
  ad$fpc1601_75 <- 0
  ad$fpc1601_75[ad$fpc1601perdet>0.75]  <- 1
  ad$fpc1601_75[is.na(ad$fpc1601perdet)] <- NA

  ad$fpc1602_25 <- 0
  ad$fpc1602_25[ad$fpc1602perdet>0.25]  <- 1
  ad$fpc1602_25[is.na(ad$fpc1602perdet)] <- NA
  
  ad$fpc1602_50 <- 0
  ad$fpc1602_50[ad$fpc1602perdet>0.50]  <- 1
  ad$fpc1602_50[is.na(ad$fpc1602perdet)] <- NA
  
  ad$fpc1602_75 <- 0
  ad$fpc1602_75[ad$fpc1602perdet>0.75]  <- 1
  ad$fpc1602_75[is.na(ad$fpc1602perdet)] <- NA

# reference group <25%

  ad$fmc1601_25ref <- 0
  ad$fmc1601_25ref[ad$fmc1601perdet>0.25]  <- 1
  ad$fmc1601_25ref[is.na(ad$fmc1601perdet)] <- NA
  
  ad$fmc1601_50ref <- NA
  ad$fmc1601_50ref[ad$fmc1601perdet>0.50]  <- 1
  ad$fmc1601_50ref[ad$fmc1601perdet<0.25]  <- 0
  ad$fmc1601_50ref[is.na(ad$fmc1601perdet)] <- NA
  
  ad$fmc1601_75ref <- NA
  ad$fmc1601_75ref[ad$fmc1601perdet>0.75]  <- 1
  ad$fmc1601_75ref[ad$fmc1601perdet<0.25]  <- 0
  ad$fmc1601_75ref[is.na(ad$fmc1601perdet)] <- NA

  ad$fmc1602_25ref <- 0
  ad$fmc1602_25ref[ad$fmc1602perdet>0.25]  <- 1
  ad$fmc1602_25ref[is.na(ad$fmc1602perdet)] <- NA
  
  ad$fmc1602_50ref <- NA
  ad$fmc1602_50ref[ad$fmc1602perdet>0.50]  <- 1
  ad$fmc1602_50ref[ad$fmc1602perdet<0.25]  <- 0
  ad$fmc1602_50ref[is.na(ad$fmc1602perdet)] <- NA
  
  ad$fmc1602_75ref <- NA
  ad$fmc1602_75ref[ad$fmc1602perdet>0.75]  <- 1
  ad$fmc1602_75ref[ad$fmc1602perdet<0.25]  <- 0
  ad$fmc1602_75ref[is.na(ad$fmc1602perdet)] <- NA

  ad$fpc1601_25ref <- 0
  ad$fpc1601_25ref[ad$fpc1601perdet>0.25]  <- 1
  ad$fpc1601_25ref[is.na(ad$fpc1601perdet)] <- NA
  
  ad$fpc1601_50ref <- NA
  ad$fpc1601_50ref[ad$fpc1601perdet>0.50]  <- 1
  ad$fpc1601_50ref[ad$fpc1601perdet<0.25]  <- 0
  ad$fpc1601_50ref[is.na(ad$fpc1601perdet)] <- NA
  
  ad$fpc1601_75ref <- NA
  ad$fpc1601_75ref[ad$fpc1601perdet>0.75]  <- 1
  ad$fpc1601_75ref[ad$fpc1601perdet<0.25]  <- 0
  ad$fpc1601_75ref[is.na(ad$fpc1601perdet)] <- NA

  ad$fpc1602_25ref <- 0
  ad$fpc1602_25ref[ad$fpc1602perdet>0.25]  <- 1
  ad$fpc1602_25ref[is.na(ad$fpc1602perdet)] <- NA
  
  ad$fpc1602_50ref <- NA
  ad$fpc1602_50ref[ad$fpc1602perdet>0.50]  <- 1
  ad$fpc1602_50ref[ad$fpc1602perdet<0.25]  <- 0
  ad$fpc1602_50ref[is.na(ad$fpc1602perdet)] <- NA
  
  ad$fpc1602_75ref <- NA
  ad$fpc1602_75ref[ad$fpc1602perdet>0.75]  <- 1
  ad$fpc1602_75ref[ad$fpc1602perdet<0.25]  <- 0
  ad$fpc1602_75ref[is.na(ad$fpc1602perdet)] <- NA

  # pooled by assay 

  ad$fmc_25 <- 0
  ad$fmc_25[ad$fmcperdet>0.25] <- 1
  ad$fmc_25[is.na(ad$fmcperdet)] <- NA
  
  ad$fmc_50 <- 0
  ad$fmc_50[ad$fmcperdet>0.50] <- 1
  ad$fmc_50[is.na(ad$fmcperdet)] <- NA
  
  ad$fmc_75 <- 0
  ad$fmc_75[ad$fmcperdet>0.75] <- 1
  ad$fmc_75[is.na(ad$fmcperdet)] <- NA
  
  ad$fpc_25 <- 0
  ad$fpc_25[ad$fpcperdet>0.25] <- 1
  ad$fpc_25[is.na(ad$fpcperdet)] <- NA
  
  ad$fpc_50 <- 0
  ad$fpc_50[ad$fpcperdet>0.50] <- 1
  ad$fpc_50[is.na(ad$fpcperdet)] <- NA
  
  ad$fpc_75 <- 0
  ad$fpc_75[ad$fpcperdet>0.75] <- 1
  ad$fpc_75[is.na(ad$fpcperdet)] <- NA
  
  ad$fmc_25ref <- 0
  ad$fmc_25ref[ad$fmcperdet>0.25] <- 1
  ad$fmc_25ref[is.na(ad$fmcperdet)] <- NA
  
  ad$fmc_50ref <- NA
  ad$fmc_50ref[ad$fmcperdet>0.50] <- 1
  ad$fmc_50ref[ad$fmcperdet<0.25] <- 0
  ad$fmc_50ref[is.na(ad$fmcperdet)] <- NA
  
  ad$fmc_75ref <- NA
  ad$fmc_75ref[ad$fmcperdet>0.75] <- 1
  ad$fmc_75ref[ad$fmcperdet<0.25] <- 0
  ad$fmc_75ref[is.na(ad$fmcperdet)] <- NA
  
  ad$fpc_25ref <- 0
  ad$fpc_25ref[ad$fpcperdet>0.25] <- 1
  ad$fpc_25ref[is.na(ad$fpcperdet)] <- NA
  
  ad$fpc_50ref <- NA
  ad$fpc_50ref[ad$fpcperdet>0.50] <- 1
  ad$fpc_50ref[ad$fpcperdet<0.25] <- 0
  ad$fpc_50ref[is.na(ad$fpcperdet)] <- NA
  
  ad$fpc_75ref <- NA
  ad$fpc_75ref[ad$fpcperdet>0.75] <- 1
  ad$fpc_75ref[ad$fpcperdet<0.25] <- 0
  ad$fpc_75ref[is.na(ad$fpcperdet)] <- NA


  # risk variable
  ad$risk=NA
  ad$risk[ad$beach=="Fairhope"| ad$beach=="Goddard"|
              (ad$beach=="Doheny" & ad$berm=="Open")|
              (ad$beach=="Avalon" & ad$groundwater=="Above median flow")]=1
  ad$risk[ad$beach=="Malibu"|ad$beach=="Mission Bay"|
              (ad$beach=="Doheny" & ad$berm=="Closed")|
              (ad$beach=="Avalon" & ad$groundwater=="Below median flow")]=0
  ad$risk=factor(ad$risk)
  levels(ad$risk)=c("Low","High")
  
  # dropping malibu 1602 results
  ad$fmc1602.pres[ad$beach=="Malibu"]=NA
  ad$fpc1602.pres[ad$beach=="Malibu"]=NA
  
  ad$fmc1602[ad$beach=="Malibu"]=NA
  ad$fpc1602[ad$beach=="Malibu"]=NA
  
  # create a 3-level categorical variable of non-swimmers, swimmers when no coliphage, and swimmers when coliphage detected
  ad$swim.fmc <- ad$fmc.pres+1
  ad$swim.fmc[ad$anycontact=="No"] <- 0
  ad$swim.fmc <- factor(ad$swim.fmc,levels=c(0,1,2),labels=c("Non-swimmers","NoFMC","FMCpresent"))
  
  ad$swim.fpc <- ad$fpc.pres+1
  ad$swim.fpc[ad$anycontact=="No"] <- 0
  ad$swim.fpc <- factor(ad$swim.fpc,levels=c(0,1,2),labels=c("Non-swimmers","NoFPC","FPCpresent"))
  
  # these aren't right: 
  
  ad$swim.fmc1601 <- as.numeric(as.character(ad$fmc1601.pres))+1
  ad$swim.fmc1601[ad$anycontact=="No"] <- 0
  ad$swim.fmc1601 <- factor(ad$swim.fmc1601,levels=c(0,1,2),labels=c("Non-swimmers","NoFMC","FMCpresent"))
  
  ad$swim.fmc1602 <- as.numeric(as.character(ad$fmc1602.pres))+1
  ad$swim.fmc1602[ad$anycontact=="No"] <- 0
  ad$swim.fmc1602 <- factor(ad$swim.fmc1602,levels=c(0,1,2),labels=c("Non-swimmers","NoFMC","FMCpresent"))
  
  ad$swim.fpc1601 <- as.numeric(as.character(ad$fpc1601.pres))+1
  ad$swim.fpc1601[ad$anycontact=="No"] <- 0
  ad$swim.fpc1601 <- factor(ad$swim.fpc1601,levels=c(0,1,2),labels=c("Non-swimmers","NoFPC","FPCpresent"))
  
  ad$swim.fpc1602 <- as.numeric(as.character(ad$fpc1602.pres))+1
  ad$swim.fpc1602[ad$anycontact=="No"] <- 0
  ad$swim.fpc1602 <- factor(ad$swim.fpc1602,levels=c(0,1,2),labels=c("Non-swimmers","NoFPC","FPCpresent"))
  
  # create a 3-level categorical variable of non-swimmers, swimmers when entero<35, and swimmers when entero>35
  ad$swim.ent <- as.numeric(as.character(ad$entero35))+1
  ad$swim.ent[ad$anycontact=="No"] <- 0
  ad$swim.ent <- factor(ad$swim.ent,levels=c(0,1,2),labels=c("Non-swimmers","Entero<35","Entero>35"))
  
  # create a 4-level categorical variable of non-swimmers, coliphage, entero>35
  ad$swim.fmc.ent <- NA
  ad$swim.fmc.ent[ad$fmc.pres==1 & ad$entero35==0] <- 1
  ad$swim.fmc.ent[ad$fmc.pres==0 & ad$entero35==1] <- 2
  ad$swim.fmc.ent[ad$fmc.pres==1 & ad$entero35==1] <- 3  
  ad$swim.fmc.ent[ad$anycontact=="No"] <- 0
  ad$swim.fmc.ent <- factor(ad$swim.fmc.ent)
  
  ad$swim.fpc.ent <- NA
  ad$swim.fpc.ent[ad$fpc.pres==1 & ad$entero35==0] <- 1
  ad$swim.fpc.ent[ad$fpc.pres==0 & ad$entero35==1] <- 2
  ad$swim.fpc.ent[ad$fpc.pres==1 & ad$entero35==1] <- 3  
  ad$swim.fpc.ent[ad$anycontact=="No"] <- 0
  ad$swim.fpc.ent <- factor(ad$swim.fpc.ent)
  
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
# Function to calculate a marginally adjusted
# risk difference curve for continuous exposures
# among swimmers 
#-------------------------------------------
marginal.pY.risk <- function(form,data,nameX,pX) {
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
  
  # risk among swimmers who were not exposed to coliphage
  A0=mean(data$gici10[data[[paste(nameX,".pres",sep="")]]==0])
  
  # for each level of the exposure, get the marginal average predicted probability
  pY <- as.numeric(rep(NA,length(pX)))
  # and get risk difference compared to risk among unexposed
  rd <- as.numeric(rep(NA,length(pX)))
  pd <- data
  for (i in 1:length(pY)) {
    pd[,c(nameX)]<- pX[i]
    pY[i] <- mean(predict(fit,newdata=pd,type="response"))
    rd[i] <- pY[i]-A0
  }
  
  return(rd)
}
# 
# #-------------------------------------------
# # Function to calculate a marginally adjusted
# # risk difference curve for continuous exposures
# # with non-swimmers as the reference group
# #-------------------------------------------
# marginal.pY.risk.nonswimmer <- function(form,data,nameX,pX) {
#   # calculates a marginally adjusted dose-response curve for
#   # an exposure (e.g., entero1600, enteroQPCR), adjusted for covariates
#   #
#   # arguments:
#   # form: glm model formula, must include "entero1600" as a covariate
#   # data: data for glm estimation, including the var specified in nameX
#   # nameX : string of the name of the exposure variable: "entero1600" or "enteroQPCR"
#   # pX    : levels of the exposure for predicted probabilities
#   #
#   # returned objects:
#   # pY : vector of predicted probabilities at pX
#   
#   # fit the model
#   fit <- glm(form,family=poisson(link="log"),data=data)
#   
#   # risk among non-swimmers
#   A0=0
#   
#   # for each level of the exposure, get the marginal average predicted probability
#   pY <- as.numeric(rep(NA,length(pX)))
#   # and get risk difference compared to risk among unexposed
#   rd <- as.numeric(rep(NA,length(pX)))
#   pd <- data
#   for (i in 1:length(pY)) {
#     pd[,c(nameX)]<- pX[i]
#     pd$
#     pY[i] <- mean(predict(fit,newdata=pd,type="response"))
#     rd[i] <- pY[i]-A0
#   }
#   
#   return(rd)
# }

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



#-------------------------------------------
# Stratified bootstrap resampling function
# Samples n_i households from stratum i,
# where stratum in this study is beach cohort
# 
# This function is tailored to store 
# exposure-outcome dose response curves 
#-------------------------------------------

boot.pY.risk <- function(fmla,data,nameX,ID,iter,dots=TRUE) {
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
  rd <- matrix(NA,ncol=length(pX),nrow=iter)
  
  # Create the bootstrap samples of size n_i, within each stratum i
  IDlist <- unique(data$ID)
  bs <- matrix(sample(IDlist, length(IDlist)*iter, replace=TRUE), ncol=iter)
  
  # calculate P(Y|X,W) for each bootstrap sample
  for (bb in 1:iter) {   
    bd <- merge(data.frame(ID=bs[,bb]),data,by="ID",all.x=TRUE)
    rd[bb,] <- tryCatch( do.call(marginal.pY.risk,args=list(form=fmla,data=bd,nameX=nameX,pX=pX)), 
                         error=function(e) rep(NA,length(pX))) 
    if(dots) cat(".",sep="")
    if(dots && bb %% 50 == 0) cat(bb,"\n")
  }
  if(dots) cat("\n Bootstrap Run Time:",round(difftime(Sys.time(),start.time,units="mins"),3)," Minutes \n\n")
  
  # if there were any bootstrap samples where the model failed to converge 
  # for any reason, caught by tryCatch(), then report that.
  # also triage any estimates that are clearly from a failed model that sort-of converged
  # (due to model instability from sparse data in a bootstrap sample)
  rd[(rd[,1]>10 | rd[,1]< -10), ] <- NA
  Nna <- sum(ifelse(is.na(rd[,1])==TRUE,1,0))
  if(Nna>0) {
    cat("\n   Warning: ",Nna," bootstrap replicates failed to converge to sensible estimates \n")
    cat("\n   Bootstrap estimates are based on the remaining ",iter-Nna," replicates\n\n")
  }
  
  # calculate the mean and percentile 95% confidence intervals for the statistics
  bootest  <- apply(rd,2,function(x) mean(x,na.rm=TRUE))
  boot95lb <- apply(rd,2,function(x) quantile(x,prob=0.025,na.rm=TRUE))
  boot95ub <- apply(rd,2,function(x) quantile(x,prob=0.975,na.rm=TRUE))
  
  list(bootest=bootest,boot95lb=boot95lb,boot95ub=boot95ub,pX=pX,rd=rd)
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




