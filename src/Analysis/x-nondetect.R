##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# Percent of non-detects among individuals
##########################################

rm(list=ls())
library(foreign)

# --------------------------------------
# load the and pre-preprocess the 
# analysis dataset
# (refer to the base functions script
# for details on the pre-processing)
# --------------------------------------
beaches13=read.csv("~/Documents/CRG/coliphage/13beaches-data/final/13beaches-analysis.csv")

# load base functions
source("~/Documents/CRG/coliphage/13beaches-coliphage/src/Analysis/0-base-functions.R")

data=preprocess.6beaches(beaches13)

# restrict to 6 beaches with coliphage data
beach.list=c("Avalon","Doheny","Malibu","Mission Bay",
             "Fairhope","Goddard")

all=data[data$beach %in% beach.list,]

# drop individuals with no water quality information
all=subset(all,nowq==0)
# subset to non-missing exposure categories
# to make the robust CI calcs work
all=subset(all,all$bodycontact=="Yes")

myfn=function(data,beachn){
  n.fmc1601=nrow(data[data$beach==beachn & !is.na(data$fmc1601),])
  nd.fmc1601=nrow(data[data$beach==beachn & data$fmc1601==-1 & !is.na(data$fmc1601),])
  n.fmc1602=nrow(data[data$beach==beachn & !is.na(data$fmc1602),])
  nd.fmc1602=nrow(data[data$beach==beachn & data$fmc1602==-1 & !is.na(data$fmc1602),])
  n.fpc1601=nrow(data[data$beach==beachn & !is.na(data$fpc1601),])
  nd.fpc1601=nrow(data[data$beach==beachn & data$fpc1601==-1 & !is.na(data$fpc1601),])
  n.fpc1602=nrow(data[data$beach==beachn & !is.na(data$fpc1602),])
  nd.fpc1602=nrow(data[data$beach==beachn & data$fpc1602==-1 & !is.na(data$fpc1602),])
  
  n.entero=nrow(data[data$beach==beachn & !is.na(data$entero),])
  nd.entero=nrow(data[data$beach==beachn & data$entero==-1 & !is.na(data$entero),])
  df=data.frame(beach=rep(beachn,5),lab=c("Somatic coliphage EPA 1601",
                      "Somatic coliphage EPA 1602",
                      "Male-specific coliphage EPA 1601",
                      "Male-specific coliphage EPA 1602","Enterococcus"),
    n=c(n.fmc1601,n.fmc1602,n.fpc1601,n.fpc1602,n.entero),
                nd=c(nd.fmc1601,nd.fmc1602,nd.fpc1601,nd.fpc1602,nd.entero))
  
  df$pernd=(df$nd/df$n)*100
  return(df)
}


res=apply(as.matrix(c("Avalon","Doheny","Malibu","Mission Bay",
                  "Fairhope","Goddard")),1,myfn,data=all)

tab=rbind(res[[1]],res[[2]],res[[3]][c(3,5),],
          res[[4]][c(1,3,5),],res[[5]][c(3,5),],
          res[[6]][c(3,5),])

tab$pernd=sprintf("%0.1f",tab$pernd)

epa1601=grep("EPA 1601",tab$lab)
epa1602=grep("EPA 1602",tab$lab)
entero=grep("Enterococcus",tab$lab)

tab$test=""
tab$test[epa1601]="EPA 1601"
tab$test[epa1602]="EPA 1602"
tab$test[entero]="EPA 1600/Enterolert"

tab$lab=as.character(tab$lab)
tab$lab[tab$lab=="Somatic coliphage EPA 1601"]="Somatic coliphage"
tab$lab[tab$lab=="Somatic coliphage EPA 1602"]="Somatic coliphage"
tab$lab[tab$lab=="Male-specific coliphage EPA 1601"]="Male-specific coliphage"
tab$lab[tab$lab=="Male-specific coliphage EPA 1602"]="Male-specific coliphage"

tab=tab[,c("beach","lab","test","n","nd","pernd")]

save(tab,file="~/Documents/CRG/coliphage/Results/Tables/Table-non-detect-indiv.RData")




