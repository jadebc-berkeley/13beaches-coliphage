##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file makes tables
##########################################


rm(list=ls())

setwd("~/Dropbox/Coliphage/")

# --------------------------------------
# paste n followed by percent in parentheses
# input numeric variables
# output string
# --------------------------------------
n.p.paren=function(n,p){
  paste(n," (",as.numeric(sprintf("%0.1f",p*100)),")",sep="")
}



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

all$beach=droplevels(all$beach)

# --------------------------------------
# table 1: demographics and water exposure by beach
# --------------------------------------
# individuals by beach
n.beach=table(all$beach)

# households by beach
all$count=1
hh.df=aggregate(all$count,list(hhid=all$hhid,beach=all$beach),sum)
hh.beach=table(hh.df$beach)

# age by beach
age.beach=table(as.character(all$agecat[all$agecat!="Missing"]),
  all$beach[all$agecat!="Missing"])
age.beach.p=prop.table(table(as.character(all$agecat[all$agecat!="Missing"]),
  all$beach[all$agecat!="Missing"]),2)
age.beach.tot=colSums(age.beach)
age.beach=rbind(age.beach,age.beach.tot)
age.beach.p=rbind(age.beach.p,rep(1,6))

# sex by beach
fem.beach=table(all$female,all$beach)
fem.beach.p=prop.table(table(all$female,all$beach),2)
fem.beach.tot=colSums(fem.beach)
fem.beach=rbind(fem.beach,fem.beach.tot)
fem.beach.p=rbind(fem.beach.p,rep(1,6))

# race by beach
white.beach=table(as.character(all$racewhite[all$racewhite!="Missing"]),
  all$beach[all$racewhite!="Missing"])
white.beach.p=prop.table(table(as.character(all$racewhite[all$racewhite!="Missing"]),
   all$beach[all$racewhite!="Missing"]),2)
white.beach.tot=colSums(white.beach)
white.beach=rbind(white.beach,white.beach.tot)
white.beach.p=rbind(white.beach.p,rep(1,6))

# water exposure by beach
any.beach=table(all$anycontact,all$beach)[2,]
any.beach.p=prop.table(table(all$anycontact,all$beach),2)[2,]

body.beach=table(all$bodycontact,all$beach)[2,]
body.beach.p=prop.table(table(all$bodycontact,all$beach),2)[2,]

swall.beach=table(all$swallwater,all$beach)[2,]
swall.beach.p=prop.table(table(all$swallwater,all$beach),2)[2,]

time=aggregate(all$watertime[all$bodycontact=="No"],
   list(all$beach[all$bodycontact=="No"]),mean,na.rm=TRUE)
time.beach=sprintf("%0.0f",time[,2]*60)
time.beach[time.beach=="NaN"]="--"

n.print=rbind(age.beach,fem.beach,white.beach,any.beach,body.beach,
              swall.beach)
p.print=rbind(age.beach.p,fem.beach.p,white.beach.p,any.beach.p,
              body.beach.p,swall.beach.p)

table.print=cbind(n.p.paren(n.print[,1],p.print[,1]),
                             n.p.paren(n.print[,2],p.print[,2]),
                             n.p.paren(n.print[,3],p.print[,3]),
                             n.p.paren(n.print[,4],p.print[,4]),
                             n.p.paren(n.print[,5],p.print[,5]),
                             n.p.paren(n.print[,6],p.print[,6]))

table.print=data.frame(rbind(as.character(n.beach),as.character(hh.beach),table.print,
                             time.beach))

table.print$label=c("Individuals","Households","0-4","5-14","15-24","25-34","35-44",
      "45-54","55-64","65-74","75+","Total","Male","Female","Total","Not white",
      "White","Total","Any contact","Body contact","Swallowed water","Minutes swam (mean)")

for(i in 3:22){
  table.print$label[i]=paste("~~~",table.print$label[i],sep="")
}
colnames(table.print)=c("Avalon","Doheny","Fairhope","Goddard","Malibu","Mission Bay",
                        "Label")

table.print=table.print[,c("Label","Avalon","Doheny","Fairhope","Goddard","Malibu",
                      "Mission Bay")]

# add category labels
agelab=c("Age (years)",rep(NA,6))
sexlab=c("Sex",rep(NA,6))
racelab=c("Race",rep(NA,6))
waterlab=c("Water exposure",rep(NA,6))
table.print=rbind(table.print[1:2,],agelab,table.print[3:12,],
                      sexlab,table.print[13:15,],racelab,table.print[16:18,],
                      waterlab,table.print[19:22,])

save(table.print,file="~/Dropbox/Coliphage/Results/Tables/desctab.RData")








