########################################################
# Tables for beaches coliphage analysis
########################################################
rm(list=ls())
library(foreign)

########################################################
# Functions
########################################################
#-------------------------------------------------------
# Function for restricting data
#-------------------------------------------------------
process.df=function(df){
  df=subset(df,df[["swimtype"]]=="bodycontact")
  df=subset(df,df[["outcome"]]=="hcgi3ci3" | df$outcome=="hcgi3ci10")
  indicator=grep("sspda",df[["indicator"]])
  df=df[indicator,]
  df=df[!is.na(df[["or"]]),]
  return(df)
}

#-------------------------------------------------------
# Function for restricting data - range of water exposures
#-------------------------------------------------------
process.sdf=function(df){
  df=subset(df,df[["outcome"]]=="hcgi3ci3" | df$outcome=="hcgi3ci10")
  df=subset(df,df[["conditions"]]=="combined")
  indicator=grep("sspda",df[["indicator"]])
  df=df[indicator,]
  df=df[!is.na(df[["or"]]),]
  return(df)
}


#-------------------------------------------------------
# Function for pasting point ests with CIs
#-------------------------------------------------------
pretty.ci=function(pt,lb,ub,digits){
  myformat=paste("%0.",digits,"f",sep="")
#   pt=as.numeric(sprintf(myformat,pt))
#   lb=as.numeric(sprintf(myformat,lb))
#   ub=as.numeric(sprintf(myformat,ub))
  pt=sprintf(myformat,pt)
  lb=sprintf(myformat,lb)
  ub=sprintf(myformat,ub)
  paste(pt," (",lb,",",ub,")",sep="")
}

#-------------------------------------------------------
# Function to organize data into table
#-------------------------------------------------------
table.prep=function(df,cond1,cond2,cond3,beach,digits){
  #name 
  #name=as.character(df[["groupname"]][1])
  #name.start=regexpr("EPA",name)[1]
  #name.print=substr(name,start=name.start,stop=nchar(name))
  name.print=beach
  myformat=paste("%0.",digits,"f",sep="")
  
  #------------
  #3-day recall
  #------------
  #point est
  or.b.o.pt.3=df[["or"]][df[["conditions"]]==cond1
                         & df[["adj"]]=="adjusted" & df[["outcome"]]=="hcgi3ci3"]
  or.b.c.pt.3=df[["or"]][df[["conditions"]]==cond2
                         & df[["adj"]]=="adjusted"& df[["outcome"]]=="hcgi3ci3"]
  or.comb.pt.3=df[["or"]][df[["conditions"]]==cond3
                          & df[["adj"]]=="adjusted"& df[["outcome"]]=="hcgi3ci3"]
  
  #95% CI
  if(cond1 %in% df[["conditions"]]==TRUE){
    or.b.o.se.3=df[["se"]][df[["conditions"]]==cond1
                           & df[["adj"]]=="adjusted" & df[["outcome"]]=="hcgi3ci3"]
    or.b.o.lb.3=or.b.o.pt.3-(qnorm(0.975)*or.b.o.se.3)
    or.b.o.ub.3=or.b.o.pt.3+(qnorm(0.975)*or.b.o.se.3)
    or.b.o.3=pretty.ci(or.b.o.pt.3,or.b.o.lb.3,or.b.o.ub.3,digits)
  }else{
    or.b.o.3="NA"
  }
  
  if(cond2 %in% df[["conditions"]]==TRUE){
    or.b.c.se.3=df[["se"]][df[["conditions"]]==cond2
                           & df[["adj"]]=="adjusted" & df[["outcome"]]=="hcgi3ci3"] 
    or.b.c.lb.3=or.b.c.pt.3-(qnorm(0.975)*or.b.c.se.3)
    or.b.c.ub.3=or.b.c.pt.3+(qnorm(0.975)*or.b.c.se.3)
    or.b.c.3=pretty.ci(or.b.c.pt.3,or.b.c.lb.3,or.b.c.ub.3,digits)
  }else{
    or.b.c.3="NA"
  }
  
  if(cond3 %in% df[["conditions"]]==TRUE){
    or.comb.se.3=df[["se"]][df[["conditions"]]==cond3
                            & df[["adj"]]=="adjusted" & df[["outcome"]]=="hcgi3ci3"] 
    or.comb.lb.3=or.comb.pt.3-(qnorm(0.975)*or.comb.se.3)
    or.comb.ub.3=or.comb.pt.3+(qnorm(0.975)*or.comb.se.3)
    or.comb.3=pretty.ci(or.comb.pt.3,or.comb.lb.3,or.comb.ub.3,digits)
  }else{
    or.comb.3="NA"
  }  
  
  
  #------------
  #10-day recall
  #------------  
  or.b.o.pt.10=df[["or"]][df[["conditions"]]==cond1
                          & df[["adj"]]=="adjusted" & df[["outcome"]]=="hcgi3ci10"]
  or.b.c.pt.10=df[["or"]][df[["conditions"]]==cond2 
                          & df[["adj"]]=="adjusted"& df[["outcome"]]=="hcgi3ci10"]
  or.comb.pt.10=df[["or"]][df[["conditions"]]==cond3
                           & df[["adj"]]=="adjusted"& df[["outcome"]]=="hcgi3ci10"]
  
  #95% CI
  if(cond1 %in% df[["conditions"]]==TRUE){
    or.b.o.se.10=df[["se"]][df[["conditions"]]==cond1 
                            & df[["adj"]]=="adjusted" & df[["outcome"]]=="hcgi3ci10"] 
    or.b.o.lb.10=or.b.o.pt.10-(qnorm(0.975)*or.b.o.se.10)
    or.b.o.ub.10=or.b.o.pt.10+(qnorm(0.975)*or.b.o.se.10)
    or.b.o.10=pretty.ci(or.b.o.pt.10,or.b.o.lb.10,or.b.o.ub.10,digits)
  }else{
    or.b.o.10="NA"
  }  
  
  if(cond2 %in% df[["conditions"]]==TRUE){
    or.b.c.se.10=df[["se"]][df[["conditions"]]==cond2
                            & df[["adj"]]=="adjusted" & df[["outcome"]]=="hcgi3ci10"] 
    or.b.c.lb.10=or.b.c.pt.10-(qnorm(0.975)*or.b.c.se.10)
    or.b.c.ub.10=or.b.c.pt.10+(qnorm(0.975)*or.b.c.se.10)
    or.b.c.10=pretty.ci(or.b.c.pt.10,or.b.c.lb.10,or.b.c.ub.10,digits)
  }else{
    or.b.c.10="NA"
  }  
  
  if(cond3 %in% df[["conditions"]]==TRUE){
    or.comb.se.10=df[["se"]][df[["conditions"]]==cond3
                             & df[["adj"]]=="adjusted" & df[["outcome"]]=="hcgi3ci10"] 
    or.comb.lb.10=or.comb.pt.10-(qnorm(0.975)*or.comb.se.10)
    or.comb.ub.10=or.comb.pt.10+(qnorm(0.975)*or.comb.se.10)
    or.comb.10=pretty.ci(or.comb.pt.10,or.comb.lb.10,or.comb.ub.10,digits)
  }else{
    or.comb.10="NA"
  }
  
  out=c(name.print,or.b.o.3,or.b.c.3,or.comb.3,or.b.o.10,or.b.c.10,or.comb.10)
  return(out)
}

#-------------------------------------------------------
# Function to organize data into table - varying levels water exposure
#-------------------------------------------------------
table.prep.s=function(df,beach){
  name.print=beach
  
  #------------
  #3-day recall
  #------------
  #point est
  or.b.pt.3=df[["or"]][df[["swimtype"]]=="bodycontact"
                         & df[["adj"]]=="adjusted" & df[["outcome"]]=="hcgi3ci3"]
  or.h.pt.3=df[["or"]][df[["swimtype"]]=="headunder"
                         & df[["adj"]]=="adjusted"& df[["outcome"]]=="hcgi3ci3"]
  or.s.pt.3=df[["or"]][df[["swimtype"]]=="swallwater"
                          & df[["adj"]]=="adjusted"& df[["outcome"]]=="hcgi3ci3"]
  
  #
  n.b.3=df$n[df$swimtype=="bodycontact" & df$adj=="adjusted" & df$outcome=="hcgi3ci3"]
  n.h.3=df$n[df$swimtype=="headunder" & df$adj=="adjusted" & df$outcome=="hcgi3ci3"]
  n.s.3=df$n[df$swimtype=="swallwater" & df$adj=="adjusted" & df$outcome=="hcgi3ci3"]

  n.b.10=df$n[df$swimtype=="bodycontact" & df$adj=="adjusted" & df$outcome=="hcgi3ci10"]
  n.h.10=df$n[df$swimtype=="headunder" & df$adj=="adjusted" & df$outcome=="hcgi3ci10"]
  n.s.10=df$n[df$swimtype=="swallwater" & df$adj=="adjusted" & df$outcome=="hcgi3ci10"]
  
  
  #95% CI
    or.b.se.3=df[["se"]][df[["swimtype"]]=="bodycontact"
                           & df[["adj"]]=="adjusted" & df[["outcome"]]=="hcgi3ci3"] 
    or.b.lb.3=or.b.pt.3-(qnorm(0.975)*or.b.se.3)
    or.b.ub.3=or.b.pt.3+(qnorm(0.975)*or.b.se.3)
    or.b.3=pretty.ci(or.b.pt.3,or.b.lb.3,or.b.ub.3)

    or.h.se.3=df[["se"]][df[["swimtype"]]=="headunder"
                           & df[["adj"]]=="adjusted" & df[["outcome"]]=="hcgi3ci3"] 
    or.h.lb.3=or.h.pt.3-(qnorm(0.975)*or.h.se.3)
    or.h.ub.3=or.h.pt.3+(qnorm(0.975)*or.h.se.3)
    or.h.3=pretty.ci(or.h.pt.3,or.h.lb.3,or.h.ub.3)
  
    or.s.se.3=df[["se"]][df[["swimtype"]]=="swallwater"
                            & df[["adj"]]=="adjusted" & df[["outcome"]]=="hcgi3ci3"] 
    or.s.lb.3=or.s.pt.3-(qnorm(0.975)*or.s.se.3)
    or.s.ub.3=or.s.pt.3+(qnorm(0.975)*or.s.se.3)
    or.s.3=pretty.ci(or.s.pt.3,or.s.lb.3,or.s.ub.3)
  
  #------------
  #10-day recall
  #------------  
  or.b.pt.10=df[["or"]][df[["swimtype"]]=="bodycontact"
                          & df[["adj"]]=="adjusted" & df[["outcome"]]=="hcgi3ci10"]
  or.h.pt.10=df[["or"]][df[["swimtype"]]=="headunder" 
                          & df[["adj"]]=="adjusted"& df[["outcome"]]=="hcgi3ci10"]
  or.s.pt.10=df[["or"]][df[["swimtype"]]=="swallwater"
                           & df[["adj"]]=="adjusted"& df[["outcome"]]=="hcgi3ci10"]
  
  #95% CI
    or.b.se.10=df[["se"]][df[["swimtype"]]=="bodycontact" 
                            & df[["adj"]]=="adjusted" & df[["outcome"]]=="hcgi3ci10"] 
    or.b.lb.10=or.b.pt.10-(qnorm(0.975)*or.b.se.10)
    or.b.ub.10=or.b.pt.10+(qnorm(0.975)*or.b.se.10)
    or.b.10=pretty.ci(or.b.pt.10,or.b.lb.10,or.b.ub.10)
  
    or.h.se.10=df[["se"]][df[["swimtype"]]=="headunder"
                            & df[["adj"]]=="adjusted" & df[["outcome"]]=="hcgi3ci10"] 
    or.h.lb.10=or.h.pt.10-(qnorm(0.975)*or.h.se.10)
    or.h.ub.10=or.h.pt.10+(qnorm(0.975)*or.h.se.10)
    or.h.10=pretty.ci(or.h.pt.10,or.h.lb.10,or.h.ub.10)

    or.s.se.10=df[["se"]][df[["swimtype"]]=="swallwater"
                             & df[["adj"]]=="adjusted" & df[["outcome"]]=="hcgi3ci10"] 
    or.s.lb.10=or.s.pt.10-(qnorm(0.975)*or.s.se.10)
    or.s.ub.10=or.s.pt.10+(qnorm(0.975)*or.s.se.10)
    or.s.10=pretty.ci(or.s.pt.10,or.s.lb.10,or.s.ub.10)
  
  out=list(recall3=c(name.print,n.b.3,or.b.3,n.h.3,or.h.3,n.s.3,or.s.3),
        recall10=c(name.print,n.b.10,or.b.10,n.h.10,or.h.10,n.s.10,or.s.10))
  return(out)
}

########################################################
# DESCRIPTIVE TABLE
########################################################
avalon=read.dta("/Users/jadebc/Dropbox/Coliphage/Data/Untouched/avalon3.dta",convert.factors=FALSE)
doheny=read.dta("/Users/jadebc/Dropbox/Coliphage/Data/Untouched/doheny3.dta",convert.factors=FALSE)
malibu=read.dta("/Users/jadebc/Dropbox/Coliphage/Data/Untouched/malibu3.dta",convert.factors=FALSE)

#n
a.indiv=nrow(avalon)
d.indiv=nrow(doheny)
m.indiv=nrow(malibu)
indiv=cbind(a.indiv,d.indiv,m.indiv)

a.hh=length(unique(avalon$id))
d.hh=length(unique(doheny$id))
m.hh=length(unique(malibu$id))
hh=cbind(a.hh,d.hh,m.hh)

#age sex race body immersion 
a.age=as.matrix(prop.table(table(avalon$agecat1)))
d.age=as.matrix(prop.table(table(doheny$agecat1)))
m.age=as.matrix(prop.table(table(malibu$agecat1)))
age=cbind(a.age,d.age,m.age)

#removing person with missing sex info from sex denom
malibu$sex[malibu$sex==7]=NA
a.sex=as.matrix(prop.table(table(avalon$sex)))
d.sex=as.matrix(prop.table(table(doheny$sex)))
m.sex=as.matrix(prop.table(table(malibu$sex)))
sex=cbind(a.sex,d.sex,m.sex)

#race
a.race=as.matrix(prop.table(table(avalon$white)))
d.race=as.matrix(prop.table(table(doheny$white)))
m.race=as.matrix(prop.table(table(malibu$white)))
race=cbind(a.race,d.race,m.race)

#water contact
a.bodycontact=as.matrix(prop.table(table(avalon$bodycontact))[2])
d.bodycontact=as.matrix(prop.table(table(doheny$bodycontact))[2])
m.bodycontact=as.matrix(prop.table(table(malibu$bodycontact))[2])
bodycontact=cbind(a.bodycontact,d.bodycontact,m.bodycontact)

a.headunder=as.matrix(prop.table(table(avalon$headunder))[2])
d.headunder=as.matrix(prop.table(table(doheny$headunder))[2])
m.headunder=as.matrix(prop.table(table(malibu$headunder))[2])
headunder=cbind(a.headunder,d.headunder,m.headunder)

a.swall=as.matrix(prop.table(table(avalon$swallwater))[2])
d.swall=as.matrix(prop.table(table(doheny$swallwater))[2])
m.swall=as.matrix(prop.table(table(malibu$swallwater))[2])
swall=cbind(a.swall,d.swall,m.swall)

desctable=rbind(indiv,hh,age,sex,race,bodycontact, headunder,swall)
lab=c("Individuals","Households","0-5",
    "5.1-10","10.1-20","20.1-30","30.1-40","40.1-50",">50",
    "Male","Female","Not white","White","Body contact","Head under",
    "Swallowed water" )
for(i in 3:16){
  lab[i]=paste("~~~",lab[i],sep="")
}
desctable=cbind(lab,desctable)
rownames(desctable)=NULL
desctable.print=desctable
desctable.print[3:16,2:4]=sprintf("%0.1f",as.numeric(desctable.print[3:16,2:4])*100)
agelab=c("Age (years)",NA,NA,NA)
sexlab=c("Sex",NA,NA,NA)
racelab=c("Race",NA,NA,NA)
waterlab=c("Water exposure",NA,NA,NA)
desctable.print=rbind(desctable.print[1:2,],agelab,desctable.print[3:9,],
  sexlab,desctable.print[10:11,],racelab,desctable.print[12:13,],
  waterlab,desctable.print[14:16,])
rownames(desctable.print)=NULL
save(desctable.print,file="/Users/jadebc/Dropbox/PhD/Dissertation/Coliphage/Data/Temp/desctab.RData")

########################################################
# HEALTH RESULTS
########################################################
#all enrollees
table(avalon$hcgi3ci3)
prop.table(table(avalon$hcgi3ci3))
table(avalon$hcgi3ci10)
prop.table(table(avalon$hcgi3ci10))

table(doheny$hcgi3ci3)
prop.table(table(doheny$hcgi3ci3))
table(doheny$hcgi3ci10)
prop.table(table(doheny$hcgi3ci10))

table(malibu$hcgi3ci3)
prop.table(table(malibu$hcgi3ci3))
table(malibu$hcgi3ci10)
prop.table(table(malibu$hcgi3ci10))

#body contact
table(avalon$hcgi3ci3[avalon$bodycontact==1])
prop.table(table(avalon$hcgi3ci3[avalon$bodycontact==1]))
table(avalon$hcgi3ci10[avalon$bodycontact==1])
prop.table(table(avalon$hcgi3ci10[avalon$bodycontact==1]))

table(doheny$hcgi3ci3[doheny$bodycontact==1])
prop.table(table(doheny$hcgi3ci3[doheny$bodycontact==1]))
table(doheny$hcgi3ci10[doheny$bodycontact==1])
prop.table(table(doheny$hcgi3ci10[doheny$bodycontact==1]))

table(malibu$hcgi3ci3[malibu$bodycontact==1])
prop.table(table(malibu$hcgi3ci3[malibu$bodycontact==1]))
table(malibu$hcgi3ci10[malibu$bodycontact==1])
prop.table(table(malibu$hcgi3ci10[malibu$bodycontact==1]))

#head immersion
table(avalon$hcgi3ci3[avalon$headunder==1])
prop.table(table(avalon$hcgi3ci3[avalon$headunder==1]))
table(avalon$hcgi3ci10[avalon$headunder==1])
prop.table(table(avalon$hcgi3ci10[avalon$headunder==1]))

table(doheny$hcgi3ci3[doheny$headunder==1])
prop.table(table(doheny$hcgi3ci3[doheny$headunder==1]))
table(doheny$hcgi3ci10[doheny$headunder==1])
prop.table(table(doheny$hcgi3ci10[doheny$headunder==1]))

table(malibu$hcgi3ci3[malibu$headunder==1])
prop.table(table(malibu$hcgi3ci3[malibu$headunder==1]))
table(malibu$hcgi3ci10[malibu$headunder==1])
prop.table(table(malibu$hcgi3ci10[malibu$headunder==1]))

#swallowed water
table(avalon$hcgi3ci3[avalon$swallwater==1])
prop.table(table(avalon$hcgi3ci3[avalon$swallwater==1]))
table(avalon$hcgi3ci10[avalon$swallwater==1])
prop.table(table(avalon$hcgi3ci10[avalon$swallwater==1]))

table(doheny$hcgi3ci3[doheny$swallwater==1])
prop.table(table(doheny$hcgi3ci3[doheny$swallwater==1]))
table(doheny$hcgi3ci10[doheny$swallwater==1])
prop.table(table(doheny$hcgi3ci10[doheny$swallwater==1]))

table(malibu$hcgi3ci3[malibu$swallwater==1])
prop.table(table(malibu$hcgi3ci3[malibu$swallwater==1]))
table(malibu$hcgi3ci10[malibu$swallwater==1])
prop.table(table(malibu$hcgi3ci10[malibu$swallwater==1]))

########################################################
# DOHENY
########################################################

#-------------------------------------------------------
# Load data
#-------------------------------------------------------
d_14FMC0511=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/doheny/Doheny_Public_Indicator_Results-June-2013/doheny_14FMC0511_ind.csv",
                        header=TRUE)
# d_14FMC0512=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/doheny/Doheny_Public_Indicator_Results-June-2013/doheny_14FMC0512_ind.csv",
#                          header=TRUE)
d_14FPC0511=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/doheny/Doheny_Public_Indicator_Results-June-2013/doheny_14FPC0511_ind.csv",
                        header=TRUE)
# d_14FPC0512=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/doheny/Doheny_Public_Indicator_Results-June-2013/doheny_14FPC0512_ind.csv",
#                          header=TRUE)
# d_14FPC0611=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/doheny/Doheny_Public_Indicator_Results-June-2013/doheny_14FPC0611_ind.csv",
#                      header=TRUE)
# d_14FPC0612=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/doheny/Doheny_Public_Indicator_Results-June-2013/doheny_14FPC0612_ind.csv",
#                      header=TRUE)
# d_14FPC0711=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/doheny/Doheny_Public_Indicator_Results-June-2013/doheny_14FPC0711_ind.csv",
#                      header=TRUE)
# d_14FPC0712=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/doheny/Doheny_Public_Indicator_Results-June-2013/doheny_14FPC0712_ind.csv",
#                      header=TRUE)
# d_14FPC3412=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/doheny/Doheny_Public_Indicator_Results-June-2013/doheny_14FPC3412_ind.csv",
#                      header=TRUE)
# d_14FPC3511=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/doheny/Doheny_Public_Indicator_Results-June-2013/doheny_14FPC3511_ind.csv",
#                      header=TRUE)
# d_14FPC3512=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/doheny/Doheny_Public_Indicator_Results-June-2013/doheny_14FPC3512_ind.csv",
#                      header=TRUE)
# d_14FPC3513=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/doheny/Doheny_Public_Indicator_Results-June-2013/doheny_14FPC3513_ind.csv",
#                      header=TRUE)
# d_14FPC3514=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/doheny/Doheny_Public_Indicator_Results-June-2013/doheny_14FPC3514_ind.csv",
#                      header=TRUE)
# d_14FPC3515=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/doheny/Doheny_Public_Indicator_Results-June-2013/doheny_14FPC3515_ind.csv",
#                      header=TRUE)
d_16FPC0811=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/doheny/Doheny_Public_Indicator_Results-June-2013/doheny_16FPC0811_ind.csv",
                     header=TRUE)
d_16FMC0811=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/doheny/Doheny_Public_Indicator_Results-June-2013/doheny_16FMC0811_ind.csv",
                     header=TRUE)

#-------------------------------------------------------
# Process data
#-------------------------------------------------------
d_14FMC0511.p=process.df(d_14FMC0511)
# d_14FMC0512.p=process.df(d_14FMC0512)
d_14FPC0511.p=process.df(d_14FPC0511)
# d_14FPC0512.p=process.df(d_14FPC0512)
# d_14FPC0611.p=process.df(d_14FPC0611)
# d_14FPC0612.p=process.df(d_14FPC0612)
# d_14FPC0711.p=process.df(d_14FPC0711)
# d_14FPC0712.p=process.df(d_14FPC0712)
# d_14FPC3412.p=process.df(d_14FPC3412)
# d_14FPC3511.p=process.df(d_14FPC3511)
# d_14FPC3512.p=process.df(d_14FPC3512)
# d_14FPC3513.p=process.df(d_14FPC3513)
# d_14FPC3514.p=process.df(d_14FPC3514)
# d_14FPC3515.p=process.df(d_14FPC3515)
d_16FPC0811.p=process.df(d_16FPC0811)
d_16FMC0811.p=process.df(d_16FMC0811)

d_14FMC0511.s=process.sdf(d_14FMC0511)
d_14FPC0511.s=process.sdf(d_14FPC0511)
d_16FPC0811.s=process.sdf(d_16FPC0811)
d_16FMC0811.s=process.sdf(d_16FMC0811)

#-------------------------------------------------------
# Preparing table rows
#-------------------------------------------------------
d_14FMC0511.p.t=table.prep(d_14FMC0511.p,"combined","sand berm closed","sand berm open","Doheny beach",digits=2)
# d_14FMC0512.p.t=table.prep(d_14FMC0512.p,"combined","sand berm closed","sand berm open","Doheny beach")
d_14FPC0511.p.t=table.prep(d_14FPC0511.p,"combined","sand berm closed","sand berm open","Doheny beach",digits=2)
# d_14FPC0512.p.t=table.prep(d_14FPC0512.p,"combined","sand berm closed","sand berm open","Doheny beach")
# d_14FPC0611.p.t=table.prep(d_14FPC0611.p,"combined","sand berm closed","sand berm open","Doheny beach")
# d_14FPC0612.p.t=table.prep(d_14FPC0612.p,"combined","sand berm closed","sand berm open","Doheny beach")
# d_14FPC0711.p.t=table.prep(d_14FPC0711.p,"combined","sand berm closed","sand berm open","Doheny beach")
# d_14FPC0712.p.t=table.prep(d_14FPC0712.p,"combined","sand berm closed","sand berm open","Doheny beach")
# d_14FPC3412.p.t=table.prep(d_14FPC3412.p,"combined","sand berm closed","sand berm open","Doheny beach")
# d_14FPC3511.p.t=table.prep(d_14FPC3511.p,"combined","sand berm closed","sand berm open","Doheny beach")
# d_14FPC3512.p.t=table.prep(d_14FPC3512.p,"combined","sand berm closed","sand berm open","Doheny beach")
# d_14FPC3513.p.t=table.prep(d_14FPC3513.p,"combined","sand berm closed","sand berm open","Doheny beach")
# d_14FPC3514.p.t=table.prep(d_14FPC3514.p,"combined","sand berm closed","sand berm open","Doheny beach")
# d_14FPC3515.p.t=table.prep(d_14FPC3515.p,"combined","sand berm closed","sand berm open","Doheny beach")
d_16FPC0811.p.t=table.prep(d_16FPC0811.p,"combined","sand berm closed","sand berm open","Doheny beach",digits=2)
d_16FMC0811.p.t=table.prep(d_16FMC0811.p,"combined","sand berm closed","sand berm open","Doheny beach",digits=2)


d_14FMC0511.s.t=table.prep.s(d_14FMC0511.s,"Doheny beach")
d_14FPC0511.s.t=table.prep.s(d_14FPC0511.s,"Doheny beach")
d_16FPC0811.s.t=table.prep.s(d_16FPC0811.s,"Doheny beach")
d_16FMC0811.s.t=table.prep.s(d_16FMC0811.s,"Doheny beach")


########################################################
# AVALON
########################################################

#-------------------------------------------------------
# Load data
#-------------------------------------------------------

a_14FMC0511=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/Avalon/Avalon_Public_Indicator_Results-June-2013/avalon_14FMC0511_ind.csv",
                     header=TRUE)
a_14FPC0511=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/Avalon/Avalon_Public_Indicator_Results-June-2013/avalon_14FPC0511_ind.csv",
                     header=TRUE)
#a_14FPC0512=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/Avalon/Avalon_Public_Indicator_Results-June-2013/avalon_14FPC0512_ind.csv",
#                     header=TRUE)
# a_14FPC0611=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/Avalon/Avalon_Public_Indicator_Results-June-2013/avalon_14FPC0611_ind.csv",
#                      header=TRUE)
# a_14FPC0711=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/Avalon/Avalon_Public_Indicator_Results-June-2013/avalon_14FPC0711_ind.csv",
#                      header=TRUE)
# a_14FPC3412=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/Avalon/Avalon_Public_Indicator_Results-June-2013/avalon_14FPC3412_ind.csv",
#                      header=TRUE)
# a_14FPC3511=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/Avalon/Avalon_Public_Indicator_Results-June-2013/avalon_14FPC3511_ind.csv",
#                      header=TRUE)
# a_14FPC3512=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/Avalon/Avalon_Public_Indicator_Results-June-2013/avalon_14FPC3512_ind.csv",
#                      header=TRUE)
# a_14FPC3513=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/Avalon/Avalon_Public_Indicator_Results-June-2013/avalon_14FPC3513_ind.csv",
#                      header=TRUE)
# a_14FPC3514=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/Avalon/Avalon_Public_Indicator_Results-June-2013/avalon_14FPC3514_ind.csv",
#                      header=TRUE)
# a_14FPC3515=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/Avalon/Avalon_Public_Indicator_Results-June-2013/avalon_14FPC3515_ind.csv",
#                      header=TRUE)
a_16FPC0811=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/Avalon/Avalon_Public_Indicator_Results-June-2013/avalon_16FPC0811_ind.csv",
                     header=TRUE)
a_16FMC0811=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/Avalon/Avalon_Public_Indicator_Results-June-2013/avalon_16FMC0811_ind.csv",
                     header=TRUE)

#-------------------------------------------------------
# Process data
#-------------------------------------------------------
a_14FMC0511.p=process.df(a_14FMC0511)
a_14FPC0511.p=process.df(a_14FPC0511)
# a_14FPC0512.p=process.df(a_14FPC0512)
# a_14FPC0611.p=process.df(a_14FPC0611)
# a_14FPC0711.p=process.df(a_14FPC0711)
# a_14FPC3412.p=process.df(a_14FPC3412)
# a_14FPC3511.p=process.df(a_14FPC3511)
# a_14FPC3512.p=process.df(a_14FPC3512)
# a_14FPC3513.p=process.df(a_14FPC3513)
# a_14FPC3514.p=process.df(a_14FPC3514)
# a_14FPC3515.p=process.df(a_14FPC3515)
a_16FPC0811.p=process.df(a_16FPC0811)
a_16FMC0811.p=process.df(a_16FMC0811)

a_14FMC0511.s=process.sdf(a_14FMC0511)
a_14FPC0511.s=process.sdf(a_14FPC0511)
a_16FPC0811.s=process.sdf(a_16FPC0811)
a_16FMC0811.s=process.sdf(a_16FMC0811)

#-------------------------------------------------------
# Preparing table rows
#-------------------------------------------------------
a_14FMC0511.p.t=table.prep(a_14FMC0511.p,"combined","groundwater below median","groundwater above median","Avalon beach",digits=2)
a_14FPC0511.p.t=table.prep(a_14FPC0511.p,"combined","groundwater below median","groundwater above median","Avalon beach",digits=2)
#a_14FPC0512.p.t=table.prep(a_14FPC0512.p,"combined","groundwater below median","groundwater above median","Avalon beach")
# a_14FPC0611.p.t=table.prep(a_14FPC0611.p,"combined","groundwater below median","groundwater above median","Avalon beach")
# a_14FPC0711.p.t=table.prep(a_14FPC0711.p,"combined","groundwater below median","groundwater above median","Avalon beach")
# a_14FPC3412.p.t=table.prep(a_14FPC3412.p,"combined","groundwater below median","groundwater above median","Avalon beach")
# a_14FPC3511.p.t=table.prep(a_14FPC3511.p,"combined","groundwater below median","groundwater above median","Avalon beach")
# a_14FPC3512.p.t=table.prep(a_14FPC3512.p,"combined","groundwater below median","groundwater above median","Avalon beach")
# a_14FPC3513.p.t=table.prep(a_14FPC3513.p,"combined","groundwater below median","groundwater above median","Avalon beach")
# a_14FPC3514.p.t=table.prep(a_14FPC3514.p,"combined","groundwater below median","groundwater above median","Avalon beach")
# a_14FPC3515.p.t=table.prep(a_14FPC3515.p,"combined","groundwater below median","groundwater above median","Avalon beach")
a_16FPC0811.p.t=table.prep(a_16FPC0811.p,"combined","groundwater below median","groundwater above median","Avalon beach",digits=2)
a_16FMC0811.p.t=table.prep(a_16FMC0811.p,"combined","groundwater below median","groundwater above median","Avalon beach",digits=2)

a_14FMC0511.s.t=table.prep.s(a_14FMC0511.s,"Avalon beach")
a_14FPC0511.s.t=table.prep.s(a_14FPC0511.s,"Avalon beach")
a_16FPC0811.s.t=table.prep.s(a_16FPC0811.s,"Avalon beach")
a_16FMC0811.s.t=table.prep.s(a_16FMC0811.s,"Avalon beach")


########################################################
# MALIBU
########################################################
#-------------------------------------------------------
# Load data
#-------------------------------------------------------
m_14FPC0511=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/malibu/Malibu_Public_Results_Indicators-June-2013/malibu_14FPC0511_ind.csv",
                    header=TRUE)
# m_14FPC0512=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/malibu/Malibu_Public_Results_Indicators-June-2013/malibu_14FPC0512_ind.csv",
#                      header=TRUE)
# m_14FPC3412=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/malibu/Malibu_Public_Results_Indicators-June-2013/malibu_14FPC3412_ind.csv",
#                      header=TRUE)
# m_14FPC3511=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/malibu/Malibu_Public_Results_Indicators-June-2013/malibu_14FPC3511_ind.csv",
#                      header=TRUE)
# m_14FPC3512=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/malibu/Malibu_Public_Results_Indicators-June-2013/malibu_14FPC3512_ind.csv",
#                      header=TRUE)
# m_14FPC3513=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/malibu/Malibu_Public_Results_Indicators-June-2013/malibu_14FPC3513_ind.csv",
#                      header=TRUE)
# m_14FPC3514=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/malibu/Malibu_Public_Results_Indicators-June-2013/malibu_14FPC3514_ind.csv",
#                      header=TRUE)
# m_14FPC3515=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/malibu/Malibu_Public_Results_Indicators-June-2013/malibu_14FPC3515_ind.csv",
#                      header=TRUE)
m_16FMC0811=read.csv("/Users/jadebc/Dropbox/Coliphage/0_From Ben/Results/Beaches-Public_Results/indicator-analysis-results/malibu/Malibu_Public_Results_Indicators-June-2013/malibu_16FMC0811_ind.csv",
                     header=TRUE)

#-------------------------------------------------------
# Process data
#-------------------------------------------------------
m_14FPC0511.p=process.df(m_14FPC0511)
# m_14FPC0512.p=process.df(m_14FPC0512)
# m_14FPC3412.p=process.df(m_14FPC3412)
# m_14FPC3511.p=process.df(m_14FPC3511)
# m_14FPC3512.p=process.df(m_14FPC3512)
# m_14FPC3513.p=process.df(m_14FPC3513)
# m_14FPC3514.p=process.df(m_14FPC3514)
# m_14FPC3515.p=process.df(m_14FPC3515)
m_16FMC0811.p=process.df(m_16FMC0811)

m_14FPC0511.s=process.sdf(m_14FPC0511)
m_16FMC0811.s=process.sdf(m_16FMC0811)


#-------------------------------------------------------
# Preparing table rows
#-------------------------------------------------------
m_14FPC0511.p.t=table.prep(m_14FPC0511.p,"combined","sand berm closed","sand berm open","Malibu beach",digits=2)
# m_14FPC0512.p.t=table.prep(m_14FPC0512.p,"combined","sand berm closed","sand berm open","Malibu beach")
# m_14FPC3412.p.t=table.prep(m_14FPC3412.p,"combined","sand berm closed","sand berm open","Malibu beach")
# m_14FPC3511.p.t=table.prep(m_14FPC3511.p,"combined","sand berm closed","sand berm open","Malibu beach")
# m_14FPC3512.p.t=table.prep(m_14FPC3512.p,"combined","sand berm closed","sand berm open","Malibu beach")
# m_14FPC3513.p.t=table.prep(m_14FPC3513.p,"combined","sand berm closed","sand berm open","Malibu beach")
# m_14FPC3514.p.t=table.prep(m_14FPC3514.p,"combined","sand berm closed","sand berm open","Malibu beach")
# m_14FPC3515.p.t=table.prep(m_14FPC3515.p,"combined","sand berm closed","sand berm open","Malibu beach")
m_16FMC0811.p.t=table.prep(m_16FMC0811.p,"combined","sand berm closed","sand berm open","Malibu beach",digits=2)

m_14FPC0511.s.t=table.prep.s(m_14FPC0511.s,"Malibu beach")
m_16FMC0811.s.t=table.prep.s(m_16FMC0811.s,"Malibu beach")

#-------------------------------------------------------
# Organizing data into table - by condition
#-------------------------------------------------------
assoc.tab.1601fmc=rbind(a_14FMC0511.p.t,d_14FMC0511.p.t)
assoc.tab.1602fmc=rbind(a_16FMC0811.p.t,d_16FMC0811.p.t,m_16FMC0811.p.t)
assoc.tab.1601fpc=rbind(a_14FPC0511.p.t,d_14FPC0511.p.t,m_14FPC0511.p.t)
assoc.tab.1602fpc=rbind(a_16FPC0811.p.t,d_16FPC0811.p.t)

rownames(assoc.tab.1601fmc)=NULL
rownames(assoc.tab.1602fmc)=NULL
rownames(assoc.tab.1601fpc)=NULL
rownames(assoc.tab.1602fpc)=NULL

save(assoc.tab.1601fmc,file="/Users/jadebc/Dropbox/Coliphage/Data/Temp/assoctab_1601fmc.RData")
save(assoc.tab.1602fmc,file="/Users/jadebc/Dropbox/Coliphage/Data/Temp/assoctab_1602fmc.RData")
save(assoc.tab.1601fpc,file="/Users/jadebc/Dropbox/Coliphage/Data/Temp/assoctab_1601fpc.RData")
save(assoc.tab.1602fpc,file="/Users/jadebc/Dropbox/Coliphage/Data/Temp/assoctab_1602fpc.RData")

#-------------------------------------------------------
# Organizing data into table - by water exposure
#-------------------------------------------------------
wat.or.1601fmc.3=rbind(a_14FMC0511.s.t$recall3,d_14FMC0511.s.t$recall3)
wat.or.1601fmc.10=rbind(a_14FMC0511.s.t$recall10,d_14FMC0511.s.t$recall10)
wat.or.1602fmc.3=rbind(a_16FMC0811.s.t$recall3,d_16FMC0811.s.t$recall3,m_16FMC0811.s.t$recall3)
wat.or.1602fmc.10=rbind(a_16FMC0811.s.t$recall10,d_16FMC0811.s.t$recall10,m_16FMC0811.s.t$recall10)
wat.or.1601fpc.3=rbind(a_14FPC0511.s.t$recall3,d_14FPC0511.s.t$recall3,m_14FPC0511.s.t$recall3)
wat.or.1601fpc.10=rbind(a_14FPC0511.s.t$recall10,d_14FPC0511.s.t$recall10,m_14FPC0511.s.t$recall10)
wat.or.1602fpc.3=rbind(a_16FPC0811.s.t$recall3,d_16FPC0811.s.t$recall3)
wat.or.1602fpc.10=rbind(a_16FPC0811.s.t$recall10,d_16FPC0811.s.t$recall10)

save(wat.or.1601fmc.3,file="/Users/jadebc/Dropbox/Coliphage/Data/Temp/wator_1601fmc3.RData")
save(wat.or.1601fmc.10,file="/Users/jadebc/Dropbox/Coliphage/Data/Temp/wator_1601fmc10.RData")
save(wat.or.1602fmc.3,file="/Users/jadebc/Dropbox/Coliphage/Data/Temp/wator_1602fmc3.RData")
save(wat.or.1602fmc.10,file="/Users/jadebc/Dropbox/Coliphage/Data/Temp/wator_1602fmc10.RData")
save(wat.or.1601fpc.3,file="/Users/jadebc/Dropbox/Coliphage/Data/Temp/wator_1601fpc3.RData")
save(wat.or.1601fpc.10,file="/Users/jadebc/Dropbox/Coliphage/Data/Temp/wator_1601fpc10.RData")
save(wat.or.1602fpc.3,file="/Users/jadebc/Dropbox/Coliphage/Data/Temp/wator_1602fpc3.RData")
save(wat.or.1602fpc.10,file="/Users/jadebc/Dropbox/Coliphage/Data/Temp/wator_1602fpc10.RData")


########################################################
# INDICATORS
########################################################
avalon.ind=read.dta("~/Dropbox/Coliphage/Data/Untouched/avalon_inddata.dta")
doheny.ind=read.dta("~/Dropbox/Coliphage/Data/Untouched/doheny_inddata.dta")
malibu.ind=read.dta("~/Dropbox/Coliphage/Data/Untouched/malibu_inddata.dta")

#number of samples per site
length(unique(avalon.ind$sampleid))
length(unique(doheny.ind$sampleid))
length(unique(malibu.ind$sampleid))

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

all.inds.dl=subset(all.inds,all.inds$log10>0)

#-------------------------------------------------------
# Summary of log10 concentration
#-------------------------------------------------------
gm_mean = function(a){
  if(length(a)==0){
    NA
  }else{
    prod(a)^(1/length(a))
  }
}

ind.table=function(x,qual){
  #min max
  min=sprintf("%0.1f",quantile(x[x>=0],probs=c(0,1),na.rm=TRUE)[1])
  max=sprintf("%0.0f",quantile(x[x>=0],probs=c(0,1),na.rm=TRUE)[2])
  #geometric mean
  gm=sprintf("%0.1f",gm_mean(x[x>=0]))
  #non-detects
  #nd=length(x[x<0]) - this is if using log10 
  nd=sprintf("%0.0f",table(qual)[1])
  #n
  #n=length(x)
  n=table(is.na(x))[1]
  if(n==nd){
   min="--"
   max="--"
   gm="--"
  }
  out=data.frame(n=paste(n),min=min,max=max,gm=gm,nd=nd)
  rownames(out)=NULL
  return(out)
}

a.1601fmc.tab=ind.table(all.inds$result[all.inds$label=="EPA 1601 F- Coliphage" & all.inds$beach=="Avalon"],
                        all.inds$qualifier[all.inds$label=="EPA 1601 F- Coliphage" & all.inds$beach=="Avalon"])
a.1602fmc.tab=ind.table(all.inds$result[all.inds$label=="EPA 1602 F- Coliphage" & all.inds$beach=="Avalon"],
                        all.inds$qualifier[all.inds$label=="EPA 1602 F- Coliphage" & all.inds$beach=="Avalon"])
a.1601fpc.tab=ind.table(all.inds$result[all.inds$label=="EPA 1601 F+ Coliphage" & all.inds$beach=="Avalon"],
                        all.inds$qualifier[all.inds$label=="EPA 1601 F+ Coliphage" & all.inds$beach=="Avalon"])
a.1602fpc.tab=ind.table(all.inds$result[all.inds$label=="EPA 1602 F+ Coliphage" & all.inds$beach=="Avalon"],
                        all.inds$qualifier[all.inds$label=="EPA 1602 F+ Coliphage" & all.inds$beach=="Avalon"])

d.1601fmc.tab=ind.table(all.inds$result[all.inds$label=="EPA 1601 F- Coliphage" & all.inds$beach=="Doheny"],
                        all.inds$qualifier[all.inds$label=="EPA 1601 F- Coliphage" & all.inds$beach=="Doheny"])
d.1602fmc.tab=ind.table(all.inds$result[all.inds$label=="EPA 1602 F- Coliphage" & all.inds$beach=="Doheny"],
                        all.inds$qualifier[all.inds$label=="EPA 1602 F- Coliphage" & all.inds$beach=="Doheny"])
d.1601fpc.tab=ind.table(all.inds$result[all.inds$label=="EPA 1601 F+ Coliphage" & all.inds$beach=="Doheny"],
                        all.inds$qualifier[all.inds$label=="EPA 1601 F+ Coliphage" & all.inds$beach=="Doheny"])
d.1602fpc.tab=ind.table(all.inds$result[all.inds$label=="EPA 1602 F+ Coliphage" & all.inds$beach=="Doheny"],
                        all.inds$qualifier[all.inds$label=="EPA 1602 F+ Coliphage" & all.inds$beach=="Doheny"])

#1601fmc not measured in malibu
m.1601fmc.tab=data.frame(n="--",min="--",max="--",gm="--",nd="--")
m.1602fmc.tab=ind.table(all.inds$result[all.inds$label=="EPA 1602 F- Coliphage" & all.inds$beach=="Malibu"],
                        all.inds$result[all.inds$label=="EPA 1602 F- Coliphage" & all.inds$beach=="Malibu"])
m.1601fpc.tab=ind.table(all.inds$result[all.inds$label=="EPA 1601 F+ Coliphage" & all.inds$beach=="Malibu"],
                        all.inds$qualifier[all.inds$label=="EPA 1601 F+ Coliphage" & all.inds$beach=="Malibu"])
m.1602fpc.tab=ind.table(all.inds$result[all.inds$label=="EPA 1602 F+ Coliphage" & all.inds$beach=="Malibu"],
                        all.inds$qualifier[all.inds$label=="EPA 1602 F+ Coliphage" & all.inds$beach=="Malibu"])

ind.table=rbind(a.1601fmc.tab,d.1601fmc.tab,m.1601fmc.tab,
                a.1602fmc.tab,d.1602fmc.tab,m.1602fmc.tab,
                a.1601fpc.tab,d.1601fpc.tab,m.1601fpc.tab,
                a.1602fpc.tab,d.1602fpc.tab,m.1602fpc.tab)
ind.table$lab=rep(c("~~~Avalon","~~~Doheny","~~~Malibu"),4)
ind.table=rbind(rep(NA,ncol(ind.table)),ind.table[1:3,],
                rep(NA,ncol(ind.table)),ind.table[4:6,], 
                rep(NA,ncol(ind.table)),ind.table[7:9,],
                rep(NA,ncol(ind.table)),ind.table[10:12,])
ind.table$lab[1]="EPA 1601 F- Coliphage"
ind.table$lab[5]="EPA 1602 F- Coliphage"
ind.table$lab[9]="EPA 1601 F+ Coliphage"
ind.table$lab[13]="EPA 1602 F+ Coliphage"
ind.table$lab[4]="~~~Malibu$\\ast$"
ind.table=ind.table[,c(6,1:5)]

save(ind.table,file="/Users/jadebc/Dropbox/Coliphage/Data/Temp/indtable.RData")



########################################################
# ASSOCIATION BETWEEN COLIPHAGE AND ENTEROCOCCUS
########################################################

#need to reorganize data so unique obs are indexed by date-time
all.inds$timeid=paste(all.inds$sampledate,all.inds$sampletime,sep=" ")

#data frame with only the info needed
all.inds.simp=data.frame(sampleid=all.inds$sampleid,
                         timeid=all.inds$timeid,log10=all.inds$log10,
                         beach=all.inds$beach,label=all.inds$label)

all.inds.simp=all.inds.simp[order(all.inds.simp$timeid,all.inds.simp$sampleid),]

#averaging over lab replicates within date, time site
all.inds.av=aggregate(list(log10=all.inds.simp$log10),list(timeid=all.inds.simp$timeid,
                                                           sampleid=all.inds.simp$sampleid,beach=all.inds.simp$beach,
                                                           label=all.inds.simp$label),gm_mean)

#reshaping long to wide
all.inds.wide=reshape(all.inds.av,direction="wide",
                      timevar="label",idvar=c("timeid","sampleid","beach"))

all.inds.wide=all.inds.wide[order(all.inds.wide$sampleid),]

colnames(all.inds.wide)[4:8]=c("ent","fmc1601","fpc1601","fmc1602","fpc1602")

#How often was coliphage detected when water exceeded standards (Enterococcus 104 CFU/100 ml)?
all.inds.wide$entcut=cut(all.inds.wide$ent,log10(104))
levels(all.inds.wide$entcut)=c("0","1")

#fmc1601 - number of non-detects when ent criteria failed/met
table(is.na(all.inds.wide$fmc1601[all.inds.wide$entcut=="0"]))
table(is.na(all.inds.wide$fmc1601[all.inds.wide$entcut=="1"]))

#fmc1602 - number of non-detects when ent criteria failed/met - CHECK DTL FOR MALIBU
table(is.na(all.inds.wide$fmc1602[all.inds.wide$entcut=="0"]))
table(is.na(all.inds.wide$fmc1602[all.inds.wide$entcut=="1"]))
prop.table(table(all.inds.wide$fmc1602[all.inds.wide$entcut=="0"]>-1))
prop.table(table(all.inds.wide$fmc1602[all.inds.wide$entcut=="1"]>-1))

#fpc1601 - number of non-detects when ent criteria failed/met - CHECK DTL FOR MALIBU
table(is.na(all.inds.wide$fpc1601[all.inds.wide$entcut=="0"]))
table(is.na(all.inds.wide$fpc1601[all.inds.wide$entcut=="1"]))
prop.table(table(all.inds.wide$fpc1601[all.inds.wide$entcut=="0"]>-1))
prop.table(table(all.inds.wide$fpc1601[all.inds.wide$entcut=="1"]>-1))

#fmc1602 - number of non-detects when ent criteria failed/met - CHECK DTL FOR MALIBU
table(is.na(all.inds.wide$fpc1602[all.inds.wide$entcut=="0"]))
table(is.na(all.inds.wide$fpc1602[all.inds.wide$entcut=="1"]))
prop.table(table(all.inds.wide$fpc1602[all.inds.wide$entcut=="0"]>-1))
prop.table(table(all.inds.wide$fpc1602[all.inds.wide$entcut=="1"]>-1))

#remove missing values and non-detects for smooth splines
ent.fmc1601=all.inds.wide[,c(1:5)]
ent.fmc1601$ent[ent.fmc1601$ent<0]=NA
ent.fmc1601$fmc1601[ent.fmc1601$fmc1601<0]=NA
ent.fmc1601=ent.fmc1601[complete.cases(ent.fmc1601),]

ent.fmc1602=all.inds.wide[,c(1:4,7)]
ent.fmc1602$ent[ent.fmc1602$ent<0]=NA
ent.fmc1602$fmc1602[ent.fmc1602$fmc1602<0]=NA
ent.fmc1602=ent.fmc1602[complete.cases(ent.fmc1602),]

ent.fpc1601=all.inds.wide[,c(1:4,6)]
ent.fpc1601$ent[ent.fpc1601$ent<0]=NA
ent.fpc1601$fpc1601[ent.fpc1601$fpc1601<0]=NA
ent.fpc1601=ent.fpc1601[complete.cases(ent.fpc1601),]

ent.fpc1602=all.inds.wide[,c(1:4,8)]
ent.fpc1602$ent[ent.fpc1602$ent<0]=NA
ent.fpc1602$fpc1602[ent.fpc1602$fpc1602<0]=NA
ent.fpc1602=ent.fpc1602[complete.cases(ent.fpc1602),]

#smooth splines
#fmc1601.ss=smooth.spline(ent.fmc1601$ent,ent.fmc1601$fmc1601,df=7)
pdf(file="/Users/jadebc/Dropbox/Coliphage/Results/Figures/scatter-ent-fmc1601.pdf",
    onefile=TRUE,width=10,height=6)
par(mai=c(1.2,1.25,0.75,0.75))
plot(ent.fmc1601$ent,ent.fmc1601$fmc1601,pch=21,
     ylab="Log10 EPA 1601 F- Coliphage \nconcentration (PFU/100 ml)",
     xlab="Log10 Enterococcus concentration (CFU/100 ml)")
lines(fmc1601.ss,col="black",lwd=2)
abline(v=log10(104),col="red",lwd=2)
legend("topleft",lty=c(0,1,1),pch=c(16,NA,NA),lwd=c(0,2,2),
       col=c("#0000FF60","black","red"),
       legend=c("Observed concentration","Smooth spline","Enterococcus cutoff"))
dev.off()

pdf(file="/Users/jadebc/Dropbox/Coliphage/Results/Figures/scatter-ent-fpc1601.pdf",
    onefile=TRUE,width=10,height=6)
par(mai=c(1.2,1.25,0.75,0.75))
fpc1601.ss=smooth.spline(ent.fpc1601$ent,ent.fpc1601$fpc1601,df=7)
plot(ent.fpc1601$ent,ent.fpc1601$fpc1601,pch=16,col="#0000FF30",
     ylab=expression(Log[10]~"EPA 1601 F+ Coliphage \nconcentration (PFU/100 ml)"),
     xlab=expression(Log[10]~"Enterococcus concentration (CFU/100 ml)"))
lines(fpc1601.ss,col="black",lwd=2)
dev.off()
bquote(Log[10] ~ "EPA 1601 F+ Coliphage concentration (PFU/100 ml)")

pdf(file="/Users/jadebc/Dropbox/Coliphage/Results/Figures/scatter-ent-fmc1602.pdf",
    onefile=TRUE,width=10,height=6)
par(mai=c(1.2,1.25,0.75,0.75))
fmc1602.ss=smooth.spline(ent.fmc1602$ent,ent.fmc1602$fmc1602,df=7)
plot(ent.fmc1602$ent,ent.fmc1602$fmc1602,pch=16,col="#0000FF30",
    ylab=expression(Log[10]~"EPA 1602 F- Coliphage \nconcentration (PFU/100 ml)"),
    xlab=expression(Log[10]~"Enterococcus concentration (CFU/100 ml)"))
lines(fmc1602.ss,col="black",lwd=2)
dev.off()

pdf(file="/Users/jadebc/Dropbox/Coliphage/Results/Figures/scatter-ent-fpc1602.pdf",
    onefile=TRUE,width=10,height=6)
par(mai=c(1.2,1.25,0.75,0.75))
fpc1602.ss=smooth.spline(ent.fpc1602$ent,ent.fpc1602$fpc1602)
plot(ent.fpc1602$ent,ent.fpc1602$fpc1602,pch=16,col="#0000FF30",
     ylab=expression(Log[10]~"EPA 1602 F+ Coliphage \nconcentration (PFU/100 ml)"),
     xlab=expression(Log[10]~"Enterococcus concentration (CFU/100 ml)"))
lines(fpc1601.ss,col="black",lwd=2)
dev.off()

# confidence bands for splines --------------------------------------------
#resample data points
sp.resampler <- function() {
  n <- nrow(sp.frame)
  resample.rows <- sample(1:n,size=n,replace=TRUE)
  return(sp.frame[resample.rows,])
}

#estimate splines on resampled points
sp.spline.estimator <- function(data,m=300) {
  # Fit spline to data, set df
  fit <- smooth.spline(x=data[,1],y=data[,2],cv=FALSE)
  #fit <- smooth.spline(x=data[,1],y=data[,2],df=df)
  # Set up a grid of m evenly-spaced points on which to evaluate the spline
  eval.grid <- seq(from=min(data[,1]),to=max(data[,1]),length.out=m)
  # Slightly inefficient to re-define the same grid every time we call this,
  # but not a big overhead
  # Do the prediction and return the predicted values
  return(predict(fit,x=eval.grid)$y) # We only want the predicted values
}

#confidence bands using quantiles
sp.spline.cis <- function(B,alpha,m=300,data) {
  spline.main <- sp.spline.estimator(data,m=m)
  # Draw B boottrap samples, fit the spline to each
  spline.boots <- replicate(B,sp.spline.estimator(sp.resampler(),m=m))
  # Result has m rows and B columns
  cis.lower <- 2*spline.main - apply(spline.boots,1,quantile,probs=1-alpha/2)
  cis.upper <- 2*spline.main - apply(spline.boots,1,quantile,probs=alpha/2)
  return(list(main.curve=spline.main,lower.ci=cis.lower,upper.ci=cis.upper,
              x=seq(from=min(data[,1]),to=max(data[,1]),length.out=m)))
}


#conf band for F+ 1601
sp.frame.fpc1601 <- data.frame(ent=ent.fpc1601$ent,col=ent.fpc1601$fpc1601)
set.seed(98)
sp.cis.fpc1601=sp.spline.cis(B=1000,alpha=0.05,data=sp.frame.fpc1601)

plot(ent.fpc1601$ent,ent.fpc1601$fpc1601,pch=16,col="#0000FF30",
     ylab=expression(Log[10]~"EPA 1601 F+ Coliphage \nconcentration (PFU/100 ml)"),
     xlab=expression(Log[10]~"Enterococcus concentration (CFU/100 ml)"),
     ylim=c(-10,10))
lines(x=sp.cis.fpc1601$x,y=sp.cis.fpc1601$main.curve)
lines(x=sp.cis.fpc1601$x,y=sp.cis.fpc1601$lower.ci,lty=2)
lines(x=sp.cis.fpc1601$x,y=sp.cis.fpc1601$upper.ci,lty=2)

#conf band for F+ 1602
sp.frame.fpc1602 <- data.frame(ent=ent.fpc1602$ent,col=ent.fpc1602$fpc1602)
set.seed(98)
sp.cis.fpc1602=sp.spline.cis(B=1000,alpha=0.05,data=sp.frame.fpc1602)

plot(ent.fpc1602$ent,ent.fpc1602$fpc1602,pch=16,col="#0000FF30",
     ylab=expression(Log[10]~"EPA 1601 F+ Coliphage \nconcentration (PFU/100 ml)"),
     xlab=expression(Log[10]~"Enterococcus concentration (CFU/100 ml)"),
     ylim=c(-10,10))
lines(x=sp.cis.fpc1602$x,y=sp.cis.fpc1602$main.curve)
lines(x=sp.cis.fpc1602$x,y=sp.cis.fpc1602$lower.ci,lty=2)
lines(x=sp.cis.fpc1602$x,y=sp.cis.fpc1602$upper.ci,lty=2)


#box plots --------------------------------------------
all.inds.wide$entcut=NA
all.inds.wide$entcut[all.inds.wide$ent>=log10(104)]=1
all.inds.wide$entcut[all.inds.wide$ent<log10(104)]=0
all.inds.wide$entcut=as.factor(all.inds.wide$entcut)

#reshape for ggplot
all.inds.bx=melt(all.inds.wide, id.vars=c("timeid","sampleid","beach","entcut"))

#remove non-detects and missings
all.inds.bx$value[all.inds.bx$value<0]=NA
all.inds.bx=subset(all.inds.bx,!is.na(all.inds.bx$value))
all.inds.bx=subset(all.inds.bx,!is.na(all.inds.bx$entcut))
all.inds.bx=subset(all.inds.bx,all.inds.bx$variable!="ent")

ggplot(all.inds.bx,aes(x=variable,y=value,fill=entcut))+geom_boxplot()+
  coord_flip()



