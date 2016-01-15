########################################################
# TABLES FOR COLIPHAGE PAPER
########################################################
rm(list=ls())
library(foreign)

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
# Function for pasting point ests with CIs
#-------------------------------------------------------
pretty.ci=function(pt,lb,ub){
  pt=as.numeric(sprintf("%0.2f",pt))
  lb=as.numeric(sprintf("%0.2f",lb))
  ub=as.numeric(sprintf("%0.2f",ub))
  paste(pt," (",lb,",",ub,")",sep="")
}

#-------------------------------------------------------
# Function to organize data into table
#-------------------------------------------------------
table.prep=function(df,cond1,cond2,cond3){
  #name 
  name=as.character(df[["groupname"]][1])
  
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
    or.b.o.3=pretty.ci(or.b.o.pt.3,or.b.o.lb.3,or.b.o.ub.3)
  }else{
    or.b.o.3="NA"
  }
  
  if(cond2 %in% df[["conditions"]]==TRUE){
    or.b.c.se.3=df[["se"]][df[["conditions"]]==cond2
                           & df[["adj"]]=="adjusted" & df[["outcome"]]=="hcgi3ci3"] 
    or.b.c.lb.3=or.b.c.pt.3-(qnorm(0.975)*or.b.c.se.3)
    or.b.c.ub.3=or.b.c.pt.3+(qnorm(0.975)*or.b.c.se.3)
    or.b.c.3=pretty.ci(or.b.c.pt.3,or.b.c.lb.3,or.b.c.ub.3)
  }else{
    or.b.c.3="NA"
  }
  
  if(cond3 %in% df[["conditions"]]==TRUE){
    or.comb.se.3=df[["se"]][df[["conditions"]]==cond3
                            & df[["adj"]]=="adjusted" & df[["outcome"]]=="hcgi3ci3"] 
    or.comb.lb.3=or.comb.pt.3-(qnorm(0.975)*or.comb.se.3)
    or.comb.ub.3=or.comb.pt.3+(qnorm(0.975)*or.comb.se.3)
    or.comb.3=pretty.ci(or.comb.pt.3,or.comb.lb.3,or.comb.ub.3)
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
    or.b.o.10=pretty.ci(or.b.o.pt.10,or.b.o.lb.10,or.b.o.ub.10)
  }else{
    or.b.o.10="NA"
  }  
  
  if(cond2 %in% df[["conditions"]]==TRUE){
    or.b.c.se.10=df[["se"]][df[["conditions"]]==cond2
                            & df[["adj"]]=="adjusted" & df[["outcome"]]=="hcgi3ci10"] 
    or.b.c.lb.pt.10=or.b.c.10-(qnorm(0.975)*or.b.c.se.10)
    or.b.c.ub.pt.10=or.b.c.10+(qnorm(0.975)*or.b.c.se.10)
    or.b.c.10=pretty.ci(or.b.c.pt.10,or.b.c.lb.10,or.b.c.ub.10)
  }else{
    or.b.c.10="NA"
  }  
  
  if(cond3 %in% df[["conditions"]]==TRUE){
    or.comb.se.10=df[["se"]][df[["conditions"]]==cond3
                             & df[["adj"]]=="adjusted" & df[["outcome"]]=="hcgi3ci10"] 
    or.comb.lb.10=or.comb.pt.10-(qnorm(0.975)*or.comb.se.10)
    or.comb.ub.10=or.comb.pt.10+(qnorm(0.975)*or.comb.se.10)
    or.comb.10=pretty.ci(or.comb.pt.10,or.comb.lb.10,or.comb.ub.10)
  }else{
    or.comb.10="NA"
  }
  
  out=c(name,or.b.o.3,or.b.c.3,or.comb.3,or.b.o.10,or.b.c.10,or.comb.10)
  return(out)
}


########################################################
# DOHENY
########################################################

#-------------------------------------------------------
# Load data
#-------------------------------------------------------
d_14FMC0511=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/doheny/Doheny_Public_Indicator_Results-June-2013/doheny_14FMC0511_ind.csv",
                         header=TRUE)
d_14FMC0512=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/doheny/Doheny_Public_Indicator_Results-June-2013/doheny_14FMC0512_ind.csv",
                         header=TRUE)
d_14FPC0511=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/doheny/Doheny_Public_Indicator_Results-June-2013/doheny_14FPC0511_ind.csv",
                         header=TRUE)
d_14FPC0512=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/doheny/Doheny_Public_Indicator_Results-June-2013/doheny_14FPC0512_ind.csv",
                         header=TRUE)
d_14FPC0611=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/doheny/Doheny_Public_Indicator_Results-June-2013/doheny_14FPC0611_ind.csv",
                     header=TRUE)
d_14FPC0612=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/doheny/Doheny_Public_Indicator_Results-June-2013/doheny_14FPC0612_ind.csv",
                     header=TRUE)
d_14FPC0711=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/doheny/Doheny_Public_Indicator_Results-June-2013/doheny_14FPC0711_ind.csv",
                     header=TRUE)
d_14FPC0712=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/doheny/Doheny_Public_Indicator_Results-June-2013/doheny_14FPC0712_ind.csv",
                     header=TRUE)
d_14FPC3412=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/doheny/Doheny_Public_Indicator_Results-June-2013/doheny_14FPC3412_ind.csv",
                     header=TRUE)
d_14FPC3511=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/doheny/Doheny_Public_Indicator_Results-June-2013/doheny_14FPC3511_ind.csv",
                     header=TRUE)
d_14FPC3512=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/doheny/Doheny_Public_Indicator_Results-June-2013/doheny_14FPC3512_ind.csv",
                     header=TRUE)
d_14FPC3513=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/doheny/Doheny_Public_Indicator_Results-June-2013/doheny_14FPC3513_ind.csv",
                     header=TRUE)
d_14FPC3514=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/doheny/Doheny_Public_Indicator_Results-June-2013/doheny_14FPC3514_ind.csv",
                     header=TRUE)
d_14FPC3515=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/doheny/Doheny_Public_Indicator_Results-June-2013/doheny_14FPC3515_ind.csv",
                     header=TRUE)

#-------------------------------------------------------
# Process data
#-------------------------------------------------------
d_14FMC0511.p=process.df(d_14FMC0511)
d_14FMC0512.p=process.df(d_14FMC0512)
d_14FPC0511.p=process.df(d_14FPC0511)
d_14FPC0512.p=process.df(d_14FPC0512)
d_14FPC0611.p=process.df(d_14FPC0611)
d_14FPC0612.p=process.df(d_14FPC0612)
d_14FPC0711.p=process.df(d_14FPC0711)
d_14FPC0712.p=process.df(d_14FPC0712)
d_14FPC3412.p=process.df(d_14FPC3412)
d_14FPC3511.p=process.df(d_14FPC3511)
d_14FPC3512.p=process.df(d_14FPC3512)
d_14FPC3513.p=process.df(d_14FPC3513)
d_14FPC3514.p=process.df(d_14FPC3514)
d_14FPC3515.p=process.df(d_14FPC3515)

#-------------------------------------------------------
# Organizing data into table
#-------------------------------------------------------
d_14FMC0511.p.t=table.prep(d_14FMC0511.p,"sand berm open","sand berm closed","combined")
d_14FMC0512.p.t=table.prep(d_14FMC0512.p,"sand berm open","sand berm closed","combined")
d_14FPC0511.p.t=table.prep(d_14FPC0511.p,"sand berm open","sand berm closed","combined")
d_14FPC0512.p.t=table.prep(d_14FPC0512.p,"sand berm open","sand berm closed","combined")
d_14FPC0611.p.t=table.prep(d_14FPC0611.p,"sand berm open","sand berm closed","combined")
d_14FPC0612.p.t=table.prep(d_14FPC0612.p,"sand berm open","sand berm closed","combined")
d_14FPC0711.p.t=table.prep(d_14FPC0711.p,"sand berm open","sand berm closed","combined")
d_14FPC0712.p.t=table.prep(d_14FPC0712.p,"sand berm open","sand berm closed","combined")
d_14FPC3412.p.t=table.prep(d_14FPC3412.p,"sand berm open","sand berm closed","combined")
d_14FPC3511.p.t=table.prep(d_14FPC3511.p,"sand berm open","sand berm closed","combined")
d_14FPC3512.p.t=table.prep(d_14FPC3512.p,"sand berm open","sand berm closed","combined")
d_14FPC3513.p.t=table.prep(d_14FPC3513.p,"sand berm open","sand berm closed","combined")
d_14FPC3514.p.t=table.prep(d_14FPC3514.p,"sand berm open","sand berm closed","combined")
d_14FPC3515.p.t=table.prep(d_14FPC3515.p,"sand berm open","sand berm closed","combined")

or.tab=rbind(d_14FMC0511.p.t,d_14FMC0512.p.t,d_14FPC0511.p.t,
            d_14FPC0512.p.t,d_14FPC0611.p.t,d_14FPC0612.p.t,
             d_14FPC0711.p.t,d_14FPC0712.p.t,d_14FPC3412.p.t,
             d_14FPC3511.p.t,d_14FPC3512.p.t,d_14FPC3513.p.t,
             d_14FPC3514.p.t,d_14FPC3515.p.t)

colnames(or.tab)=c("Coliphage indicator","3day - berm open", "3day - berm closed","3day - combined",
                   "10day - berm open", "10day - berm closed","10day - combined")

rownames(or.tab)=NULL
save(or.tab,file="")

########################################################
# AVALON
########################################################

#-------------------------------------------------------
# Load data
#-------------------------------------------------------

a_14FMC0511=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/Avalon/Avalon_Public_Indicator_Results-June-2013/avalon_14FMC0511_ind.csv",
                     header=TRUE)
a_14FPC0511=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/Avalon/Avalon_Public_Indicator_Results-June-2013/avalon_14FPC0511_ind.csv",
                     header=TRUE)
a_14FPC0512=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/Avalon/Avalon_Public_Indicator_Results-June-2013/avalon_14FPC0512_ind.csv",
                     header=TRUE)
a_14FPC0611=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/Avalon/Avalon_Public_Indicator_Results-June-2013/avalon_14FPC0611_ind.csv",
                     header=TRUE)
a_14FPC0711=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/Avalon/Avalon_Public_Indicator_Results-June-2013/avalon_14FPC0711_ind.csv",
                     header=TRUE)
a_14FPC3412=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/Avalon/Avalon_Public_Indicator_Results-June-2013/avalon_14FPC3412_ind.csv",
                     header=TRUE)
a_14FPC3511=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/Avalon/Avalon_Public_Indicator_Results-June-2013/avalon_14FPC3511_ind.csv",
                     header=TRUE)
a_14FPC3512=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/Avalon/Avalon_Public_Indicator_Results-June-2013/avalon_14FPC3512_ind.csv",
                     header=TRUE)
a_14FPC3513=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/Avalon/Avalon_Public_Indicator_Results-June-2013/avalon_14FPC3513_ind.csv",
                     header=TRUE)
a_14FPC3514=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/Avalon/Avalon_Public_Indicator_Results-June-2013/avalon_14FPC3514_ind.csv",
                     header=TRUE)
a_14FPC3515=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/Avalon/Avalon_Public_Indicator_Results-June-2013/avalon_14FPC3515_ind.csv",
                     header=TRUE)

#-------------------------------------------------------
# Process data
#-------------------------------------------------------
a_14FMC0511.p=process.df(a_14FMC0511)
a_14FPC0511.p=process.df(a_14FPC0511)
a_14FPC0512.p=process.df(a_14FPC0512)
a_14FPC0611.p=process.df(a_14FPC0611)
a_14FPC0711.p=process.df(a_14FPC0711)
a_14FPC3412.p=process.df(a_14FPC3412)
a_14FPC3511.p=process.df(a_14FPC3511)
a_14FPC3512.p=process.df(a_14FPC3512)
a_14FPC3513.p=process.df(a_14FPC3513)
a_14FPC3514.p=process.df(a_14FPC3514)
a_14FPC3515.p=process.df(a_14FPC3515)

#-------------------------------------------------------
# Organizing data into table
#-------------------------------------------------------
a_14FMC0511.p.t=table.prep(a_14FMC0511.p,"groundwater above median","groundwater below median","combined")
a_14FPC0511.p.t=table.prep(a_14FPC0511.p,"groundwater above median","groundwater below median","combined")
a_14FPC0512.p.t=table.prep(a_14FPC0512.p,"groundwater above median","groundwater below median","combined")
a_14FPC0611.p.t=table.prep(a_14FPC0611.p,"groundwater above median","groundwater below median","combined")
a_14FPC0711.p.t=table.prep(a_14FPC0711.p,"groundwater above median","groundwater below median","combined")
a_14FPC3412.p.t=table.prep(a_14FPC3412.p,"groundwater above median","groundwater below median","combined")
a_14FPC3511.p.t=table.prep(a_14FPC3511.p,"groundwater above median","groundwater below median","combined")
a_14FPC3512.p.t=table.prep(a_14FPC3512.p,"groundwater above median","groundwater below median","combined")
a_14FPC3513.p.t=table.prep(a_14FPC3513.p,"groundwater above median","groundwater below median","combined")
a_14FPC3514.p.t=table.prep(a_14FPC3514.p,"groundwater above median","groundwater below median","combined")
a_14FPC3515.p.t=table.prep(a_14FPC3515.p,"groundwater above median","groundwater below median","combined")

or.tab=rbind(a_14FMC0511.p.t,a_14FPC0511.p.t,
             a_14FPC0512.p.t,a_14FPC0611.p.t,a_14FPC3412.p.t,
            a_14FPC3512.p.t,a_14FPC3513.p.t)

colnames(or.tab)=c("Coliphage indicator","3day - gw above median", "3day - below median","3day - combined",
                   "10day - gw above median", "10day - gw below median","10day - combined")

rownames(or.tab)=NULL


########################################################
# MALIBU
########################################################
#-------------------------------------------------------
# Load data
#-------------------------------------------------------
m_14FPC0511=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/malibu/Malibu_Public_Results_Indicators-June-2013/malibu_14FPC0511_ind.csv",
                     header=TRUE)
m_14FPC0512=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/malibu/Malibu_Public_Results_Indicators-June-2013/malibu_14FPC0512_ind.csv",
                     header=TRUE)
m_14FPC3412=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/malibu/Malibu_Public_Results_Indicators-June-2013/malibu_14FPC3412_ind.csv",
                     header=TRUE)
m_14FPC3511=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/malibu/Malibu_Public_Results_Indicators-June-2013/malibu_14FPC3511_ind.csv",
                     header=TRUE)
m_14FPC3512=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/malibu/Malibu_Public_Results_Indicators-June-2013/malibu_14FPC3512_ind.csv",
                     header=TRUE)
m_14FPC3513=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/malibu/Malibu_Public_Results_Indicators-June-2013/malibu_14FPC3513_ind.csv",
                     header=TRUE)
m_14FPC3514=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/malibu/Malibu_Public_Results_Indicators-June-2013/malibu_14FPC3514_ind.csv",
                     header=TRUE)
m_14FPC3515=read.csv("/Users/jadebc/Dropbox/Beaches/Results/Beaches-Public_Results/indicator-analysis-results/malibu/Malibu_Public_Results_Indicators-June-2013/malibu_14FPC3515_ind.csv",
                     header=TRUE)


#-------------------------------------------------------
# Process data
#-------------------------------------------------------
m_14FPC0511.p=process.df(m_14FPC0511)
m_14FPC0512.p=process.df(m_14FPC0512)
m_14FPC3412.p=process.df(m_14FPC3412)
m_14FPC3511.p=process.df(m_14FPC3511)
m_14FPC3512.p=process.df(m_14FPC3512)
m_14FPC3513.p=process.df(m_14FPC3513)
m_14FPC3514.p=process.df(m_14FPC3514)
m_14FPC3515.p=process.df(m_14FPC3515)

#-------------------------------------------------------
# Organizing data into table
#-------------------------------------------------------
m_14FPC0511.p.t=table.prep(m_14FPC0511.p,"sand berm open","sand berm closed","combined")
m_14FPC0512.p.t=table.prep(m_14FPC0512.p,"sand berm open","sand berm closed","combined")
m_14FPC3412.p.t=table.prep(m_14FPC3412.p,"sand berm open","sand berm closed","combined")
m_14FPC3511.p.t=table.prep(m_14FPC3511.p,"sand berm open","sand berm closed","combined")
m_14FPC3512.p.t=table.prep(m_14FPC3512.p,"sand berm open","sand berm closed","combined")
m_14FPC3513.p.t=table.prep(m_14FPC3513.p,"sand berm open","sand berm closed","combined")
m_14FPC3514.p.t=table.prep(m_14FPC3514.p,"sand berm open","sand berm closed","combined")
m_14FPC3515.p.t=table.prep(m_14FPC3515.p,"sand berm open","sand berm closed","combined")

or.tab=rbind(d_14FPC0511.p.t,d_14FPC0512.p.t,d_14FPC3412.p.t,
             d_14FPC3512.p.t)

colnames(or.tab)=c("Coliphage indicator","3day - berm open", "3day - berm closed","3day - combined",
                   "10day - berm open", "10day - berm closed","10day - combined")

rownames(or.tab)=NULL

