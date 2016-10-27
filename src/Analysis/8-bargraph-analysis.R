##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 11/10/15

# Analyses for bar graph comparing risk of illness
# under different conditions
##########################################

rm(list=ls())

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

write.csv(all,file="~/Documents/CRG/coliphage/13beaches-data/final/13beaches-analysis-processed.csv",row.names=FALSE)

#-----------------------------------------
# Non-swimmers
#-----------------------------------------
ns=prop.table(table(all$anycontact,all$gici10),1)[1,2]
se.ns=sqrt(ns*(1-ns)/sum(table(all$anycontact,all$gici10)))
lb.ns=ns-(qnorm(0.975)*se.ns)
ub.ns=ns+(qnorm(0.975)*se.ns)

#-----------------------------------------
# Swimmers not exposed to coliphage
#-----------------------------------------
swim=prop.table(table(all$bodycontact[all$fpc.pres==0 & all$fmc.pres==0],
                      all$gici10[all$fpc.pres==0 & all$fmc.pres==0]),1)[2,2]
se.swim=sqrt(swim*(1-swim)/sum(table(all$bodycontact[all$fpc.pres==0 & all$fmc.pres==0],
                                 all$gici10[all$fpc.pres==0 & all$fmc.pres==0])))
lb.swim=swim-(qnorm(0.975)*se.swim)
ub.swim=swim+(qnorm(0.975)*se.swim)

#-----------------------------------------
# Swimmers exposed to somatic coliphage
#-----------------------------------------
s.som=prop.table(table(all$fmc.pres[all$bodycontact=="Yes"],
                       all$gici10[all$bodycontact=="Yes"]),1)[2,2]
se.som=sqrt(swim*(1-swim)/sum(table(all$fmc.pres[all$bodycontact=="Yes"],
                                     all$gici10[all$bodycontact=="Yes"])))
lb.som=s.som-(qnorm(0.975)*se.som)
ub.som=s.som+(qnorm(0.975)*se.som)

#-----------------------------------------
# Swimmers exposed to male-specific coliphage
#-----------------------------------------
s.male=prop.table(table(all$fpc.pres[all$bodycontact=="Yes"],
                       all$gici10[all$bodycontact=="Yes"]),1)[2,2]
se.male=sqrt(swim*(1-swim)/sum(table(all$fpc.pres[all$bodycontact=="Yes"],
                                    all$gici10[all$bodycontact=="Yes"])))
lb.male=s.male-(qnorm(0.975)*se.male)
ub.male=s.male+(qnorm(0.975)*se.male)

#-----------------------------------------
# Swimmers exposed to somatic coliphage + 
# enterococcus > 35
#-----------------------------------------
s.som.ent=prop.table(table(all$fmc.pres[all$bodycontact=="Yes" & all$entero35==1],
                       all$gici10[all$bodycontact=="Yes" & all$entero35==1]),1)[2,2]
se.som.ent=sqrt(swim*(1-swim)/sum(table(all$fmc.pres[all$bodycontact=="Yes" & all$entero35==1],
                                    all$gici10[all$bodycontact=="Yes" & all$entero35==1])))
lb.som.ent=s.som.ent-(qnorm(0.975)*se.som.ent)
ub.som.ent=s.som.ent+(qnorm(0.975)*se.som.ent)

#-----------------------------------------
# Swimmers exposed to male-specific coliphage + 
# enterococcus > 35
#-----------------------------------------
s.male.ent=prop.table(table(all$fpc.pres[all$bodycontact=="Yes" & all$entero35==1],
        all$gici10[all$bodycontact=="Yes" & all$entero35==1]),1)[2,2]
se.male.ent=sqrt(swim*(1-swim)/sum(table(all$fpc.pres[all$bodycontact=="Yes" & all$entero35==1],
        all$gici10[all$bodycontact=="Yes" & all$entero35==1])))
lb.male.ent=s.male.ent-(qnorm(0.975)*se.male.ent)
ub.male.ent=s.male.ent+(qnorm(0.975)*se.male.ent)

#-----------------------------------------
# Test of trend
#-----------------------------------------
all$exp.fmc=NA
all$exp.fmc[all$bodycontact=="No"]=0
all$exp.fmc[all$bodycontact=="Yes" & all$fmc.pres==0]=1
all$exp.fmc[all$bodycontact=="Yes" & all$fmc.pres==1]=2
all$exp.fmc[all$bodycontact=="Yes" & all$fmc.pres==1 & all$entero35==1]=3

temp=glm(gici10~exp.fmc,family="binomial"(link="log"),data=all)

all$exp.fpc=NA
all$exp.fpc[all$bodycontact=="No"]=0
all$exp.fpc[all$bodycontact=="Yes" & all$fpc.pres==0]=1
all$exp.fpc[all$bodycontact=="Yes" & all$fpc.pres==1]=2
all$exp.fpc[all$bodycontact=="Yes" & all$fpc.pres==1 & all$entero35==1]=3

temp=glm(gici10~exp.fpc,family="binomial"(link="log"),data=all)

save(ns, se.ns, lb.ns, ub.ns,
     swim, se.swim, lb.swim, ub.swim,
     s.som, se.som, lb.som, ub.som,
     s.male, se.male, lb.male, ub.male,
     s.som.ent, se.som.ent, lb.som.ent, ub.som.ent,
     s.male.ent, se.male.ent, lb.male.ent, ub.male.ent, 
     
     file="~/Documents/CRG/coliphage/results/rawoutput/bargraph-output.Rdata")



