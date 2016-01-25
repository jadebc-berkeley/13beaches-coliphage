##########################################
# Coliphage analysis - 6 beaches
# v1 by Jade 7/13/15

# This file makes a table with regression output
# For a table with the combined and interaction
# results for body immersion

# Single & joint indicators
##########################################


rm(list=ls())

load("~/Documents/CRG/coliphage/results/rawoutput/regress-3day-body.Rdata")
load("~/Documents/CRG/coliphage/results/rawoutput/regress-3day-body-entero.Rdata")
load("~/Documents/CRG/coliphage/results/rawoutput/regress-3day-body-joint.Rdata")

# ------------------------------------------------
# combine output
# ------------------------------------------------

# n's --------------------------------------------
gici3.n.body.pool=c(all.n3.fmc1601,all.n3.fmc1602,
                     all.n3.fpc1601,all.n3.fpc1602)

gici3.n.body.high=c(all.n3.fmc1602.high,
                     all.n3.fpc1601.high,all.n3.fpc1602.high)

gici3.n.body.low=c(all.n3.fmc1602.low,
                    all.n3.fpc1601.low,all.n3.fpc1602.low)

# add comma
gici3.n.body.pool=format(gici3.n.body.pool,scientific=FALSE,big.mark=",")
gici3.n.body.high=format(gici3.n.body.high,scientific=FALSE,big.mark=",")
gici3.n.body.low=format(gici3.n.body.low,scientific=FALSE,big.mark=",")

# entero n's  --------------------------------------------

entero.n.body.pool=c(n3.entero35.fmc1601,n3.entero35.fmc1602,
                     n3.entero35.fpc1601,n3.entero35.fpc1602)

entero.n.body.high=c(n3.entero35.fmc1602.high,n3.entero35.fpc1601.high,
                     n3.entero35.fpc1602.high)

entero.n.body.low=c(n3.entero35.fmc1602.low,n3.entero35.fpc1601.low,
                    n3.entero35.fpc1602.low)

# add comma
entero.n.body.pool=format(entero.n.body.pool,scientific=FALSE,big.mark=",")
entero.n.body.high=format(entero.n.body.high,scientific=FALSE,big.mark=",")
entero.n.body.low=format(entero.n.body.low,scientific=FALSE,big.mark=",")

# joint n's --------------------------------------------
gici3.n.joint=c(all.n3.fmc1601.joint,all.n3.fmc1602.joint,
                 all.n3.fpc1601.joint,all.n3.fpc1602.joint)

gici3.n.joint.high=c(all.n3.fmc1602.high.joint,
                      all.n3.fpc1601.high.joint,all.n3.fpc1602.high.joint)

gici3.n.joint.low=c(all.n3.fmc1602.low.joint,
                     all.n3.fpc1601.low.joint,all.n3.fpc1602.low.joint)

# add comma
gici3.n.joint=format(gici3.n.joint,scientific=FALSE,big.mark=",")
gici3.n.joint.high=format(gici3.n.joint.high,scientific=FALSE,big.mark=",")
gici3.n.joint.low=format(gici3.n.joint.low,scientific=FALSE,big.mark=",")

# single indicator estimates --------------------------------------------
gici3.body.pool=list(overall.fit3.fmc1601, overall.fit3.fmc1602,
                      overall.fit3.fpc1601, overall.fit3.fpc1602)

gici3.body.pool.high=list(overall.fit3.fmc1602.high,
                           overall.fit3.fpc1601.high, overall.fit3.fpc1602.high)

gici3.body.pool.low=list(overall.fit3.fmc1602.low,
                          overall.fit3.fpc1601.low, overall.fit3.fpc1602.low)

# entero estimates  --------------------------------------------

gici3.entero.pool=list(overall.fit3.entero.fmc1601,overall.fit3.entero.fmc1602,
                        overall.fit3.entero.fpc1601,overall.fit3.entero.fpc1602)
gici3.entero.high=list(overall.fit3.entero.high.fmc1602,overall.fit3.entero.high.fpc1601,
                        overall.fit3.entero.high.fpc1602)
gici3.entero.low=list(overall.fit3.entero.low.fmc1602,overall.fit3.entero.low.fpc1601,
                       overall.fit3.entero.low.fpc1602)


# joint indicator estimates  --------------------------------------------
gici3.joint.int=list(overall.fit3.fmc1601.int, overall.fit3.fmc1602.int,
                      overall.fit3.fpc1601.int, overall.fit3.fpc1602.int)

gici3.joint.high.int=list(overall.fit3.fmc1602.high.int,
                           overall.fit3.fpc1601.high.int, overall.fit3.fpc1602.high.int)

gici3.joint.low.int=list(overall.fit3.fmc1602.low.int,
                          overall.fit3.fpc1601.low.int, overall.fit3.fpc1602.low.int)


# ------------------------------------------------
# function to make table row with exponentiated point estimate
# and 95% ci in parentheses
# ------------------------------------------------
mkrow.pool=function(out){
  pt.est=out[2,1]
  lb=pt.est-qnorm(.975)*out[2,2]
  ub=pt.est+qnorm(.975)*out[2,2]
  paste(sprintf("%0.2f",exp(pt.est))," (",
        sprintf("%0.2f",exp(lb)),",",
        sprintf("%0.2f",exp(ub)), ")",sep="")
}

mkrow.joint=function(out){
  row=grep("ent4",rownames(out))
  pt.est=out[row,1]
  lb=pt.est-qnorm(.975)*out[4,2]
  ub=pt.est+qnorm(.975)*out[4,2]
  paste(sprintf("%0.2f",exp(pt.est))," (",
        sprintf("%0.2f",exp(lb)),",",
        sprintf("%0.2f",exp(ub)), ")",sep="")
}

# ------------------------------------------------
# convert results into table format
# pooled results
# ------------------------------------------------
# single indicator results
gici3.tab.body.pool=data.frame(combined=unlist(lapply(gici3.body.pool,mkrow.pool)))
gici3.tab.body.pool.high=data.frame(combined=unlist(lapply(gici3.body.pool.high,mkrow.pool)))
gici3.tab.body.pool.low=data.frame(combined=unlist(lapply(gici3.body.pool.low,mkrow.pool)))

gici3.tab.body.pool.high$combined=as.character(gici3.tab.body.pool.high$combined)
gici3.tab.body.pool.high=rbind("",gici3.tab.body.pool.high)
gici3.tab.body.pool.low$combined=as.character(gici3.tab.body.pool.low$combined)
gici3.tab.body.pool.low=rbind("",gici3.tab.body.pool.low)

gici3.tab=data.frame(cbind(gici3.n.body.pool,gici3.tab.body.pool,
                            c("",gici3.n.body.low),gici3.tab.body.pool.low,
                            c("",gici3.n.body.high),gici3.tab.body.pool.high))
colnames(gici3.tab)=NULL
gici3.tab[,1]=as.character(gici3.tab[,1])
gici3.tab[,2]=as.character(gici3.tab[,2])
gici3.tab[,3]=as.character(gici3.tab[,3])
gici3.tab[,4]=as.character(gici3.tab[,4])
gici3.tab[,5]=as.character(gici3.tab[,5])
gici3.tab[,6]=as.character(gici3.tab[,6])

# entero results
entero.tab.pool=data.frame(combined=unlist(lapply(gici3.entero.pool,mkrow.pool)))
entero.tab.high=data.frame(combined=unlist(lapply(gici3.entero.high,mkrow.pool)))
entero.tab.high$combined=as.character(entero.tab.high$combined)
entero.tab.high=rbind("",entero.tab.high)
entero.tab.low=data.frame(combined=unlist(lapply(gici3.entero.low,mkrow.pool)))
entero.tab.low$combined=as.character(entero.tab.low$combined)
entero.tab.low=rbind("",entero.tab.low)

entero.tab=data.frame(cbind(entero.n.body.pool,entero.tab.pool,
                            c("",entero.n.body.low),entero.tab.low,
                            c("",entero.n.body.high),entero.tab.high))
colnames(entero.tab)=NULL
entero.tab[,1]=as.character(entero.tab[,1])
entero.tab[,2]=as.character(entero.tab[,2])
entero.tab[,3]=as.character(entero.tab[,3])
entero.tab[,4]=as.character(entero.tab[,4])
entero.tab[,5]=as.character(entero.tab[,5])
entero.tab[,6]=as.character(entero.tab[,6])

# joint results
joint.tab.pool=data.frame(combined=unlist(lapply(gici3.joint.int,mkrow.joint)))
joint.tab.high=data.frame(combined=unlist(lapply(gici3.joint.high.int,mkrow.joint)))
joint.tab.high$combined=as.character(joint.tab.high$combined)
joint.tab.high=rbind("",joint.tab.high)
joint.tab.low=data.frame(combined=unlist(lapply(gici3.joint.low.int,mkrow.joint)))
joint.tab.low$combined=as.character(joint.tab.low$combined)
joint.tab.low=rbind("",joint.tab.low)

joint.tab=data.frame(cbind(gici3.n.joint,joint.tab.pool,
                            c("",gici3.n.joint.low),joint.tab.low,
                            c("",gici3.n.joint.high),joint.tab.high))
colnames(joint.tab)=NULL
joint.tab[,1]=as.character(joint.tab[,1])
joint.tab[,2]=as.character(joint.tab[,2])
joint.tab[,3]=as.character(joint.tab[,3])
joint.tab[,4]=as.character(joint.tab[,4])
joint.tab[,5]=as.character(joint.tab[,5])
joint.tab[,6]=as.character(joint.tab[,6])

colnames(gici3.tab)=c("npool","pool","nlow","low","nhigh","high")
colnames(entero.tab)=c("npool","pool","nlow","low","nhigh","high")
colnames(joint.tab)=c("npool","pool","nlow","low","nhigh","high")

gici3.tab.out=rbind(gici3.tab[1,],entero.tab[1,],joint.tab[1,],
                     gici3.tab[2,],entero.tab[2,],joint.tab[2,],
                     gici3.tab[3,],entero.tab[3,],joint.tab[3,],
                     gici3.tab[4,],entero.tab[4,],joint.tab[4,])

lab=c(rep(c("Coliphage detected","Entero > 35 CFU/100 ml",
          "Coliphage detected $\\&$ Entero > 35 CFU/100 ml"),4))
gici3.tab.out=cbind(lab,gici3.tab.out)

rownames(gici3.tab.out)=NULL

coli=c(rep("Somatic coliphage (EPA 1601)",3),
       rep("Somatic coliphage (EPA 1602)",3),
       rep("Male-specific coliphage (EPA 1601)",3),
       rep("Male-specific coliphage (EPA 1602)",3))

gici3.tab.out=cbind(coli,gici3.tab.out)

save(gici3.tab.out,file="~/Documents/CRG/coliphage/Results/Tables/CIR-3.RData")




