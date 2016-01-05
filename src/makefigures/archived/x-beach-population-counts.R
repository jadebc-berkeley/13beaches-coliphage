


rm(list=ls())
library(xtable)

d <- read.csv("~/Dropbox/13beaches/data/final/13beaches-analysis.csv")




pop <- table(d$beach)
pop0to4 <- table(d$beach[d$agestrat=="(0, 4]"])
pop5to10 <- table(d$beach[d$agestrat=="(4, 10]"])

swim <- table(d$beach[d$bodycontact=="Yes"])
swim0to4 <- table(d$beach[d$bodycontact=="Yes" & d$agestrat=="(0, 4]"])
swim5to10 <- table(d$beach[d$bodycontact=="Yes" & d$agestrat=="(4, 10]"])

Ntab <- cbind(pop,pop0to4,pop5to10,swim,swim0to4,swim5to10)

print(xtable(Ntab),type="html",file="~/dropbox/13beaches/aim1-results/rawoutput/aim1-beach-population-counts.xls")