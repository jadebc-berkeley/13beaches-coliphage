capture log close
set more off
clear all

log using "~/Documents/CRG/coliphage/13beaches-coliphage/src/analysis/9-test-trend.log", text replace

******************************************
* Coliphage analysis

* Test of trend
******************************************


insheet using "~/Documents/CRG/coliphage/13beaches-data/final/13beaches-analysis-processed.csv", clear

replace fmcpres="." if fmcpres=="NA"
replace fpcpres="." if fpcpres=="NA"
destring fmcpres fpcpres, replace

gen expfmc=.
replace expfmc = 1 if bodycontact=="No"
replace expfmc = 2 if bodycontact=="Yes" & fmcpres==0
replace expfmc = 3 if bodycontact=="Yes" & fmcpres==1
replace expfmc = 4 if bodycontact=="Yes" & fmcpres==1 & entero35==1

gen expfpc=.
replace expfpc = 1 if bodycontact=="No"
replace expfpc = 2 if bodycontact=="Yes" & fpcpres==0
replace expfpc = 3 if bodycontact=="Yes" & fpcpres==1
replace expfpc = 4 if bodycontact=="Yes" & fpcpres==1 & entero35==1

glm gici10 ibn.expfmc, noconstant family(binomial) link(log)
test -3*1.expfmc - 2.expfmc + 3.expfmc + 3*4.expfmc = 0

glm gici10 ibn.expfpc, noconstant family(binomial) link(log)
test -3*1.expfpc - 2.expfpc + 3.expfpc + 3*4.expfpc = 0




log close
