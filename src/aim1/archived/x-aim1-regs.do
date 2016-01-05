capture log close
set more off
clear all

log using "~/dropbox/13beaches/programs/aim1/1-aim1-regs.log", text replace

*----------------------------------------
* 1-aim1-regs.do
* ben arnold
*
* estimate the risk of GI illness
* associated with a 1-log increase in
* water quality indicator concentrations
*----------------------------------------

*----------------------------------------
* input files:
*	13beaches-analysis.dta
*
* output files:
*	xxx
*
*----------------------------------------


*----------------------------------------
* Load the data and restrict to
* the analysis population
*----------------------------------------
use "~/dropbox/13beaches/data/final/13beaches-analysis.dta", clear


keep if swimmer==1

drop if venfest==1

drop if nowq==1

drop if gibase==1




*----------------------------------------
* prelim regressions to check against
* results in R
*----------------------------------------

xi: glm gici10 i.berm*avgdyenterocfu if beach=="Doheny", family(binomial) link(log)
xi: glm gici10 i.berm*avgdyenterocfu if beach=="Doheny", family(binomial) link(log) cluster(hhid) robust 

glm, eform

lincom avgdyenterocfu+_IberXavgdy_1, eform


glm gici10 avgdyenterocfu if beach=="Malibu", family(binomial) link(log) eform cluster(hhid) robust

* xi: logit gici10 i.berm*avgdyenterocfu if beach=="Doheny", cluster(hhid) robust
* xi: poisson gici10 i.berm*avgdyenterocfu if beach=="Doheny", cluster(hhid) robust







log close
exit




