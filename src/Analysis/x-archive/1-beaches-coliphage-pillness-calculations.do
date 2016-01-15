capture log close
set more off
clear all

log using "/Users/jadebc/Dropbox/Coliphage/Programs/1-beaches-coliphage-pillness-calculations.log", text replace

*-----------------------------------------
* 1-beaches-coliphage-pillness-calculations.do
*
* Ben Arnold
*
* Calculate the marginal predicted probability
* of diarrhea and GI illness at observed
* levels of different Coliphage indicators
* for Avalon, Doheny, and Malibu studies
*
* Use site-specific daily average values
* Repeat the calculations under all conditions
* and under "high risk" conditions, where
* high risk means berm open at Doheny and Malibu
* and above median groundwater flow for Avalon
*
* run calculations for swimmers with body immersion
* and under different environmental conditions
* Avalon = combined, groundwater above/below median
* Doheny/Malibu = combined, berm open/closed
* 
* version 1 (17 sep 2013)
*
*-----------------------------------------

*-----------------------------------------
* input files:
*   avalon2.dta
*   doheny2.dta
*   malibu2.dta
*
* output files:
*   avalon-coliphage-pillness.dta
*	doheny-coliphage-pillness.dta
*	malibu-coliphage-pillness.dta
*-----------------------------------------

*-----------------------------------------
* source the backward selection algorithm
*-----------------------------------------
qui do "~/Dropbox/Coliphage/Programs/2sup-bs_reg.do"


*------------------------------------------
* Function to calculate the probability of
* illness at each level of the indicator
* (exposure) of interest, adjusted at the
* mean of covariates.
*
* fit models on all body immersion swimmers
*------------------------------------------

capture program drop pillness
program define pillness, rclass
	version 10
	syntax [if], outcome(varlist max=1 numeric) exposure(varlist max=1 numeric) Xvars(string) id(varlist max=1 numeric)
	* outcome = health outcome of interest
	* exposure = indicator of interest (log scale)
	* Xvars = covariates to consider in adjusted models
	* id = id for independent units (for robust SE calculation)
	
	marksample touse

	* run the backward selection algorithm
	* info returned from bs_reg: 
	* covariates selected in r(covarlist), sparse data flag in r(sparse)
	qui bs_reg if `touse', model(logit) outcome(`outcome') exposure(`exposure') covar("`xvars'") keepvar() unit(1) criteria(0.05)
  		local covarlist = r(covarlist)
  		* use continuous age rather than age category
		* to get adjusted probabilities at the mean
		* value
		local covarlist = subinstr("`covarlist'","i.agecat1","ageyrx",.)
		local covarlist = subinstr("`covarlist'","i.agecat2","ageyrx",.)
		if("`covarlist'"==".") local covarlist = ""
	
	* estimate a logistic model
	* use (v10 for backward compatability with
	* the -adjust- function)
	version 10
	logistic `outcome' `exposure' `covarlist' if `touse', cluster(`id') robust

	* store the estimated OR and its SE and CI
	gen orn = e(N)
	gen or = exp(_b[`exposure'])
	gen se = exp(_b[`exposure'])*_se[`exposure']
	gen orlb = exp(_b[`exposure'] - invnormal(0.975)*_se[`exposure'])
	gen orub = exp(_b[`exposure'] + invnormal(0.975)*_se[`exposure'])
		label var orn "N indivs used to estimate OR"
		label var or "Odds ratio"
		label var se "SE of odds ratio"
		label var orlb "OR lower 95% CI"
		label var orub "OR upper 95% CI"

	* calculate the predicted probability of
	* illness at each level of the indicator
	* adjusted at the means of the covariates
	adjust `covarlist', se gen(xb xbse) nokey noheader
	gen xblb = xb - 1.96*xbse
	gen xbub = xb + 1.96*xbse
	gen pr = exp(xb)/(1+exp(xb))
	gen prse = pr*(1-pr)*xbse
	gen prlb = exp(xblb)/(1+exp(xblb))
	gen prub = exp(xbub)/(1+exp(xbub))
	drop xb xbse xblb xbub
	label var pr "predicted probability of illness"
	label var prse "Pr(illness) SE"
	label var prlb "Pr(illness) lower 95% CI"
	label var prub "Pr(illness) upper 95% CI"

	* collapse to unique values of the indicator
	* drop individuals who were not exposed
	keep if `touse'
	bysort `exposure': gen n = _N
	label var n "N indivs exposed"
	bysort `exposure': keep if _n==1
	drop if `exposure'==.
	rename `exposure' log10lev
		label var log10lev "Log10 indicator concentration"
	gen lev = 10^log10lev
		label var lev "Indicator concentration"
	keep lev log10lev orn or se orlb orub n pr prse prlb prub

	gen outcome = "`outcome'"
		label var outcome "Health outcome"
	gen indicator = "`exposure'"
		label var indicator "Water quality indicator code"
	order outcome indicator n orn or se orlb orub pr prse prlb prub

end pillness





*-----------------------------------------
* AVALON
*-----------------------------------------

* load data and drop incomplete interviews
* and individuals with GI illness at baseline
use "~/Dropbox/Coliphage/Data/Untouched/avalon2.dta", clear
keep if pout==1
drop if gibase==1

*-----------------------------------------
* restrict to key covariates and the 
* coliphage indicators of interest
* restrict to the quantiative phage indicators
* EPA 1601 and 1602 + / -
*-----------------------------------------

keep id diarrheaci10 hcgi3ci10 groundwater anycontact bodycontact headunder swallwater ageyrx agecat1 sex white animy sick come dig rawfood allergy multiswim year sspda14FMC0511 sspda14FPC0511 sspda16FMC0811 sspda16FPC0811 


*------------------------------------------
* set adjustment covariates for selection
*------------------------------------------
local Xvars = "i.agecat1 sex white animy sick come dig rawfood allergy multiswim year"

*------------------------------------------
* F- Coliphage EPA 1601
*------------------------------------------
*tempfile fmc1601a fmc1601b fmc1601c
*preserve
pillness if bodycontact==1, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda14FMC0511)
gen cond = "Combined"
save `fmc1601a'
restore
preserve

pillness if bodycontact==1 & groundwater==0, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda14FMC0511)
gen cond = "Groundwater below median"
save `fmc1601b'
restore

preserve
pillness if bodycontact==1 & groundwater==1, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda14FMC0511)
gen cond = "Groundwater above median"
save `fmc1601c'
restore

*------------------------------------------
* F+ Coliphage EPA 1601
*------------------------------------------
tempfile fpc1601a fpc1601b fpc1601c
preserve
pillness if bodycontact==1, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda14FPC0511)
gen cond = "Combined"
save `fpc1601a'
restore

preserve
pillness if bodycontact==1 & groundwater==0, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda14FPC0511)
gen cond = "Groundwater below median"
save `fpc1601b'
restore

preserve
pillness if bodycontact==1 & groundwater==1, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda14FPC0511)
gen cond = "Groundwater above median"
save `fpc1601c'
restore


*------------------------------------------
* F- Coliphage EPA 1602
*------------------------------------------
tempfile fmc1602a fmc1602b fmc1602c
preserve
pillness if bodycontact==1, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda16FMC0811) 
gen cond = "Combined"
save `fmc1602a'
restore

preserve
pillness if bodycontact==1 & groundwater==0, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda16FMC0811) 
gen cond = "Groundwater below median"
save `fmc1602b'
restore

preserve
pillness if bodycontact==1 & groundwater==1, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda16FMC0811) 
gen cond = "Groundwater above median"
save `fmc1602c'
restore

*------------------------------------------
* F+ Coliphage EPA 1602
*------------------------------------------
tempfile fpc1602a fpc1602b fpc1602c
preserve
pillness if bodycontact==1, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda16FPC0811) 
gen cond = "Combined"
save `fpc1602a'
restore

/* data are too sparse to estimate this association
preserve
pillness if bodycontact==1 & groundwater==0, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda16FPC0811) 
gen cond = "Groundwater below median"
save `fpc1602b'
restore
*/

pillness if bodycontact==1 & groundwater==1, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda16FPC0811) 
gen cond = "Groundwater above median"

*------------------------------------------
* append results into a single output file
*------------------------------------------
append using `fpc1602a'
* append using `fpc1602b'

append using `fpc1601a'
append using `fpc1601b'
append using `fpc1601c'

append using `fmc1602a'
append using `fmc1602b'
append using `fmc1602c'

append using `fmc1601a'
append using `fmc1601b'
append using `fmc1601c'

order outcome indicator cond lev log10lev
sort outcome indicator cond lev log10lev
label var cond "Environmental conditions / population"

label data "Avalon coliphage P(illness), created by 1-beaches-coliphage-pillness-calculations.do"
save "/Users/jadebc/Dropbox/Coliphage/Data/Temp/avalon-coliphage-pillness.dta", replace



*-----------------------------------------
* DOHENY
*-----------------------------------------

* load data and drop incomplete interviews
* and individuals with GI illness at baseline
use "/Users/jadebc/Dropbox/Coliphage/Data/Untouched/doheny2.dta", clear
keep if pout==1
drop if gibase==1

*-----------------------------------------
* restrict to key covariates and the 
* coliphage indicators of interest
* restrict to the quantiative phage indicators
* EPA 1601 and 1602 + / -
*-----------------------------------------

keep id diarrheaci10 hcgi3ci10 berm anycontact bodycontact headunder swallwater ageyrx agecat1 sex white animy sick come dig rawfood allergy multiswim year sspda14FMC0511 sspda14FPC0511 sspda16FMC0811 sspda16FPC0811 



*------------------------------------------
* set adjustment covariates for selection
*------------------------------------------
local Xvars = "i.agecat1 sex white animy sick come dig rawfood allergy multiswim year"

*------------------------------------------
* F- Coliphage EPA 1601
*------------------------------------------
tempfile fmc1601a fmc1601b fmc1601c
preserve
pillness if bodycontact==1, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda14FMC0511)
gen cond = "Combined"
save `fmc1601a'
restore
preserve

pillness if bodycontact==1 & berm==0, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda14FMC0511)
gen cond = "Berm closed"
save `fmc1601b'
restore

/* data are too sparse to estimate this association
preserve
pillness if bodycontact==1 & berm==1, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda14FMC0511)
gen cond = "Berm open"
save `fmc1601c'
restore
*/

*------------------------------------------
* F+ Coliphage EPA 1601
*------------------------------------------
tempfile fpc1601a fpc1601b fpc1601c
preserve
pillness if bodycontact==1, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda14FPC0511)
gen cond = "Combined"
save `fpc1601a'
restore

preserve
pillness if bodycontact==1 & berm==0, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda14FPC0511)
gen cond = "Berm closed"
save `fpc1601b'
restore

preserve
pillness if bodycontact==1 & berm==1, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda14FPC0511)
gen cond = "Berm open"
save `fpc1601c'
restore


*------------------------------------------
* F- Coliphage EPA 1602
*------------------------------------------
tempfile fmc1602a fmc1602b fmc1602c
preserve
pillness if bodycontact==1, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda16FMC0811) 
gen cond = "Combined"
save `fmc1602a'
restore

preserve
pillness if bodycontact==1 & berm==0, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda16FMC0811) 
gen cond = "Berm closed"
save `fmc1602b'
restore

preserve
pillness if bodycontact==1 & berm==1, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda16FMC0811) 
gen cond = "Berm open"
save `fmc1602c'
restore

*------------------------------------------
* F+ Coliphage EPA 1602
*------------------------------------------
tempfile fpc1602a fpc1602b fpc1602c
preserve
pillness if bodycontact==1, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda16FPC0811) 
gen cond = "Combined"
save `fpc1602a'
restore

preserve
pillness if bodycontact==1 & berm==0, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda16FPC0811) 
gen cond = "Berm closed"
save `fpc1602b'
restore

pillness if bodycontact==1 & berm==1, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda16FPC0811) 
gen cond = "Berm open"

*------------------------------------------
* append results into a single output file
*------------------------------------------
append using `fpc1602a'
append using `fpc1602b'

append using `fpc1601a'
append using `fpc1601b'
append using `fpc1601c'

append using `fmc1602a'
append using `fmc1602b'
append using `fmc1602c'

append using `fmc1601a'
append using `fmc1601b'
* append using `fmc1601c'

order outcome indicator cond lev log10lev
sort outcome indicator cond lev log10lev
label var cond "Environmental conditions / population"

label data "Doheny coliphage P(illness), created by 1-beaches-coliphage-pillness-calculations.do"
save "/Users/jadebc/Dropbox/Coliphage/Data/Temp/doheny-coliphage-pillness.dta", replace


*-----------------------------------------
* MALIBU
*-----------------------------------------

* load data and drop incomplete interviews
* and individuals with GI illness at baseline
use "/Users/jadebc/Dropbox/Coliphage/Data/Untouched/malibu2.dta", clear
keep if pout==1
drop if gibase==1

*-----------------------------------------
* restrict to key covariates and the 
* coliphage indicators of interest
* restrict to the quantiative phage indicators
* EPA 1601 and 1602 + / -
*-----------------------------------------

keep id diarrheaci10 hcgi3ci10 berm anycontact bodycontact headunder swallwater ageyrx agecat2 female white animy sick come dig rawfood allergy multiswim fu12 sspda14FPC0511 sspda16FMC0811 sspda16FPC0811 



*------------------------------------------
* set adjustment covariates for selection
*------------------------------------------
local Xvars = "i.agecat2 female white animy sick come dig rawfood allergy multiswim fu12"

*------------------------------------------
* F- Coliphage EPA 1601
*------------------------------------------

* NOT MEASURED AT MALIBU

*------------------------------------------
* F+ Coliphage EPA 1601
*------------------------------------------
tempfile fpc1601a fpc1601b fpc1601c
preserve
pillness if bodycontact==1, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda14FPC0511)
gen cond = "Combined"
save `fpc1601a'
restore

preserve
pillness if bodycontact==1 & berm==0, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda14FPC0511)
gen cond = "Berm closed"
save `fpc1601b'
restore

preserve
pillness if bodycontact==1 & berm==1, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda14FPC0511)
gen cond = "Berm open"
save `fpc1601c'
restore


*------------------------------------------
* F- Coliphage EPA 1602
*------------------------------------------
tempfile fmc1602a fmc1602b fmc1602c
preserve
pillness if bodycontact==1, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda16FMC0811) 
gen cond = "Combined"
save `fmc1602a'
restore

preserve
pillness if bodycontact==1 & berm==0, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda16FMC0811) 
gen cond = "Berm closed"
save `fmc1602b'
restore


pillness if bodycontact==1 & berm==1, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda16FMC0811) 
gen cond = "Berm open"
save `fmc1602c'


*------------------------------------------
* F+ Coliphage EPA 1602
*------------------------------------------

/* Data are too sparse for any estimation
tempfile fpc1602a fpc1602b fpc1602c
preserve
pillness if bodycontact==1, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda16FPC0811) 
gen cond = "Combined"
save `fpc1602a'
restore

preserve
pillness if bodycontact==1 & berm==0, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda16FPC0811) 
gen cond = "Berm closed"
save `fpc1602b'
restore

pillness if bodycontact==1 & berm==1, outcome(hcgi3ci10) xvars("`Xvars'") id(id) exposure(sspda16FPC0811) 
gen cond = "Berm open"

*/ 

*------------------------------------------
* append results into a single output file
*------------------------------------------
* append using `fpc1602a'
* append using `fpc1602b'

append using `fpc1601a'
append using `fpc1601b'
append using `fpc1601c'

append using `fmc1602a'
append using `fmc1602b'
* append using `fmc1602c'


order outcome indicator cond lev log10lev
sort outcome indicator cond lev log10lev
label var cond "Environmental conditions / population"

label data "Doheny coliphage P(illness), created by 1-beaches-coliphage-pillness-calculations.do"
save "/Users/jadebc/Dropbox/Coliphage/Data/Temp/malibu-coliphage-pillness.dta", replace







log close
exit
