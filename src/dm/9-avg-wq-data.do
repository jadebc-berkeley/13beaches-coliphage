capture log close
set more off
clear all

log using "~/dropbox/13beaches/src/dm/9-avg-wq-data.log", text replace

*----------------------------------------
* 9-avg-wq-data.do
* ben arnold
*
* impute non-detect values in datasets
* calculate daily average log10 values
* 
* append the water quality datasets into
* a single file and calculate quantiles
* of indicator distributions
*
*
*----------------------------------------

*----------------------------------------
* input files:
*	13beaches-wq-samples.dta
*
* output files:
*	13beaches-wq.dta / .csv
*
*----------------------------------------


*----------------------------------------
* Do the averaging separately by broad
* study type because of small differences
* in data collection
*----------------------------------------

*----------------------------------------
* EPA NEEAR water quality data
* Calculate daily average values for
* Entero 1600
* Entero QPCR
* F-plus Coliphage 1601
*----------------------------------------
use "~/dropbox/13beaches/data/final/13beaches-wq-samples.dta", clear
drop if inlist(beach,"Avalon","Doheny","Malibu","Mission Bay")

* impute non-detects with 0.1
local vars "entero1600cfu enteroQPCRcce fpc1601mpn"
foreach var of local vars {
	di as res _n "Replacements of `var'"
	replace `var' = 0.1 if `var'_nd==1
}

gen logentero1600 = log10(entero1600cfu)
	label var logentero1600 "log10 Entero EPA 1600"
gen logenteropcr = log10(enteroQPCRcce)
	label var logenteropcr "log10 Entero qPCR EPA 1611"
gen logfpc1601 = log10(fpc1601mpn)
	label var logfpc1601 "log10 F+ Coliphage EPA 1601"
	
local stubs "entero1600 enteropcr fpc1601"
foreach stub of local stubs {
	sort beach coldate
	by beach coldate: egen _dy = mean(log`stub')
	by beach coldate: egen _am = mean(log`stub') if coltime==8
	by beach coldate: egen _md = mean(log`stub') if coltime==11
	by beach coldate: egen _pm = mean(log`stub') if coltime==15
	by beach coldate: egen _sh = mean(log`stub') if depth==1
	by beach coldate: egen _wt = mean(log`stub') if depth==2
	by beach coldate: egen avgdy`stub' = max(_dy)
	by beach coldate: egen avgam`stub' = max(_am)
	by beach coldate: egen avgmd`stub' = max(_md)
	by beach coldate: egen avgpm`stub' = max(_pm)
	by beach coldate: egen avgsh`stub' = max(_sh)
	by beach coldate: egen avgwt`stub' = max(_wt)
	drop _dy - _wt
	local stublab : var label log`stub'
	label var avgdy`stub' "Daily Avg `stublab'"
	label var avgam`stub' "AM Avg `stublab'"
	label var avgmd`stub' "Mid-Day Avg `stublab'"
	label var avgpm`stub' "PM Avg `stublab'"
	label var avgsh`stub' "shin Avg `stublab'"
	label var avgwt`stub' "waist Avg `stublab'"
}


by beach coldate: keep if _n==1
keep beach coldate avg*

gen byte marine = !inlist(beach,"West","Washington Park","Silver","Huntington")
	label var marine "Marine beach"

order beach coldate marine
sort beach coldate

* store a daily average dataset
tempfile neearwq
save `neearwq'



*----------------------------------------
* Avalon / Doheny / Malibu (ADM) water quality data
* Calculate daily average values for
* Entero 1600
* Entero QPCR 1611 ddct
* F-plus Coliphage 1601
* F-plus Coliphage 1602
* F-minus Coliphage 1601
* F-minus Coliphage 1602
*----------------------------------------
use "~/dropbox/13beaches/data/final/13beaches-wq-samples.dta", clear
keep if inlist(beach,"Avalon","Doheny","Malibu")

* impute non-detects with 0.1
local vars "entero1600cfu enteroQPCRcce fpc1601mpn fpc1602mpn fmc1601mpn fmc1602mpn"
foreach var of local vars {
	di as res _n "Replacements of `var'"
	replace `var' = 0.1 if `var'_nd==1
}

gen logentero1600 = log10(entero1600cfu)
	label var logentero1600 "log10 Entero EPA 1600"
gen logenteropcr = log10(enteroQPCRcce)
	label var logenteropcr "log10 Entero qPCR EPA 1611"
gen logfpc1601 = log10(fpc1601mpn)
	label var logfpc1601 "log10 F+ Coliphage EPA 1601"
gen logfpc1602 = log10(fpc1602mpn)
	label var logfpc1602 "log10 F+ Coliphage EPA 1602"
gen logfmc1601 = log10(fmc1601mpn)
	label var logfmc1601 "log10 F- Coliphage EPA 1601"
gen logfmc1602 = log10(fmc1602mpn)
	label var logfmc1602 "log10 F- Coliphage EPA 1602"
	
local stubs "entero1600 enteropcr fpc1601 fpc1602 fmc1601 fmc1602"
foreach stub of local stubs {
	sort beach beachcode coldate
	by beach beachcode coldate: egen _dy = mean(log`stub')
	by beach beachcode coldate: egen _am = mean(log`stub') if coltime==8
	by beach beachcode coldate: egen _md = mean(log`stub') if coltime==11
	by beach beachcode coldate: egen _pm = mean(log`stub') if coltime==15
	* not doing shin + waist summaries for these beaches - waist only collected
	* for a few samples from Avalon
	* by beach beachcode coldate: egen _sh = mean(log`stub') if depth==1
	* by beach beachcode coldate: egen _wt = mean(log`stub') if depth==2
	by beach beachcode coldate: egen avgdy`stub' = max(_dy)
	by beach beachcode coldate: egen avgam`stub' = max(_am)
	by beach beachcode coldate: egen avgmd`stub' = max(_md)
	by beach beachcode coldate: egen avgpm`stub' = max(_pm)
	* by beach beachcode coldate: egen avgsh`stub' = max(_sh)
	* by beach beachcode coldate: egen avgwt`stub' = max(_wt)
	drop _dy - _pm
	local stublab : var label log`stub'
	label var avgdy`stub' "Daily Avg `stublab'"
	label var avgam`stub' "AM Avg `stublab'"
	label var avgmd`stub' "Mid-Day Avg `stublab'"
	label var avgpm`stub' "PM Avg `stublab'"
	* label var avgsh`stub' "shin Avg `stublab'"
	* label var avgwt`stub' "waist Avg `stublab'"
}


by beach beachcode coldate: keep if _n==1
keep beach beachcode coldate avg*

gen byte marine = 1
	label var marine "Marine beach"

order beach beachcode coldate marine
sort beach beachcode coldate

* store a daily average dataset
tempfile admwq
save `admwq'



*----------------------------------------
* Mission Bay water quality data
* Calculate daily average values for
* Entero 1600
* Entero enterolert
* Entero QPCR
* F-plus Coliphage 1601
* F-minus Coliphage 1601
*----------------------------------------
use "~/dropbox/13beaches/data/final/13beaches-wq-samples.dta", clear
keep if inlist(beach,"Mission Bay")

* impute non-detects as 0.1 (consistent w/ other datasets)
local vlist "entero1600cfu enteroQPCRcce enteroELTmpn fpc1601mpn fmc1601mpn"
foreach var of local vlist {
	di as res _n "Imputing NDs for `var'"
	replace `var' = 0.1 if `var'_nd==1
}

* calculate log10 values for each indicator
gen logentero1600 = log10(entero1600cfu)
	label var logentero1600 "log10 Entero EPA 1600"
gen logenteropcr = log10(enteroQPCRcce)
	label var logenteropcr "log10 Entero qPCR EPA 1611"
gen logenteroELT = log10(enteroELTmpn)
	label var logenteroELT "log10 Entero Enterolert"
gen logfpc1601 = log10(fpc1601mpn)
	label var logfpc1601 "log10 F+ Coliphage EPA 1601"
gen logfmc1601 = log10(fmc1601mpn)
	label var logfmc1601 "log10 F- Coliphage EPA 1601"
	
* calculate daily average log10 values
local stubs "entero1600 enteropcr enteroELT fpc1601 fmc1601"
foreach stub of local stubs {
	sort beach beachcode coldate
	by beach beachcode coldate: egen _dy = mean(log`stub')
	by beach beachcode coldate: egen avgdy`stub' = max(_dy)
	drop _dy
	local stublab : var label log`stub'
	label var avgdy`stub' "Daily Avg `stublab'"
}

* keep 1 obs per beach and day, and subset to relevant variables
bysort beach beachcode coldate: keep if _n==1

* drop 3 days (may 24, 25, 26) with completely empty data
drop if coldate < mdy(5,27,2003)

keep beach beachcode coldate avgdy*

gen marine = 1
	label var marine "Marine beach"
order beach beachcode marine coldate 


* store a daily average dataset
tempfile mbwq
save `mbwq'

*----------------------------------------
* append averaged indicator data
*----------------------------------------

use `admwq', clear
append using `mbwq'
order beach beachcode marine coldate
append using `neearwq'


*----------------------------------------
* calculate quartiles of the indicator
* distribution
*----------------------------------------

* Enterococcus EPA 1600 quartiles
* exclude site C from Malibu and Doheny from the
* calculation because they are very different
*
* use Enterolert data from Mission Bay since that is used in the
* final calculation
*
gen entero = avgdyentero1600
replace entero = avgdyenteroELT if beach=="Mission Bay"
sum entero
replace entero = . if inlist(beachcode,"Doheny-C","Malibu-C")
sum entero, d
	local p25 = r(p25)
	local p50 = r(p50)
	local p75 = r(p75)
replace entero = avgdyentero1600 if inlist(beachcode,"Doheny-C","Malibu-C")
gen int qavgdyentero1600 = 1
	replace qavgdyentero1600 = 2 if (entero> `p25') & (entero<= `p50')
	replace qavgdyentero1600 = 3 if (entero> `p50') & (entero<= `p75')
	replace qavgdyentero1600 = 4 if (entero> `p75') 
	replace qavgdyentero1600 = . if entero==.
	label var qavgdyentero1600 "Quartiles of Enterococcus EPA 1600"
	tab qavgdyentero1600 if !inlist(beachcode,"Doheny-C","Malibu-C")
	tab qavgdyentero1600
	drop entero 

* Enterococcus QPCR quartiles
* exclude site C from Malibu and Doheny from the
* calculation because they are very different
sum avgdyenteropcr, d
sum avgdyenteropcr if !inlist(beachcode,"Doheny-C","Malibu-C"), d
	local p25 = r(p25)
	local p50 = r(p50)
	local p75 = r(p75)
gen int qavgdyenteropcr = 1
	replace qavgdyenteropcr = 2 if (avgdyenteropcr> `p25') & (avgdyenteropcr <= `p50')
	replace qavgdyenteropcr = 3 if (avgdyenteropcr> `p50') & (avgdyenteropcr <= `p75')
	replace qavgdyenteropcr = 4 if (avgdyenteropcr> `p75') 
	replace qavgdyenteropcr = . if avgdyenteropcr ==.
	label var qavgdyenteropcr "Quartiles of Enterococcus EPA qPCR 1611"
	tab qavgdyenteropcr if !inlist(beachcode,"Doheny-C","Malibu-C")
	tab qavgdyenteropcr


*----------------------------------------
* save file
*----------------------------------------
note: 13-Beaches Averaged Water Quality Data, created by 9-avg-wq-data.do
note: samples below detection were imputed with 0.1 before calculating log10 averages
compress
note: All indicators are means of log10 values, summarized daily over all samples, or separately for AM, Mid-Day, PM, shin-depth, or waist-depth samples 
note: The enteropcr2-4 variables are alternative Enterococcus qPCR calculations at Avalon/Doheny/Malibu beaches only
note: The Mission Bay water quality data are separate for 6 beaches within the area and need to merge to individuals separately by beachcode
note: Avalon water quality data are summarized separately for sites A/B/C and D; site D is very different. See Yau et al 2014 Figs 1 and 2.
note: Doheny and Malibu water quality data are summarized separately for sites C; at both beaches site C was located in a lagoon. At Doheny, site E is also summarized separately, since it was nearly 1 mile from the other sampling sites (See Fig 1 of Arnold et al. 2013 and SI Fig 1 of Colford et al. 2012)
note: Quartile categories for Entero 1600 and qPCR exclude water quality measurements from Doheny and Malibu sites C in the lagoons because they are not representative of the broader swimmer exposure. The Entero 1600 quartiles were calculated using Enterolert data from Mission Bay
note: AM/Mid-Day/PM and shin/waist disaggregations were not available for Mission Bay due to different sampling methods at that beach


compress
label data "13 beaches water quality data, created by 9-avg-wq-data.do"
saveold "~/dropbox/13beaches/data/final/13beaches-wq.dta", replace version(12)
outsheet using "~/dropbox/13beaches/data/final/13beaches-wq.csv", comma replace

desc


* write a codebook to separate file
log close
log using "~/dropbox/13beaches/data/final/13beaches-wq-codebook.txt", text replace
desc, s
notes
codebook
log close

exit





