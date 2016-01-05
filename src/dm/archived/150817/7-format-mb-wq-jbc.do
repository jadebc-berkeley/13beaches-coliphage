capture log close
set more off
clear all

log using "~/dropbox/13beaches/src/dm/7-format-mb-wq-jbc.log", text replace

*----------------------------------------
* 7-format-mb-wq.do
* ben arnold
*
* format the Mission Bay water quality data
* water quality data.
* 
* calculate log averages using a standardized
* imputation approach for values at 
* the detection limits
* for values below LOD, use LOD/2
* for values at upper LOD, use upper LOD
*
* version 3 (15 jul 2015) by Jade
* modify script to output a dataset with 
* a row for each sample
*
* version 2 (21 feb 2015)
* updated the beach code
*
* version 1 (3 feb 2015)
*
*----------------------------------------

*----------------------------------------
* input files:
*	mb_analysis_final.dta
*
* output files:
*	mb-wq.dta
*
*----------------------------------------


*----------------------------------------
* read in the Mission Bay analysis dataset
* and subset it to the water quality
* indicators 
*----------------------------------------
use "~/dropbox/13beaches/data/untouched/missionbay/mb_analysis_final.dta", clear

*----------------------------------------
* Create dataset with one row for each sample
* Subset to final variables and save data
*----------------------------------------

*----------------------------------------
* subset to the water quality indicators
* of interest, and reduce the data to
* 1 obs per day of sampling for each beach
*----------------------------------------

* note: for the enterococcus measurements, its unclear what samples 1-4 indicate
* so, we are carrying forward with the daily geometric mean for the entire beach

keep sdate station g_mean_beach_entero g_mean_beach_entero_pcr result_phage_m qualifier_phage_m result_phage_s qualifier_phage_s 

gen beachcode = "Mission Bay " + substr(station,1,1)
	label var beachcode "Mission bay beach code"
	drop station
	
bysort sdate beachcode: keep if _n==1
order sdate beachcode
rename sdate coldate
label var coldate "Sample collection date"

gen beach = "Mission Bay"
	label var beach "Study Beach"
	
gen marine = 1
	label var marine "Marine Beach"
	
ren result_phage_m fpc1601 
ren result_phage_s fmc1601
ren qualifier_phage_m qualfpc1601
ren qualifier_phage_s qualfmc1601

order beach coldate beachcode marine

* drop averaged variables
drop g_mean*

* drop days/beach sites with no water quality information
drop if (fpc1601==.) & (fmc1601==.)

* convert from wide to long format
foreach var of varlist fpc1601 fmc1601{
	ren `var' result`var'
}

reshape long result qual, i(beach coldate beachcode) j(groupindex) string

drop if result==.

* generate variables to match other datasets
gen indicator=""
replace indicator="Fminus coliphage" if groupindex=="fmc1601"
replace indicator="Fplus coliphage" if groupindex=="fpc1601"

gen analysismethod=""
replace analysismethod="EPA 1601" if groupindex=="fmc1601" | groupindex=="fpc1601"

gen log10 = log10(result)
	
drop groupindex

compress
label data "Mission Bay beach-level water quality sample data, created by 7-format-mb-wq-sample.do"
save "~/dropbox/13beaches/data/final/mb-wq-samples.dta", replace

desc

*----------------------------------------
* subset to the water quality indicators
* of interest, and reduce the data to
* 1 obs per day of sampling for each beach
*----------------------------------------
use "~/dropbox/13beaches/data/untouched/missionbay/mb_analysis_final.dta", clear

* note: for the enterococcus measurements, its unclear what samples 1-4 indicate
* so, we are carrying forward with the daily geometric mean for the entire beach

keep sdate station g_mean_beach_entero g_mean_beach_entero_pcr result_phage_m qualifier_phage_m result_phage_s qualifier_phage_s 

gen beachcode = "Mission Bay " + substr(station,1,1)
	label var beachcode "Mission bay beach code"
	drop station
	
bysort sdate beachcode: keep if _n==1
order sdate beachcode



*----------------------------------------
* Enterococcus EPA 1600 and QPCR
*----------------------------------------
gen avgdyentero1600 = log10(g_mean_beach_entero)
	label var avgdyentero1600 "Daily Avg log10 Entero EPA 1600"

gen avgdyenteropcr = log10(g_mean_beach_entero_pcr)
	label var avgdyenteropcr "Daily Avg log10 Entero EPA qPCR"


*----------------------------------------
* Coliphage EPA 1601
*----------------------------------------
* recode the phage results as 0.1 for non-detects (consistent with NEEAR)
replace result_phage_m = 0.1 if qualifier_phage_m=="<"
replace result_phage_s = 0.1 if qualifier_phage_s=="<"

gen avgdyfpc1601 = log10(result_phage_m)
	label var avgdyfpc1601 "Daily Avg log10 F-plus Coliphage EPA 1601"
gen avgdyfmc1601 = log10(result_phage_s)
	label var avgdyfmc1601 "Daily Avg log10 F-minus Coliphage EPA 1601"


*----------------------------------------
* Subset to final variables and save data
*----------------------------------------
rename sdate coldate
label var coldate "Sample collection date"

gen beach = "Mission Bay"
	label var beach "Study Beach"
	
gen marine = 1
	label var marine "Marine Beach"

keep beach coldate beachcode marine avgdy*
order beach coldate beachcode marine

* drop days/beach sites with no water quality information
drop if (avgdyentero1600==.) & (avgdyenteropcr==.) & (avgdyfpc1601==.) & (avgdyfmc1601==.)



compress
label data "Mission Bay beach-level water quality data, created by 7-format-mb-wq.do"
save "~/dropbox/13beaches/data/final/mb-wq.dta", replace

desc

log close
exit




