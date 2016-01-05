capture log close
set more off
clear all

log using "~/dropbox/13beaches/src/dm/7-format-mb-wq.log", text replace

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
* version 3 (17 aug 2015)
* major re-write to generate a sample-specific dataset
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
*   mb-wq-samples.dta / .csv
*	mb-wq.dta / .csv
*
*----------------------------------------


*----------------------------------------
* read in the Mission Bay analysis dataset
* and subset it to the water quality
* indicators 
*----------------------------------------
use "~/dropbox/13beaches/data/untouched/missionbay/mb_analysis_final.dta", clear




*----------------------------------------
* get 1 obs per day / beach / station
* restrict to water quality variables
*----------------------------------------
bysort sdate station: keep if _n==1

keep sdate station sampletime_ecoli_comp - max_beach_bact_pcr

*----------------------------------------
* Need to restrict to samples slightly
* Differently because there were more
* samples for EPA qPCR than for 1600 or 
* for phage (see Colford 2007 Appendix)
*
* EPA qCPR   - 2 samples per day at each site
* Enterolert - 4 samples per station at each measuring point
* EPA 1600   - 1 composite sample per day per beach
* F- phage   - 1 composite sample per day per beach
* F+ phage   - 1 composite sample per day per beach
*----------------------------------------

* create a beach code
gen beachcode = "Mission Bay " + substr(station,1,1)
	label var beachcode "Mission bay beach code"
order beachcode station sdate


* all of the samples were collected around mid-day, so 
* but not standardized with the NEEAR or Avalon/Doheny/Malibu
* sampling times.  So just go with daily averages at this beach
sum sampletime_entero* sampletime_pcr, d


* break off the Entero qPCR data and the enterolert data
preserve
keep sdate station beachcode result_entero_pcr  *entero1* *entero2* *entero3* *entero4*

* get beach / location averages for qPCR (only 1 sample per beach/location)
gen enteroQPCRcce = result_entero_pcr
	label var enteroQPCRcce "Entero qPCR EPA 1611, CCE/100ml"
gen byte enteroQPCRcce_nd = enteroQPCRcce==0
	replace enteroQPCRcce_nd = . if enteroQPCRcce==.
	label define ND 0 "Detected" 1 "Below detection"
	label values enteroQPCRcce_nd ND
	label var enteroQPCRcce_nd "Entero qPCR EPA 1611, below detection"

tempfile enteropcr
save `enteropcr'



* now, reshape the Enterolert data to long format
drop result_entero_pcr
keep beachcode station sdate sampletime_entero* result_entero* qualifier_entero*
reshape long sampletime_entero result_entero qualifier_entero, i(beachcode station sdate) j(samplenum)

* create a sampleid variable
gen sampleid = station+"-"+string(samplenum)
	label var sampleid "Water sample ID"

* format the variable
rename result_entero enteroELTmpn
	label var enteroELTmpn "Entero Enterolert, MPN/100ml"
gen byte enteroELTmpn_nd = qualifier_entero=="<"
	label var enteroELTmpn_nd "Entero Enterolert, below detection"
	label values enteroELTmpn_nd ND

* drop empty observations and set ND values to 0 (consistent w/ other assays)
drop qualifier_entero
drop if enteroELTmpn==.
replace enteroELTmpn=0 if enteroELTmpn_nd==1
rename sampletime_entero coltime
label var samplenum "Mission Bay sample number replicate (Enterolert)"

* append the qPCR data
append using `enteropcr'
sort sdate beachcode station sampleid


keep sdate station beachcode sampleid enteroQPCRcce* enteroELTmpn*
keep if enteroQPCRcce!=. | enteroELTmpn!=.
tempfile enterodata
save `enterodata'
restore

* downsample data to just one obs per date and beach
bysort sdate beachcode: keep if _n==1

* rename the stations for these composite samples
replace station = substr(station,1,1)+"0"

* Entero 1600 processing
gen entero1600cfu = result_entero_comp
	label var entero1600cfu "Entero EPA 1600, CFU/100ml"
gen byte entero1600cfu_nd = 1 if qualifier_entero_comp=="<"
	replace entero1600cfu_nd = 0 if entero1600cfu_nd==. & entero1600cfu!=.
	label values entero1600cfu_nd ND
	label var entero1600cfu_nd "Entero EPA 1600, below detection"
	replace entero1600cfu=0 if entero1600cfu_nd==1

* F+ Coliphage 1601
gen fpc1601mpn = result_phage_m
	label var fpc1601mpn "F+ Coliphage EPA 1601, MPN/100ml"
gen byte fpc1601mpn_nd = 1 if qualifier_phage_m=="<"
	replace fpc1601mpn_nd = 0 if fpc1601mpn_nd==. & fpc1601mpn!=.
	label values fpc1601mpn_nd ND
	label var fpc1601mpn_nd "F+ Coliphage EPA 1601, below detection"
	replace fpc1601mpn=0 if fpc1601mpn_nd==1
	
* F- Coliphage 1601
gen fmc1601mpn = result_phage_s
	label var fmc1601mpn "F- Coliphage EPA 1601, MPN/100ml"
gen byte fmc1601mpn_nd = 1 if qualifier_phage_s=="<"
	replace fmc1601mpn_nd = 0 if fmc1601mpn_nd==. & fmc1601mpn!=.
	label values fmc1601mpn_nd ND
	label var fmc1601mpn_nd "F- Coliphage EPA 1601, below detection"
	replace fmc1601mpn=0 if fmc1601mpn_nd==1
	
* restrict to variables of interest
keep sdate beachcode station entero1600* fpc1601* fmc1601*

* append the Entero qPCR data
append using `enterodata'
sort sdate beachcode station sampleid



*----------------------------------------
* final variable labeling for water sample
* dataset
*----------------------------------------
gen beach = "Mission Bay"
	label var beach "Study beach"
order beach beachcode sdate station sampleid
rename sdate coldate
	label var coldate "Sample collection date"
rename station stationid

notes: Careful: this is a weird dataset because it is in long format. Read all of the notes!
notes: Values below the detection limit for all indictors are set to 0
notes: Values are missing if a sample was not tested for a particular indicator
notes stationid: Water sample station - a 0 in Mission Bay station indicates a composite sample, whereas a letter (A, B,...) indicates station location
notes stationid: Some of the stationids are replicated within day for enterolert and qPCR samples. The unique identifier in the data is coldate/stationid/samplenum
notes sampleid: Values 1-4 for enterolert assay only.
compress
label data "Mission Bay water sample data, formatted by 7-format-mb-wq.do"
save "~/dropbox/13beaches/data/final/mb-wq-samples.dta", replace
outsheet using "~/dropbox/13beaches/data/final/mb-wq-samples.csv", comma replace

desc
notes


*----------------------------------------
* calculate mean log10 values by day/beach
*----------------------------------------
use "~/dropbox/13beaches/data/final/mb-wq-samples.dta", clear


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

*----------------------------------------
*  save a daily log10 average dataset
*----------------------------------------

compress
label data "Mission Bay beach-level water quality data, created by 7-format-mb-wq.do"
save "~/dropbox/13beaches/data/final/mb-wq.dta", replace
outsheet using "~/dropbox/13beaches/data/final/mb-wq.csv", comma replace

desc

log close
exit




