capture log close
set more off
clear all

log using "~/dropbox/13beaches/src/dm/6-format-adm-wq.log", text replace

*----------------------------------------
* 6-format-adm-wq.do
* ben arnold
*
* format the Avalon, Doheny, Malibu
* water quality data.
* 
* calculate log averages using a standardized
* imputation approach for values at 
* the detection limits
* for values below LOD, use 0.1, consistent with NEEAR data
* for values at upper LOD, use upper LOD
*
* Only include Entero 1600 and QPCR
* for comparability to the NEEAR data
*
* version 5 (7 aug 2015)
* major re-write to output an individual sample dataset
*
* version 4 (17 apr 2015)
* Further split out site E at Doheny as a separate water quality site
*
* version 3 (21 feb 2015)
* created a separate site-specific match for each beach. At Avalon
* since site D is so different + geographically separate from A/B/C 
* (see Figs 1 and 2 of Yau et al. 2014)
* At Doheny and Malibu, Site C was in the lagoon.
*
* version 2 (13 feb 2015)
* add Coliphage indicator variables for EPA 1601 and 1602
*
* version 1 (25 feb 2014)
*
*----------------------------------------

*----------------------------------------
* input files:
*	avalon_inddata.dta
*   doheny_inddata.dta
*   malibu_inddata.dta
*
* output files:
*	adm-wq.dta
*
*----------------------------------------


*----------------------------------------
* read in the beaches indicator data
* create an indicator code that 
* excludes beach and year info
*----------------------------------------
use "~/dropbox/beaches/indicators/avalon_inddata.dta", clear
gen indcode = substr(groupindex,1,8) + substr(groupindex,13,1)
keep if inlist(indcode,"02ENT2411","02ENT2412","02ENT2413","02ENT2414","12ENT0411","14FMC0511","16FMC0811","14FPC0511","16FPC0811")
tempfile avalon
save `avalon'

use "~/dropbox/beaches/indicators/doheny_inddata.dta", clear
gen indcode = substr(groupindex,1,8) + substr(groupindex,13,1)
keep if inlist(indcode,"02ENT2411","02ENT2412","02ENT2413","02ENT2414","15ENT0411","14FMC0511","16FMC0811","14FPC0511","16FPC0811")
tempfile doheny
save `doheny'

use "~/dropbox/beaches/indicators/malibu_inddata.dta", clear
gen indcode = substr(groupindex,1,8) + substr(groupindex,13,1)
keep if inlist(indcode,"02ENT2411","02ENT2412","02ENT2413","02ENT2414","12ENT0411","14FMC0511","16FMC0811","14FPC0511","16FPC0811")

append using `avalon'
append using `doheny'

* standardize the EPA 1600 indicator code for Doheny
replace indcode = "12ENT0411" if indcode=="15ENT0411"

tab indcode beach

*----------------------------------------
* Create a standard beach code for Avalon
* to differentiate sites A/B/C from site D
* and in Doheny/Malibu to differentiate site C
* (in the lagoon) from the other sites
*----------------------------------------
gen str beachcode = ""
	replace beachcode = "Avalon-ABC" if beach=="Avalon" & stationid!="A_D"
	replace beachcode = "Avalon-D" if beach=="Avalon" & stationid=="A_D"
	replace beachcode= "Doheny-ABD" if beach=="Doheny" & (stationid!="D_C" & stationid!="D_E")
	replace beachcode= "Doheny-C" if beach=="Doheny" & stationid=="D_C"
	replace beachcode= "Doheny-E" if beach=="Doheny" & stationid=="D_E"
	replace beachcode= "Malibu-ABDE" if beach=="Malibu" & stationid!="M_C"
	replace beachcode= "Malibu-C" if beach=="Malibu" & stationid=="M_C"
	assert beachcode!= ""
	label var beachcode "Water quality sampling location beach code"

	
	
*----------------------------------------
* create some standard variables
* consistent with the NEEAR water quality
* dataset
*----------------------------------------
rename sampledate coldate
	
gen coltime = 8 if inlist(sampletime,"07:00:00","08:00:00")
	replace coltime = 11 if inlist(sampletime,"12:00:00","13:00:00")
	replace coltime = 15 if inlist(sampletime,"15:00:00")
	label var coltime "Sample collection time (approx)"
	

gen depth = 1 if sampledepth==0.5
	replace depth = 2 if sampledepth==1
	label var depth "Sample depth"
	capture drop label depth
	label define depth 1 "Shin depth" 2 "Waist depth"
	label values depth depth

	
*----------------------------------------
* There are 35 duplicate samples in
* Malibu that are lab replicates
* Identify them, drop any that are NDs
* and then average the others
*----------------------------------------
duplicates report beach stationid beachcode coldate coltime depth indcode
duplicates tag beach stationid beachcode coldate coltime depth indcode, gen(dups)
tab beach dups
drop if dups==1 & qualifier=="<"
bysort beach stationid beachcode coldate coltime depth indcode: egen _x = mean(result)
replace result = _x if dups==1
bysort beach stationid beachcode coldate coltime depth indcode: keep if _n==1
drop _x
duplicates report beach stationid beachcode coldate coltime depth indcode


*----------------------------------------
* reshape the data from long to wide
*----------------------------------------
rename qualifier nd
keep beach stationid beachcode sampleid coldate coltime depth nd result indcode

reshape wide result nd, i(beach stationid beachcode sampleid coldate coltime depth) j(indcode, string)


*----------------------------------------
* rename the variables so they are no longer encoded
*----------------------------------------

local codenames "02ENT2411 02ENT2412 02ENT2413 02ENT2414 12ENT0411 14FMC0511 16FMC0811 14FPC0511 16FPC0811"
local newnames "enteroQPCRcce enteroQPCRcce2 enteroQPCRcce3 enteroQPCRcce4 entero1600cfu fmc1601mpn fmc1602mpn fpc1601mpn fpc1602mpn"
capture label drop ND
label define ND 0 "Detected" 1 "Below detection"
local i = 1
foreach var of local codenames {
	local newname = word("`newnames'",`i')
	gen `newname' =  result`var'
	gen byte `newname'_nd = 1 if (nd`var'!="" & nd`var'!=">") | (`newname'==0)
		replace `newname'_nd = 0 if `newname'_nd==. & `newname' !=.
		label values `newname'_nd ND
	replace `newname' = 0 if `newname'_nd==1
	local i = `i'+1
}

drop nd* result*


label var enteroQPCRcce  "Entero EPA 1611 qPCR ddct CCE/100ml"
label var enteroQPCRcce2 "Entero EPA 1611 qPCR ddct CCE/100ml (w/inhib)"
label var enteroQPCRcce3 "Entero EPA 1611 qPCR dct CCE/100ml"
label var enteroQPCRcce4 "Entero EPA 1611 qPCR dct CCE/100ml (w/inhib)"
label var entero1600cfu  "Entero EPA 1600 CFU/100ml"
label var fpc1601mpn "F+ Coliphage EPA 1601 MPN/100ml"
label var fpc1602mpn "F+ Coliphage EPA 1602 MPN/100ml"
label var fmc1601mpn "F- Coliphage EPA 1601 MPN/100ml"
label var fmc1602mpn "F- Coliphage EPA 1602 MPN/100ml"
foreach var of local newnames {
	local vlab : var lab `var'
	label var `var'_nd "`vlab', below detection"
}


*----------------------------------------
* Final variable cleanup / reorder
*----------------------------------------

* for one value at Doheny, replace a 0.09 F+ Coliphage 1601 value with ND
replace fpc1601mpn_nd = 1 if (fpc1601mpn>0 & fpc1601mpn<0.1) & beach=="Doheny"
replace fpc1601mpn    = 0 if (fpc1601mpn>0 & fpc1601mpn<0.1) & beach=="Doheny"


order beach beachcode stationid sampleid coldate coltime depth
notes: Values below the detection limit for all indictors are set to 0
notes: Values are missing if a sample was not tested for a particular indicator
label data "Avalon/Doheny/Malibu water sample data, formatted by 6-format-adm-wq.do"
save "~/dropbox/13beaches/data/final/adm-wq-samples.dta", replace
outsheet using "~/dropbox/13beaches/data/final/adm-wq-samples.csv", comma replace


*----------------------------------------
* Calculate daily average values for
* Entero 1600
* Entero QPCR 1611 ddct
* F-plus Coliphage 1601
* F-plus Coliphage 1602
* F-minus Coliphage 1601
* F-minus Coliphage 1602
*----------------------------------------

use "~/dropbox/13beaches/data/final/adm-wq-samples.dta", clear


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

*----------------------------------------
* save a daily average dataset
*----------------------------------------
note: samples below detection were imputed with 0.1 before calculating log10 averages
compress
label data "Avalon/Doheny/Malibu daily average water quality data, created by 6-format-adm-wq.do"
save "~/dropbox/13beaches/data/final/adm-wq.dta", replace
outsheet using "~/dropbox/13beaches/data/final/adm-wq.csv", comma replace

desc

log close
exit






