capture log close
set more off
clear all

log using "~/dropbox/13beaches/src/dm/6-format-adm-wq-sample.log", text replace

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
* version 5 (15 jul 2015) by Jade
* modify script to output a dataset with 
* a row for each sample
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
use "~/dropbox/13beaches/data/untouched/adm/avalon_inddata.dta", clear
gen indcode = substr(groupindex,1,8) + substr(groupindex,13,1)
keep if inlist(indcode,"02ENT2411","02ENT2412","02ENT2413","02ENT2414","12ENT0411","14FMC0511","16FMC0811","14FPC0511","16FPC0811")
tempfile avalon
save `avalon'

use "~/dropbox/13beaches/data/untouched/adm/doheny_inddata.dta", clear
gen indcode = substr(groupindex,1,8) + substr(groupindex,13,1)
keep if inlist(indcode,"02ENT2411","02ENT2412","02ENT2413","02ENT2414","15ENT0411","14FMC0511","16FMC0811","14FPC0511","16FPC0811")
tempfile doheny
save `doheny'

use "~/dropbox/13beaches/data/untouched/adm/malibu_inddata.dta", clear
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
* and in Doheny to differentiate site C
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

* creating temporary file for use in calculating daily averages below
tempfile adm
save `adm'
	
*----------------------------------------
* Create dataset with one row for each sample
* Subset to final variables and save data
*----------------------------------------
* calculate log 10 values for the indicator counts
capture drop log10
gen log10 = log10(result)

* impute values below detection limit at 0.1
replace result = 0.1 if (qualifier=="ND" | qualifier=="<" | result==0 ) & inlist(indcode,"12ENT0411","02ENT2411","02ENT2412","02ENT2413","02ENT2414")
replace result = 0.1 if (qualifier=="<") & inlist(indcode,"14FMC0511","16FMC0811","14FPC0511","16FPC0811")

replace log10 = log10(0.1) if (qualifier=="ND" | qualifier=="<" | result==0 ) & inlist(indcode,"12ENT0411","02ENT2411","02ENT2412","02ENT2413","02ENT2414")
replace log10 = log10(0.1) if (qualifier=="<") & inlist(indcode,"14FMC0511","16FMC0811","14FPC0511","16FPC0811")

* for one value at Doheny, replace a 0.09 F+ Coliphage 1601 value with 0.1
replace log10 = log10(0.1) if indcode=="14FPC0511" & result< 0.1


rename sampledate coldate
gen marine = 1
	label var marine "Marine Beach"
order beach coldate marine

* drop columns with no water quality information
drop f26-f34 comments posneg

* drop columns that are not needed
drop year stationid

compress
label data "ADM water quality sample data, created by 6-format-adm-wq-sample.do"
save "~/dropbox/13beaches/data/final/adm-wq-samples.dta", replace

desc


*----------------------------------------
* impute values for non-detects
*----------------------------------------
use `adm', clear

* calculate log 10 values for the indicator counts
capture drop log10
gen log10 = log10(result)

* impute values below detection limit at 0.1
replace log10 = log10(0.1) if (qualifier=="ND" | qualifier=="<" | result==0 ) & inlist(indcode,"12ENT0411","02ENT2411","02ENT2412","02ENT2413","02ENT2414")

replace log10 = log10(0.1) if (qualifier=="<") & inlist(indcode,"14FMC0511","16FMC0811","14FPC0511","16FPC0811")

* for one value at Doheny, replace a 0.09 F+ Coliphage 1601 value with 0.1
replace log10 = log10(0.1) if indcode=="14FPC0511" & result< 0.1


*----------------------------------------
* calculate daily averages
*----------------------------------------

local inds "02ENT2411 02ENT2412 02ENT2413 02ENT2414 12ENT0411 14FMC0511 16FMC0811 14FPC0511 16FPC0811"
sort beach beachcode sampledate
foreach ind of local inds {
	di as res "indicator: `ind'"
	* daily average
	by beach beachcode sampledate: egen _aa = mean(log10) if indcode=="`ind'"
	by beach beachcode sampledate: egen avgdy`ind' = max(_aa)

	* am average
	by beach beachcode sampledate: egen _bb = mean(log10) if indcode=="`ind'" & sampletime2==1
	by beach beachcode sampledate: egen avgam`ind' = max(_bb)
	
	* mid-day average
	by beach beachcode sampledate: egen _cc = mean(log10) if indcode=="`ind'" & sampletime2==2
	by beach beachcode sampledate: egen avgmd`ind' = max(_cc)
	
	* pm average
	by beach beachcode sampledate: egen _dd = mean(log10) if indcode=="`ind'" & sampletime2==3
	by beach beachcode sampledate: egen avgpm`ind' = max(_dd)
	
	drop _aa _bb _cc _dd
}

keep beach beachcode sampledate avg*
by beach beachcode sampledate: keep if _n==1


*----------------------------------------
* rename and label variables to be 
* consistent with NEEAR
*----------------------------------------

local inds "02ENT2411 02ENT2412 02ENT2413 02ENT2414 12ENT0411"
local names "pcr pcr2 pcr3 pcr4 cfu"
local vars "avgdy avgam avgmd avgpm"
local labs "Daily AM Mid-Day PM"
local i = 1
foreach var of local vars {
	local lab = word("`labs'",`i')
	rename `var'02ENT2411 `var'enteropcr
		label var `var'enteropcr "`lab' Avg log10 Entero QPCR ddct"
	rename `var'02ENT2412 `var'enteropcr2
		label var `var'enteropcr2 "`lab' Avg log10 Entero QPCR ddct w/i"
	rename `var'02ENT2413 `var'enteropcr3
		label var `var'enteropcr3 "`lab' Avg log10 Entero QPCR dct"
	rename `var'02ENT2414 `var'enteropcr4
		label var `var'enteropcr4 "`lab' Avg log10 Entero QPCR dct w/i"
	rename `var'12ENT0411 `var'entero1600
		label var `var'entero1600 "`lab' Avg log10 Entero EPA 1600"
	
	rename `var'14FMC0511 `var'fmc1601
		label var `var'fmc1601 "`lab' Avg log10 F-minus Coliphage EPA 1601"
	rename `var'16FMC0811 `var'fmc1602
		label var `var'fmc1602 "`lab' Avg log10 F-minus Coliphage EPA 1602"
	rename `var'14FPC0511 `var'fpc1601
		label var `var'fpc1601 "`lab' Avg log10 F-plus Coliphage EPA 1601"
	rename `var'16FPC0811 `var'fpc1602
		label var `var'fpc1602 "`lab' Avg log10 F-plus Coliphage EPA 1602"
			
	local i = `i'+1
}


*----------------------------------------
* save file
*----------------------------------------

rename sampledate coldate
gen marine = 1
	label var marine "Marine Beach"
order beach coldate marine

compress
label data "ADM water quality data, created by 6-format-adm-wq.do"
save "~/dropbox/13beaches/data/final/adm-wq.dta", replace

desc

log close
exit





