capture log close
set more off
clear all

log using "~/13beaches-coliphage/src/dm/5-format-neear-wq.log", text replace

*----------------------------------------
* 5-format-neear-wq.do
* ben arnold
*
* format the NEEAR water quality data from
* Tim Wade.  label variables
*
*
*----------------------------------------

*----------------------------------------
* input files:
*
*  Entero 1600 datasets:
*	cfudata.txt
*	wqentero1600marine.txt
*	wqentero16002009.txt
*	
*  Entero qPCR datasets:
*	wqepcrddctfresh.txt
*	wqepcrddctmarine.txt
*	wqepcrddct2009.txt
*
*  Coliphage dataset:
*	wqphagemarine.dta
*
* output files:
*	neear-wq-samples.dta
*
*----------------------------------------


*----------------------------------------
* read in the enterococcus 1600 sample
* datasets
*----------------------------------------

*------------------------------
* Freshwater Entero 1600 data
insheet using "~/13beaches-data/untouched/neear/cfudata.txt", clear

* reformat collection date
gen coldate = date(collectiondate,"MDY")
	format coldate %d
	label var coldate "Sample Collection Date"

* reformat collection time
gen coltime = 8 if collectiontime=="8:00:00 AM"
	replace coltime = 11 if collectiontime=="11:00:00 AM"
	replace coltime = 15 if collectiontime=="3:00:00 PM"
	label var coltime "Sample Collection Time"

rename count_cfu entero1600cfu
	
keep  beach coldate coltime sampleid sampleloc depth entero1600cfu
order beach coldate coltime sampleid sampleloc depth entero1600cfu
sort  beach coldate coltime sampleid sampleloc

tempfile fresh1600
save `fresh1600'

*------------------------------
* Marine Entero 1600 data
insheet using "~/13beaches-data/untouched/neear/wqentero1600marine.txt", clear

* reformat collection date
gen coldate = date(collection_date,"MDY")
	format coldate %d
	label var coldate "Sample Collection Date"

rename collection_time coltime
rename sample_id sampleid
rename count entero1600cfu
	
keep  beach coldate coltime sampleid sampleloc swimloc depth entero1600cfu
order beach coldate coltime sampleid sampleloc swimloc depth entero1600cfu
sort  beach coldate coltime sampleid sampleloc

tempfile marine1600
save `marine1600'

*------------------------------
* 2009 studies Entero 1600 data
insheet using "~/13beaches-data/untouched/neear/wqentero16002009.txt", clear

* reformat collection date
gen coldate = date(collection_date,"MDY")
	format coldate %d
	label var coldate "Sample Collection Date"

* reformat collection time
gen coltime = 8 if collection_time=="08:00:00"
	replace coltime = 11 if collection_time=="11:00:00"
	replace coltime = 15 if collection_time=="15:00:00"
	label var coltime "Sample Collection Time"

rename sample_id sampleid
rename count entero1600cfu
	
keep  beach coldate coltime sampleid sampleloc swimloc depth entero1600cfu
order beach coldate coltime sampleid sampleloc swimloc depth entero1600cfu
sort  beach coldate coltime sampleid sampleloc


*------------------------------
* Append the 3 datasets
* and final format

append using `marine1600'
append using `fresh1600'

* check for duplicates and fix 2 IDs
duplicates list beach sampleid
	replace sampleid = "71150115SS" if (sampleid=="71050115SS") & beach=="SB" & coldate==mdy(7,11,2004)
	replace sampleid = "71150115W" if (sampleid=="71050115W") & beach=="WP" & coldate==mdy(7,11,2004)
	duplicates report beach sampleid

* label depth
label var depth "Sample depth"
	capture label drop depth
label define depth 1 "Shin depth" 2 "Waist depth"
label values depth depth

label var beach "Beach"
label var sampleid "Water sample ID"
	note sampleid: for NEEAR beaches this is date (3) loc (1) type(2) time (02) beach(1)
label var sampleloc "Water sample location"
label var swimloc "Swim area location"
label var entero1600cfu "Enterococcus 1600 CFU/100ml"

tempfile neearentero
save `neearentero'


*----------------------------------------
* read in the enterococcus qPCR sample
* datasets (all identical format, so easier)
* merge to the existing dataset
*----------------------------------------

insheet using "~/13beaches-data/untouched/neear/wqepcrddctfresh.txt", clear
tempfile freshqpcr
save `freshqpcr'
insheet using "~/13beaches-data/untouched/neear/wqepcrddctmarine.txt", clear
tempfile marineqpcr
save `marineqpcr'
insheet using "~/13beaches-data/untouched/neear/wqepcrddct2009.txt", clear
append using `freshqpcr'
append using `marineqpcr'


* reformat collection date
gen coldate = date(collection_date,"YMD")
	format coldate %d
	label var coldate "Sample Collection Date"

* reformat collection time
gen coltime = 8 if collection_time=="08:00:00"
	replace coltime = 11 if collection_time=="11:00:00"
	replace coltime = 15 if collection_time=="15:00:00"
	label var coltime "Sample Collection Time"


rename sample_id sampleid_neearqpcr
	label var sampleid_neearqpcr "Water Sample ID, NEEAR Entero QPCR data"

rename count enteroQPCRcce
	label var enteroQPCRcce "Entero 1611 qPCR ddct CCE/100ml"

label var beach "Beach"
label var sampleid "Water Sample ID"
label var sampleloc "Sample location"
label var depth "Sample depth"
	capture label drop depth
label define depth 1 "Shin depth" 2 "Waist depth"
label values depth depth

label var dl "Below detection limit"
	capture label drop dl
	label define dl 0 "Detected" 1 "Below detection"
	label values dl dl
label var passes_qc "Passes Salmon QC criterion"

order beach sampleid coldate coltime sampleloc depth dl passes_qc enteroQPCRcce
keep beach sampleid coldate coltime sampleloc depth dl passes_qc enteroQPCRcce

sort beach coldate coltime sampleloc depth
tempfile neearqpcr
save `neearqpcr'


* Merge
use `neearentero', clear
sort beach coldate coltime sampleloc depth
merge 1:1 beach coldate coltime sampleloc depth using `neearqpcr'
list beach coldate coltime sampleloc sampleid sampleid_neearqpcr if _merge !=3
replace sampleid = sampleid_neearqpcr if _merge !=3 & sampleid==""
drop _merge

save `neearentero', replace


*----------------------------------------
* read in the Coliphage dataset
* merge to the existing dataset
*----------------------------------------

use "~/13beaches-data/untouched/neear/wqphagemarine.dta", clear

rename sample_id sampleid_neearphage
	label var sampleid_neearphage "Water Sample ID, NEEAR Coliphage Data"
	

rename collection_date coldate
rename collection_time coltime
destring sampleloc, replace force

* restrict to non-redundant variables
keep  beach coldate coltime sampleloc depth sampleid_neearphage beach mpn_* frozen
sort beach coldate coltime sampleloc depth
tempfile neearphage
save `neearphage'

* merge to the main dataset.
use `neearentero', clear
sort beach coldate coltime sampleloc depth
merge 1:1  beach coldate coltime sampleloc depth using `neearphage'
	assert _merge !=2
	drop _merge


*----------------------------------------
* Variable cleanup and renaming in the
* combined dataset
*----------------------------------------

rename sampleid sampleid_neear1600
	label var sampleid_neear1600 "Water Sample ID, NEEAR Entero 1600 data"
order beach sampleid* coldate coltime sampleloc depth


label var frozen "Sample frozen (coliphage assays only)"
label var mpn_rna_clat5 "MPN RNA Coliphage, 5-hour Clat assay"
label var mpn_dna_clat5 "MPN DNA Coliphage, 5-hour Clat assay"
label var mpn_group1_clat5 "MPN Group 1 RNA Coliphage, 5-hour Clat assay"
label var mpn_group2_clat5 "MPN Group 2 RNA Coliphage, 5-hour Clat assay"
label var mpn_group3_clat5 "MPN Group 3 RNA Coliphage, 5-hour Clat assay"
label var mpn_group4_clat5 "MPN Group 4 RNA Coliphage, 5-hour Clat assay"
label var mpn_spot24 "F-plus Coliphage EPA 1601"
	
	rename mpn_rna_clat5 rnaclat5mpn
	rename mpn_dna_clat5 dnaclat5mpn
	rename mpn_group1_clat5 group1clat5mpn
	rename mpn_group2_clat5 group2clat5mpn
	rename mpn_group3_clat5 group3clat5mpn
	rename mpn_group4_clat5 group4clat5mpn
	rename mpn_spot24 fpc1601mpn

gen stationid = string(sampleloc)
	label var stationid "Water sample station/location ID"
	drop sampleloc
	
* make beach labels more informative
replace beach = "Avalon" if beach=="AV"
replace beach = "Boqueron" if beach=="BB"
replace beach = "Doheny" if beach=="DO"
replace beach = "Edgewater" if beach=="EB"
replace beach = "Fairhope" if beach=="FB"
replace beach = "Goddard" if beach=="GB"
replace beach = "Huntington" if beach=="HB"
replace beach = "Malibu" if beach=="MA"
replace beach = "Surfside" if beach=="MB"
replace beach = "Silver" if beach=="SB"
replace beach = "West" if beach=="WB"
replace beach = "Washington Park" if beach=="WP"



*----------------------------------------
* Flag samples that fell below
* the detection limit
*----------------------------------------

local vars "entero1600cfu fpc1601mpn rnaclat5mpn dnaclat5mpn group1clat5mpn group2clat5mpn group3clat5mpn group4clat5mpn"
capture label drop ND
label define ND 0 "Detected" 1 "Below detection"
foreach var of local vars {
	local vlab : var lab `var'
	gen `var'_nd = (`var'==0)
	replace `var'_nd = . if `var'==.
	label values `var'_nd ND
	label var `var'_nd "`vlab', below detection"
}

rename dl enteroQPCRcce_nd
rename passes_qc enteroQPCRcce_qc


order beach coldate coltime sampleid_neear1600 sampleid_neearqpcr sampleid_neearphage stationid swimloc depth frozen entero1600* enteroQPCRcce enteroQPCRcce_nd enteroQPCRcce_qc fpc1601* rnaclat5* dnaclat5* group1* group2* group3* group4*


*----------------------------------------
* save an individual water-sample dataset
*----------------------------------------
note: for Entero QPCR samples that failed QC criteron, the value was imputed based on the average of the other two samples at the same depth and time
compress
label data "NEEAR water sample data, formatted by 5-format-neear-wq.do"
save "~/13beaches-data/final/neear-wq-samples.dta", replace
outsheet using "~/13beaches-data/final/neear-wq-samples.csv", comma replace

desc
notes

log close
exit





	
