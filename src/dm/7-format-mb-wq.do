capture log close
set more off
clear all

log using "~/Documents/CRG/coliphage/13beaches-coliphage/src/dm/7-format-mb-wq.log", text replace


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
* version 4 (17 sep 2015)
* caught a small number of 0 qPCR samples
* that should be flagged as NDs but were not
*
* version 3 (19 aug 2015)
* major re-write to generate a sample-specific dataset
* 
* version 2 (21 feb 2015)
* updated the beach code
*
* version 1 (3 feb 2015)
*
*----------------------------------------



*------------------------------------------
* load the culture method sample data
* and prune out empty rows
*------------------------------------------

insheet using "~/Documents/CRG/coliphage/13beaches-data/untouched/missionbay/wq/tblMicrobiologyData.csv", comma clear

ren station stationid

* confirm Entero EPA 1600 and E.coli IDEXX and Total Coliforms IDEXX 
* were only run on composite samples, then drop non-composite records
assert result == -99 if analysismethod=="EPA 1600" & strpos(stationid,"COMP")==0
assert result == -99 if parametercode=="E. Coli" & analysismethod=="Colilert (96 Well Tray)" & strpos(stationid,"COMP")==0
assert result == -99 if parametercode=="Total Coliforms" & analysismethod=="Colilert (96 Well Tray)" & strpos(stationid,"COMP")==0

drop if analysismethod=="EPA 1600" & strpos(stationid,"COMP")==0
drop if parametercode=="E. Coli" & analysismethod=="Colilert (96 Well Tray)" & strpos(stationid,"COMP")==0
drop if parametercode=="Total Coliforms" & analysismethod=="Colilert (96 Well Tray)" & strpos(stationid,"COMP")==0

* confirm Entero Enterolert IDEXX, MF-Fecal Colif and MF-Total COlif 
* were not run on composite samples, then drop those composite records
assert result == -99 if parametercode=="Enterococcus" & analysismethod=="Colilert (96 Well Tray)" & strpos(stationid,"COMP")>0
assert result == -99 if parametercode=="Fecal Coliforms" & analysismethod=="MF-Fecal Coliform (APHA 9" & strpos(stationid,"COMP")>0
assert result == -99 if parametercode=="Total Coliforms" & analysismethod=="MF-Total Coliform (APHA 9" & strpos(stationid,"COMP")>0

drop if parametercode=="Enterococcus" & analysismethod=="Colilert (96 Well Tray)" & strpos(stationid,"COMP")>0
drop if parametercode=="Fecal Coliforms" & analysismethod=="MF-Fecal Coliform (APHA 9" & strpos(stationid,"COMP")>0
drop if parametercode=="Total Coliforms" & analysismethod=="MF-Total Coliform (APHA 9" & strpos(stationid,"COMP")>0

* extract sample date
gen coldate = date(substr(sampledate,1,8),"MDY",2020)
	format coldate %d
	label var coldate "Sample collection date"

	
*------------------------------------------
* Create indicator stub names and
* then reshape the data to wide format
*------------------------------------------	

keep stationid coldate sampletime parametercode qualifier result units analysismethod
order stationid coldate sampletime parametercode qualifier result units analysismethod

* compress indicator and analysis method into a single field for reshaping
gen ind = "entero1600cfu" if analysismethod=="EPA 1600"
	replace ind = "enteroELTmpn" if parametercode=="Enterococcus" & analysismethod=="Colilert (96 Well Tray)"
	replace ind = "ecoliCLTmpn" if parametercode=="E. Coli" & analysismethod=="Colilert (96 Well Tray)"
	replace ind = "tcolCLTmpn" if parametercode=="Total Coliforms" & analysismethod=="Colilert (96 Well Tray)"
	replace ind = "fcolifcfu" if parametercode=="Fecal Coliforms" & analysismethod=="MF-Fecal Coliform (APHA 9"
	replace ind = "tcolifcfu" if parametercode=="Total Coliforms" & analysismethod=="MF-Total Coliform (APHA 9"
	assert ind!=""
	
keep stationid coldate sampletime ind result qualifier
order stationid coldate sampletime ind result qualifier

sort stationid coldate sampletime 
reshape wide result qualifier, i(stationid coldate sampletime) j(ind, string)

*------------------------------------------
* Flag non-detect values and set them to 0
* to be consistent with other datasets
*------------------------------------------	

label define ND 0 "Detected" 1 "Below detection"
local vars "entero1600cfu enteroELTmpn ecoliCLTmpn tcolCLTmpn fcolifcfu tcolifcfu"
foreach var of local vars  {
	rename result`var' `var'
	gen `var'_nd = 1 if inlist(qualifier`var',"<","ND")
		replace `var'_nd = 0 if (`var'!=.) & (`var'_nd==.)
		label values `var'_nd ND
	replace `var' = 0 if `var'_nd==1
		
	
}
drop qualifier*
order stationid coldate sampletime enteroELT* fcolif* tcolif* entero1600* ecoliCLT*

*------------------------------------------
* Label variables and add notes
*------------------------------------------

* data cleaning: there is a single sample with a missing collection time
* the time is 1245 based on the qPCR data (identified when merging)
replace sampletime = 1245 if sampletime==. & stationid=="LL3C" & coldate==mdy(6,29,2003)

gen sampletype = "Hourly sample"
	replace sampletype = "Daily composite sample" if strpos(stationid,"COMP")>0
	assert sampletype!=""
	label var sampletype "Water sample type"
	tab stationid sampletype

rename sampletime coltime
label var coltime "Sample collection time"
label var enteroELTmpn "Enterococcus, Enterolert MPN/100ml"
label var enteroELTmpn_nd "Enterococcus, Enterolert below detection"
label var fcolifcfu "Fecal coliform, APHA 9222D CFU/100ml"
label var fcolifcfu_nd "Fecal coliform, APHA 9222D below detection"
label var tcolifcfu "Total coliform, APHA 9222B CFU/100ml"
label var tcolifcfu_nd "Total coliform, APHA 9222B below detection"
label var entero1600cfu "Enterococcus, EPA 1600 CFU/100ml"
label var ecoliCLTmpn "E. coli, Colilert MPN/100ml"
label var ecoliCLTmpn_nd "E. coli, Colilert below detection"
label var tcolCLTmpn "Total coliform, Colilert MPN/100ml"
label var tcolCLTmpn_nd "Total coliform, Colilert below detection"

tempfile culture
save `culture'

*------------------------------------------
* load the qPCR sample data
*------------------------------------------

insheet using "~/Documents/CRG/coliphage/13beaches-data/untouched/missionbay/wq/PCR Data All.csv", comma clear

* extract sample date
gen coldate = date(substr(sampledate,1,8),"MDY",2020)
	format coldate %d
	label var coldate "Sample collection date"

	
*------------------------------------------
* There are 34 sets of duplicates in the data
* Aug 19, 2015
* Discussed this issue with John Griffith
* at SCCWRP, and he said that if there are
* duplicates that they are replicates and
* should be retained in the dataset
* (not dropped).  He suggested taking the
* average of the replicates, so that's
* what we are doing here, imputing at 0.1
* for non-detects before averaging.
*------------------------------------------	

* Identify the duplicates
ren station stationid
sort stationid sampledate sampletime parametercode
duplicates report stationid sampledate sampletime parametercode
duplicates tag stationid sampledate sampletime parametercode, gen(dups)

* calculate the average for each stationid / date / time / assay
gen impres = result
	replace impres = 0.1 if inlist(qualifier,"<","ND")
gen logimpres = log(impres)
by stationid sampledate sampletime parametercode: egen mulogimpres = mean(logimpres)
replace result = exp(mulogimpres) if dups==1
by stationid sampledate sampletime parametercode: keep if _n==1

* flag only values as non-detect if both samples were non-detect
replace qualifier = "" if dups==1 & (result <0.09 | result>0.11)
drop *impres* dups


duplicates report stationid sampledate sampletime parametercode


*------------------------------------------
* Create indicator stub names and
* then reshape the data to wide format
*------------------------------------------	
gen ind = "enteroQPCRcce" if parametercode=="Enterococcus faecalis"
	replace ind = "bactQPCRcce" if parametercode=="Bacteroides"
	assert ind!=""
	
keep stationid coldate sampletime ind result qualifier
order stationid coldate sampletime ind result qualifier
reshape wide result qualifier, i(stationid coldate sampletime) j(ind, string)


*------------------------------------------
* Flag non-detect values and set them to 0
* to be consistent with other datasets
* assume samples with values of 0 are NDs
*------------------------------------------	

local vars "enteroQPCRcce bactQPCRcce"
foreach var of local vars  {
	rename result`var' `var'
	gen `var'_nd = 1 if inlist(qualifier`var',"<","ND") | (`var'==0)
		replace `var'_nd = 0 if (`var'!=.) & (`var'_nd==.)
		label values `var'_nd ND
	replace `var' = 0 if `var'_nd==1
		
	
}
drop qualifier*
order stationid coldate sampletime enteroQPCRcce* bactQPCRcce*

*------------------------------------------
* Label variables and add notes
*------------------------------------------	

gen sampletype = "qPCR sample"
	label var sampletype "Water sample type"


rename sampletime coltime
label var coltime "Sample collection time"
label var enteroQPCRcce "Enterococcus, qPCR EPA 1611 CCE/100ml"
label var enteroQPCRcce_nd "Enterococcus, qPCR EPA 1611 below detection"
label var bactQPCRcce "Bacteroides, qPCR CCE/100ml"
label var bactQPCRcce_nd "Bacteroides, qPCR below detection"

sort stationid coldate coltime
tempfile qpcr
save `qpcr'


*------------------------------------------
* load the coliphage sample data
*------------------------------------------

insheet using "~/Documents/CRG/coliphage/13beaches-data/untouched/missionbay/wq/tblPhageData.csv", comma clear

* extract sample date
gen coldate = date(substr(sampledate,1,8),"MDY",2020)
	format coldate %d
	label var coldate "Sample collection date"

* code stationid to be consistent with the other data
gen stationid = ""
	replace stationid = "COMP1X" if samplesite=="Visitor's"
	replace stationid = "COMP2X" if samplesite=="Tecolote"
	replace stationid = "COMP3X" if samplesite=="L. Lagoon"
	replace stationid = "COMP4X" if samplesite=="DeAnza"
	replace stationid = "COMP5X" if samplesite=="Crown Point"
	replace stationid = "COMP6X" if samplesite=="Bonita Cove"
	assert stationid!=""
	label var stationid "stationid"

* ensure no duplicate obs
duplicates report stationid sampledate coliphag


*------------------------------------------
* Create indicator stub names and
* then reshape the data to wide format
*------------------------------------------	
gen ind = "fpc1601mpn" if coliphageenrichment == "male-specific"
	replace ind = "fmc1601mpn" if coliphageenrichment == "Somatic"
	assert ind!=""
	
keep stationid coldate ind result qualifier
order stationid coldate ind result qualifier
reshape wide result qualifier, i(stationid coldate) j(ind, string)


*------------------------------------------
* Flag non-detect values and set them to 0
* to be consistent with other datasets
*------------------------------------------	

local vars "fpc1601mpn fmc1601mpn"
foreach var of local vars  {
	rename result`var' `var'
	gen `var'_nd = 1 if inlist(qualifier`var',"<","ND")
		replace `var'_nd = 0 if (`var'!=.) & (`var'_nd==.)
		label values `var'_nd ND
	replace `var' = 0 if `var'_nd==1
		
	
}
drop qualifier*
order stationid coldate fpc1601mpn* fmc1601mpn*

*------------------------------------------
* Label variables and add notes
*------------------------------------------	

gen sampletype = "Daily composite sample"
	label var sampletype "Water sample type"

label var fpc1601mpn "F+ Coliphage, EPA 1601 MPN/100ml"
label var fpc1601mpn_nd "F+ Coliphage, EPA 1601 below detection"
label var fmc1601mpn "F- Coliphage, EPA 1601 MPN/100ml"
label var fmc1601mpn_nd "F- Coliphage, EPA 1601 below detection"

sort stationid coldate
tempfile phage
save `phage'

*------------------------------------------
* SAMPLE-LEVEL DATASET SYNTHESIS
*------------------------------------------

*------------------------------------------
* Merge the Coliphage data onto the
* culture data for daily composite samples
*------------------------------------------	
use `culture', clear
keep if sampletype =="Daily composite sample"
sort stationid coldate
merge 1:1 stationid coldate using `phage'
	assert _merge !=2
	drop _merge

tempfile comps
save `comps'


*------------------------------------------
* Merge the qPCR data onto the
* culture data for hourly samples
*------------------------------------------	
use `culture', clear
drop if sampletype =="Daily composite sample"
sort stationid coldate coltime
merge 1:1 stationid coldate coltime using `qpcr'
	assert _merge !=2
	drop _merge
	
* append the composite samples
append using `comps'


*------------------------------------------
* Final variable organization and labeling
*------------------------------------------	


* create a beach number for consistency with
* the epidemiology data
gen beach = "Mission Bay"
	label var beach "Beach"
gen beachcode = "Mission Bay " + substr(stationid,strlen(stationid)-1,1)
	label var beachcode "Mission Bay Beach Code"

	
* label depth
gen byte depth = 1
label var depth "Sample depth"
	capture label drop depth
label define depth 1 "Shin depth" 2 "Waist depth"
label values depth depth

order sampletype beachcode stationid coldate coltime depth enteroELTmpn - tcolifcfu_nd enteroQPCRcce - bactQPCRcce_nd 
sort sampletype beachcode stationid coldate coltime depth


* add dataset and variable notes
notes: All non-detect values are set to 0 and flagged with _nd indicators
notes enteroQPCRcce: only analyzed in 2 samples per day
notes bactQPCRcce: only analyzed in 2 samples per day
notes entero1600cfu: only analyzed in daily composite samples
notes ecoliCLTmpn: only analyzed in daily composite samples
notes tcolCLTmpn: only analyzed in daily composite samples
notes fpc1601mpn: only analyzed in daily composite samples
notes fmc1601mpn: only analyzed in daily composite samples

*------------------------------------------
* Save a sample-level dataset
*------------------------------------------

compress
label data "Mission Bay water sample data, formatted by 7-format-mb-wq.do"
save "~/Documents/CRG/coliphage/13beaches-data/final/mb-wq-samples.dta", replace
outsheet using "~/Documents/CRG/coliphage/13beaches-data/final/mb-wq-samples.csv", comma replace

desc
codebook, c


log close
exit

