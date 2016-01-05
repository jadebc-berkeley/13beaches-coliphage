capture log close
set more off
clear all

log using "~/dropbox/13beaches/src/dm/5-format-neear-wq-jbc.log", text replace

*----------------------------------------
* 5-format-neear-wq.do
* ben arnold
*
* format the NEEAR water quality data from
* Tim Wade.  label variables
*
* version 3 (15 jul 2015) by Jade
* modify script to output a dataset with 
* a row for each sample
*
* version 2 (4 feb 2015)
* added coliphage data provide by Tim Wade 3 Feb 2015
*
* version 1 (25 feb 2014)
*
*----------------------------------------

*----------------------------------------
* input files:
*	epcrmeans.txt
*   cfumeans.txt
*
* output files:
*	neear-wq.dta
*
*----------------------------------------

*********************************************
* create dataset with one row for each sample
*********************************************

/*----------------------------------------
* read in the enterococcus dataset
*----------------------------------------
insheet using "~/dropbox/13beaches/data/untouched/neear/cfudata.txt", clear

* reformat collection date
gen coldate = date(collectiondate,"YMD")
	format coldate %d
	label var coldate "Sample Collection Date"
	
sort beach coldate

rename count_cfu entero1600
label var entero1600 "log10 Entero EPA 1600"

tempfile cfu_sample
save `cfu_sample'



*/
*----------------------------------------
* read in the Coliphage dataset
* merge to the existing dataset
*----------------------------------------

use "~/dropbox/13beaches/data/untouched/neear/wqphagemarine.dta", clear

*----------------------------------------
* make beach labels more informative
* add labels to a few other variables
*----------------------------------------
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

gen marine=1
	label var marine "Marine Beach"

ren collection_date coldate
ren sample_id sampleid

label var beach "Study beach"
label var marine "Marine beach"
label var coldate "Sample Collection Date"

order beach coldate
sort beach coldate

order beach coldate marine

* replace non detects with 0.1
local vlist "rna_clat5 group1_clat5 group2_clat5 group3_clat5 group4_clat5 dna_clat5 spot24"
foreach var of local vlist {
	replace mpn_`var' = 0.1 if mpn_`var'==0
}

ren mpn_spot24 fpc1601
	
* convert from wide to long format
foreach var of varlist mpn_rna_clat5 mpn_group1_clat5-fpc1601 {
	ren `var' result`var'
}

reshape long result, i(sampleid) j(groupindex) string

* generate variables to match other datasets
gen indicator=""
replace indicator="Fminus coliphage" if groupindex=="fmc1601"
replace indicator="Fplus coliphage" if groupindex=="fpc1601"
#delimit;
replace indicator="Fplus RNA coliphage" if groupindex=="mpn_group1_clat5"|
groupindex=="mpn_group2_clat5"|groupindex=="mpn_group3_clat5"|groupindex=="mpn_group4_clat5"|
groupindex=="mpn_rna_clat5";
replace indicator="Fplus DNA coliphage" if groupindex=="mpn_dna_clat5";
#delimit cr

gen analysismethod=""
replace analysismethod="EPA 1601" if groupindex=="fmc1601" | groupindex=="fpc1601"
#delimit;
replace analysismethod="CLAT5" if groupindex=="mpn_dna_clat5" | groupindex=="mpn_group1_clat5"|
groupindex=="mpn_group2_clat5"|groupindex=="mpn_group3_clat5"|groupindex=="mpn_group4_clat5"|
groupindex=="mpn_rna_clat5";
#delimit cr

gen qualifier=""
replace qualifier = "<" if result==0.1

gen log10 = log10(result)

drop groupindex

compress
label data "NEEAR water quality sample data, created by 5-format-neear-wq.do"
save "~/dropbox/13beaches/data/final/neear-wq-samples.dta", replace

*********************************************
* create dataset with daily averages
*********************************************
*----------------------------------------
* read in the CFU dataset
* rename variables
*----------------------------------------

insheet using "~/dropbox/13beaches/data/untouched/neear/cfumeans.txt", clear
sort beach collection_date

rename logcountdy avgdyentero1600
rename logcount8  avgamentero1600
rename logcount11 avgmdentero1600
rename logcount15 avgpmentero1600
rename logcount1  avgshentero1600
rename logcount2  avgwtentero1600

label var avgdyentero1600 "Daily Avg log10 Entero EPA 1600"
label var avgamentero1600 "AM Avg log10 Entero EPA 1600"
label var avgmdentero1600 "Mid-Day Avg log10 Entero EPA 1600"
label var avgpmentero1600 "PM Avg log10 Entero EPA 1600"
label var avgshentero1600 "shin Avg log10 Entero EPA 1600"
label var avgwtentero1600 "waist Avg log10 Entero EPA 1600"

tempfile cfu
save `cfu'

*----------------------------------------
* read in the QPCR dataset
* merge to the CFU dataset
*----------------------------------------
insheet using "~/dropbox/13beaches/data/untouched/neear/epcrmeans.txt", clear

rename logcountdy avgdyenteropcr
rename logcount8  avgamenteropcr
rename logcount11 avgmdenteropcr
rename logcount15 avgpmenteropcr
rename logcount1  avgshenteropcr
rename logcount2  avgwtenteropcr

label var avgdyenteropcr "Daily Avg log10 Entero qPCR"
label var avgamenteropcr "AM Avg log10 Entero qPCR"
label var avgmdenteropcr "Mid-Day Avg log10 Entero qPCR"
label var avgpmenteropcr "PM Avg log10 Entero qPCR"
label var avgshenteropcr "shin Avg log10 Entero qPCR"
label var avgwtenteropcr "waist Avg log10 Entero qPCR"

sort beach collection_date
merge 1:1 beach collection_date using `cfu'
list if _merge!=3
drop _merge


* reformat collection date
gen coldate = date(collection_date,"YMD")
	format coldate %d
	label var coldate "Sample Collection Date"

sort beach coldate
save `cfu', replace

*----------------------------------------
* calculate daily averages across all times and locations within each beach
* calculate am, mid-day, pm averages
* calculate shin and waist depth averages
*----------------------------------------
use "~/dropbox/13beaches/data/untouched/neear/wqphagemarine.dta", clear

local vlist "rna_clat5 group1_clat5 group2_clat5 group3_clat5 group4_clat5 dna_clat5 spot24"
foreach var of local vlist {

	* replace non detects with 0.1
	replace mpn_`var' = 0.1 if mpn_`var'==0

	gen log10`var' = log10(mpn_`var')
	
	sort beach collection_date
	by beach collection_date: egen avgdyphage_`var' = mean(log10`var')
	
	by beach collection_date: egen avgamphage_`var' = mean(log10`var') if collection_time==8
	by beach collection_date: egen avgmdphage_`var' = mean(log10`var') if collection_time==11
	by beach collection_date: egen avgpmphage_`var' = mean(log10`var') if collection_time==15
	
	by beach collection_date: egen avgshphage_`var' = mean(log10`var') if depth==1
	by beach collection_date: egen avgwtphage_`var' = mean(log10`var') if depth==2

}

* collapse the data by beach / date
collapse (max) avg*, by(beach collection_date)
rename collection_date coldate

* label variables
local vlist "avgdy avgam avgmd avgpm avgsh avgwt"
foreach v of local vlist {
	rename `v'phage_spot24 `v'fpc1601
}
label var avgdyfpc1601 "Daily Avg log10 F-plus Coliphage EPA 1601"
label var avgamfpc1601 "AM Avg log10 F-plus Coliphage EPA 1601"
label var avgmdfpc1601 "Mid-Day Avg log10 F-plus Coliphage EPA 1601"
label var avgpmfpc1601 "PM Avg log10 F-plus Coliphage EPA 1601"
label var avgshfpc1601 "shin Avg log10 F-plus Coliphage EPA 1601"
label var avgwtfpc1601 "waist Avg log10 F-plus Coliphage EPA 1601"

forvalues i = 1/4 {
	label var avgdyphage_group`i'_clat5 "Daily Avg log10 Group `i' RNA Coliphage, 5-hour Clat assay"
	label var avgamphage_group`i'_clat5 "AM Avg log10 Group `i' RNA Coliphage, 5-hour Clat assay"
	label var avgmdphage_group`i'_clat5 "Mid-Day Avg log10 Group `i' RNA Coliphage, 5-hour Clat assay"
	label var avgpmphage_group`i'_clat5 "PM Avg log10 Group `i' RNA Coliphage, 5-hour Clat assay"
	label var avgshphage_group`i'_clat5 "shin Avg log10 Group `i' RNA Coliphage, 5-hour Clat assay"
	label var avgwtphage_group`i'_clat5 "waist Avg log10 Group `i' RNA Coliphage, 5-hour Clat assay"
}

label var avgdyphage_rna_clat5 "Daily Avg log10 RNA Coliphage, 5-hour Clat assay"
label var avgamphage_rna_clat5 "AM Avg log10 RNA Coliphage, 5-hour Clat assay"
label var avgmdphage_rna_clat5 "Mid-Day Avg log10 RNA Coliphage, 5-hour Clat assay"
label var avgpmphage_rna_clat5 "PM Avg log10 RNA Coliphage, 5-hour Clat assay"
label var avgshphage_rna_clat5 "shin Avg log10 RNA Coliphage, 5-hour Clat assay"
label var avgwtphage_rna_clat5 "waist Avg log10 RNA Coliphage, 5-hour Clat assay"

label var avgdyphage_dna_clat5 "Daily Avg log10 DNA Coliphage, 5-hour Clat assay"
label var avgamphage_dna_clat5 "AM Avg log10 DNA Coliphage, 5-hour Clat assay"
label var avgmdphage_dna_clat5 "Mid-Day Avg log10 DNA Coliphage, 5-hour Clat assay"
label var avgpmphage_dna_clat5 "PM Avg log10 DNA Coliphage, 5-hour Clat assay"
label var avgshphage_dna_clat5 "shin Avg log10 DNA Coliphage, 5-hour Clat assay"
label var avgwtphage_dna_clat5 "waist Avg log10 DNA Coliphage, 5-hour Clat assay"

order beach coldate *1601

sort beach coldate
tempfile phage
save `phage'

use `cfu', clear
merge 1:1 beach coldate using `phage'
assert _merge != 2
tab beach _merge
drop _merge

*----------------------------------------
* make beach labels more informative
* add labels to a few other variables
*----------------------------------------
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

label var beach "Study beach"
label var marine "Marine beach"
label var collection_date "Sample Collection Date"


order beach collection_date coldate
sort beach coldate

drop collection_date
order beach coldate marine

compress
label data "NEEAR water quality data, created by 5-format-neear-wq.do"
save "~/dropbox/13beaches/data/final/neear-wq.dta", replace

desc

log close
exit





