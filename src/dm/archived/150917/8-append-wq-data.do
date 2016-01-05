capture log close
set more off
clear all

log using "~/dropbox/13beaches/src/dm/8-append-wq-data.log", text replace

*----------------------------------------
* 8-append-wq-data.do
* ben arnold
*
* append the water quality datasets into
* a single file and calculate quantiles
* of indicator distributions
*
* version 7 (19 aug 2015)
* updated the quartile calculations to include
* Enterolert data from Mission Bay
*
* version 6 (12 jun 2015)
* output a codebook for the 13beaches-wq.dta/.csv datasets
*
* version 5 (24 apr 2015)
* revised the indicator quartile calculations
* to exclude the samples from lagoons in Doheny and Malibu (sites C)
* since very, very few individuals were exposed at those sites
* and the measurements are not representative of the distributions
* experienced in the population
*
* version 4 (20 feb 2015)
* added indicator quantile calculations
* 
* version 3 (12 feb 2015)
* added .csv export
*
* version 2 (4 feb 2015)
* added Mission Bay data
*
* version 1 (28 feb 2014)
*
*----------------------------------------

*----------------------------------------
* input files:
*	neear-wq.dta
*   adm-wq.dta
*   mb-wq.dta
*
* output files:
*	13beaches-wq.dta / .csv
*
*----------------------------------------


*----------------------------------------
* read in the beaches indicator data
* append 
*----------------------------------------

use "~/dropbox/13beaches/data/final/adm-wq.dta", clear

append using "~/dropbox/13beaches/data/final/mb-wq.dta"

order beach beachcode marine coldate

append using  "~/dropbox/13beaches/data/final/neear-wq.dta"

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
note: TS 13-Beaches Water Quality Data, created by 8-append-wq-data.do
note: All indicators are means of log10 values, summarized daily over all samples, or separately for AM, Mid-Day, PM, shin-depth, or waist-depth samples 
note: The enteropcr2-4 variables are alternative Enterococcus qPCR calculations at Avalon/Doheny/Malibu beaches only
note: The clat5 series of coliphage indicators are unique to 2 NEEAR beaches (Fairhope Goddard) 
note: The Mission Bay water quality data are separate for 6 beaches within the area and need to merge to individuals separately by beachcode
note: Avalon water quality data are summarized separately for sites A/B/C and D; site D is very different. See Yau et al 2014 Figs 1 and 2.
note: Doheny and Malibu water quality data are summarized separately for sites C; at both beaches site C was located in a lagoon. At Doheny, site E is also summarized separately, since it was nearly 1 mile from the other sampling sites (See Fig 1 of Arnold et al. 2013 and SI Fig 1 of Colford et al. 2012)
note: Quartile categories for Entero 1600 and qPCR exclude water quality measurements from Doheny and Malibu sites C in the lagoons because they are not representative of the broader swimmer exposure. The Entero 1600 quartiles were calculated using Enterolert data from Mission Bay
note: AM/Mid-Day/PM and shin/waist disaggregations were not available for Mission Bay due to different sampling methods at that beach


compress
label data "13 beaches water quality data, created by 8-append-wq-data.do"
save "~/dropbox/13beaches/data/final/13beaches-wq.dta", replace
outsheet using "~/dropbox/13beaches/data/final/13beaches-wq.csv", comma replace

desc


* write a codebook to separate file
log close
log using "~/dropbox/13beaches/data/final/13beaches-wq-codebook.txt", text replace
desc, s
*aorder
codebook
log close

exit





