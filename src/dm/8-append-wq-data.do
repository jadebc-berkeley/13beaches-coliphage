
capture log close
set more off
clear all

log using "~/Documents/CRG/coliphage/13beaches-coliphage/src/dm/8-append-wq-data.log", text replace

*----------------------------------------
* 8-append-wq-data.do
* ben arnold
*
* impute non-detect values in datasets
* calculate daily average log10 values
* 
* append the water quality samples datasets 
* into a single file 
*
*
*----------------------------------------

*----------------------------------------
* input files:
*	neear-wq-samples.dta
*   adm-wq-samples.dta
*   mb-wq-samples.dta
*
* output files:
*	13beaches-wq-samples.dta / .csv
*
*----------------------------------------



*----------------------------------------
* read in the beaches sample indicator data
* standardize
* append 
*----------------------------------------

use "~/Documents/CRG/coliphage/13beaches-data/final/neear-wq-samples.dta", clear

append using "~/Documents/CRG/coliphage/13beaches-data/final/adm-wq-samples.dta"

append using  "~/Documents/CRG/coliphage/13beaches-data/final/mb-wq-samples.dta"


* restrict the data to sample labels and common WQ indicators
* Entero 1600
* Entero Enterolert (Mission Bay only)
* Entero qPCR 1611
* Coliphage (F+/- 1601/1602)

keep beach beachcode coldate coltime sampleid* stationid swimloc depth frozen sampletype entero1600cfu* enteroELTmpn* enteroQPCRcce enteroQPCRcce_nd enteroQPCRcce_qc fpc1601mpn* fpc1602mpn* fmc1601mpn* fmc1602mpn*
order beach beachcode coldate coltime sampleid* stationid swimloc depth frozen sampletype entero1600cfu* enteroELTmpn* enteroQPCRcce enteroQPCRcce_nd enteroQPCRcce_qc fpc1601mpn* fpc1602mpn* fmc1601mpn* fmc1602mpn*

notes drop *
notes drop _dta *
notes : 13-Beaches Individual Sample Water Quality Data, created by 8-append-wq-data.do
notes : Values below the detection limit for all indictors are set to 0
notes : Values are missing if a sample was not tested for a particular indicator
notes : All non-detect values are set to 0 and flagged with _nd indicators
notes : sampleid numbers are separate by type of indicator in the NEEAR studies and that is why there are multiple sampleid variables in this dataset
notes coltime: Mission Bay has different sample collection times than other beaches
notes sampletype: Only relevant for Mission Bay Samples
notes frozen: Only recorded for NEEAR coliphage samples
notes entero1600cfu: At Mission Bay, only analyzed in daily composite samples
notes enteroQPCRcce: At Mission Bay, only analyzed in 2 samples per day
notes enteroQPCRcce_qc: Only recorded for NEEAR Entero qPCR samples.  For NEEAR Entero qPCR samples that failed QC criteron, the value was imputed based on the average of the other two samples at the same depth and time


compress
label data "13 beaches water quality data, created by 8-append-wq-data.do"
saveold "~/Documents/CRG/coliphage/13beaches-data/final/13beaches-wq-samples.dta", replace version(12)
outsheet using "~/Documents/CRG/coliphage/13beaches-data/final/13beaches-wq-samples.csv", comma replace

desc

* write a codebook to separate file
log close
log using "~/Documents/CRG/coliphage/13beaches-data/final/13beaches-wq-samples-codebook.txt", text replace
desc, s
notes
codebook
log close




