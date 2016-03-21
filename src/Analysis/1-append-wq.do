

use "~/Documents/CRG/coliphage/13beaches-data/final/13beaches-wq-samples.dta", clear

preserve
use "~/Documents/CRG/coliphage/13beaches-data/final/13beaches-epi.dta", clear
keep beach beachcode intdate berm groundwater
ren intdate coldate
keep if beach=="Avalon" | beach=="Doheny" | beach=="Malibu" | beach=="Fairhope" | beach=="Goddard" | beach=="Mission Bay"

collapse (mean) berm groundwater , by(beach beachcode coldate)
tempfile epi
save `epi'
restore

* merge epi and wq datasets 
merge m:1 beach beachcode coldate using `epi'

* drop wq observations with no matching swimmers
drop if _m==1

* drop if there are no wq obs for a swimmer
drop if _m==2

drop _m

* drop sites with unusual data per 13 beaches protocol
drop if beachcode=="Avalon-D"
drop if beachcode=="Doheny-E"
drop if beachcode=="Doheny-C"
drop if beachcode=="Malibu-C"

* make risk variable
gen risk=.
replace risk = 1 if berm==1 | groundwater==1
replace risk = 1 if beach=="Fairhope" | beach=="Goddard" 
replace risk = 0 if berm== 0 | groundwater==0
replace risk = 0 if beach=="Mission Bay" 
replace risk = 0 if beach=="Malibu" 

label define riskl 1 "High" 0 "Low"
label values risk riskl

* rename variables
ren fmc1601mpn fmc1601
ren fmc1602mpn fmc1602
ren fpc1601mpn fpc1601
ren fpc1602mpn fpc1602

ren fmc1601mpn_nd fmc1601_nd
ren fmc1602mpn_nd fmc1602_nd
ren fpc1601mpn_nd fpc1601_nd
ren fpc1602mpn_nd fpc1602_nd


* drop 1602 values for Malibu
replace fmc1602 = . if beach=="Malibu"
replace fpc1602 = . if beach=="Malibu"
replace fmc1602_nd = . if beach=="Malibu"
replace fpc1602_nd = . if beach=="Malibu"

outsheet using "~/Documents/CRG/coliphage/13beaches-data/temp/beaches-coli-ent-wq.csv", comma replace 
