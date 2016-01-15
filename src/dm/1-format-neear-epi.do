capture log close
set more off
clear all

log using "~/Documents/CRG/coliphage/13beaches-coliphage/src/dm/1-format-neear-epi.log", text replace

*----------------------------------------
* 1-format-neear-epi.do
* ben arnold
*
* format the NEEAR epi data from
* Tim Wade.
*
* subset the data to the most 
* relevant variables.
*
* label variables
*----------------------------------------

*----------------------------------------
* input files:
*	healthfinalUCB.dta
*
* output files:
*	neear-epi.dta
*
*----------------------------------------


*----------------------------------------
* read in the dataset
*----------------------------------------
use "~/Documents/CRG/coliphage/13beaches-data/untouched/neear/healthfinalUCB.dta", clear

*----------------------------------------
* rename and if necessary label 
* variables of interest
*----------------------------------------

label var indid "Individual ID"
rename hh hhid
label var hhid "Household ID"

label var beach "Beach"
label var race "Racial category"


* rename some variables to be consistent with ADM

label define sex 1 "Male" 2 "Female"
label values sex sex

gen watertime = water
	label var watertime "Time in water (mins)"

local symps "stomach diarrhea nausea vomiting urinarytractinfection fever headache sorethroat cough cold runnynose earache wateryeyes eyeinfection cut rash"
foreach symp of local symps {
	rename `symp'startdate `symp'stdt
	label var `symp'stdt "Start date for `symp'"
}

* rename some of the impact variables from "uri" to "uti"
local vlist "stayhome stayhomedays stopdaily stopdailydays othersmiss othersmissdays phonedoc visitdoc visitdoctimes diagnosis emergencyroom visitertimes hospitalized hospitaldays anyprescrdrug ownmoneyfordrugs anyotcmeds ownmoneyforotc"
foreach var of local vlist {
	rename `var'_uri `var'_uti
}


label var swimloc "Swim location"

rename awash algaewash
rename amouth algaemouth
rename aller allergy
rename nose noseplugs
rename plug earplugs
rename swetdry sanddry
rename smouth sandmouth
rename swash sandwash

rename dsun sunmin

rename eyeinfect eyeinfection
rename stomachache stomach
rename vomit vomiting
rename uti urinarytractinfection
rename cuts cut
drop earachelist

* make rawfood variable to be consistent w/ ADM
gen byte rawfood = (eggs_any==1)|(fish_any==1)|(rawmeat_any==1)
	label var rawfood "Ate undercooked eggs/meat/fish"


* label sex variable
label values sex sex

label var venfest "Silver Beach festival day with 1000s of food vendors"

label define beachtype 1 "Lake" 2 "River" 3 "Ocean" 6 "Other" 7 "Refused" 8 "Don't know" 9 "N/A"
label values beachtype beachtype


*----------------------------------------
* fix a few coding inconsistencies
* identified during variable reconciliation
* after appending all of the datasets together
*----------------------------------------
* affects huntington + west beaches
recode racewhite raceblack raceasian raceindian racehaw hisp (2=0)

* recode blockface
recode blockface (2=0)

* affects edgewater, fairhope, goddard beaches
recode pool (2=0) (8=0) (9=0)

* affects Huntington and west beachs
recode coldallergy (2=0) (9=0)
recode othersmiss_* (2=0)


* label define yesno 0 "No" 1 "Yes"
local vlist "sandmouth seatdr sandwash algae algaemouth blockface wave *still *allergy stayhome_* stopdaily_* phonedoc_* visitdoc_* emergencyroom_* hospitalized_* anyprescrdrug_* anyotcmeds_* prot sunburn sunburn1-sunburn7 blockface"
foreach var of varlist `vlist' {
	label values `var' yesno
}

* add-back some code labels that are in the NEEAR codebook 
* but for some reason didn't come through with the dataset


label var tanning "What happens when in sun repeatedly without sunscreen during the summer?"
label define tanning 1 "Dark Tan" 2 "Some Tanning" 3 "No Tan, Some Freckles" 4 "Repeated Sunburns" 5 "Other (specify)" 6 "Never Go Out in the Sun"
label values tanning tanning

* label define beachtype 1 "Lake" 2 "River" 3 "Ocean" 6 "Other"
label values beachtype beachtype

* standardize race categories to be consistent with the ADM data
gen byte race2 = .
replace race2 = 1 if race==1 & hisp!=1
replace race2 = 2 if race==1 & hisp==1
replace race2 = 3 if race!=1 & hisp==1
replace race2 = 4 if race==2 & hisp!=1
replace race2 = 5 if race==3 & hisp!=1
replace race2 = 6 if race==4 & hisp!=1
replace race2 = 7 if race==6 & hisp!=1
replace race2 = 8 if race==7 & hisp!=1
replace race2 = 9 if race==.
drop race
label drop race
label define race 1 "white" 2 "white, hispanic" 3 "non-white, hispanic" 4 "black" 5 "asian" 6 "american indian" 7 "multiple races" 8 "other" 9 "missing"
rename race2 race
label values race race
label var race "Racial category"



*----------------------------------------
* restrict to variables of interest
*----------------------------------------


* variables to check on:
* did NEEAR collect info on income?

# delimit ;
drop racecat agecat1 milescat body3 swall3 wave3 
e1 e2 hcgi hcresp eye swim1face chronany gag3
comecat flag head3 mouth3 kindan* 
racewhite raceblack raceasian raceindian racehaw

/* variables dropped after cross-check with ADM -- not collected in ADM */
aeatdr airmat
animunk*
anycuts
arr_time
bodyph
bsurf canoe surf wrun wski wsurf wtube jski waveb
ksurf pboat raft sail scuba snorkel
*startday
cond
dmach drink drinks dvend
enrolled
faceph gagph gagwater getfacewet
missingperson
mouthph
privatepool
shell specifybeachtype
swallph swim1base swim2base swimface
totalph
total water timeinwaterin
wadeph wadingpool waterpark waterph
hadcontact
mach vend
otherswimlocation anyotherswimming
tanos
statusc
transect
;
# delimit cr


*---------------------------------------------
* output variable information for harmonization
*---------------------------------------------
preserve
desc, replace clear
label data "variables for neear-epi.dta"
save "~/Documents/CRG/coliphage/13beaches-data/temp/neear-epi-vars.dta", replace
restore


*---------------------------------------------
* Save a dataset to combine with other beaches
*---------------------------------------------
order beach indid hhid intdate teledate

label data "NEEAR epi data, created by 1-format-neear-epi.do"
save "~/Documents/CRG/coliphage/13beaches-data/final/neear-epi.dta", replace

log close



* write a codebook
log using "~/Documents/CRG/coliphage/13beaches-data/final/neear-epi-codebook.txt", text replace
*aorder
codebook
log close
exit






