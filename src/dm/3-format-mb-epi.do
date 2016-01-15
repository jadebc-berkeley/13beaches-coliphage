capture log close
set more off
clear all

log using "~/13beaches-coliphage/src/dm/3-format-mb-epi.log", text replace

*----------------------------------------
* 3-format-mb-epi.do
* ben arnold
*
* format the Mission Bay 
* dataset so that it can be 
* appended to the other beaches datasets 
* (NEEAR, Avalon, Doheny, Malibu). 
*
*
* use the NEEAR data as the name template
* 
* version 3 (21 feb 2015)
* retained beachcode to merge to water quality data
*
* version 2 (13 nov 2014)
* final data cleaning
*
* version 1 (17 july 2014)
*
*----------------------------------------

*----------------------------------------
* input files:
*	mb_analysis_final.dta
*
* output files:
*	mb-epi.dta
*
*----------------------------------------


*----------------------------------------
* read in the Mission Bay Analysis data
*----------------------------------------


use "~/13beaches-data/untouched/missionbay/mb_analysis_final.dta", clear


* drop individuals with no follow-up data
keep if catiyn==1


*----------------------------------------
* rename variables and subset to be
* consistent with the NEEAR data
*----------------------------------------




* make unique household and individual IDs
capture drop beach
gen beach = "MissionBay"
	label var beach "Beach"
	
rename hhid hhid2
gen hhid = string(hhid2)
	label var hhid "Household ID"

gen indid = string(fullid,"%9.0f")
	label var indid "Individual ID"

order beach hhid indid

rename sdate intdate
	label var intdate "Beach interview date"

gen _m = floor(catidate/1000000)
gen _d = floor(catidate/10000)- _m*100
gen _y = catidate-(_m*1000000)-(_d*10000)
gen teledate = mdy(_m,_d,_y)
	format teledate %d
	label var teledate "Phone interview date"
	drop _m _d _y


label var age "age in years"

* Not Measured: algaemouth
* Not Measured: algaewash

rename allergies allergy
	label var allergy "Has allergies"

rename saq16a anim_base
rename oact_a anim_int
	label var anim_int "Person: contact with animals"
gen anim_any = (anim_base==1)|(anim_int==1)
	label var anim_any "Any animal contact"
	
rename water anycontact
	label var anycontact "Any contact with water"

* Note: for Mission bay, cannot separate out medication use by OTC vs. prescription
* any medication use (either OTC or prescription) is stored in the common OTC variables
rename gh8_h anyotcmeds_ear
	recode anyotcmeds_ear (5=0)
	label var anyotcmeds_ear "Used OTC medications for syndrome (Ear_infec)"
rename gh7_h anyotcmeds_eye
	recode anyotcmeds_eye (5=0)
	label var anyotcmeds_eye "Used OTC medications for syndrome (Eye_infec)"

gen anyotcmeds_gas = (gh12_h==1) | (gh13_h==1) | (gh14_h==1) | (gh15_h==1)
	replace anyotcmeds_gas = . if (gh12_h==.) & (gh13_h==.) & (gh14_h==.) & (gh15_h==.)
	label var anyotcmeds_gas "Used OTC medications for syndrome (GI_illness)"
	
gen anyotcmeds_res = (gh16_h==1) | (gh17_h==1) | (gh18_h==1) | (gh19_h==1)
	replace anyotcmeds_res = . if (gh16_h==.) & (gh17_h==.) & (gh18_h==.) & (gh19_h==.)
	label var anyotcmeds_res "Used OTC medications for syndrome (Respiratory)"
rename gh10_h anyotcmeds_skn
	recode anyotcmeds_skn (5=0)
	label var anyotcmeds_skn "Used OTC medications for syndrome (Skin_Rash)"
	
* Not Measured: aanyotcmeds_uti
* Not Measured: aanyprescrdrug_ear
* Not Measured: aanyprescrdrug_eye
* Not Measured: aanyprescrdrug_gas
* Not Measured: aanyprescrdrug_res
* Not Measured: aanyprescrdrug_skn
* Not Measured: aanyprescrdrug_uti

rename chronresp asth
	label var asth "chronic respiratory problems"

* Not Measured: berm

drop block
rename saq8b block
	label var block "Wear sunscreen-block today (R)"

rename water_shoulders bodycontact
	label var bodycontact "Body contact with water"

rename buried bsand
	label var bsand "Had body buried in sand (R)"

* Not Measured: cold
* Not Measured: coldallergy
* Not Measured: colddays
* Not Measured: coldstdt
* Not Measured: coldstill

* Not Measured: come

* cough is measured 2 ways in MB. combine into one, using the min(start date)
rename cough cough2
gen cough = (cough2==1) | (cough_phlegm==1)
	replace cough = . if (cough2==.) & (cough_phlegm==.)
	label var cough "Person had cough"

	* clean 3 values with impossible months
	replace gh16_b1 = 5 if gh16_b1<5 & gh16_b1!=.
gen coughdt1 = mdy(gh16_b1,gh16_b2,2003)
gen coughdt2 = mdy(gh17_b1,gh17_b2,2003)
gen coughstdt = min(coughdt1,coughdt2)
	drop coughdt1 coughdt2
	format coughstdt %d
	label var coughstdt "cough start date"

* Not Measured: coughallergy
* Not Measured: coughdays
* Not Measured: coughstill

rename scrapes cut
rename saq12 cutbase
	label var cutbase "Cut self or have open cut (R)"
gen cutstdt = mdy(gh11_b1,gh11_b2,2003)
	format cutstdt %d
	label var cutstdt "cut start date"
* Not Measured: cutdays
* Not Measured: cutstill


gen diarrheastdt = mdy(gh14_b1,gh14_b2,2003)
	format diarrheastdt %d
	label var diarrheastdt "diarrhea start date"
* Not Measured: diarrheadays
* Not Measured: diarrheanumber
* Not Measured: diarrheastill

* clean 1 apparent data entry error in the diarrhea start date
replace diarrheastdt = mdy(6,1,2003) if indid=="40004601"

gen differentbeach = (sw1c_c>0 & sw1c_c!=. & sw1c_c!=98) | (sw1d_c>0 & sw1d_c!=. & sw1d_c!=98)
	label var differentbeach "Summary variable: swam at other beach"
	
rename saq11 drinksy
	recode drinksy (2=0)

* note: sunhr and sunmin not calculable for Mission Bay -- only collected as categorical data
* create sunhr_cat rather than sunhr
rename saq9b sunhr_cat
	label var sunhr_cat "Total hours in direct sunlight (category)"
	
* Not Measured: sunmin	


* earache and ear discharge need to be combined.  combine into one, using the min(start date)
rename earache earache2
gen earache = (earache2==1) | (ear_discharge==1)
	replace earache = . if (earache2==.) & (ear_discharge==.)
	label var earache "Person had earache/infection/runny ears since beach interview"

gen eardt1 = mdy(gh8_b1,gh9_b2,2003)
gen eardt2 = mdy(gh9_b1,gh9_b2,2003)
gen earachestdt = min(eardt1,eardt2)
	drop eardt1 eardt2
	format earachestdt %d
	label var earachestdt "earache start date"
	
* Not Measured: earacheallergy
* Not Measured: earachedays
* Not Measured: earachestill

rename saq6a earplugs

rename saq10 eatfood


rename saq15b eggs_base
rename food_b eggs_int
	label var eggs_int "Person:eaten runny/raw eggs since interview"
drop eggs_any
gen eggs_any = (eggs_base==1)|(eggs_int==1)
	label var eggs_any "Any undercooked eggs"

* Not Measured: emergencyroom_ear
* Not Measured: emergencyroom_eye
* Not Measured: emergencyroom_gas
* Not Measured: emergencyroom_res
* Not Measured: emergencyroom_skn
* Not Measured: emergencyroom_uti

rename eye_infection_past eyebase
	label var eyebase "Eye infection at baseline"

rename eye_irritation eyeinfection
	* recode 1 implausible value
	replace gh7_b2 = 26 if gh7_b2==1 & gh7_b1==1
	replace gh7_b1 = 5 if gh7_b1==1
gen eyeinfectionstdt = mdy(gh7_b1,gh7_b2,2003)
	format eyeinfectionstdt %d
	label var eyeinfectionstdt "eye infection start date"

* Not Measured: eyeinfectiondays
* Not Measured: eyeinfectionstill

gen feverstdt = mdy(gh5_b1,gh5_b2,2003)
	format feverstdt %d
	label var feverstdt "fever start date"
	
* Not Measured: feverdays
* Not Measured: feverstill
* Not Measured: fevertemp
* Not Measured: fevertemptaken

rename food_e fish_int
	recode fish_int (5=0)
rename saq15g fish_base
gen fish_any = (fish_base==1)|(fish_int==1)
	label var fish_any "Any raw fish"


rename saq10a food

gen gibase = (diarrhea_past==1) | (vomiting_past==1)
	label var gibase "GI illness at baseline"
	* flag 1 individual who reported GI symptoms before the interview date
	replace gibase=1 if indid=="20002502"


rename chrongi gichron
	label var gichron "Chronic GI problems/Crohn's/IBS"
	* assume that the large number of individuals with missing values (3703) do not have chronic GI problems
	replace gichron = 0 if gichron==.

drop gicontact
rename saq16b gicontact_base

* Not Measured: groundwater

* Not Measured: headache
* Not Measured: headachedays
* Not Measured: headachestdt
* Not Measured: headachestill

rename water_face_under headunder
	label var headunder "Head under water"


* household income at MB uses different categories
* create hhinc_mb
rename isum hhinc_mb
	label var hhinc_mb "Household Pre-tax income, 2002 (Mission Bay categories)"

* Not Measured: hospitaldays_ear
* Not Measured: hospitaldays_eye
* Not Measured: hospitaldays_gas
* Not Measured: hospitaldays_res
* Not Measured: hospitaldays_skn
* Not Measured: hospitaldays_uti

* Not Measured: hospitalized_ear
* Not Measured: hospitalized_eye
* Not Measured: hospitalized_gas
* Not Measured: hospitalized_res
* Not Measured: hospitalized_skn
* Not Measured: hospitalized_uti

rename saq6d mask

rename water_mouth mouthwater

gen nauseastdt = mdy(gh12_b1,gh12_b2,2003)
	format nauseastdt %d
	label var nauseastdt "nausea start date"
	
* Not Measured: nauseadays
* Not Measured: nauseastill

rename saq6b noseplugs

* Not Measured: othersmiss_ear
* Not Measured: othersmiss_eye
* Not Measured: othersmiss_gas
* Not Measured: othersmiss_res
* Not Measured: othersmiss_skn
* Not Measured: othersmiss_uti
* Not Measured: othersmissdays_ear
* Not Measured: othersmissdays_eye
* Not Measured: othersmissdays_gas
* Not Measured: othersmissdays_res
* Not Measured: othersmissdays_skn
* Not Measured: othersmissdays_uti	

* Not Measured: ownmoneyfordrugs_ear
* Not Measured: ownmoneyfordrugs_eye
* Not Measured: ownmoneyfordrugs_gas
* Not Measured: ownmoneyfordrugs_res
* Not Measured: ownmoneyfordrugs_skn
* Not Measured: ownmoneyfordrugs_uti
* Not Measured: ownmoneyforotc_ear
* Not Measured: ownmoneyforotc_eye
* Not Measured: ownmoneyforotc_gas
* Not Measured: ownmoneyforotc_res
* Not Measured: ownmoneyforotc_skn
* Not Measured: ownmoneyforotc_uti

* Not Measured: phonedoc_ear
* Not Measured: phonedoc_eye
* Not Measured: phonedoc_gas
* Not Measured: phonedoc_res
* Not Measured: phonedoc_skn
* Not Measured: phonedoc_uti

* Not Measured: pool

rename gh4_a preg
recode preg (5=0)

rename saq8c prothat
rename saq9a prot

rename sw2a_a publicpool
recode publicpool (5=0)


* create a standardize race category consistent with other beaches
gen byte race2 = . 
replace race2 = 1 if race==1 & inlist(hisp,5,.)
replace race2 = 2 if race==1 & inlist(hisp,1,2,3,4)
replace race2 = 3 if race!=1 & (race==7 | inlist(hisp,1,2,3,4) )
replace race2 = 4 if race==2 & inlist(hisp,5,.)
replace race2 = 5 if race==4 & inlist(hisp,5,.)
replace race2 = 6 if race==3 & inlist(hisp,5,.)
replace race2 = 7 if race==6 & inlist(hisp,5,.)
replace race2 = 8 if race==5 & inlist(hisp,5,.)
replace race2 = 9 if race==.
drop race
label drop race
label define race 1 "white" 2 "white, hispanic" 3 "non-white, hispanic" 4 "black" 5 "asian" 6 "american indian" 7 "multiple races" 8 "other" 9 "missing"
rename race2 race
label values race race
label var race "Racial category"

* standardize hisp var
gen byte hisp2 = inlist(hisp,1,2,3,4) | (race==7 & hisp==5)
drop hisp
rename hisp2 hisp
label var hisp "Hispanic"

	
rename skin_rash rash
rename rash_past rashbase
gen rashstdt = mdy(gh10_b1,gh10_b2,2003)
	format rashstdt %d
	label var rashstdt "skin rash start date"

* Not Measured: rashdays
* Not Measured: rashstill

rename food_d rawmeat_int
	recode rawmeat_int (5=0)
rename saq15d rawmeat_base
gen rawmeat_any = (rawmeat_int==1) | (rawmeat_base==1)
	label var rawmeat_any "Any raw meat"

gen byte rawfood = (eggs_any==1)|(fish_any==1)|(rawmeat_any==1)
	label var rawfood "Ate undercooked eggs/meat/fish"

rename saq8a repel

rename nasal_congestion runnynose

gen runnynosestdt = mdy(gh18_b1,gh18_b2,2003)
	format runnynosestdt %d
	label var runnynosestdt "congestion/runny nose start date"

* Not Measured: runnynoseallergy
* Not Measured: runnynosedays
* Not Measured: runnynosestill

gen samebeach = (sw1a_a==1) & (sw1a_b==1)
	label var samebeach "Swim/wade: same beach"
	
* Not Measured: sanddry
* Not Measured: sandmouth
* Not Measured: sandwash

gen shade = prot
	label var shade "Use protective equipment (R)"

rename gh1_b skinchron
recode skinchron (5=0)

rename sore_throat_past sorebase

rename sore_throat sorethroat
gen sorethroatstdt = mdy(gh19_b1,gh19_b2,2003)
	format sorethroatstdt %d
	label var sorethroatstdt "sore throat start date"

* Not Measured: sorethroatallergy
* Not Measured: sorethroatdays
* Not Measured: sorethroatstill

gen byte stayhome_ear = (gh8_d==1) | (gh9_d==1)
	replace stayhome_ear = . if (gh8_d==.) & (gh9_d==.)
	label var stayhome_ear "Missed work/school for syndrome (Ear_infec)"

gen byte stayhome_eye = (gh7_d==1)
	replace stayhome_eye = . if (gh7_d==.)
	label var stayhome_eye "Missed work/school for syndrome (Eye_infec)"
	
gen byte stayhome_gas = (gh12_d==1) | (gh13_d==1) | (gh14_d==1) | (gh15_d==1)
	replace stayhome_gas = . if (gh12_d==.) & (gh13_d==.) & (gh14_d==.) & (gh15_d==.)
	label var stayhome_gas "Missed work/school for syndrome (GI_illness)"
gen byte stayhome_res = (gh16_d==1) | (gh17_d==1) | (gh18_d==1) | (gh19_d==1)
	replace stayhome_res = . if (gh16_d==.) & (gh17_d==.) & (gh18_d==.) & (gh19_d==.)
gen byte stayhome_skn = (gh10_d==1)
	replace stayhome_skn = . if (gh10_d==.)
	label var stayhome_skn "Missed work/school for syndrome (Skin_rash)"

* Not Measured: stayhome_uti

gen stayhomedays_ear = max(gh8_e,gh9_e)
	label var stayhomedays_ear "Number of days missed from work/school (Ear_infec)"
gen stayhomedays_eye = gh7_e
	label var stayhomedays_eye "Number of days missed from work/school (Eye_infec)"
gen stayhomedays_gas = max(gh12_e,gh13_e,gh14_e,gh15_e)
	label var stayhomedays_gas "Number of days missed from work/school (GI_illness)"
gen stayhomedays_res = max(gh16_e,gh17_e,gh18_e,gh19_e)
	label var stayhomedays_res "Number of days missed from work/school (Respiratory)"
gen stayhomedays_skn = gh10_e
	label var stayhomedays_skn "Number of days missed from work/school (Skin_rash)"
	
* Not Measured: stayhomedays_uti
* Not Measured: stmatch

rename cramps stomach
gen stomachstdt = mdy(gh15_b1,gh15_b2,2003)
	format stomachstdt %d
	label var stomachstdt "stomach pain or cramps start date"

* Not Measured: stomachdays
* Not Measured: stomachstill

* NOTE: In Mission Bay, these questions asked about whether an
* individual -Limited- their daily activities due to the symptoms
* but for the other beaches the questions asked about whether they
* -missed- their daily activities due to the symptoms

gen byte stopdaily_ear = (gh8_f==1) | (gh9_f==1)
	replace stopdaily_ear = . if (gh8_f==.) & (gh9_f==.)
	label var stopdaily_ear "Missed other activities (Ear_infec)"
gen byte stopdaily_eye = (gh7_f==1)
	replace stopdaily_eye = . if (gh7_f==.)
	label var stopdaily_eye "Missed other activities (Eye_infec)"
gen byte stopdaily_gas = (gh12_f==1) | (gh13_f==1) | (gh14_f==1) | (gh15_f==1)
	replace stopdaily_gas = . if (gh12_f==.) & (gh13_f==.) & (gh14_f==.) & (gh15_f==.)
	label var stopdaily_gas "Missed other activities (GI_illness)"
gen byte stopdaily_res = (gh16_f==1) | (gh17_f==1) | (gh18_f==1) | (gh19_f==1)
	replace stopdaily_res = . if (gh16_f==.) & (gh17_f==.) & (gh18_f==.) & (gh19_f==.)
	label var stopdaily_res "Missed other activities (Respiratory)"
gen byte stopdaily_skn = gh10_f
	label var stopdaily_skn "Missed other activities (Skin_rash)"

* Not Measured: stopdaily_uti

gen stopdailydays_ear = max(gh8_g,gh9_g)
	label var stopdailydays_ear "Days missed from other activities (Ear_infec)"
gen stopdailydays_eye = gh9_g
	label var stopdailydays_eye "Days missed from other activities (Eye_infec)"
gen stopdailydays_gas = max(gh12_g,gh13_g,gh14_g,gh15_g)
	label var stopdailydays_gas "Days missed from other activities (GI_illness)"
gen stopdailydays_res  = max(gh16_g,gh17_g,gh18_g,gh19_g)
	label var stopdailydays_res "Days missed from other activities (Respiratory)"
gen stopdailydays_skn = gh10_g
	label var stopdailydays_skn "Days missed from other activities (Skin_rash)"

* Not Measured: stopdailydays_uti

rename sunburn_past sunbase

* Not Measured: sunburn
* Not Measured: sunburndays
* Not Measured: sunburnstdt
* Not Measured: sunburnstill

rename water_swallow swallwater

gen swam = anycontact
	label var swam "Person went swimming/wading"

* Not Measured: totdays

* Not Measured: urinarytractinfection
* Not Measured: urinarytractinfectiondays
* Not Measured: urinarytractinfectionstdt
* Not Measured: urinarytractinfectionstill
* Not Measured: utibase


gen byte visitdoc_ear = (gh8_c==1) | (gh9_c==1)
	replace visitdoc_ear = . if (gh8_c==.) & (gh9_c==.)
	label var visitdoc_ear "Visited doctor/nurse/clinic (Ear_infec)"
gen byte visitdoc_eye = (gh7_c==1)
	replace visitdoc_eye = . if (gh7_c==.)
	label var visitdoc_eye "Visited doctor/nurse/clinic (Eye_infec)"
gen byte visitdoc_gas =  (gh12_c==1) | (gh13_c==1) | (gh14_c==1) | (gh15_c==1)
	replace visitdoc_gas = . if (gh12_c==.) & (gh13_c==.) & (gh14_c==.) & (gh15_c==.)
	label var visitdoc_gas "Visited doctor/nurse/clinic (GI_illness)"
gen byte visitdoc_res = (gh16_c==1) | (gh17_c==1) | (gh18_c==1) | (gh19_c==1)
	replace visitdoc_res =. if (gh16_c==.) & (gh17_c==.) & (gh18_c==.) & (gh19_c==.)
	label var visitdoc_res "Visited doctor/nurse/clinic (Respiratory)"
gen byte visitdoc_skn = (gh10_c==1)
	replace visitdoc_skn = . if (gh10_c==.)
	label var visitdoc_skn "Visited doctor/nurse/clinic (Skin_rash)"

* Not Measured: visitdoc_uti

* Not Measured: visitdoctimes_ear
* Not Measured: visitdoctimes_eye
* Not Measured: visitdoctimes_gas
* Not Measured: visitdoctimes_res
* Not Measured: visitdoctimes_skn
* Not Measured: visitdoctimes_uti

* Not Measured: visitertimes_ear
* Not Measured: visitertimes_eye
* Not Measured: visitertimes_gas
* Not Measured: visitertimes_res
* Not Measured: visitertimes_skn
* Not Measured: visitertimes_uti

rename vomiting_past vomitingbase
gen vomitingstdt = mdy(gh13_b1,gh13_b2,2003)
	format vomitingstdt %d
	label var vomitingstdt "vomiting start date"

* Not Measured: vomitingdays
* Not Measured: vomitingnumber
* Not Measured: vomitingstill

gen watertime = water_minutes
	label var watertime "Time in water (mins)"
	replace watertime = 0 if anycontact==0

* Not Measured: wateryeyes
* Not Measured: wateryeyesallergy
* Not Measured: wateryeyesdays
* Not Measured: wateryeyesstdt
* Not Measured: wateryeyesstill

* NOTE: in MB the wave variable only includes windsurfing and "boarding"
gen byte wave = (saq5c==1) | (saq5f==1)
	label var wave "boog board/wave rid/surf/windsurf/wake-kite board"

* NOTE: in MB the working variable includes persons who work or go to school
rename work working
	label var working "Person working for pay"

* Not Measured: beachtype
* Not Measured: blockface
* Not Measured: diagnosis_ear
* Not Measured: diagnosis_eye
* Not Measured: diagnosis_gas
* Not Measured: diagnosis_res
* Not Measured: diagnosis_skn
* Not Measured: diagnosis_uti
* Not Measured: gicontact_any
* Not Measured: gicontact_int
* Not Measured: lat
* Not Measured: longi
* Not Measured: meanbathers
* Not Measured: miles

* Not Measured: protbrim
* Not Measured: protslv
* Not Measured: reapp
* Not Measured: seatdr
* Not Measured: spf
* Not Measured: statusc

* Not Measured: sunburn1
* Not Measured: sunburn2
* Not Measured: sunburn3
* Not Measured: sunburn4
* Not Measured: sunburn5
* Not Measured: sunburn6
* Not Measured: sunburn7
* Not Measured: sunburnother

* Not Measured: swimloc
* Not Measured: tanning
* Not Measured: tanos
* Not Measured: transect
* Not Measured: venfest

recode sex 0=2
label define sex 1 "Male" 2 "Female"
label values sex sex

recode anim_int eggs_int stopdaily_skn (5=0)

* recode variables from 2=0
local vlist "cutbase eggs_base rawmeat_base fish_base anim_base gicontact_base skinchron block prothat shade repel eatfood food earplugs noseplugs mask preg publicpool fish_int rawmeat_int"
recode `vlist' (2=0)


*----------------------------------------
* retain the beach code to merge these
* data to the water quality data
*----------------------------------------
gen beachcode = "Mission Bay " + substr(station,1,1)
	label var beachcode "Water quality sampling location beach code"
order beach beachcode

*----------------------------------------
* subset the dataset to relevant variables
*----------------------------------------
# delimit ;
keep beach beachcode indid hhid intdate teledate age sex 
watertime 
earplugs noseplugs mask
sunhr dig bsand algae
block prothat shade repel  
eatfood food drinksy publicpool
anim_base anim_int anim_any 
gicontact_base 
fish_base fish_int fish_any
rawmeat_base rawmeat_int rawmeat_any
eggs_base eggs_int eggs_any
rawfood
asth gichron skinchron allergy preg 
anycontact bodycontact headunder mouthwater swallwater wave
diarrhea nausea vomiting fever sorethroat cough runnynose earache eyeinfection cut rash
*allergy
working 
race hisp
hhinc
*_gas *_eye *_res *_ear  *_skn
*stdt
stopdaily_gas
gibase cutbase sorebase eyebase rashbase
watertime
swam samebeach differentbeach

;
# delimit cr

* Retained in other datasets, but missing in Mission Bay
* hhmem wetsuit dsun sanddry sandmouth sandwash algaemouth algaewash come pool diarrheanumber vomitingnumber urinarytractinfection * *_uti *still *days fevertemp* headache cold wateryeyes sunburn racehaw utibase  vomitbase earbase siteid stmatch groundwater berm



*---------------------------------------------
* output variable information for harmonization
*---------------------------------------------
preserve
desc, replace clear
label data "variables for mb-epi.dta"
save "~/13beaches-data/temp/mb-epi-vars.dta", replace
restore



*----------------------------------------
* save the data for combining it with
* other beaches
*----------------------------------------
label data "Mission Bay epi data, created by 3-format-mb-epi.do"
save "~/13beaches-data/final/mb-epi.dta", replace

codebook, c

log close


* write a codebook
log using "~/13beaches-data/final/mb-epi-codebook.txt", text replace
desc, s
*aorder
codebook
log close




exit





