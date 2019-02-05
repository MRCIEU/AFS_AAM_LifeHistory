********************************************************************************
******************************* DATA *******************************************
********************************************************************************

use "\dataset_women500.dta", clear
set more off

********************************************************************************
******************************** EXPOSURES *************************************
********************************************************************************

* Age at menarche
gen age_menarche = n_2714_0_0
replace age_menarche=. if age_menarche <0
replace age_menarche=. if age_menarche <5
tab age_menarche

* Age at first sex
gen age_firstsex = n_2139_0_0
replace age_firstsex=. if age_firstsex <0
replace age_firstsex=. if age_firstsex<age_menarche 

* log of z scores too
zscore age_firstsex
gen age_firstsex_log = log(z_age_firstsex)

********************************************************************************
**************************        OUTCOMES         *****************************
********************************************************************************

* Childless
gen childless = n_2734_0_0
replace childless=. if n_2734_0_0 <0
replace childless=100 if n_2734_0_0==0
replace childless=0 if childless!=100 & childless!=.
replace childless=1 if childless==100
**recode childless (0 = 1) (1 = 2), gen(childlessplink)

* Age at first birth
gen age_fbirth_multi = n_2754_0_0
replace age_fbirth_multi=. if n_2754_0_0 <0
*replace age_fbirth_multi=. if age_fbirth_multi>45

* Age at last birth
gen age_lbirth = n_2764_0_0
replace age_lbirth=. if age_lbirth <0
* replace age_lbirth=. if age_lbirth>48

* Birthrange
gen birthrange_m = age_lbirth - age_fbirth_multi

* Number of children 
gen num_kids = n_2734_0_0
replace num_kids=. if num_kids <0
browse if num_kids >11 & num_kids!=.
replace num_kids=. if num_kids>11 // 99.99 percentile

* Qualification...
replace n_6138_0_0=. if n_6138_0_0==-3 
replace n_6138_0_0=0 if n_6138_0_0==-7 
replace n_6138_0_1=. if n_6138_0_1==-3
replace n_6138_0_1=0 if n_6138_0_1==-7
replace n_6138_0_2=. if n_6138_0_2==-3
replace n_6138_0_2=0 if n_6138_0_2==-7
replace n_6138_0_3=. if n_6138_0_3==-3
replace n_6138_0_3=0 if n_6138_0_3==-7
replace n_6138_0_4=. if n_6138_0_4==-3
replace n_6138_0_4=0 if n_6138_0_4==-7
replace n_6138_0_5=. if n_6138_0_5==-3
replace n_6138_0_5=0 if n_6138_0_5==-7
* Subjects who selected several response categories were assigned the higher qualification.
egen qualification=rowmax(n_6138_0_0 n_6138_0_1 n_6138_0_2 n_6138_0_3 n_6138_0_4 n_6138_0_5)
* iscd conversion to years
gen edu_years = .
replace edu_years=20 if qualification==1
replace edu_years=13 if qualification==2
replace edu_years=10 if qualification==3
replace edu_years=10 if qualification==4
replace edu_years=19 if qualification==5
replace edu_years=15 if qualification==6
replace edu_years=7 if qualification==0
summ edu_years
* iscd conversion to college binary
gen college=.
replace college=1 if qualification==1
replace college=0 if qualification==2
replace college=0 if qualification==3
replace college=0 if qualification==4
replace college=1 if qualification==5
replace college=0 if qualification==6
replace college=0 if qualification==0
**recode college(1=2) (0=1), gen(collegeplink)

* Age left education 
gen age_edu = n_845_0_0
replace age_edu=. if age_edu <0

* Risk-taking
gen risk = n_2040_0_0
replace risk=. if risk <0
**recode risk (0 = 1) (1 = 2), gen(riskplink)

* Ever smoked 
replace n_1239_0_0=. if n_1239_0_0 <0 //current smoking
replace n_1249_0_0=. if n_1249_0_0 <0 // past smoking

gen ever_smok = .
replace ever_smok =1 if n_1239_0_0==1 | n_1239_0_0==2 | n_1249_0_0==1 | n_1249_0_0==2 | n_1249_0_0==3
replace ever_smok=0 if n_1239_0_0==0 & n_1249_0_0==4 
**recode ever_smok (0 = 1) (1 = 2), gen(ever_smokplink)
 
* Alcohol intake
gen alc_intake = n_1558_0_0
replace alc_intake=. if alc_intake <0

* Number of sexual partners
gen num_partners = n_2149_0_0
replace num_partners=. if n_2149_0_0 <0
replace num_partners=. if num_partners>200 // 99.99th percentile

********************************************************************************
**********************        OTHER VARIABLES         **************************
********************************************************************************

* Age at assessment centre
tab n_21003_0_0
gen age_centre = n_21003_0_0

* birthyear
tab n_34_0_0
gen birthyear=n_34_0_0

*******************************************************************************
*******************                  SAVE            ***************************
********************************************************************************
keep n_eid split sex age_menarche age_firstsex z_age_firstsex age_firstsex_log ///
childless age_fbirth_multi age_lbirth birthrange_m num_kids qualification edu_years ///
college age_edu risk ever_smok alc_intake num_partners ///
age_centre birthyear
 
save "\dataset_women500_cleaned.dta", replace
outsheet using "women_cleaned500.csv", comma replace

