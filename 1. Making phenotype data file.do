********************************************************************************
*******************         MAKING DATASET WOMEN       *************************
********************************************************************************

set more off
use "\Women.dta", clear

keep n_eid n_19_0_0 n_21_0_0 n_31_0_0 n_34_0_0 n_50_0_0 n_52_0_0 n_54_0_0 n_189_0_0 ///
n_738_0_0 n_845_0_0 n_1239_0_0 n_1249_0_0 n_1259_0_0 n_1558_0_0 n_1647_0_0 n_2040_0_0  ///
 n_2129_0_0 n_2139_0_0 n_2149_0_0 n_2159_0_0 n_2178_0_0 n_2463_0_0 n_2644_0_0 n_2714_0_0 n_2724_0_0 n_2734_0_0 n_2744_0_0 n_2754_0_0 ///
 n_2764_0_0 n_2774_0_0 n_2784_0_0 n_2794_0_0 n_2804_0_0 n_2814_0_0 n_2824_0_0 n_2867_0_0 ///
 n_2887_0_0 n_2897_0_0 n_2907_0_0 n_3140_0_0 n_3160_0_0 n_3436_0_0 n_3456_0_0 n_3466_0_0 ///
 n_3486_0_0 n_3581_0_0 n_3591_0_0 n_3659_0_0 n_3669_0_0 n_3700_0_0 n_3710_0_0 n_3720_0_0 ///
 n_3731_0_0 n_3829_0_0 n_3839_0_0 n_3872_0_0 n_4407_0_0 n_4418_0_0 n_4429_0_0 n_4440_0_0 ///
 n_4451_0_0 n_4462_0_0  n_5057_0_0 n_6138_0_0 n_6138_0_1 n_6138_0_2 n_6138_0_3 ///
 n_6138_0_4 n_6138_0_5 n_6142_0_0 n_6142_0_1 n_6142_0_2 n_6142_0_3 n_6142_0_4 n_6142_0_5 ///
 n_6142_0_6 n_6183_0_0 n_6194_0_0 n_6218_0_0 n_20016_0_0 n_20047_0_0 n_20048_0_0 n_20115_0_0 ///
 n_20116_0_0 n_20117_0_0 n_21000_0_0 n_21001_0_0 n_21002_0_0 n_21003_0_0 ///
 n_23104_0_0 n_23106_0_0 n_40007_0_0
 
save "\dataset_women500.dta", replace

********************************************************************************
*************             MERGING WITHDRAWALS/SPLIT         ********************
********************************************************************************

merge 1:1 n_eid using "\withdrawals.dta"
drop if _merge ==2
drop if _merge ==3
drop _merge

* merge in interim genetic release data so can have a var to identify who is in 150 or 350k (interim vs. full genetic data releases)
merge m:1 n_eid using "\Genetic parameters.dta", keepus(n_22003_0_0 n_22001_0_0)
drop if _merge==1 
drop if _merge==2
drop _merge
gen split = .
replace split = 1 if n_22001_0_0!=.
replace split = 2 if split!=1

********************************************************************************
*******************           SAVING WOMEN DATASET         *********************
********************************************************************************

tab n_31_0_0
gen sex = n_31_0_0
tab sex

save "\dataset_women500.dta", replace












