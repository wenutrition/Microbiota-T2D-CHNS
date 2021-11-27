
**--Glucose--

use "D:\CNHS\data_pre\gphe_prospect_northshouth.dta",clear 
drop if glu_v2==. 
tab district
sort district
qui foreach diet of varlist rice-others {
     by district: egen `diet'_mean= mean(`diet')
	 by district: replace `diet'=`diet'_mean if `diet'==.
}
inspect rice wheat
tab District
drop rice_mean-others_mean

sort District
qui foreach diet of varlist rice-others {
     by District: egen `diet'_mean= mean(`diet')
	 by District: replace `diet'=`diet'_mean if `diet'==.
}

gen veg=dveg+lveg

save "D:\CNHS\data_pre\gphe_prospect_northshouth_glu.dta",replace

use "D:\CNHS\data_pre\gphe_prospect_northshouth_glu.dta",clear
keep if District=="North" 

preserve
keep g125-g1830 glu_v2
rename glu_v2 glu
save "D:\CNHS\data_pre\gnorth_glu.dta",replace 
restore

tempname coef
tempfile res
postfile `coef' str200(model micro district) float(n rr lul uul p) str200(cov)  using "`res'", replace
global cov1   energy MET age   sex Glu_field BMI i.alcohol i.smoke i.education i.marrige  income i.city  cityscore 
global cov2   energy MET age   sex Glu_field BMI i.alcohol i.smoke i.education i.marrige  income i.city  cityscore rice wheat fruit nuts pork poultry milk egg fish OIL_VEG OIL_ANI veg


forvalue m=1/2{
foreach var of varlist std_ln_g125-std_ln_g1830 { 
    mixed std_glu_v2 `var' ${cov`m'} ||district:, covariance(unstructured) 
test `var'
 post `coef'   ("`m'") ("`var'") ("`i'") (e(N)) (_b[`var'])  (_b[`var']-1.96*_se[`var'])  (_b[`var']+1.96*_se[`var'])  (chi2tail(1,(_b[`var']/_se[`var'])^2))  ("${cov`m'}")
 }
 }

  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\CNHS\results_final_final\G_glu_update.xlsx",  firstrow(variables) sheet("CNHS_north_abundance") sheetreplace  
restore


use "D:\CNHS\data_pre\gphe_prospect_northshouth_glu.dta",clear
keep if District=="South"
tempname coef
tempfile res
postfile `coef' str200(model micro district) float(n rr lul uul p) str200(cov)  using "`res'", replace
global cov1   energy MET age   sex Glu_field BMI i.alcohol i.smoke i.education i.marrige  income i.city  cityscore 
global cov2   energy MET age   sex Glu_field BMI i.alcohol i.smoke i.education i.marrige  income i.city  cityscore rice wheat fruit nuts pork poultry milk egg fish OIL_VEG OIL_ANI


forvalue m=1/2{
foreach var of varlist std_ln_g125-std_ln_g1830 std_shannon-std_pielou_e{ 
    mixed std_glu_v2 `var' ${cov`m'} ||district:, covariance(unstructured) 
test `var'
 post `coef'   ("`m'") ("`var'") ("`i'") (e(N)) (_b[`var'])  (_b[`var']-1.96*_se[`var'])  (_b[`var']+1.96*_se[`var'])  (chi2tail(1,(_b[`var']/_se[`var'])^2))  ("${cov`m'}")
 }
 }

  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\CNHS\results_final_final\G_glu_update.xlsx",  firstrow(variables) sheet("CNHS_south_abundance") sheetreplace  
restore


 *--meta analysis-
 
import excel "D:\CNHS\results_final_final\G_glu_update.xlsx", sheet("CNHS_south") firstrow clear
gen District="South"
append using "D:\CNHS\data_pre\G_level_glucose_north.dta"
gen g=substr(micro,8,.)
keep if model=="1"
merge m:m g using "D:\CNHS\data_pre\g_taxo.dta" 
keep if _merge==3
drop _merge
merge m:m g using "D:\CNHS\data_pre\model1_glucose_metap.dta"
keep if p_meta<0.05
sort District taxo
metan rr lul uul, label(namevar=District) nooverall  nowt nobox  dp(4) ///
by(taxo)fixed effect(rr)  ///
         boxopt( mcolor(navy8) msymbol(square) ) ///
		  pointopt( msymbol(square) mcolor(navy8) msize(small) ///
            mlabposition(1) ) ///
         ciopt( lcolor(navy8) lwidth(medium) ) ///
///
 graphregion(fcolor(white)lcolor(white))  subtitle("")
graph save Graph "D:\CNHS\results_final_final\micro_glucose_CNHS_model1.gph", replace
graph export "D:\CNHS\results_final_final\micro_glucose_CNHS_model1.pdf", as(pdf) replace


**---------------------------Describe the T2D related microbiome distribution---------------

use "D:\CNHS\data_pre\gphe_prospect.dta",clear // 2772
gen District=""
replace District="North" if district==11 | district==21  | district==23  | district==37 | district==41  | district==61 //652
replace District="South" if district==31 | district==32  | district==33  | district==42  | district==43 | district==55| district==53 | district==45 | district==52 // 1189
egen g_all=rowtotal(g125-g1830)
 foreach micro of varlist g125-g1830{
   replace `micro'=`micro'/g_all
   } 
keep District g389	g932	g956	g964	g980	g1030	g1056	g1071	g1354	g1519	g1705	g414	g889	g908	g958	g965	g1024	g1038	g1674	g766	g849	g845	g900	g1792	g875
egen g_all=rowtotal(g389-g1792)

table1,by(District) vars(g389 conts \	g932 conts\ g956 conts\ g964 conts\ g980 conts\ g1030 conts\ g1056 conts\ g1071 conts\ g1354	conts\ g1519	conts\ g1705	conts\ g414 conts\g889	conts\ g908 conts\	g958 conts\ g965	conts\ g1024	conts\ g1038	conts\ g1674	conts\ g766 conts\ g849 conts\ g845 conts\	g900 conts\ g1792	conts\ g875  conts)format(%30.25f) onecol test pdp(2) saving (D:\CNHS\results_paper\selectmicroabundance_des, replace)
table1,vars(g389 conts \	g932 conts\ g956 conts\ g964 conts\ g980 conts\ g1030 conts\ g1056 conts\ g1071 conts\ g1354	conts\ g1519	conts\ g1705	conts\ g414 conts\g889	conts\ g908 conts\	g958 conts\ g965	conts\ g1024	conts\ g1038	conts\ g1674	conts\ g766 conts\ g849 conts\ g845 conts\	g900 conts\ g1792	conts\ g875  conts)format(%30.25f) onecol test pdp(2) saving (D:\CNHS\results_paper\selectmicroabundanceoverall_des, replace)


*---------------Sensisentive analysis (Multiple imputation of missing values)------------

*--multiple inputation--

use "D:\CNHS\data_pre\gphe_prospect_northshouth.dta",clear 
egen g_all=rowtotal(g125-g1830)

qui foreach micro of varlist g125-g1830{
	gen ln_`micro'=ln(`micro'+1) 
}
sort District
qui foreach micro of varlist ln_g125-ln_g1830 {
     by District: egen `micro'_mean= mean(`micro')
     by District: egen `micro'_sd  = sd(`micro')
     by District: gen  std_`micro' = (`micro'-`micro'_mean)/`micro'_sd 
	 drop `micro'_sd `micro'_mean
}
keep eid std_ln_g125-std_ln_g1830 
save "D:\CNHS\data_pre\gphe_prospect_glevel.dta",replace


use "D:\CNHS\data_pre\gphe_prospect.dta",clear // 2772
keep eid Glu_field glu_v2 HbA1c hba1c_v2 insulin ins_v2
merge 1:1 eid using "D:\CNHS\data_pre\phenotype_raw1.dta"
keep if _merge==3
drop _merge 
merge 1:1 eid using "D:\CNHS\data_pre\gphe_prospect_glevel.dta"
keep if _merge==3

tab smoke 
tab alcohol
replace smoke=. if smoke==9
replace alcohol=. if alcohol==9

gen veg=dveg+lveg
gen District=.
replace District=0 if district==11 | district==21  | district==23  | district==37 | district==41  | district==61 //North
replace District=1 if district==31 | district==32  | district==33  | district==42  | district==43 | district==55| district==53 | district==45 | district==52 // South

replace marrige=0 if marrige!=2
replace marrige=1 if marrige==2
tab marrige
inspect rice wheat fruit nuts pork poultry milk egg fish OIL_VEG OIL_ANI veg energy MET age  sex BMI alcohol smoke education marrige  income  city  cityscore District 

mi set wide
set seed 29390
mi register imputed Glu_field  HbA1c  insulin rice wheat fruit nuts pork poultry milk egg fish OIL_VEG OIL_ANI veg  energy  MET BMI cityscore income alcohol smoke //需填充的变量
mi register regular  age  sex  marrige  education    city   District //非缺失变量
mi describe

mi impute chained (regress) rice wheat fruit nuts pork poultry milk egg fish OIL_VEG OIL_ANI veg  energy  MET BMI cityscore income (logit) alcohol smoke = age  sex  marrige  education    city   District, noisily add(5)

save "D:\CNHS\data_pre\gphe_prospect_northshouth_missinginpute.dta",replace
