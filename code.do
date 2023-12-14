net install grc1leg, from("http://www.stata.com/users/vwiggins") replace
net install gr0075, from("http://www.stata-journal.com/software/sj18-4") replace
ssc install labutil, replace
ssc install sencode, replace
ssc install panelview, all
*net install panelview, all replace from("https://yiqingxu.org/packages/panelview_stata")

ssc install bacondecomp, replace
ssc install twowayfeweights, replace
ssc install did_multiplegt 
ssc install drdid, replace
ssc install csdid, replace
ssc install hdfe, replace
ssc install avar, replace
ssc install eventstudyweights, replace
ssc install eventstudyinteract, replace
ssc install did_imputation, replace
ssc install event_plot, replace
ssc install reghdfe, replace
ssc install estout, replace
search svmat2 //dm79

**-------Part1 Data Preprocessing-------**

// Load data, control the maximum number of variables
cd "C:\Users\HUAWEI\Desktop\Summer School\DID reading list for undergraduate students\project"
use "C:\Users\HUAWEI\Desktop\Images and Code\HRS_long.dta", clear
set matsize 800
global outcomes "oop_spend riearnsemp "

drop if wave < 7 // Keep data for the last 5 years of calendar time
bys hhidpn: gen N = _N // N is the number of observations for each ID
keep if N == 5 // Keep individuals with exactly 5 periods of data

bys hhidpn: egen flag = min(evt_time) // flag is the min(wave - E)
drop if flag >= 0 & flag != . // Remove IDs where the first hospitalization is before or in the 7th period
drop if flag == . // Remove NA values
drop flag

bys hhidpn: egen wave_hosp_copy = min(wave_hosp)
replace wave_hosp = wave_hosp_copy
drop wave_hosp_copy // Set hospitalization waves to the minimum hospitalization wave for each ID

tab wave_hosp

keep if ever_hospitalized // Keep samples with hospitalization in waves 8-11
gen post = (wave >= wave_hosp)

xi i.wave // Convert the categorical variable wave into dummy variables, construct a cohort dummy variable
tab wave_hosp, gen(wave_hosp_)

keep if age_hosp <= 59
// Count the number of observations for each cohort {0, 1, 2, 3}
count if wave_hosp == 7 & wave == 7
global c_7 = r(N) // Global macro
count if wave_hosp == 8 & wave == 7
global c_8 = r(N)
count if wave_hosp == 9 & wave == 7
global c_9 = r(N)
count if wave_hosp == 10 & wave == 7
global c_10 = r(N)

gen F4event = evt_time == -4
tab F4event
forvalues l = 3(-1)1 {
    gen F`l'event = evt_time == -`l'
}

forvalues l = 0/3 {
    gen L`l'event = evt_time == `l'
}
save HRS_long.dta, replace

*-------Panel Data Visualization-------*
tab wave_hosp
panelview oop_spend post, i(hhidpn) t(wave) type(treat) xtitle("Wave") ytitle("Patient ID") ///
title("Treatment Status") bytiming prepost ylabel(none)

**-------Part2 Identify Negative Weight Issues--------**

*-------De Chaisemartin and d'Haultfoeuille to examine the extent of negative weights-------*
use HRS_long.dta, clear
twowayfeweights oop_spend hhidpn wave post, type(feTR)

/*Results:
Under the common trends assumption, beta estimates a weighted sum of 1248 ATTs. 
979 ATTs receive a positive weight, and 269 receive a negative weight.
The sum of the positive weights is equal to 1.2307007.
The sum of the negative weights is equal to -.23070082.

beta is compatible with a DGP where the average of those ATTs is equal to 0,
while their standard deviation is equal to 2263.8308.

beta is compatible with a DGP where those ATTs all are of a different sign than beta,
while their standard deviation is equal to 5267.3159.
*/

/*Explanation:
It corresponds to the minimal value of the standard deviation of the treatment effect across the treated groups and time periods under which beta and the average treatment effect on the treated (ATT) could be of opposite signs. 
When that number is large, this means that beta and the ATT can only be of opposite signs if there is a lot of treatment effect heterogeneity across groups and time periods. 
When that number is low, this means that beta and the ATT can be of opposite signs even if there is not a lot of treatment effect heterogeneity across groups and time periods.
*/

*-------Goodman-Bacon-------*
xtset hhidpn wave
bacondecomp riearnsemp post, ddetail robust // Dependent variable is income

mat list e(sumdd)
mat A = e(sumdd)
clear
svmat2 A, names(col) rnames(group) full

tab group
sum
g temp = Beta * TotalWeight

// Summing for each group separately
preserve
collapse (sum) temp TotalWeight, by(group)
g BetaGroup = temp / TotalWeight
list group BetaGroup TotalWeight

/*Results: Negative weights exist in the Late_v_Early group
     +------------------------------------+
     |        group   BetaGr~p   TotalW~t |
     |------------------------------------|
  1. | Early_v_Late   -3351.67   .4668037 |
  2. | Late_v_Early   2605.522   .5331963 |
     +------------------------------------+
*/

*-------Event Study Method-------*
// Using eventstudyweights function to calculate weights
use HRS_long.dta, clear
drop F1event
drop F4event

eventstudyweights F*event L*event, absorb(hhidpn wave) cohort(wave_hosp) rel_time(evt_time) saveweights("weights")

import excel "weights.xlsx", clear firstrow

keep F2event wave_hosp evt_time
tab evt_time

// Relative time, displaying the F2 cohort, i.e., wave_hosp in {8,9,10,11}, weights for each CATT in the cohort. Positive weights indicate that the heterogeneity of treatment effects in other periods affects F2.
preserve

reshape wide F2event, i(evt_time) j(wave_hosp)
graph twoway line F2event* evt_time, xtitle("relative time") ///
ytitle("weight in TWFE F2event coefficient") graphregion(fcolor(white)) scheme(sj)

restore
// Check the sum of weights using SUN's weight decomposition
tab evt_time
recode evt_time (-1 = -1) (-4 = -1) ///
(-2 = 1) (nonm = 0), g(evt_time_group)
tab evt_time_group
collapse (sum) F2event, by(evt_time_group)
list
/* -1 represents excluded relative time points (including -1 and -4)
    1 represents the relative time points to be tested (here, F2)
	0 represents other relative time points
     +----------------------+
     | evt_ti~p     F2event |
     |----------------------|
  1. |       -1          -1 |
  2. |        0   3.539e-16 |
  3. |        1           1 |
     +----------------------+
*/

**-------Part3 Estimation-------

*-------TWFE Estimation-------
use HRS_long.dta, clear
reghdfe oop_spend post, a(hhidpn wave) cluster(wave_hosp)
/* Regression results:
HDFE Linear regression                                  Number of obs =  2,624
Absorbing 2 HDFE groups                                 F(1,655) =       30.66
Statistics robust to heteroskedasticity                 Prob > F =      0.0000
                                                        R-squared =     0.4053
                                                        Adj R-squared = 0.2054
                                                        Within R-sq.=   0.0184
Number of clusters (hhidpn)=656                   Root MSE =   6992.7668

                               (Std. err. adjusted for 656 clusters in
>  hhidpn)
------------------------------------------------------------------------------
             |               Robust
   oop_spend | Coefficient  std. err.      t    P>|t|     [95% conf. interval]
-------------+----------------------------------------------------------------
        post |   3041.337   549.2324     5.54   0.000     1962.868    4119.806
       _cons |   2006.972   266.0345     7.54   0.000     1484.589    2529.355
------------------------------------------------------------------------------

Absorbed degrees of freedom:
-----------------------------------------------------+
 Absorbed FE | Categories  - Redundant  = Num. Coefs |
-------------+---------------------------------------|
      hhidpn |       656         656           0    *|
        wave |         4           0           4     |
-----------------------------------------------------+
* = FE nested within cluster; treated as redundant for DoF computation
*/

*-------Dynamic DID-------
use HRS_long.dta, clear
drop F4event
drop F1event
reghdfe oop_spend F*event L*event, a(hhidpn wave) cluster(hhidpn)

event_plot, default_look stub_lag(L#event) stub_lead(F#event) together ///
graph_opt(xtitle("Relative time") ytitle("Coefficient") xlabel(-4(1)3) ///
	title("Hospitalization and Out-of-pocket spend"))
graph export fig_ols.pdf, replace

estimates store ols_d 

/* Results:
HDFE Linear regression                            Number of obs   =      3,280
Absorbing 2 HDFE groups                           F(   6,    655) =       9.50
Statistics robust to heteroskedasticity           Prob > F        =     0.0000
                                                  R-squared       =     0.3470
                                                  Adj R-squared   =     0.1805
                                                  Within R-sq.    =     0.0218
Number of clusters (hhidpn)  =        656         Root MSE        =  7036.6646

                               (Std. err. adjusted for 656 clusters in hhidpn)
------------------------------------------------------------------------------
             |               Robust
   oop_spend | Coefficient  std. err.      t    P>|t|     [95% conf. interval]
-------------+----------------------------------------------------------------
     F3event |   148.9752   792.1174     0.19   0.851     -1406.42    1704.371
     F2event |   202.6305   480.2132     0.42   0.673    -740.3125    1145.574
     L0event |   3012.641   511.1341     5.89   0.000     2008.982      4016.3
     L1event |   888.0903    664.332     1.34   0.182     -416.387    2192.568
     L2event |    1171.52   983.4892     1.19   0.234    -759.6522    3102.692
     L3event |   1914.034   1425.849     1.34   0.180    -885.7519     4713.82
       _cons |   2422.494   387.7893     6.25   0.000     1661.034    3183.954
------------------------------------------------------------------------------

Absorbed degrees of freedom:
-----------------------------------------------------+
 Absorbed FE | Categories  - Redundant  = Num. Coefs |
-------------+---------------------------------------|
      hhidpn |       656         656           0    *|
        wave |         5           0           5     |
-----------------------------------------------------+
* = FE nested within cluster; treated as redundant for DoF computation
*/

*-------Inverse-Weighted Estimation-------
use HRS_long.dta, clear
gen lastcohort = wave_hosp == . & post == 0 // dummy for never-treated cohort
tab wave_hosp lastcohort, m
tab lastcohort

preserve 
drop if wave_hosp == . & post == 1 // remove the always-treated group
// Estimate the true model with dynamic treatment effects
gen last_cohort = (wave_hosp == 11)
eventstudyinteract oop_spend F3-F2 L0-L2 ///
	if ever_hospitalized & wave < 11 , cohort(wave_hosp) absorb(hhidpn wave) ///
	control_cohort(last_cohort) vce(cluster hhidpn)
	
/*result：
IW estimates for dynamic effects                     Number of obs =     2,624
Absorbing 2 HDFE groups                              F(9, 655)     =      6.42
                                                     Prob > F      =    0.0000
                                                     R-squared     =    0.4113
                                                     Adj R-squared =    0.2101
                                                     Root MSE      = 6971.9298
                               (Std. err. adjusted for 656 clusters in hhidpn)
------------------------------------------------------------------------------
             |               Robust
   oop_spend | Coefficient  std. err.      t    P>|t|     [95% conf. interval]
-------------+----------------------------------------------------------------
     F3event |   591.0464   1272.833     0.46   0.643    -1908.278    3090.371
     F2event |   352.6393   697.6528     0.51   0.613    -1017.266    1722.545
     L0event |   2960.045   540.9075     5.47   0.000     1897.923    4022.167
     L1event |   529.7669   586.9677     0.90   0.367    -622.7985    1682.332
     L2event |   800.1065   1010.619     0.79   0.429    -1184.337     2784.55
------------------------------------------------------------------------------
*/

//short model without dynamic treatment effect
eventstudyinteract oop_spend F3-F2 L0-L2 ///
	if ever_hospitalized & wave < 11 , cohort(wave_hosp) absorb(hhidpn wave) ///
	control_cohort(last_cohort)

/*result：
IW estimates for dynamic effects                     Number of obs =     2,624
Absorbing 2 HDFE groups                              F(9, 1956)    =      6.31
                                                     Prob > F      =    0.0000
                                                     R-squared     =    0.4113
                                                     Adj R-squared =    0.2105
                                                     Root MSE      = 6970.1474
------------------------------------------------------------------------------
   oop_spend | Coefficient  Std. err.      t    P>|t|     [95% conf. interval]
-------------+----------------------------------------------------------------
     F3event |   591.0464   1446.019     0.41   0.683    -2244.853    3426.946
     F2event |   352.6393   813.6326     0.43   0.665    -1243.039    1948.317
     L0event |   2960.045   608.1964     4.87   0.000     1767.264    4152.826
     L1event |   529.7669   999.1559     0.53   0.596    -1429.755    2489.289
     L2event |   800.1065   1371.291     0.58   0.560    -1889.239    3489.452
------------------------------------------------------------------------------
*/