********************************************************************************
* UIP Application: Duration Diff-in-Diff
* Replication of UIP_application_revision.m
* Ben Deaner and Hyejin Ku, UCL
********************************************************************************

clear all
set more off

* Parameters
local burn_in = 152
local t_star = 202
local cohort_end = 273
local extrap_end = 365
local initial_PBD_end = 210
local bootstrap_reps = 10000
local level = 0.95
local parallel_level = 0.6

* Load data
import delimited "restud_data.csv", clear

* Durations in days
gen double duration_days = dur * 7

* Reform date: earliest begin date among post-reform individuals
summarize begin if after == 1, meanonly
local reform_day = r(min)
di "Reform date: `reform_day'"

********************************************************************************
* Sample selection
********************************************************************************

* Eligible: PBD extension (t39==1) but NOT replacement rate change (tr==0)
gen byte eligible = (tr == 0) & (t39 == 1)

* Untreated cohort: became unemployed between
*   (reform - PBD_end - 365) and (reform - cohort_end)
gen byte untreated = (begin > `reform_day' - `initial_PBD_end' - 365) & ///
                     (begin <= `reform_day' - `cohort_end')

* Treated cohort: became unemployed between
*   (reform - PBD_end) and (reform + 365 - cohort_end)
gen byte treated = (begin <= `reform_day' + 365 - `cohort_end') & ///
                   (begin > `reform_day' - `initial_PBD_end')

* Keep eligible individuals in one of the two cohorts
keep if eligible == 1 & (treated == 1 | untreated == 1)

* Treatment indicator
gen byte D = treated

* Absorbed time = spell duration in days
gen double absorbed_time = duration_days

* Covariate: start day of year, aligned across cohorts
gen double start_day = D * begin + (1 - D) * (begin + 365)

* Report sample sizes
count if D == 1
di "Treated: " r(N)
count if D == 0
di "Untreated: " r(N)

********************************************************************************
* Figure: Raw mean outcomes
********************************************************************************

preserve
    local T = `extrap_end'
    clear
    quietly set obs `T'
    gen int t = _n
    gen double y1 = .
    gen double y2 = .
    save "_temp_periods.dta", replace
restore

forvalues t = 1/`extrap_end' {
    quietly count if D == 1
    local n1 = r(N)
    quietly count if D == 1 & absorbed_time <= `t'
    local y1_`t' = r(N) / `n1'

    quietly count if D == 0
    local n0 = r(N)
    quietly count if D == 0 & absorbed_time <= `t'
    local y2_`t' = r(N) / `n0'
}

preserve
    use "_temp_periods.dta", clear
    forvalues t = 1/`extrap_end' {
        quietly replace y1 = `y1_`t'' in `t'
        quietly replace y2 = `y2_`t'' in `t'
    }

    twoway (line y1 t, lcolor(blue) lwidth(medthick)) ///
           (line y2 t, lcolor(red) lwidth(medthick)), ///
           xline(210, lpattern(dash_dot)) xline(273) ///
           legend(order(1 "Treated Cohort" 2 "Untreated") pos(5) ring(0)) ///
           xtitle("t") ytitle("Mean Outcome") ///
           scheme(s2color)
    graph export "UIP_fig_levels.png", replace width(960) height(720)
restore
erase "_temp_periods.dta"


********************************************************************************
* Common Dynamics specification
********************************************************************************

di _newline(2) "=============================================="
di "COMMON DYNAMICS SPECIFICATION"
di "=============================================="

durationdid absorbed_time D, ///
    tstar(`t_star') ///
    covariates(start_day) ///
    burnin(`burn_in') ///
    extrapend(`extrap_end') ///
    spec(cd) ///
    breps(`bootstrap_reps') ///
    level(`level') ///
    plevel(`parallel_level') ///
    prefix("1UIP") ///
    seed(1)


********************************************************************************
* Proportional Hazards specification
********************************************************************************

di _newline(2) "=============================================="
di "PROPORTIONAL HAZARDS SPECIFICATION"
di "=============================================="

durationdid absorbed_time D, ///
    tstar(`t_star') ///
    covariates(start_day) ///
    burnin(`burn_in') ///
    extrapend(`extrap_end') ///
    spec(ph) ///
    breps(`bootstrap_reps') ///
    level(`level') ///
    plevel(`parallel_level') ///
    prefix("2UIP") ///
    seed(1)


di _newline "Application complete."
