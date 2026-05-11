This repo contains replication files for the paper "Causal Duration Analysis with  Diff-in-Diff" by Ben Deaner and Hyejin Ku.

Replication files for the paper are provided in MATLAB (note that we have not included the dataset)

We also include a STATA command in the .ado file and STATA .do file that replicates the results in the paper. The .ado file is written in the Mata language and is designed to be general purpose.


INSTRUCTIONS FOR USING THE STATA COMMAND are below. Note that the .ado file has (at time of writing) been tested only on the data used in the paper. 
Please email b.deaner@ucl.ac.uk with any queries or feedback on the command.



# durationdid — Duration Difference-in-Differences

Stata command implementing the estimator from Deaner and Ku (2025), "Causal Duration Analysis with Diff-in-Diff."

---

## Getting Started

`durationdid` is a standard Stata command. Although it uses Mata internally for performance, no knowledge of Mata is required — you call it from the Stata command window (or a `.do` file) just like any other Stata command:

```stata
durationdid absorbed_time treatment, tstar(202) breps(1000)
```

### Preparing your data

The command requires a cross-sectional dataset (one row per individual) with a variable `absorbed_time` recording the time at which each individual first enters the absorbing state. If your data are already in this form — for example, a variable recording the duration of an unemployment spell in days — you can use it directly.

If instead your data consist of a **panel of binary outcomes** $Y_{i,1}, Y_{i,2}, \ldots$ indicating whether individual $i$ has reached the absorbing state by period $t$, you will need to convert them to the `absorbed_time` format. The required variable is the first period in which $Y_{i,t} = 1$.

**Panel data in wide format** (separate variables `y1`, `y2`, ..., `yT` in one row per individual):

```stata
* Find the first period where Y = 1
gen absorbed_time = .
forvalues t = 1/100 {
    replace absorbed_time = `t' if y`t' == 1 & absorbed_time == .
}
* Individuals never absorbed: set to a value beyond the extrapolation end
replace absorbed_time = 101 if absorbed_time == .
```

**Panel data in long format** (variables `id`, `t`, `Y` with multiple rows per individual):

```stata
* Keep the first period where Y = 1 for each individual
keep if Y == 1
bysort id (t): keep if _n == 1
rename t absorbed_time

* Merge back to get individuals who were never absorbed
* (they will have absorbed_time == . after the merge)
```

Alternatively, in long format you can compute `absorbed_time` directly:

```stata
bysort id (t): egen absorbed_time = min(cond(Y == 1, t, .))
* Keep one row per individual
bysort id: keep if _n == 1
* Set never-absorbed individuals beyond the extrapolation end
replace absorbed_time = 101 if absorbed_time == .
```

In all cases, individuals who are never observed in the absorbing state should have `absorbed_time` set to a value exceeding the `extrapend` parameter (or the maximum of `absorbed_time` if `extrapend` is not specified). These individuals will be treated as right-censored — they contribute to the survival estimates but are never counted as having reached the absorbing state.

---

## Syntax

```
durationdid absorbed_time treatment [if] [in], tstar(#) [options]
```

### Required arguments

| Argument | Description |
|---|---|
| `absorbed_time` | Numeric variable containing the time at which each individual first enters the absorbing state (e.g., exits unemployment). Individuals who are never observed in the absorbing state should have values exceeding `extrapend`. |
| `treatment` | Binary variable equal to 1 for the treated group and 0 for the untreated group. |
| `tstar(#)` | Final pre-treatment period. Treatment is assumed to occur strictly after this time. |

### Options

| Option | Default | Description |
|---|---|---|
| `covariates(`*varlist*`)` | *(none)* | Variables used for covariate balancing. Weights are computed as density ratios of covariate cell frequencies between groups, among individuals not yet in the absorbing state at time 1. Cells with fewer than 2 untreated individuals are dropped. |
| `burnin(#)` | `1` | First period used in estimation of *c*. Periods before `burnin` receive zero weight in the estimation of the level difference (or ratio). |
| `extrapend(#)` | `max(absorbed_time)` | Final period for extrapolation. Results are computed for periods 1 through `extrapend`. |
| `spec(`*string*`)` | `cd` | Specification for the identifying restriction. `cd` = common dynamics (fixed level difference in hazard rates). `ph` = proportional hazards (fixed ratio of hazard rates). |
| `breps(#)` | `1000` | Number of block bootstrap replications for inference. |
| `level(#)` | `0.95` | Confidence level for pointwise and uniform confidence bands. |
| `plevel(#)` | `0.6` | Level for the parallel/proportional trends specification test. |
| `prefix(`*string*`)` | `durationdid` | Prefix for saved output files and exported figures. |
| `nograph` | *(off)* | Suppress all graph output. |
| `seed(#)` | `12345` | Random number seed for the bootstrap. |

---

## Description

`durationdid` estimates average treatment effects in settings where the outcome is a binary indicator that an individual has entered an absorbing state (e.g., exited unemployment, passed an exam, left a marriage). Standard difference-in-differences is generally inconsistent in such settings because mean outcomes converge mechanically over time. This command replaces the parallel trends assumption on mean outcomes with restrictions on group-specific hazard rates.

The command computes the negative log survival R_{k,t} = −ln(1 − E[Y_{i,t} | G_i = k]) for each group and period, then estimates the parameter *c* that governs the relationship between group-specific time-average hazards. Under the **common dynamics** specification (`spec(cd)`), *c* is the constant difference in hazard rates between groups. Under **proportional hazards** (`spec(ph)`), *c* is the constant ratio.

Counterfactual mean outcomes for the treated group are imputed by extrapolating the estimated relationship to the post-treatment period and inverting the log transformation. Treatment effects are the difference between observed and imputed counterfactual outcomes.

If `covariates()` is specified, covariate balancing weights are computed using discrete cell frequencies, following the approach of Abadie (2005). For each covariate cell, the weight for untreated individuals is the ratio of the cell's frequency among treated survivors to its frequency among untreated survivors. Cells with fewer than 2 untreated individuals are dropped and remaining weights are rescaled.

### Inference

Block bootstrap inference is used, resampling individuals (with their complete outcome histories) with replacement. Both pointwise and uniform confidence bands are computed for treatment effects and imputed counterfactual outcomes. The uniform bands provide simultaneous coverage over all post-treatment periods.

### Specification test

A test for pre-treatment parallel (or proportional) trends is computed. Under the null, the difference (or ratio) of time-average hazards is constant over the pre-treatment period. The test statistic is the maximum absolute standardized deviation from constancy, and the p-value is computed from the bootstrap distribution.

---

## Stored results

`durationdid` stores the following in `r()`:

| Result | Description |
|---|---|
| `r(c_hat)` | Estimated *c* (level difference or ratio) |
| `r(pvalue)` | P-value from the parallel/proportional trends test |
| `r(N_treated)` | Number of treated observations |
| `r(N_untreated)` | Number of untreated observations |
| `r(tstar)` | Final pre-treatment period |
| `r(spec)` | Specification used (`cd` or `ph`) |

### Saved datasets

Two datasets are saved to the working directory:

**`<prefix>_results.dta`** — Period-level results (one observation per period, *t* = 1, ..., `extrapend`):

| Variable | Description |
|---|---|
| `t` | Period |
| `H1` | Estimated time-average hazard for the treated group |
| `H2` | Estimated time-average hazard for the untreated group |
| `H1_imp` | Imputed counterfactual time-average hazard for the treated group |
| `EY1_imp` | Imputed counterfactual mean outcome for the treated group |
| `tau` | Estimated average treatment effect |
| `mean_Y1` | Observed mean outcome for the treated group |
| `mean_Y2` | Observed mean outcome for the untreated group |
| `ci_tau_u_lo`, `ci_tau_u_hi` | Uniform confidence band for the ATT |
| `ci_tau_p_lo`, `ci_tau_p_hi` | Pointwise confidence band for the ATT |
| `ci_EY1_u_lo`, `ci_EY1_u_hi` | Uniform confidence band for the imputed counterfactual |
| `ci_EY1_p_lo`, `ci_EY1_p_hi` | Pointwise confidence band for the imputed counterfactual |

**`<prefix>_delta.dta`** — Pre-treatment specification test (one observation per pre-treatment period from `burnin` to `tstar`):

| Variable | Description |
|---|---|
| `t` | Period |
| `delta` | Deviation of the time-average hazard difference (or ratio) from its value at *t*\* |
| `ci_lo`, `ci_hi` | Uniform confidence band for `delta` |

### Exported figures

If `nograph` is not specified, five figures are exported as PNG files:

| File | Description |
|---|---|
| `<prefix>_Hs.png` | Time-average hazards for both groups (full sample period) |
| `<prefix>_H0.png` | Time-average hazards with imputed counterfactual (from `burnin` onward) |
| `<prefix>_point.png` | Observed and imputed mean outcomes with confidence bands |
| `<prefix>_taus.png` | Treatment effect estimates with confidence bands |
| `<prefix>_test.png` | Pre-treatment specification test with confidence bands |

---

## Examples

### Basic usage

Estimate treatment effects under common dynamics with 1000 bootstrap replications:

```stata
durationdid duration treatment, tstar(202)
```

### With covariates and proportional hazards

```stata
durationdid duration treatment, tstar(202) covariates(start_day) spec(ph) breps(5000)
```

### Suppress graphs and use a custom burn-in

```stata
durationdid duration treatment, tstar(202) burnin(150) nograph
```

### Full application example

```stata
* Load and prepare data
import delimited "restud_data.csv", clear
gen double absorbed_time = dur * 7
gen byte D = treated_cohort
gen double start_day = D * begin + (1 - D) * (begin + 365)

* Common dynamics with covariate balancing
durationdid absorbed_time D, tstar(202) covariates(start_day) ///
    burnin(152) extrapend(365) spec(cd) breps(10000) ///
    level(0.95) plevel(0.6) prefix("1UIP") seed(1)

* Display results
di "Estimated c: " r(c_hat)
di "P-value: " r(pvalue)

* Load and inspect period-level results
use "1UIP_results.dta", clear
list t tau ci_tau_u_lo ci_tau_u_hi if t > 202 & t <= 220
```

---

## Installation

Copy `durationdid.ado` to your personal ado directory. To find this directory, type in Stata:

```stata
sysdir
```

The file should be placed in the directory listed next to `PERSONAL`. Alternatively, place it in the current working directory.

---

## Requirements

- Stata 16.0 or later (uses frames and Mata features from version 16).
- The core computation is implemented in Mata for performance.

---

## References

Deaner, B. and H. Ku (2025). "Causal Duration Analysis with Diff-in-Diff."

Abadie, A. (2005). "Semiparametric Difference-in-Differences Estimators." *Review of Economic Studies*, 72(1), 1–19.

Wooldridge, J.M. (2023). "Simple Approaches to Nonlinear Difference-in-Differences with Panel Data." *Econometrics Journal*, 26(3), C31–C66.

---

## Authors

Ben Deaner and Hyejin Ku, University College London.
