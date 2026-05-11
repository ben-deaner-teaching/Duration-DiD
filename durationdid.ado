*! durationdid v1.0 - Duration Difference-in-Differences
*! Ben Deaner and Hyejin Ku, UCL
*! Implements the estimator from "Causal Duration Analysis with Diff-in-Diff"

capture program drop durationdid
program define durationdid, rclass
    version 16.0
    syntax varlist(min=2 max=2) [if] [in], ///
        Tstar(integer)                      /// Final pre-treatment period
        [                                   ///
        COVariates(varlist)                 /// Covariate variable(s) for balancing
        SWeights(varname)                   /// Sampling weights variable
        BURNin(integer 1)                   /// First period used in estimation of c
        EXTRAPend(integer 0)                /// Last period to extrapolate to (0 = max)
        SPEC(string)                        /// "cd" for common dynamics (default), "ph" for proportional hazards
        BREPS(integer 1000)                 /// Bootstrap replications
        Level(real 0.95)                    /// Confidence level
        PLevel(real 0.6)                    /// Parallel trends test level
        PREfix(string)                      /// Prefix for saved figures
        NOGRaph                             /// Suppress graphs
        SEED(integer 12345)                 /// Random seed
        ]

    * Parse inputs
    tokenize `varlist'
    local absorbed_time `1'
    local D `2'
    if "`spec'" == "" local spec "cd"
    if !inlist("`spec'", "cd", "ph") {
        di as error "spec() must be cd (common dynamics) or ph (proportional hazards)"
        exit 198
    }
    if "`prefix'" == "" local prefix "durationdid"

    marksample touse
    if "`covariates'" != "" {
        markout `touse' `covariates'
    }

    * Determine extrapolation end
    if `extrapend' == 0 {
        quietly summarize `absorbed_time' if `touse'
        local extrapend = floor(r(max))
    }

    local T = `extrapend'

    * Report sample sizes
    quietly count if `touse' & `D' == 1
    local n_treated = r(N)
    quietly count if `touse' & `D' == 0
    local n_untreated = r(N)
    di as text "Treated observations: `n_treated'"
    di as text "Untreated observations: `n_untreated'"

    * Load data into Mata and run estimation
    preserve
    quietly keep if `touse'

    local sw_var ""
    if "`sweights'" != "" {
        local sw_var "`sweights'"
    }

    if "`covariates'" != "" {
        mata: _dd_main("`absorbed_time'", "`D'", "`covariates'", "`sw_var'", ///
            `tstar', `burnin', `T', "`spec'", `breps', `level', `plevel', `seed')
    }
    else {
        mata: _dd_main("`absorbed_time'", "`D'", "", "`sw_var'", ///
            `tstar', `burnin', `T', "`spec'", `breps', `level', `plevel', `seed')
    }
    restore

    * Retrieve scalar results
    local c_hat = r(c_hat)
    local pvalue = r(pvalue)

    di as text _newline "Estimated c: " %9.6f `c_hat'
    di as text "Parallel trends test p-value: " %6.4f `pvalue'

    * Produce graphs
    if "`nograph'" == "" {
        quietly {
            preserve
            use "`prefix'_results.dta", clear

            * Figure 1: Time-average hazards (full sample)
            twoway (line H1 t, lcolor(blue) lwidth(medthick)) ///
                   (line H2 t, lcolor(red) lwidth(medthick)), ///
                   xline(`tstar', lpattern(dash_dot)) ///
                   legend(order(1 "Treated" 2 "Untreated") pos(5) ring(0)) ///
                   xtitle("t") ytitle("Time-Average Hazard") ///
                   scheme(s2color) name(Hs, replace)
            graph export "`prefix'_Hs.png", replace width(960) height(720)

            * Figure 2: Imputed time-average hazards
            twoway (line H1 t if t >= `burnin', lcolor(blue) lwidth(medthick)) ///
                   (line H1_imp t if t >= `burnin', lcolor(blue) lwidth(medthick) lpattern(dash)) ///
                   (line H2 t if t >= `burnin', lcolor(red) lwidth(medthick)), ///
                   xline(`tstar', lpattern(dash_dot)) ///
                   legend(order(1 "Treated" 2 "Imputed" 3 "Untreated") pos(7) ring(0)) ///
                   xtitle("t") ytitle("Time-Average Hazard") ///
                   scheme(s2color) name(H0, replace)
            graph export "`prefix'_H0.png", replace width(960) height(720)

            * Figure 3: Imputed mean outcomes with CIs
            twoway (line mean_Y1 t if t >= `burnin', lcolor(blue) lwidth(medthick)) ///
                   (line EY1_imp t if t > `burnin', lcolor(blue) lwidth(medthick) lpattern(dash)) ///
                   (line mean_Y2 t if t >= `burnin', lcolor(red) lwidth(medthick)) ///
                   (rarea ci_EY1_u_lo ci_EY1_u_hi t if t > `tstar', color(gs10%40) lwidth(none)) ///
                   (rarea ci_EY1_p_lo ci_EY1_p_hi t if t > `tstar', color(gs6%40) lwidth(none)), ///
                   xline(`tstar', lpattern(dash_dot)) ///
                   legend(order(1 "Treated" 2 "Imputed" 3 "Untreated" ///
                          4 "`=`level'*100'% Uniform CIs" 5 "`=`level'*100'% Pointwise CIs") pos(5) ring(0)) ///
                   xtitle("t") ytitle("Mean Outcome") ///
                   scheme(s2color) name(point, replace)
            graph export "`prefix'_point.png", replace width(960) height(720)

            * Figure 4: Treatment effects with CIs
            twoway (line tau t if t > `tstar', lcolor(black) lwidth(medthick)) ///
                   (rarea ci_tau_u_lo ci_tau_u_hi t if t > `tstar', color(gs10%40) lwidth(none)) ///
                   (rarea ci_tau_p_lo ci_tau_p_hi t if t > `tstar', color(gs6%40) lwidth(none)), ///
                   yline(0, lpattern(dash) lcolor(black)) ///
                   xline(`tstar', lpattern(dash_dot)) ///
                   legend(order(1 "Estimated Effect" ///
                          2 "`=`level'*100'% Uniform CIs" 3 "`=`level'*100'% Pointwise CIs") pos(5) ring(0)) ///
                   xtitle("t") ytitle("ATT") ///
                   scheme(s2color) name(taus, replace)
            graph export "`prefix'_taus.png", replace width(960) height(720)

            restore

            * Figure 5: Parallel trends test
            preserve
            use "`prefix'_delta.dta", clear
            twoway (line delta t, lcolor(green) lwidth(medthick)) ///
                   (rarea ci_lo ci_hi t, color(gs10%40) lwidth(none)), ///
                   yline(0, lpattern(dash) lcolor(black)) ///
                   xline(`tstar', lpattern(dash_dot)) ///
                   legend(order(1 "delta" 2 "`=`plevel'*100'% Uniform CIs") pos(1) ring(0)) ///
                   xtitle("t") ///
                   scheme(s2color) name(test, replace)
            graph export "`prefix'_test.png", replace width(960) height(720)
            restore
        }
    }

    * Return scalars
    return scalar c_hat = `c_hat'
    return scalar pvalue = `pvalue'
    return scalar N_treated = `n_treated'
    return scalar N_untreated = `n_untreated'
    return scalar tstar = `tstar'
    return local spec "`spec'"
end


********************************************************************************
* MATA IMPLEMENTATION
********************************************************************************
mata:
mata set matastrict off

// =================================================================
// _dd_make_groups: Create integer group indices from covariate matrix
// Input:  X (n x p matrix)
// Output: grp (n x 1 vector of group IDs 1,...,G)
// =================================================================
real colvector _dd_make_groups(real matrix X)
{
    real scalar n, p, i, g
    real matrix X_s
    real colvector ord, grp_s, grp

    n = rows(X)
    p = cols(X)

    ord = order(X, (1..p))
    X_s = X[ord, .]

    grp_s = J(n, 1, 1)
    g = 1
    for (i = 2; i <= n; i++) {
        if (X_s[i, .] != X_s[i-1, .]) g++
        grp_s[i] = g
    }

    grp = J(n, 1, .)
    grp[ord] = grp_s
    return(grp)
}


// =================================================================
// _dd_compute_weights: Covariate balancing weights
// Matches covariate_balancing() in durationDiD.m
// sw: sampling weights (pre-normalized to mean 1)
// =================================================================
real colvector _dd_compute_weights(real colvector at,
                                   real matrix X,
                                   real colvector D,
                                   real colvector sw)
{
    real scalar n, ngroups, g
    real scalar wn_treat, wn_untreat
    real scalar frac_kept_treat, frac_kept_untr
    real colvector grp, surv
    real colvector wcnt_treat, wcnt_untreat, cnt_untreat
    real colvector wcnt_treat_i, wcnt_untreat_i, cnt_untreat_i
    real colvector weights, censor
    real colvector treat_idx, untreat_idx

    n = rows(at)
    grp = _dd_make_groups(X)
    ngroups = max(grp)

    surv = (at :>= 1)

    // Sampling-weighted cell frequencies for density ratios
    wcnt_treat = J(ngroups, 1, 0)
    wcnt_untreat = J(ngroups, 1, 0)
    // Unweighted counts for cell censoring
    cnt_untreat = J(ngroups, 1, 0)
    for (g = 1; g <= ngroups; g++) {
        wcnt_treat[g] = sum(sw :* ((grp :== g) :& (D :== 1) :& surv))
        wcnt_untreat[g] = sum(sw :* ((grp :== g) :& (D :== 0) :& surv))
        cnt_untreat[g] = sum((grp :== g) :& (D :== 0) :& surv)
    }

    // Weighted totals
    wn_treat = sum(sw :* ((D :== 1) :& surv))
    wn_untreat = sum(sw :* ((D :== 0) :& surv))

    // Map to individual level
    wcnt_treat_i = wcnt_treat[grp]
    wcnt_untreat_i = wcnt_untreat[grp]
    cnt_untreat_i = cnt_untreat[grp]

    // Weights: weighted density ratio for untreated, 1 for treated
    weights = J(n, 1, 1)
    untreat_idx = selectindex(D :== 0)
    if (length(untreat_idx) > 0 & wn_treat > 0 & wn_untreat > 0) {
        weights[untreat_idx] = (wcnt_treat_i[untreat_idx] :/ wn_treat) :/ ///
                               (wcnt_untreat_i[untreat_idx] :/ wn_untreat)
    }

    // Censor cells with < 2 untreated (based on unweighted counts)
    censor = (cnt_untreat_i :< 2)
    weights = weights :* (1 :- censor)

    // Rescale
    treat_idx = selectindex(D :== 1)
    if (length(treat_idx) > 0) {
        frac_kept_treat = mean(1 :- censor[treat_idx])
        if (frac_kept_treat > 0) {
            weights[treat_idx] = weights[treat_idx] :/ frac_kept_treat
        }
    }
    if (length(untreat_idx) > 0) {
        frac_kept_untr = mean(1 :- censor[untreat_idx])
        if (frac_kept_untr > 0) {
            weights[untreat_idx] = weights[untreat_idx] :/ frac_kept_untr
        }
    }

    return(weights)
}


// =================================================================
// _dd_weighted_means: Compute weighted mean survivals for each group
// Vectorized: builds n x T survival matrix via broadcasting
//
// E1_Yk[t] = (1/n_k) * sum_i w_i * I(at_i > t)
// =================================================================
void _dd_weighted_means(real colvector at,
                        real colvector D,
                        real colvector weights,
                        real scalar T_end,
                        real rowvector E1_Y1,
                        real rowvector E1_Y2)
{
    real colvector treat_idx, untreat_idx
    real colvector at_treat, at_untreat, w_treat, w_untreat
    real scalar n1, n2
    real rowvector periods
    real matrix surv_treat, surv_untreat

    treat_idx = selectindex(D :== 1)
    untreat_idx = selectindex(D :== 0)

    at_treat = at[treat_idx]
    at_untreat = at[untreat_idx]
    w_treat = weights[treat_idx]
    w_untreat = weights[untreat_idx]

    n1 = length(treat_idx)
    n2 = length(untreat_idx)

    periods = (1..T_end)

    // Survival matrices: surv[i,t] = I(at[i] > t)
    // Broadcasting: at (n x 1) vs periods (1 x T)
    surv_treat = (at_treat * J(1, T_end, 1)) :> (J(n1, 1, 1) * periods)
    surv_untreat = (at_untreat * J(1, T_end, 1)) :> (J(n2, 1, 1) * periods)

    // Weighted mean survival
    E1_Y1 = mean(w_treat :* surv_treat)
    E1_Y2 = mean(w_untreat :* surv_untreat)
}


// =================================================================
// _dd_estimate: Point estimates from weighted mean survivals
// Computes H, c, imputation, ATT, and delta
// =================================================================
void _dd_estimate(real rowvector E1_Y1,
                  real rowvector E1_Y2,
                  real scalar burn_in,
                  real scalar t_star,
                  real scalar T_end,
                  string scalar spec,
                  real rowvector H1,
                  real rowvector H2,
                  real rowvector H1_imp,
                  real rowvector EY1_imp,
                  real rowvector tau,
                  real rowvector delta,
                  real scalar c_hat)
{
    real rowvector periods, R1_imp, diffs, vals

    periods = (1..T_end)

    // Time-average hazards
    H1 = -(ln(E1_Y1) :- ln(E1_Y1[1])) :/ (periods :- 1)
    H2 = -(ln(E1_Y2) :- ln(E1_Y2[1])) :/ (periods :- 1)
    H1[1] = 0
    H2[1] = 0

    // Estimate c and impute
    if (spec == "cd") {
        diffs = H1[1, burn_in..t_star] - H2[1, burn_in..t_star]
        vals = select(diffs, diffs :< .)
        c_hat = (length(vals) > 0 ? mean(vals') : .)
        H1_imp = c_hat :+ H2
        delta = diffs :- (H1[t_star] - H2[t_star])
    }
    else {
        diffs = H1[1, burn_in..t_star] :/ H2[1, burn_in..t_star]
        vals = select(diffs, diffs :< .)
        c_hat = (length(vals) > 0 ? mean(vals') : .)
        H1_imp = c_hat :* H2
        delta = diffs :- (H1[t_star] / H2[t_star])
    }

    // Imputed R_{1,t} and counterfactual mean outcomes
    R1_imp = H1_imp :* (periods :- 1) :- ln(E1_Y1[1])

    EY1_imp = 1 :- E1_Y1
    if (t_star < T_end) {
        EY1_imp[1, (t_star+1)..T_end] = 1 :- exp(-R1_imp[1, (t_star+1)..T_end])
    }

    // ATT
    tau = (1 :- E1_Y1) :- EY1_imp
}


// =================================================================
// _dd_main: Orchestrates estimation and bootstrap
// Called from the Stata program via mata: _dd_main(...)
// =================================================================
void _dd_main(string scalar at_var,
              string scalar D_var,
              string scalar X_vars,
              string scalar sw_var,
              real scalar t_star,
              real scalar burn_in,
              real scalar T_end,
              string scalar spec,
              real scalar breps,
              real scalar level,
              real scalar plevel,
              real scalar seed)
{
    // ---- Load data ----
    real colvector at, D, weights, sw
    real matrix X
    real scalar n, has_covs, has_sw, b, ndelta
    string scalar prefix

    at = st_data(., at_var)
    D = st_data(., D_var)
    n = rows(at)

    has_covs = (X_vars != "")
    if (has_covs) {
        X = st_data(., X_vars)
    }

    has_sw = (sw_var != "")
    if (has_sw) {
        sw = st_data(., sw_var)
        sw = sw :/ mean(sw)
    }
    else {
        sw = J(n, 1, 1)
    }

    prefix = st_local("prefix")

    // ---- Covariate balancing weights ----
    if (has_covs) {
        weights = _dd_compute_weights(at, X, D, sw) :* sw
        printf("Treated after cell drops: %g\n", sum((D :== 1) :& (weights :> 0)))
        printf("Untreated after cell drops: %g\n", sum((D :== 0) :& (weights :> 0)))
    }
    else {
        weights = sw
    }

    // ---- Point estimates ----
    real rowvector E1_Y1, E1_Y2, H1, H2, H1_imp, EY1_imp, tau, delta
    real scalar c_hat

    _dd_weighted_means(at, D, weights, T_end, E1_Y1, E1_Y2)
    _dd_estimate(E1_Y1, E1_Y2, burn_in, t_star, T_end, spec,
                 H1, H2, H1_imp, EY1_imp, tau, delta, c_hat)

    ndelta = t_star - burn_in + 1
    printf("Estimated c = %12.8f\n", c_hat)

    // ---- Bootstrap ----
    printf("Running %g bootstrap replications...\n", breps)
    rseed(seed)

    real matrix boot_tau, boot_EY1, boot_delta, boot_idx
    real rowvector E1_Y1_b, E1_Y2_b, H1_b, H2_b, H1_imp_b
    real rowvector EY1_imp_b, tau_b, delta_b
    real colvector idx, at_b, D_b, weights_b, sw_b
    real matrix X_b
    real scalar c_b

    boot_tau = J(breps, T_end, .)
    boot_EY1 = J(breps, T_end, .)
    boot_delta = J(breps, ndelta, .)
    boot_idx = ceil(n * runiform(n, breps))

    for (b = 1; b <= breps; b++) {
        if (mod(b, 500) == 0) {
            printf("  Bootstrap %g / %g\n", b, breps)
        }

        idx = boot_idx[., b]
        at_b = at[idx]
        D_b = D[idx]
        sw_b = sw[idx]
        sw_b = sw_b :/ mean(sw_b)

        if (has_covs) {
            X_b = X[idx, .]
            weights_b = _dd_compute_weights(at_b, X_b, D_b, sw_b) :* sw_b
        }
        else {
            weights_b = sw_b
        }

        _dd_weighted_means(at_b, D_b, weights_b, T_end, E1_Y1_b, E1_Y2_b)
        _dd_estimate(E1_Y1_b, E1_Y2_b, burn_in, t_star, T_end, spec,
                     H1_b, H2_b, H1_imp_b, EY1_imp_b, tau_b, delta_b, c_b)

        boot_tau[b, .] = tau_b
        boot_EY1[b, .] = EY1_imp_b
        boot_delta[b, .] = delta_b
    }

    printf("Bootstrap complete.\n")

    // ---- Confidence bands ----
    real rowvector sd_tau, sd_EY1, sd_delta
    real matrix adev_tau, adev_EY1, adev_delta
    real colvector max_adev_tau, max_adev_EY1, max_adev_delta
    real scalar cv_tau_u, cv_EY1_u, cv_delta_u, idx_level, idx_plevel
    real rowvector cv_tau_p, cv_EY1_p
    real scalar pvalue, t
    real colvector col_sorted

    // Standard deviations
    sd_tau = sqrt(diagonal(variance(boot_tau)))'
    sd_EY1 = sqrt(diagonal(variance(boot_EY1)))'
    sd_delta = sqrt(diagonal(variance(boot_delta)))'

    // Guard against zero SDs
    sd_tau = sd_tau + (sd_tau :== 0) :* 1e-100
    sd_EY1 = sd_EY1 + (sd_EY1 :== 0) :* 1e-100
    sd_delta = sd_delta + (sd_delta :== 0) :* 1e-100

    // Standardized absolute deviations
    adev_tau = abs(boot_tau :- tau) :/ sd_tau
    adev_EY1 = abs(boot_EY1 :- EY1_imp) :/ sd_EY1
    adev_delta = abs(boot_delta :- delta) :/ sd_delta

    // Uniform critical values
    max_adev_tau = rowmax(adev_tau[., (t_star+1)..T_end])
    max_adev_EY1 = rowmax(adev_EY1[., (t_star+1)..T_end])
    max_adev_delta = rowmax(adev_delta)

    idx_level = ceil(level * breps)
    idx_plevel = ceil(plevel * breps)

    cv_tau_u = sort(max_adev_tau, 1)[idx_level]
    cv_EY1_u = sort(max_adev_EY1, 1)[idx_level]
    cv_delta_u = sort(max_adev_delta, 1)[idx_plevel]

    // Pointwise critical values
    cv_tau_p = J(1, T_end, .)
    cv_EY1_p = J(1, T_end, .)
    for (t = 1; t <= T_end; t++) {
        col_sorted = sort(adev_tau[., t], 1)
        cv_tau_p[t] = col_sorted[idx_level]
        col_sorted = sort(adev_EY1[., t], 1)
        cv_EY1_p[t] = col_sorted[idx_level]
    }

    // P-value (matches MATLAB exactly)
    pvalue = mean(max_adev_delta :> abs(max(delta :/ sd_delta)))
    printf("Parallel trends test p-value = %6.4f\n", pvalue)

    // ---- Construct CI bounds ----
    real rowvector ci_tau_u_lo, ci_tau_u_hi, ci_tau_p_lo, ci_tau_p_hi
    real rowvector ci_EY1_u_lo, ci_EY1_u_hi, ci_EY1_p_lo, ci_EY1_p_hi
    real rowvector ci_delta_lo, ci_delta_hi

    ci_tau_u_lo = tau - cv_tau_u :* sd_tau
    ci_tau_u_hi = tau + cv_tau_u :* sd_tau
    ci_tau_p_lo = tau - cv_tau_p :* sd_tau
    ci_tau_p_hi = tau + cv_tau_p :* sd_tau

    ci_EY1_u_lo = EY1_imp - cv_EY1_u :* sd_EY1
    ci_EY1_u_hi = EY1_imp + cv_EY1_u :* sd_EY1
    ci_EY1_p_lo = EY1_imp - cv_EY1_p :* sd_EY1
    ci_EY1_p_hi = EY1_imp + cv_EY1_p :* sd_EY1

    ci_delta_lo = delta - cv_delta_u :* sd_delta
    ci_delta_hi = delta + cv_delta_u :* sd_delta

    // ---- Store results in Stata ----
    st_numscalar("r(c_hat)", c_hat)
    st_numscalar("r(pvalue)", pvalue)

    // --- Main results dataset ---
    stata("clear")
    st_addobs(T_end)

    (void) st_addvar("int", "t")
    (void) st_addvar("double", "H1")
    (void) st_addvar("double", "H2")
    (void) st_addvar("double", "H1_imp")
    (void) st_addvar("double", "EY1_imp")
    (void) st_addvar("double", "tau")
    (void) st_addvar("double", "mean_Y1")
    (void) st_addvar("double", "mean_Y2")
    (void) st_addvar("double", "ci_tau_u_lo")
    (void) st_addvar("double", "ci_tau_u_hi")
    (void) st_addvar("double", "ci_tau_p_lo")
    (void) st_addvar("double", "ci_tau_p_hi")
    (void) st_addvar("double", "ci_EY1_u_lo")
    (void) st_addvar("double", "ci_EY1_u_hi")
    (void) st_addvar("double", "ci_EY1_p_lo")
    (void) st_addvar("double", "ci_EY1_p_hi")

    st_store(., "t", (1::T_end))
    st_store(., "H1", H1')
    st_store(., "H2", H2')
    st_store(., "H1_imp", H1_imp')
    st_store(., "EY1_imp", EY1_imp')
    st_store(., "tau", tau')
    st_store(., "mean_Y1", (1 :- E1_Y1)')
    st_store(., "mean_Y2", (1 :- E1_Y2)')
    st_store(., "ci_tau_u_lo", ci_tau_u_lo')
    st_store(., "ci_tau_u_hi", ci_tau_u_hi')
    st_store(., "ci_tau_p_lo", ci_tau_p_lo')
    st_store(., "ci_tau_p_hi", ci_tau_p_hi')
    st_store(., "ci_EY1_u_lo", ci_EY1_u_lo')
    st_store(., "ci_EY1_u_hi", ci_EY1_u_hi')
    st_store(., "ci_EY1_p_lo", ci_EY1_p_lo')
    st_store(., "ci_EY1_p_hi", ci_EY1_p_hi')

    stata("save " + char(34) + prefix + "_results.dta" + char(34) + ", replace")

    // --- Delta dataset ---
    stata("clear")
    st_addobs(ndelta)

    (void) st_addvar("int", "t")
    (void) st_addvar("double", "delta")
    (void) st_addvar("double", "ci_lo")
    (void) st_addvar("double", "ci_hi")

    st_store(., "t", (burn_in::t_star))
    st_store(., "delta", delta')
    st_store(., "ci_lo", ci_delta_lo')
    st_store(., "ci_hi", ci_delta_hi')

    stata("save " + char(34) + prefix + "_delta.dta" + char(34) + ", replace")
}

end
*! durationdid v1.0 - Duration Difference-in-Differences
*! Ben Deaner and Hyejin Ku, UCL
*! Implements the estimator from "Causal Duration Analysis with Diff-in-Diff"

capture program drop durationdid
program define durationdid, rclass
    version 16.0
    syntax varlist(min=2 max=2) [if] [in], ///
        Tstar(integer)                      /// Final pre-treatment period
        [                                   ///
        COVariates(varlist)                 /// Covariate variable(s) for balancing
        BURNin(integer 1)                   /// First period used in estimation of c
        EXTRAPend(integer 0)                /// Last period to extrapolate to (0 = max)
        SPEC(string)                        /// "cd" for common dynamics (default), "ph" for proportional hazards
        BREPS(integer 1000)                 /// Bootstrap replications
        Level(real 0.95)                    /// Confidence level
        PLevel(real 0.6)                    /// Parallel trends test level
        PREfix(string)                      /// Prefix for saved figures
        NOGRaph                             /// Suppress graphs
        SEED(integer 12345)                 /// Random seed
        ]

    * Parse inputs
    tokenize `varlist'
    local absorbed_time `1'
    local D `2'
    if "`spec'" == "" local spec "cd"
    if !inlist("`spec'", "cd", "ph") {
        di as error "spec() must be cd (common dynamics) or ph (proportional hazards)"
        exit 198
    }
    if "`prefix'" == "" local prefix "durationdid"

    marksample touse
    if "`covariates'" != "" {
        markout `touse' `covariates'
    }

    * Determine extrapolation end
    if `extrapend' == 0 {
        quietly summarize `absorbed_time' if `touse'
        local extrapend = floor(r(max))
    }

    local T = `extrapend'

    * Report sample sizes
    quietly count if `touse' & `D' == 1
    local n_treated = r(N)
    quietly count if `touse' & `D' == 0
    local n_untreated = r(N)
    di as text "Treated observations: `n_treated'"
    di as text "Untreated observations: `n_untreated'"

    * Load data into Mata and run estimation
    preserve
    quietly keep if `touse'

    if "`covariates'" != "" {
        mata: _dd_main("`absorbed_time'", "`D'", "`covariates'", ///
            `tstar', `burnin', `T', "`spec'", `breps', `level', `plevel', `seed')
    }
    else {
        mata: _dd_main("`absorbed_time'", "`D'", "", ///
            `tstar', `burnin', `T', "`spec'", `breps', `level', `plevel', `seed')
    }
    restore

    * Retrieve scalar results
    local c_hat = r(c_hat)
    local pvalue = r(pvalue)

    di as text _newline "Estimated c: " %9.6f `c_hat'
    di as text "Parallel trends test p-value: " %6.4f `pvalue'

    * Produce graphs
    if "`nograph'" == "" {
        quietly {
            preserve
            use "`prefix'_results.dta", clear

            * Figure 1: Time-average hazards (full sample)
            twoway (line H1 t, lcolor(blue) lwidth(medthick)) ///
                   (line H2 t, lcolor(red) lwidth(medthick)), ///
                   xline(`tstar', lpattern(dash_dot)) ///
                   legend(order(1 "Treated" 2 "Untreated") pos(5) ring(0)) ///
                   xtitle("t") ytitle("Time-Average Hazard") ///
                   scheme(s2color) name(Hs, replace)
            graph export "`prefix'_Hs.png", replace width(960) height(720)

            * Figure 2: Imputed time-average hazards
            twoway (line H1 t if t >= `burnin', lcolor(blue) lwidth(medthick)) ///
                   (line H1_imp t if t >= `burnin', lcolor(blue) lwidth(medthick) lpattern(dash)) ///
                   (line H2 t if t >= `burnin', lcolor(red) lwidth(medthick)), ///
                   xline(`tstar', lpattern(dash_dot)) ///
                   legend(order(1 "Treated" 2 "Imputed" 3 "Untreated") pos(7) ring(0)) ///
                   xtitle("t") ytitle("Time-Average Hazard") ///
                   scheme(s2color) name(H0, replace)
            graph export "`prefix'_H0.png", replace width(960) height(720)

            * Figure 3: Imputed mean outcomes with CIs
            twoway (line mean_Y1 t if t >= `burnin', lcolor(blue) lwidth(medthick)) ///
                   (line EY1_imp t if t > `burnin', lcolor(blue) lwidth(medthick) lpattern(dash)) ///
                   (line mean_Y2 t if t >= `burnin', lcolor(red) lwidth(medthick)) ///
                   (rarea ci_EY1_u_lo ci_EY1_u_hi t if t > `tstar', color(gs10%40) lwidth(none)) ///
                   (rarea ci_EY1_p_lo ci_EY1_p_hi t if t > `tstar', color(gs6%40) lwidth(none)), ///
                   xline(`tstar', lpattern(dash_dot)) ///
                   legend(order(1 "Treated" 2 "Imputed" 3 "Untreated" ///
                          4 "`=`level'*100'% Uniform CIs" 5 "`=`level'*100'% Pointwise CIs") pos(5) ring(0)) ///
                   xtitle("t") ytitle("Mean Outcome") ///
                   scheme(s2color) name(point, replace)
            graph export "`prefix'_point.png", replace width(960) height(720)

            * Figure 4: Treatment effects with CIs
            twoway (line tau t if t > `tstar', lcolor(black) lwidth(medthick)) ///
                   (rarea ci_tau_u_lo ci_tau_u_hi t if t > `tstar', color(gs10%40) lwidth(none)) ///
                   (rarea ci_tau_p_lo ci_tau_p_hi t if t > `tstar', color(gs6%40) lwidth(none)), ///
                   yline(0, lpattern(dash) lcolor(black)) ///
                   xline(`tstar', lpattern(dash_dot)) ///
                   legend(order(1 "Estimated Effect" ///
                          2 "`=`level'*100'% Uniform CIs" 3 "`=`level'*100'% Pointwise CIs") pos(5) ring(0)) ///
                   xtitle("t") ytitle("ATT") ///
                   scheme(s2color) name(taus, replace)
            graph export "`prefix'_taus.png", replace width(960) height(720)

            restore

            * Figure 5: Parallel trends test
            preserve
            use "`prefix'_delta.dta", clear
            twoway (line delta t, lcolor(green) lwidth(medthick)) ///
                   (rarea ci_lo ci_hi t, color(gs10%40) lwidth(none)), ///
                   yline(0, lpattern(dash) lcolor(black)) ///
                   xline(`tstar', lpattern(dash_dot)) ///
                   legend(order(1 "delta" 2 "`=`plevel'*100'% Uniform CIs") pos(1) ring(0)) ///
                   xtitle("t") ///
                   scheme(s2color) name(test, replace)
            graph export "`prefix'_test.png", replace width(960) height(720)
            restore
        }
    }

    * Return scalars
    return scalar c_hat = `c_hat'
    return scalar pvalue = `pvalue'
    return scalar N_treated = `n_treated'
    return scalar N_untreated = `n_untreated'
    return scalar tstar = `tstar'
    return local spec "`spec'"
end


********************************************************************************
* MATA IMPLEMENTATION
********************************************************************************
mata:
mata set matastrict off

// =================================================================
// _dd_make_groups: Create integer group indices from covariate matrix
// Input:  X (n x p matrix)
// Output: grp (n x 1 vector of group IDs 1,...,G)
// =================================================================
real colvector _dd_make_groups(real matrix X)
{
    real scalar n, p, i, g
    real matrix X_s
    real colvector ord, grp_s, grp

    n = rows(X)
    p = cols(X)

    ord = order(X, (1..p))
    X_s = X[ord, .]

    grp_s = J(n, 1, 1)
    g = 1
    for (i = 2; i <= n; i++) {
        if (X_s[i, .] != X_s[i-1, .]) g++
        grp_s[i] = g
    }

    grp = J(n, 1, .)
    grp[ord] = grp_s
    return(grp)
}


// =================================================================
// _dd_compute_weights: Covariate balancing weights
// Matches covariate_balancing() in durationDiD.m
// =================================================================
real colvector _dd_compute_weights(real colvector at,
                                   real matrix X,
                                   real colvector D)
{
    real scalar n, ngroups, n_treat, n_untreat, g
    real scalar frac_kept_treat, frac_kept_untr
    real colvector grp, surv, cnt_treat, cnt_untreat
    real colvector cnt_treat_i, cnt_untreat_i
    real colvector weights, censor
    real colvector treat_idx, untreat_idx

    n = rows(at)
    grp = _dd_make_groups(X)
    ngroups = max(grp)

    surv = (at :>= 1)

    cnt_treat = J(ngroups, 1, 0)
    cnt_untreat = J(ngroups, 1, 0)
    for (g = 1; g <= ngroups; g++) {
        cnt_treat[g] = sum((grp :== g) :& (D :== 1) :& surv)
        cnt_untreat[g] = sum((grp :== g) :& (D :== 0) :& surv)
    }

    n_treat = sum((D :== 1) :& surv)
    n_untreat = sum((D :== 0) :& surv)

    // Map counts to individual level
    cnt_treat_i = cnt_treat[grp]
    cnt_untreat_i = cnt_untreat[grp]

    // Weights: density ratio for untreated, 1 for treated
    weights = J(n, 1, 1)
    untreat_idx = selectindex(D :== 0)
    if (length(untreat_idx) > 0 & n_treat > 0 & n_untreat > 0) {
        weights[untreat_idx] = (cnt_treat_i[untreat_idx] :/ n_treat) :/ ///
                               (cnt_untreat_i[untreat_idx] :/ n_untreat)
    }

    // Censor cells with < 2 untreated; zero out weights
    censor = (cnt_untreat_i :< 2)
    weights = weights :* (1 :- censor)

    // Rescale
    treat_idx = selectindex(D :== 1)
    if (length(treat_idx) > 0) {
        frac_kept_treat = mean(1 :- censor[treat_idx])
        if (frac_kept_treat > 0) {
            weights[treat_idx] = weights[treat_idx] :/ frac_kept_treat
        }
    }
    if (length(untreat_idx) > 0) {
        frac_kept_untr = mean(1 :- censor[untreat_idx])
        if (frac_kept_untr > 0) {
            weights[untreat_idx] = weights[untreat_idx] :/ frac_kept_untr
        }
    }

    return(weights)
}


// =================================================================
// _dd_weighted_means: Compute weighted mean survivals for each group
// Vectorized: builds n x T survival matrix via broadcasting
//
// E1_Yk[t] = (1/n_k) * sum_i w_i * I(at_i > t)
// =================================================================
void _dd_weighted_means(real colvector at,
                        real colvector D,
                        real colvector weights,
                        real scalar T_end,
                        real rowvector E1_Y1,
                        real rowvector E1_Y2)
{
    real colvector treat_idx, untreat_idx
    real colvector at_treat, at_untreat, w_treat, w_untreat
    real scalar n1, n2
    real rowvector periods
    real matrix surv_treat, surv_untreat

    treat_idx = selectindex(D :== 1)
    untreat_idx = selectindex(D :== 0)

    at_treat = at[treat_idx]
    at_untreat = at[untreat_idx]
    w_treat = weights[treat_idx]
    w_untreat = weights[untreat_idx]

    n1 = length(treat_idx)
    n2 = length(untreat_idx)

    periods = (1..T_end)

    // Survival matrices: surv[i,t] = I(at[i] > t)
    // Broadcasting: at (n x 1) vs periods (1 x T)
    surv_treat = (at_treat * J(1, T_end, 1)) :> (J(n1, 1, 1) * periods)
    surv_untreat = (at_untreat * J(1, T_end, 1)) :> (J(n2, 1, 1) * periods)

    // Weighted mean survival
    E1_Y1 = mean(w_treat :* surv_treat)
    E1_Y2 = mean(w_untreat :* surv_untreat)
}


// =================================================================
// _dd_estimate: Point estimates from weighted mean survivals
// Computes H, c, imputation, ATT, and delta
// =================================================================
void _dd_estimate(real rowvector E1_Y1,
                  real rowvector E1_Y2,
                  real scalar burn_in,
                  real scalar t_star,
                  real scalar T_end,
                  string scalar spec,
                  real rowvector H1,
                  real rowvector H2,
                  real rowvector H1_imp,
                  real rowvector EY1_imp,
                  real rowvector tau,
                  real rowvector delta,
                  real scalar c_hat)
{
    real rowvector periods, R1_imp, diffs, vals

    periods = (1..T_end)

    // Time-average hazards
    H1 = -(ln(E1_Y1) :- ln(E1_Y1[1])) :/ (periods :- 1)
    H2 = -(ln(E1_Y2) :- ln(E1_Y2[1])) :/ (periods :- 1)
    H1[1] = 0
    H2[1] = 0

    // Estimate c and impute
    if (spec == "cd") {
        diffs = H1[1, burn_in..t_star] - H2[1, burn_in..t_star]
        vals = select(diffs, diffs :< .)
        c_hat = (length(vals) > 0 ? mean(vals') : .)
        H1_imp = c_hat :+ H2
        delta = diffs :- (H1[t_star] - H2[t_star])
    }
    else {
        diffs = H1[1, burn_in..t_star] :/ H2[1, burn_in..t_star]
        vals = select(diffs, diffs :< .)
        c_hat = (length(vals) > 0 ? mean(vals') : .)
        H1_imp = c_hat :* H2
        delta = diffs :- (H1[t_star] / H2[t_star])
    }

    // Imputed R_{1,t} and counterfactual mean outcomes
    R1_imp = H1_imp :* (periods :- 1) :- ln(E1_Y1[1])

    EY1_imp = 1 :- E1_Y1
    if (t_star < T_end) {
        EY1_imp[1, (t_star+1)..T_end] = 1 :- exp(-R1_imp[1, (t_star+1)..T_end])
    }

    // ATT
    tau = (1 :- E1_Y1) :- EY1_imp
}


// =================================================================
// _dd_main: Orchestrates estimation and bootstrap
// Called from the Stata program via mata: _dd_main(...)
// =================================================================
void _dd_main(string scalar at_var,
              string scalar D_var,
              string scalar X_vars,
              real scalar t_star,
              real scalar burn_in,
              real scalar T_end,
              string scalar spec,
              real scalar breps,
              real scalar level,
              real scalar plevel,
              real scalar seed)
{
    // ---- Load data ----
    real colvector at, D, weights
    real matrix X
    real scalar n, has_covs, b, ndelta
    string scalar prefix

    at = st_data(., at_var)
    D = st_data(., D_var)
    n = rows(at)

    has_covs = (X_vars != "")
    if (has_covs) {
        X = st_data(., X_vars)
    }

    prefix = st_local("prefix")

    // ---- Covariate balancing weights ----
    if (has_covs) {
        weights = _dd_compute_weights(at, X, D)
        printf("Treated after cell drops: %g\n", sum((D :== 1) :& (weights :> 0)))
        printf("Untreated after cell drops: %g\n", sum((D :== 0) :& (weights :> 0)))
    }
    else {
        weights = J(n, 1, 1)
    }

    // ---- Point estimates ----
    real rowvector E1_Y1, E1_Y2, H1, H2, H1_imp, EY1_imp, tau, delta
    real scalar c_hat

    _dd_weighted_means(at, D, weights, T_end, E1_Y1, E1_Y2)
    _dd_estimate(E1_Y1, E1_Y2, burn_in, t_star, T_end, spec,
                 H1, H2, H1_imp, EY1_imp, tau, delta, c_hat)

    ndelta = t_star - burn_in + 1
    printf("Estimated c = %12.8f\n", c_hat)

    // ---- Bootstrap ----
    printf("Running %g bootstrap replications...\n", breps)
    rseed(seed)

    real matrix boot_tau, boot_EY1, boot_delta, boot_idx
    real rowvector E1_Y1_b, E1_Y2_b, H1_b, H2_b, H1_imp_b
    real rowvector EY1_imp_b, tau_b, delta_b
    real colvector idx, at_b, D_b, weights_b
    real matrix X_b
    real scalar c_b

    boot_tau = J(breps, T_end, .)
    boot_EY1 = J(breps, T_end, .)
    boot_delta = J(breps, ndelta, .)
    boot_idx = ceil(n * runiform(n, breps))

    for (b = 1; b <= breps; b++) {
        if (mod(b, 500) == 0) {
            printf("  Bootstrap %g / %g\n", b, breps)
        }

        idx = boot_idx[., b]
        at_b = at[idx]
        D_b = D[idx]

        if (has_covs) {
            X_b = X[idx, .]
            weights_b = _dd_compute_weights(at_b, X_b, D_b)
        }
        else {
            weights_b = J(n, 1, 1)
        }

        _dd_weighted_means(at_b, D_b, weights_b, T_end, E1_Y1_b, E1_Y2_b)
        _dd_estimate(E1_Y1_b, E1_Y2_b, burn_in, t_star, T_end, spec,
                     H1_b, H2_b, H1_imp_b, EY1_imp_b, tau_b, delta_b, c_b)

        boot_tau[b, .] = tau_b
        boot_EY1[b, .] = EY1_imp_b
        boot_delta[b, .] = delta_b
    }

    printf("Bootstrap complete.\n")

    // ---- Confidence bands ----
    real rowvector sd_tau, sd_EY1, sd_delta
    real matrix adev_tau, adev_EY1, adev_delta
    real colvector max_adev_tau, max_adev_EY1, max_adev_delta
    real scalar cv_tau_u, cv_EY1_u, cv_delta_u, idx_level, idx_plevel
    real rowvector cv_tau_p, cv_EY1_p
    real scalar pvalue, t
    real colvector col_sorted

    // Standard deviations
    sd_tau = sqrt(diagonal(variance(boot_tau)))'
    sd_EY1 = sqrt(diagonal(variance(boot_EY1)))'
    sd_delta = sqrt(diagonal(variance(boot_delta)))'

    // Guard against zero SDs
    sd_tau = sd_tau + (sd_tau :== 0) :* 1e-100
    sd_EY1 = sd_EY1 + (sd_EY1 :== 0) :* 1e-100
    sd_delta = sd_delta + (sd_delta :== 0) :* 1e-100

    // Standardized absolute deviations
    adev_tau = abs(boot_tau :- tau) :/ sd_tau
    adev_EY1 = abs(boot_EY1 :- EY1_imp) :/ sd_EY1
    adev_delta = abs(boot_delta :- delta) :/ sd_delta

    // Uniform critical values
    max_adev_tau = rowmax(adev_tau[., (t_star+1)..T_end])
    max_adev_EY1 = rowmax(adev_EY1[., (t_star+1)..T_end])
    max_adev_delta = rowmax(adev_delta)

    idx_level = ceil(level * breps)
    idx_plevel = ceil(plevel * breps)

    cv_tau_u = sort(max_adev_tau, 1)[idx_level]
    cv_EY1_u = sort(max_adev_EY1, 1)[idx_level]
    cv_delta_u = sort(max_adev_delta, 1)[idx_plevel]

    // Pointwise critical values
    cv_tau_p = J(1, T_end, .)
    cv_EY1_p = J(1, T_end, .)
    for (t = 1; t <= T_end; t++) {
        col_sorted = sort(adev_tau[., t], 1)
        cv_tau_p[t] = col_sorted[idx_level]
        col_sorted = sort(adev_EY1[., t], 1)
        cv_EY1_p[t] = col_sorted[idx_level]
    }

    // P-value (matches MATLAB exactly)
    pvalue = mean(max_adev_delta :> abs(max(delta :/ sd_delta)))
    printf("Parallel trends test p-value = %6.4f\n", pvalue)

    // ---- Construct CI bounds ----
    real rowvector ci_tau_u_lo, ci_tau_u_hi, ci_tau_p_lo, ci_tau_p_hi
    real rowvector ci_EY1_u_lo, ci_EY1_u_hi, ci_EY1_p_lo, ci_EY1_p_hi
    real rowvector ci_delta_lo, ci_delta_hi

    ci_tau_u_lo = tau - cv_tau_u :* sd_tau
    ci_tau_u_hi = tau + cv_tau_u :* sd_tau
    ci_tau_p_lo = tau - cv_tau_p :* sd_tau
    ci_tau_p_hi = tau + cv_tau_p :* sd_tau

    ci_EY1_u_lo = EY1_imp - cv_EY1_u :* sd_EY1
    ci_EY1_u_hi = EY1_imp + cv_EY1_u :* sd_EY1
    ci_EY1_p_lo = EY1_imp - cv_EY1_p :* sd_EY1
    ci_EY1_p_hi = EY1_imp + cv_EY1_p :* sd_EY1

    ci_delta_lo = delta - cv_delta_u :* sd_delta
    ci_delta_hi = delta + cv_delta_u :* sd_delta

    // ---- Store results in Stata ----
    st_numscalar("r(c_hat)", c_hat)
    st_numscalar("r(pvalue)", pvalue)

    // --- Main results dataset ---
    stata("clear")
    st_addobs(T_end)

    (void) st_addvar("int", "t")
    (void) st_addvar("double", "H1")
    (void) st_addvar("double", "H2")
    (void) st_addvar("double", "H1_imp")
    (void) st_addvar("double", "EY1_imp")
    (void) st_addvar("double", "tau")
    (void) st_addvar("double", "mean_Y1")
    (void) st_addvar("double", "mean_Y2")
    (void) st_addvar("double", "ci_tau_u_lo")
    (void) st_addvar("double", "ci_tau_u_hi")
    (void) st_addvar("double", "ci_tau_p_lo")
    (void) st_addvar("double", "ci_tau_p_hi")
    (void) st_addvar("double", "ci_EY1_u_lo")
    (void) st_addvar("double", "ci_EY1_u_hi")
    (void) st_addvar("double", "ci_EY1_p_lo")
    (void) st_addvar("double", "ci_EY1_p_hi")

    st_store(., "t", (1::T_end))
    st_store(., "H1", H1')
    st_store(., "H2", H2')
    st_store(., "H1_imp", H1_imp')
    st_store(., "EY1_imp", EY1_imp')
    st_store(., "tau", tau')
    st_store(., "mean_Y1", (1 :- E1_Y1)')
    st_store(., "mean_Y2", (1 :- E1_Y2)')
    st_store(., "ci_tau_u_lo", ci_tau_u_lo')
    st_store(., "ci_tau_u_hi", ci_tau_u_hi')
    st_store(., "ci_tau_p_lo", ci_tau_p_lo')
    st_store(., "ci_tau_p_hi", ci_tau_p_hi')
    st_store(., "ci_EY1_u_lo", ci_EY1_u_lo')
    st_store(., "ci_EY1_u_hi", ci_EY1_u_hi')
    st_store(., "ci_EY1_p_lo", ci_EY1_p_lo')
    st_store(., "ci_EY1_p_hi", ci_EY1_p_hi')

    stata("save " + char(34) + prefix + "_results.dta" + char(34) + ", replace")

    // --- Delta dataset ---
    stata("clear")
    st_addobs(ndelta)

    (void) st_addvar("int", "t")
    (void) st_addvar("double", "delta")
    (void) st_addvar("double", "ci_lo")
    (void) st_addvar("double", "ci_hi")

    st_store(., "t", (burn_in::t_star))
    st_store(., "delta", delta')
    st_store(., "ci_lo", ci_delta_lo')
    st_store(., "ci_hi", ci_delta_hi')

    stata("save " + char(34) + prefix + "_delta.dta" + char(34) + ", replace")
}

end
