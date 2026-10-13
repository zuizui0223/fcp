# Conditional field-colour allocation derivations

The data are observational field quadrat counts and the p values are conditional randomized **label-allocation benchmarks**, not evidence of randomized recruitment or causal interference. Site-year populations and source transect labels are treated as fixed observed strata.

## Fixed M species inventory and co-occurring species exposure

Within each original line h and occupied quadrat j let `n_j=v_j+w_j`, where v/w are **only M species**. The other white-flowered L/S species have counts l/s, and the quadrat exposure is `a_j=(l_j+s_j)/(n_j+l_j+s_j)`.

The line totals are N_h=sum n_j, V_h=sum v_j, W_h=N_h-V_h and p_h=V_h/N_h. We use only lines with >=2 M-occupied quadrats, N_h>=10, 0<V_h<N_h, and differing exposures. D_h=N_h p_h(1-p_h)=V_h W_h/N_h. The within-line score

` t_h= Σ_j(v_j-n_j p_h) a_j / D_h`

is exactly the difference between the observed mean exposure encountered by violet-M individuals and the corresponding mean encountered by white-M individuals in line h. The original source does not provide position per individual beyond quadrat, so this describes **quadrats inhabited by each morph**, not direct plant-to-plant contact.

The siteyear score is `T_site=Σ_h D_h t_h / Σ_h D_h` over source-eligible informative lines, and the study score is the equal-siteyear mean. Site eligibility requires >=2 informative lines, >=30 M individuals, >=5 each morph in those lines. All observational holds are reported, not treated as 0.

Under conditional exchangeability with line-specific V_h and quadrat-specific n_j fixed, the exact sampling distribution per line is multivariate hypergeometric. Every score has null expected 0. We draw 9999 independent line-conditional label allocations per endpoint and report plus-one Monte Carlo one-sided p, using pre-frozen seeds. The two preplanned primary family endpoints (association and residual quadrat clumping) are Holm-adjusted; distance sensitivity remains exploratory.

## Between-quadrat Gini–Simpson within-site and within-line null

At siteyear level with M N and violet count V, p=V/N, the observed **phenotypic** between-quadrat fraction is

`B= Σ_j n_j(p_j-p)^2 / [N*p*(1-p)]`.

Let m_h be the number of occupied quadrats in line h and p_h=V_h/N_h. Then exact between-line fraction:

` B_line = [Σ_h N_h*(p_h-p)^2] / [N*p*(1-p)] `.

At **fixed V_h** and fixed n_j within each h, the conditional mean within-line numerator is

`E[Σ_{j∈h} n_j(p_j-p_h)^2 | V_h,n_j ] = N_h*p_h*(1-p_h)*(m_h-1)/(N_h-1)`

when N_h>1 (zero for one-quadrat/monomorph lines). This is a consequence of multivariate hypergeometric second moments, not fitted to the outcome. Therefore

` E[B | all line totals] = B_line + Σ_h N_h*p_h*(1-p_h)*(m_h-1)/(N_h-1) / [N*p*(1-p)] `.

The residual observed minus exact conditional expectation, Delta_line, tests fine-scale allocation after the line-to-line morph pattern is **held unchanged** in every randomized realization. This is narrower than PR182's whole-site exchangeability null (which also randomizes M colours **between lines**). It cannot tell whether local aggregation is driven by local selection, dispersal/kinship, recruitment, soil, or observer classification errors.

## Linear distance-projection sensitivity

For any line with >=3 M-occupied quadrats and non-constant recorded `distance`, fit `a_j = alpha_h + beta_h*d_j + residual_j` by **n_j-weighted** linear least squares. Use these fixed residual values in the exact same line-conditioned random-label statistic. Across sites, only siteyears with >=2 distance-informative lines, >=30 M, >=5 each morph are counted. This is a diagnostic of local linear distance co-alignment, **not a causal adjustment identifying a mechanism**, and **not** a correction for arbitrary nonlinear environmental clines. Distance in source may have different coverage/meaning in different study designs; 

Report also the unadjusted T on the **identical residualization-eligible lines** to avoid confounding statistical attenuation with a change in sample support. The matched-support attenuation diagnostic was added after reading results, so is explicitly labeled posthoc. The 2026-10-13 outcome p=0.0378 fails a leave-one-siteyear-out (Banpaku) support check and must **not** be promoted as robust evidence that L/S affects M morph even after distance.

## Multiplicity and species ecology

Two planned endpoint p values are tested in a two-member Holm family. A 9999-draw tail probability of 0.0001 is **simulation resolution**, not p literally equal to probability of zero. Siteyear correlations across time and named populations are not independent site evolutionary replicates; latest-year-only sensitivity addresses repeated named population voting but not environmental confounding.

Of all 4492 source-recorded other white species individuals, **L contributes 4463 (99.3544%)** and S only 29. The covariate is thus overwhelmingly *L* exposure in this raw survey, not equal inference about two congeneric species.

Authors already reported broad congener–M-purple association and collinearity with distance in the original publication; the new audit's useful contribution is determining what persists **within the recorded line** and what *does not* survive a local distance/leave-one-siteyear stress test. It cannot prove pollen transfer, mating outcomes or the original paper's selection hypothesis.
