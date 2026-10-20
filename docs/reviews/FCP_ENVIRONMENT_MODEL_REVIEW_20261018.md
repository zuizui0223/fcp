# Independent method review: FCP within-species environmental slopes

Date: 2026-10-18. Reviewed the existing design document, climate sampler and source metadata, the implementation `scripts/analysis/run_fcp_447500_environment_model_20261018.py`, and outputs under `results/fcp_real_447500_environment_model_20261018/`. This review is scoped to the environmental analysis; it does not certify unexecuted molecular experiments, taxonomy joins, or external prediction performance.

## Assessment

The implementation is a defensible exploratory two-stage analysis of conditional environmental associations in the selected photographed species. It fits species-specific Firth bias-reduced binomial logistic regressions, then axiswise normal-normal REML meta-analysis with modified HKSJ mean intervals. It is accurately labeled an approximate two-stage model, not a fitted one-stage binomial GLMM. The mathematics and output consistency checks below passed. No blocking implementation defect remains in the reviewed version.

The supported finding is a modest positive average BIO5 association with photographed white status, conditional on BIO14 and SRADmean. The analysis does not establish that most species have positive true slopes, a universal thermal law, local adaptation, pigment cost savings, or molecular tradeoffs. Lack of resolved heterogeneity is particularly important: BIO5 and BIO14 REML heterogeneity estimates lie on zero, and all three primary heterogeneity profile intervals include zero. That is weak resolution of heterogeneity under the working sampling model, not proof of identical biology.

## Target and selection

The primary panel is exactly 1,574 species and 86,114 photos, selected for seven-climate-variable complete cases, at least 40 classified photos and at least five white and five chromatic photos. Independent reproduction of the three-predictor panel gave the same counts. All 86,114 photos have distinct observation IDs; observer IDs are complete. Thus duplicated photos from one observation do not explain primary pseudoreplication, although shared observers and spatial sampling still can.

Selection conditions on observed colour variation and classification success. Results describe this selected panel, not all 4,475 measured species or natural plant populations. The seven-variable complete-case rule preserves the earlier panel even though three variables enter the regression.

The primary arithmetic mean gives equal weight to fitted species. Its finite-panel standard error is sqrt(sum of approximate sampling variances)/K under independent species and independent photographs within each fit. The REML hypermean is precision weighted and is a different estimand. Empirical-Bayes estimates borrow strength under a working normal distribution and depend on estimated hyperparameters. These three outputs are distinguished correctly.

Each predictor is centered within species, so the fitted slopes concern within-species climate variation. The scale is one primary-panel SD for that species, not a common physical environmental increment. Thus standardized slope comparisons cannot identify the largest causal environmental effect. Slopes per native unit are saved for other clearly specified targets.

## Model and numerical checks

The Firth objective is the ordinary binomial log likelihood plus one half the log determinant of its information matrix. Its score correction uses the appropriate leverage term. Inverse ordinary Fisher information at the Firth estimate is explicitly labeled an approximate sampling covariance; it is not an exact penalized-profile confidence interval.

Two independent analytical checks were run:

1. For a saturated two-group binomial logistic model, Firth estimation has a closed-form solution obtained by adding one half to successes and failures in each group. The implementation agreed within 4.56e-9 in coefficient magnitude. Its analytic gradient agreed with numerical differentiation within 2.32e-9.
2. When every first-stage variance is the same known value v, REML has tau-squared = max(0, sample variance of effects minus v), with the arithmetic mean as pooled mean. The implementation agreed in tau-squared within 1.67e-9 and in mean exactly to displayed precision; the modified HKSJ standard error agreed within 3.80e-10 for the chosen heterogeneous example.

The modified HKSJ variance uses max(1,q), avoiding spurious narrowing when q falls below one. Profile heterogeneity intervals use an approximate chi-square cutoff with an explicit boundary caveat. Marginal meta-analysis uses each joint slope's sampling variance; cross-axis covariance is saved but between-axis covariance is not modeled. Consequently, no joint multivariate inference or posterior ranking of the three axes follows from these summaries.

Independent primary design audit found condition numbers median 2.60, 95th percentile 5.48, maximum 12.74. Twenty-three species had a predictor pair with absolute within-species correlation exceeding .95. All primary fits converged. All saved 3-by-3 slope covariance blocks were positive definite; the minimum eigenvalue was .02210 and the largest symmetry/diagonal consistency discrepancy was 1.14e-12.

The reviewed refreshed execution receipt matched the current analysis script. The corrected outputs have no reported conditional sign probabilities when tau-squared is zero; these formerly degenerate 0/1 values are now undefined. This was an important reporting correction.

## Findings and sensitivity interpretation

Primary BIO5 equal-species mean is 0.07756 log odds per within-species SD, with independence-based interval [0.04809, 0.10704]. The precision-weighted REML mean is 0.05467 [0.03149, 0.07785]. BIO14 and SRADmean primary mean intervals include zero. All five cohorts have positive BIO5 point estimates, but the equal-species discovery and validation intervals include zero; this is not five individually decisive replications.

Primary tau-squared estimates are 0 for BIO5, 0 for BIO14, and .004902 for SRADmean. Upper profile limits for the corresponding tau values are approximately .0736, .0573, and .1447. The apparent precision of a plug-in new-species prediction interval at tau=0 is conditional on the estimated heterogeneity and can be misleading: it approximately becomes a mean interval. Do not display such an interval as independently validated generality. No held-out species, genus, or location prediction test was executed.

For SRADmean, the normal latent-positive fraction is about .439, but remains conditional on the normal exchangeability model and estimated heterogeneity. The observed positive coefficient fractions are .551, .493, and .479 for BIO5, BIO14, and SRADmean. Neither observed signs nor conditional empirical-Bayes probabilities are direct biological responder fractions. For BIO5 and BIO14, the boundary model does not support a useful direction-prevalence conclusion.

Sensitivity selection is material:

| Fit | Successfully fitted species | Photos used |
|---|---:|---:|
| Primary | 1,574 | 86,114 |
| One photo per species-observer | 1,486 | 79,063 |
| One photo per species and 1-degree latitude-longitude cell | 633 | 32,127 |
| Added linear latitude and circular-unwrapped longitude | 1,134 | 62,216 |

Thinning is response-blind via a fixed SHA256 ordering. Primary SDs remain fixed while predictors are recentered for each fitted subset. Eligibility is reapplied and matched full-data baselines are supplied, which is essential because restricted panels differ. Observer thinning leaves 433 retained species unchanged; spatial thinning leaves 210 unchanged. A one-degree cell does not enforce a minimum distance, has varying area by latitude, and does not remove spatial confounding.

The geographic fit rejects 440 species whose longitudes do not meet its predefined range/geometry rules. It is a limited linear geographic sensitivity rather than a spatial process adjustment. High leverage is more prominent there: 21 geographic fits have maximum leverage above .99, compared with one primary fit; 164 geographic fits exceed .8. Geographic SRADmean arithmetic mean is particularly sensitive and uncertain. This is a reason to keep this fit a labeled sensitivity, not select it because it produces a preferred result.

The BIO5 genus-cluster bootstrap interval [.043999, .109067] is positive; other primary axes cross zero. The bootstrap resamples 727 named genus clusters and addresses stability under one taxonomy-based grouping. It is not a full phylogenetic correction or a fixed-panel sampling confidence interval. Shared observers and other dependencies can cross genus boundaries.

## Remaining inferential limits

Photograph-level information matrices assume independent outcomes. Thinning reduces repeated observation patterns but does not estimate residual spatial or observer covariance. Date, season, cultivation, exposure, near clipping, background, and other image-quality adjustments are absent. The absence of an all-cohort near-clip table is reported. A photographed white label is not direct pigment chemistry or a confirmed heritable morph.

Firth correction stabilizes separated fits but does not guarantee accurate normal sampling approximations with small counts or influential environmental observations. First-stage estimates and their estimated variances may be associated. REML normality and independent/exchangeable species are working approximations, not established biological facts. Selection for both response states further limits population prevalence inference.

Long-term climate normals describe location, not weather during development or at photography. BIO5 is maximum temperature of the warmest month, BIO14 precipitation of the driest month, and annual mean solar radiation is not direct UV exposure. All causal, genetic, ecological fitness, and molecular interpretations require independent evidence.

The analysis is therefore suitable for the statement that the selected photo panel shows a small average conditional BIO5-white association that persists in the specified limited sensitivities. It supplies real within-species slopes and uncertainty, while leaving direction prevalence, geographic prediction, and molecular mechanism unresolved.
