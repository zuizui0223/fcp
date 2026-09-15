# P500 coupling implementation specification

**Chronology correction (2026-09-15):** this implementation was developed on
synthetic data without inspecting P500 responses in this task. However, later
inspection of GitHub Actions showed that another branch had already opened
P500 pixels before this code was authored. It is not a repository-wide
pre-opening freeze. See the [Actions audit](FCP_H123_P500_ACTIONS_CORRECTION_20260915.md).
Historical wording below describes the intended local route, not current
global pixel-access status. Do not retrofit this method as a preregistered
control for the running measurement.

Pre-response implementation work, 2026-09-15. No biological response is opened.
This specification does not replace the original protocol or authorize P500.

## Fixed scientific requirements

The original measurement-control protocol fixes a species-conditioned binary
white/nonwhite response, near-clip fraction standardized within species,
95% odds-ratio intervals per within-species SD, and the existing 0.80–1.25
interpretation tolerance. Primary H2 and sensitivity rules remain unchanged.

## Input convention now implemented

`prepare_coupling_inputs` requires a complete, externally frozen photo-to-species
mapping and one valid row per identity. It rejects missing observations,
duplicates, invalid responses, nonfinite fractions and species relabeling.
Standardization uses the arithmetic mean and sample SD (ddof=1) within each
species. A species is informative only with both response states and nonzero
computed predictor variance. Every species receives an inclusion/exclusion
audit. Excluding non-informative strata from this model does not redefine the
primary H2 denominator or permit silent image attrition.

The sample-SD convention is an explicit implementation choice made on artificial
inputs, not a historical rule attributed to the earlier protocol. No future
change may be justified by obtaining a more favorable OR.

## Conditional estimator and numeric conventions now implemented

Species-conditioned logistic likelihood is the intended estimator, with no
global or species intercept column in a conditional fit. Model-based Wald
intervals must not be described as observer-cluster robust. `fit_coupling`
validates the full cohort through the input helper before any model construction.
It uses statsmodels **0.14.6** conditional likelihood/score/Hessian directly and
SciPy **1.15.3** scalar root solving, not the conditional `fit` wrapper.
Other versions return INDETERMINATE until separately qualified.
Do not use a successful library return value as the only convergence evidence.
The [official-source audit](P500_CONDITIONAL_MODEL_REFERENCE.md) records why.

The following are new pre-response implementation choices, not historical
protocol claims:

1. For each species with m white responses, enumerate the *extremes*, not every
   arrangement: sum its m smallest standardized values and its m largest.
   Sum these over species to obtain L and U; T is the observed sum of white
   standardized values. Require T-L and U-T to exceed
   `1e-10 * max(1, abs(L), abs(U))`. Otherwise return INDETERMINATE.
   This guard includes near-boundary numerical ambiguity; it must not label all
   rejected cases as proven exact separation.
2. Seek a strictly sign-changing score bracket at symmetric bounds
   1, 2, 4, 8, 16, 20 in that order. Any numerical warning/error stops execution;
   there is no numerical rescue. A finite MLE outside these bounds is not a
   scientific negative: it is unqualified under this implementation.
3. Use Brent root solving with absolute tolerance 1e-10, relative tolerance
   1e-12, at most 100 iterations and retained convergence metadata. Require
   convergence and absolute score <=1e-7.
   [SciPy 1.15.3 documented return and tolerance contract](https://docs.scipy.org/doc/scipy-1.15.3/reference/generated/scipy.optimize.brentq.html).
4. Require finite likelihood, coefficient and positive observed information
   greater than 1e-8. Cross-check the package Hessian against a central score
   difference at +/-1e-4; require agreement at relative 1e-3 / absolute 1e-8.
5. Use SE=1/sqrt(information); exponentiate beta +/-1.959963984540054*SE.
   Nonfinite or zero endpoints fail. All warnings, numerical exceptions and
   failed checks yield no OR or interval, and no fitted-success flag.

### Why the scalar support guard applies

This is a project derivation for this **one-predictor** conditional model,
not a general multi-predictor separation detector. Conditional on each species'
white count, the joint likelihood is proportional to exp(beta*T) divided by
the sum of exp(beta*t) over the finite conditional response arrangements.
Its score is T-E_beta[t] and negative second derivative is Var_beta[t].
With nonconstant conditional support, E_beta[t] increases continuously from
L to U as beta goes from negative to positive infinity. Thus a unique finite
root exists for L<T<U; endpoints have only a limiting maximum. Ties do not
change the extrema construction. Floating-point near-endpoints are conservatively
rejected, and the numerical checks remain necessary even for interior T.

### Artificial-data evidence and remaining qualification

Tests compare analytic matched-pair coefficients, information and Wald intervals
under positive, zero and negative association, and an independently enumerated
multi-photo conditional likelihood (all arrangements, with analytic support
moments). They also cover both separated directions, tied extrema, excluded
strata, invalid cohorts, order invariance, runtime drift, injected numerical
warnings and invalid curvature. This is not simulation-based interval coverage
or proof of valid uncertainty under repeated observers. Full-size performance,
observer dependence and an authenticated execution recorder remain unresolved.
No penalty or alternate model may rescue an unqualified primary fit.

## Remaining external and scientific limits

The artificial-data fitter neither verifies authenticity of the
supplied species mapping nor supplies pre-opening chronology. The model cannot
make failed ROI validation disappear. P500 remains closed until the original
chronology and execution requirements are satisfied independently.
