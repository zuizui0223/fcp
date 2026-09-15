# P500 coupling implementation specification

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

## Model qualification still required

Species-conditioned logistic likelihood is the intended estimator, with no
global or species intercept column in a conditional fit. Model-based Wald
intervals must not be described as observer-cluster robust. Before implementing
and using the fit, specify and test finite-MLE existence/separation, optimizer
termination, score and curvature checks, warning handling, and covariance.
Do not use a successful library return value as the only convergence evidence.
The [official-source audit](P500_CONDITIONAL_MODEL_REFERENCE.md) records why.

Required synthetic benchmarks include analytic matched pairs, zero association,
separation, single-state groups, zero predictor variance, invalid inputs and
order invariance. Record all failures. No penalty or alternate model may rescue
an unqualified primary fit.

## Remaining external and scientific limits

The input helper is not a fitted model. It neither verifies authenticity of the
supplied species mapping nor supplies pre-opening chronology. The model cannot
make failed ROI validation disappear. P500 remains closed until the original
chronology and execution requirements are satisfied independently.
