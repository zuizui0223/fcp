# P500 conditional-model implementation reference

Status: pre-response reference only; accessed 2026-09-15. This note neither
authorizes image opening nor supplies missing historical records. No P500
images, responses, or fitted models were accessed for this research.

## Model meaning

`ConditionalLogit` fits grouped binary responses, conditioning away each
group's intercept. Covariates must not include an intercept. If groups are
species, the slope is a within-species conditional association, not a
between-species relationship or a causal effect. Missing-value handling must
be explicit; the default performs no missing-value checking.
[Official model documentation](https://www.statsmodels.org/stable/generated/statsmodels.discrete.conditional_models.ConditionalLogit.html).

## Source-level hazards to qualify before any real response

The inspected stable source identifies itself as statsmodels 0.15.0:

- Groups with constant responses are removed during construction, with a
  warning counting groups and observations. Track these exclusions explicitly;
  they are not missing records and do not establish that an entire cohort is
  informative.
- The Hessian is a numerical derivative of the score.
- `fit` calls the parent optimizer but reconstructs `ConditionalResults` from
  parameters and covariance. It does not copy `mle_retvals` or `mle_settings`.
- Although the signature accepts additional arguments, the inspected call to
  the parent does not forward `kwargs`, `callback`, `fargs`, or `retall`.
  Do not assume that passing a tolerance changes the optimizer, or that a
  returned conditional result exposes convergence diagnostics.
- The denominator uses exponentiation and recursive sums; numerical
  overflow/underflow needs explicit checks.

These are observations of the inspected version, not assertions about every
release. Pin and inspect the actual runtime version; an adapter needs synthetic
qualification before use.
[Conditional-model implementation](https://www.statsmodels.org/stable/_modules/statsmodels/discrete/conditional_models.html).

### Installed-runtime check

Read-only `inspect.getsource(ConditionalLogit.fit)` in the local checkout's
Python environment reports **statsmodels 0.14.6**, not the stable page's 0.15.0.
The installed method likewise neither forwards `kwargs`/`callback`/`fargs`/
`retall` nor copies optimizer diagnostics into the reconstructed result.
Unlike the inspected 0.15.0 source, 0.14.6 unconditionally calls
`rslt.cov_params()` when reconstructing results, even if `skip_hessian=True`
was passed. Do not use that option as a presumed diagnostic-only solution.
This was source inspection only, not a model fit or a runtime behavior test.
The inspected interpreter was `C:\Program Files\Python310\python.exe`;
the module was under the user's Python310 site-packages. SHA-256 of the
UTF-8 `inspect.getsource(ConditionalLogit.fit)` string was
`77f871b3770e5dc8056d9fecc43a74fdee2137a5ddb1eeaf33d959641c2f38d9`.

## Uncertainty and convergence

In the inspected parent fit, the non-Newton path obtains covariance by
inverting the negative Hessian after checking finite values and positive
eigenvalues. Failure emits a Hessian-inversion warning. Parent results carry
optimizer diagnostics and may emit a convergence warning; that does not imply
the conditional wrapper retains those diagnostics. Standard errors use the
square roots of covariance diagonal entries. Confidence limits are parameter
plus/minus a reference quantile times its standard error: normal when `use_t`
is false, Student's t otherwise. For a declared predictor increment, transform
both coefficient limits to odds-ratio limits using that same increment.
These are model-based Wald intervals, not an established observer-cluster
robust procedure. Species intercept conditioning alone does not validate
dependence assumptions for repeated observers.
[Parent likelihood and result implementation](https://www.statsmodels.org/stable/_modules/statsmodels/base/model.html).

## Conservative failure handling: proposed qualification requirements

Statsmodels warns that apparently converged fits may still have infinite or
nonunique estimates, and that binary models can suffer separation or
quasi-separation. Its general pitfalls discussion is not a guarantee that
`ConditionalLogit` detects these situations. Insufficient variation and
ill-conditioning also require attention.
[Official pitfalls](https://www.statsmodels.org/stable/pitfalls.html).

The following are project-facing recommendations inferred from those risks,
not a documented complete separation detector:

1. Validate finite, typed inputs and the frozen join first. Record both total
   and response-informative groups and rows; never silently replace excluded
   groups or redefine retention denominators after inspecting outcomes.
2. Check within-group predictor variation and design rank. Rank sufficiency
   alone does not prove existence of a finite maximum-likelihood estimate.
3. Require independently retained optimizer termination information, finite
   objective/score/parameters, and a finite positive-definite information
   matrix under predeclared numerical tolerances.
4. Qualify known separation, quasi-separation, constant-response, singular,
   and extreme-numerical synthetic cases. Warnings or finite coefficient
   values alone cannot establish estimability. If detection remains
   unqualified, report the real-data fit as not evaluable rather than clear.
5. Do not introduce penalties, drop extra groups, change predictors, or retune
   optimizer settings in response to observed support. Any unresolved
   estimator or covariance choice must be settled before responses are read.

## Access and evidence limits

The four linked official documentation/source pages and the installed 0.14.6
`ConditionalLogit.fit` method were inspected. The `stable` URLs are mutable.
The installed check is limited to that method, not a complete package audit.
No executable adapter,
separation proof, simulation-based coverage assessment, cluster-robust
qualification, or real-data result is supplied by this note. It is a source
reference for the pending method specification, not a completed validation.
