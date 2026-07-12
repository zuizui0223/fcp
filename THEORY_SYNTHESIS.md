# FCP theory synthesis — recovered core

## Central question

Why are some flowering-plant species polymorphic in floral colour while others are monomorphic?

The framework separates four logically distinct filters:

\[
\text{origin}\rightarrow\text{invasion}\rightarrow\text{persistence}\rightarrow\text{observation}.
\]

A species can appear monomorphic because a second morph never originated, originated but failed to establish, established but was lost, or persists but was not detected.

## Deterministic phase core

For reciprocal invasion exponents \(\gamma_A\) and \(\gamma_B\), define

\[
B_{\rm eff}=\frac{\gamma_A+\gamma_B}{2},\qquad
D_{\rm eff}=\frac{\gamma_A-\gamma_B}{2}.
\]

Then

\[
\gamma_A=B_{\rm eff}+D_{\rm eff},\qquad
\gamma_B=B_{\rm eff}-D_{\rm eff}.
\]

Therefore protected polymorphism exists exactly when

\[
\boxed{B_{\rm eff}>|D_{\rm eff}|}.
\]

This identity is the general phase-boundary compression. Its scientific content comes from deriving \(\gamma_A\) and \(\gamma_B\) from explicit ecological mechanisms.

## Nested mechanistic realizations

### Static, one population

\[
B_{\rm eff}=b,\qquad D_{\rm eff}=d.
\]

### Static, two patches

\[
B_{\rm eff}=b+\sqrt{h^2+m^2}-m,\qquad D_{\rm eff}=d.
\]

Thus heterogeneity promotes reciprocal invasibility while migration attenuates that contribution.

### Static, arbitrary network

The two reciprocal invasion exponents are principal eigenvalues of rare-growth-plus-migration operators. Environmental variance alone is insufficient because placement of local selective environments relative to the dispersal network matters.

### Time-varying environment

The reciprocal invasion exponents are top Lyapunov exponents, or Floquet exponents under periodic forcing. Temporal variance alone is not a universal balancing term; temporal order and covariance structure can matter.

## Finite-population and observation bridge

Deterministic positive invasion does not imply guaranteed observed polymorphism.

For a rare morph with invasion exponent \(r>0\), a weak-selection branching approximation gives establishment probability

\[
P_{\rm est}\approx 1-e^{-2r}.
\]

If new copies arise as a Poisson process with mean \(N\mu T\), then the probability of at least one established origin is approximately

\[
1-\exp\{-N\mu T P_{\rm est}\}.
\]

Conditional on true frequency \(p\), observing both morphs in a random sample of size \(n\) has probability

\[
P(\text{observed polymorphism}\mid p,n)
=1-p^n-(1-p)^n.
\]

Hence observed monomorphism is not equivalent to evolutionary monomorphism.

## Claims that are currently recoverable

1. **Phase-boundary claim**: floral-colour polymorphism can be framed by reciprocal rare-morph invasibility rather than by a list of candidate mechanisms.
2. **Spatial claim**: environmental heterogeneity contributes to polymorphism only after being filtered by connectivity; the relevant object is a spectral invasion quantity, not raw heterogeneity alone.
3. **Temporal claim**: temporal variability is not generically balancing. In the minimal additive model, only the temporal mean matters. Temporal variability requires additional ecological structure to promote polymorphism.
4. **Spatiotemporal claim**: in the general case, polymorphism is governed by reciprocal top Lyapunov/Floquet exponents.
5. **Observation claim**: cross-species comparisons must model ascertainment because sampling effort and morph frequency jointly determine whether true polymorphism is recorded.

## What should be tested by comparative synthesis

The meta-analysis should not attempt to estimate inaccessible theoretical parameters directly. It should test qualitative and semi-mechanistic predictions:

- greater environmental heterogeneity should be associated with polymorphism when connectivity is not too high;
- the heterogeneity effect should depend on dispersal/connectivity proxies;
- temporal variability alone should have weak or inconsistent effects unless paired with buffering/storage proxies;
- species with broader ranges and larger sampled populations may show more recorded polymorphism partly because of origin and detection opportunities;
- reported polymorphism should increase with sampling effort even after ecological predictors are included;
- phylogenetic non-independence must be modeled.

## Main warning

The general identity

\[
B_{\rm eff}>|D_{\rm eff}|
\]

is mathematically exact but, by itself, is not novel enough to carry the paper. Novelty must come from the explicit mechanistic reductions, the distinction among origin/invasion/persistence/observation, and a comparative test designed from those reductions.
