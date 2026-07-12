# FCP — Floral Colour Polymorphism phase theory

A mathematical and comparative framework for the question:

> Why are some plant species polymorphic in floral colour while others are monomorphic?

The project starts from phase boundaries rather than from a single candidate mechanism, then uses the recovered theory to generate explicit comparative predictions for meta-analysis.

## Recovered theoretical core

For reciprocal rare-morph invasion exponents \(\gamma_A\) and \(\gamma_B\), define

\[
B_{\rm eff}=\frac{\gamma_A+\gamma_B}{2},\qquad
D_{\rm eff}=\frac{\gamma_A-\gamma_B}{2}.
\]

Protected polymorphism exists exactly when

\[
\boxed{B_{\rm eff}>|D_{\rm eff}|}.
\]

This identity is the general phase-boundary compression. Its biological content comes from deriving the two invasion exponents from explicit ecological mechanisms.

## Stage 1 — one population

For

\[
\Delta(p)=d+b(1-2p),\qquad
\dot p=p(1-p)\Delta(p),
\]

protected polymorphism exists when

\[
\boxed{b>|d|}.
\]

See [`THEORY.md`](THEORY.md).

## Stage 2 — two patches

For opposite local selection contrasts \(+h\) and \(-h\) coupled by symmetric migration \(m\),

\[
\boxed{b+\sqrt{h^2+m^2}-m>|d|}.
\]

Thus spatial heterogeneity can promote reciprocal invasibility while migration attenuates that contribution.

See [`THEORY_2PATCH.md`](THEORY_2PATCH.md).

## Stage 3 — arbitrary landscape networks

For a weighted migration network with graph Laplacian \(L\), local selection deviations enter the rare-growth operators through principal eigenvalues. The exact phase boundary can be compressed to

\[
\boxed{b+H_{\rm bal}>|d+D_{\rm asym}|}.
\]

Environmental variance alone is therefore insufficient: the placement of selective environments relative to the connectivity network matters.

See [`THEORY_MULTIPATCH.md`](THEORY_MULTIPATCH.md).

## Stage 4 — temporal and spatiotemporal environments

In a well-mixed population with additive temporal forcing \(d(t)\),

\[
\gamma_{A|B}=b+\overline d,\qquad
\gamma_{B|A}=b-\overline d.
\]

Hence temporal variance alone does **not** generate a universal balancing contribution in the minimal model; only the temporal mean appears.

For a general time-varying landscape,

\[
\dot x=A(t)x,
\]

and invasion is governed by the top Lyapunov exponent, or by a Floquet exponent under periodic forcing. Temporal order, autocorrelation, synchrony, and covariance with movement can therefore matter even when marginal means and variances are identical.

See [`THEORY_TEMPORAL.md`](THEORY_TEMPORAL.md).

## Stage 5 — origin, persistence, and observation

The observed cross-species pattern is filtered through

\[
\text{origin}\rightarrow\text{invasion}\rightarrow\text{persistence}\rightarrow\text{observation}.
\]

Observed monomorphism is not automatically evolutionary monomorphism. If a true morph has frequency \(p\), the probability of observing both morphs in a sample of size \(n\) is

\[
1-p^n-(1-p)^n.
\]

This makes sampling effort and ascertainment part of the model rather than nuisance afterthoughts.

See [`finite_population.py`](finite_population.py), [`THEORY_SYNTHESIS.md`](THEORY_SYNTHESIS.md), and [`META_ANALYSIS_PROTOCOL.md`](META_ANALYSIS_PROTOCOL.md).

## Comparative predictions

The theory-guided meta-analysis should test, rather than freely search for, the following broad patterns:

1. spatial environmental heterogeneity is associated with more floral-colour polymorphism after observation effort is controlled;
2. that association depends on connectivity/dispersal structure;
3. temporal variability alone has no universal positive effect and should matter most when coupled to buffering, storage-like life histories, or spatial structure;
4. range size mixes ecological heterogeneity, origin opportunity, and detection opportunity and should be decomposed rather than interpreted as a single mechanism;
5. sampling effort must predict recorded polymorphism and should be modeled explicitly.

## Run

No third-party dependencies are required for the core code.

```bash
python -m unittest discover -s tests -v
python scan_phase_grid.py --n 101 --out results/phase_grid.csv
```

## Repository structure

```text
model.py                              # one-population analytical classifier
 two_patch.py                         # exact two-patch spatial extension
multi_patch.py                        # arbitrary-network spectral boundary
temporal.py                           # temporal forcing result
spatiotemporal.py                     # Lyapunov/Floquet numerical framework
finite_population.py                  # origin, establishment, and detection bridge
scan_phase_grid.py                    # phase-grid generator
tests/                                # analytical and numerical consistency tests
THEORY.md                             # one-population derivation
THEORY_2PATCH.md                      # two-patch derivation
THEORY_MULTIPATCH.md                  # arbitrary network derivation
THEORY_TEMPORAL.md                    # temporal/spatiotemporal derivation
THEORY_SYNTHESIS.md                   # recovered claims and limits
META_ANALYSIS_PLAN.md                 # initial comparative design
META_ANALYSIS_PROTOCOL.md             # operational coding and model protocol
```

## Current status

The broad theoretical scaffold is now recovered far enough to constrain the empirical synthesis. The next research step is not to add mechanisms indefinitely. It is to build the comparative dataset, test the preregistered pattern, and return to theory only where the data expose a specific failure or non-identifiability.

The main novelty should not be claimed as a new universal theory of polymorphism; general spatiotemporal polymorphism theory already exists. The safer contribution is a flower-colour-specific phase framework that connects reciprocal invasibility to cross-species occurrence while explicitly separating origin, persistence, and observation.
