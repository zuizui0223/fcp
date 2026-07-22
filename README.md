# FCP — Floral Colour Polymorphism phase theory

A mathematical and comparative framework for the question:

> Why are some plant species polymorphic in floral colour while others are monomorphic?

The project starts from phase boundaries rather than from a single candidate mechanism, then uses the recovered theory to generate explicit comparative predictions for meta-analysis.

## Current empirical paper

The current submission-oriented analysis uses **intraspecific flower-colour variation** as its umbrella concept. This is broader than strict flower-colour polymorphism and separates two spatial expressions:

- **within-population flower-colour polymorphism**: discrete colour morphs coexist within at least one population;
- **geographically structured flower-colour variation**: colour states are differentiated primarily among populations or regions.

This distinction is essential. Species-wide persistence of multiple colour states does not imply local coexistence, and geographically structured variation should not be labelled simply as FCP without qualification.

The empirical paper asks:

> Do species differ in occupied climatic niche according to whether documented intraspecific flower-colour variation occurs as within-population polymorphism or as geographic differentiation?

The dataset is a global literature-derived comparative sample of documented cases. It is not a random sample of angiosperms and does not estimate worldwide prevalence.

The current result hierarchy is:

1. no clear evidence that documented colour-variable species generally occupy broader climatic niches than matched controls;
2. a negative association between geographically structured colour variation and species-level realised moisture-niche breadth in the baseline-unambiguous subset;
3. weaker evidence after less certain enriched classifications are added;
4. no clear explanation of that association from coarse GBIF fragmentation or generic climatic turnover among unlabelled occurrence clusters.

These results identify a comparative signal, not a demonstrated causal mechanism. Morph-labelled locality data would be required to test environmental sorting, local adaptation, or physiological differentiation directly.

See [`docs/current_data_research_results.md`](docs/current_data_research_results.md) and [`docs/jbi_revision_plan.md`](docs/jbi_revision_plan.md).

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

## Stage 6 — opportunity is not maintenance; persistence is not local coexistence

Two interpretation rules are now explicit.

### Range size is an opportunity axis

Large geographic range can increase environmental heterogeneity, mutation supply, population number, persistence, and sampling effort. But range size is not itself a balancing term in the phase boundary. It should therefore be decomposed through mediating variables rather than interpreted as a direct maintenance mechanism.

### Species-wide persistence is not local coexistence

Mutual invasibility guarantees persistence at the modeled population or metapopulation scale. It does **not** guarantee that both morphs are mixed within each local population.

Realized spatial outcomes are separated into:

- `monomorphic`
- `geographic_mosaic`
- `mixed_spatial_polymorphism`
- `local_coexistence`

This distinction is implemented in [`outcome_structure.py`](outcome_structure.py). See [`THEORY_ROBUSTNESS.md`](THEORY_ROBUSTNESS.md).

The full conceptual chain is therefore

\[
\text{opportunity}
\rightarrow
\text{origin}
\rightarrow
\text{invasion/persistence}
\rightarrow
\text{spatial organization}
\rightarrow
\text{detection}.
\]

## Comparative predictions

The theory-guided meta-analysis should test, rather than freely search for, the following broad patterns:

1. spatial environmental heterogeneity is associated with more within-population floral-colour polymorphism after observation effort is controlled;
2. that association depends on connectivity and dispersal structure;
3. temporal variability alone has no universal positive effect and should matter most when coupled to buffering, storage-like life histories, or spatial structure;
4. range size mixes ecological heterogeneity, origin opportunity, persistence, and detection opportunity and should be decomposed rather than interpreted as a single mechanism;
5. geographically structured colour variation and within-population coexistence should be analysed as distinct outcomes;
6. sampling effort must predict recorded variation and should be modelled explicitly.

## Run

No third-party dependencies are required for the core code.

```bash
python -m unittest discover -s tests -v
python scan_phase_grid.py --n 101 --out results/phase_grid.csv
```

## Repository structure

```text
model.py                              # one-population analytical classifier
two_patch.py                          # exact two-patch spatial extension
multi_patch.py                        # arbitrary-network spectral boundary
temporal.py                           # temporal forcing result
spatiotemporal.py                     # Lyapunov/Floquet numerical framework
finite_population.py                  # origin, establishment, and detection bridge
outcome_structure.py                  # mosaic versus local coexistence decomposition
scan_phase_grid.py                    # phase-grid generator
tests/                                # analytical and numerical consistency tests
THEORY.md                             # one-population derivation
THEORY_2PATCH.md                      # two-patch derivation
THEORY_MULTIPATCH.md                  # arbitrary network derivation
THEORY_TEMPORAL.md                    # temporal/spatiotemporal derivation
THEORY_SYNTHESIS.md                   # recovered claims and limits
THEORY_ROBUSTNESS.md                  # range opportunity and spatial outcome distinctions
META_ANALYSIS_PLAN.md                 # initial comparative design
META_ANALYSIS_PROTOCOL.md             # operational coding and model protocol
```

## Current status

The theoretical scaffold now constrains the empirical synthesis. The immediate priority is not to add mechanisms indefinitely, but to finish the comparative paper, preserve evidence-quality sensitivity, and add only analyses that resolve a specific reviewer-facing weakness.

The main novelty should not be claimed as a new universal theory of polymorphism. The defensible contribution is a flower-colour-specific framework connecting reciprocal invasibility to cross-species evidence while explicitly separating opportunity, origin, persistence, spatial organisation, and observation.