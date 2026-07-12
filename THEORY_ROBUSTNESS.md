# Robustness note: range opportunity and spatial outcome are distinct from persistence

This note hardens the interpretation of the FCP framework around two points.

## 1. Large geographic range is an opportunity axis, not a balancing mechanism

A large range may correlate with floral-colour variation because it can increase:

- the number of environments encountered;
- the number of semi-independent populations;
- total mutation supply and opportunities for repeated origin;
- opportunities for divergent local selection;
- sampling effort and therefore the probability that rare morphs are detected.

But range size does not enter the core invasion condition as a balancing force by itself.
A species can be widespread and monomorphic, or narrowly distributed and locally polymorphic.

The causal structure should therefore be represented as

\[
R \rightarrow (H,\;N_e,\;\mu N,\;S_{\rm effort}),
\]

where `R` is range size, `H` is realized environmental/selective heterogeneity,
`N_e` is effective population size or persistence capacity, `mu N` is mutation
supply, and `S_effort` is observation effort.

Range size is consequently a predictor or opportunity proxy, not a substitute for
`b`, `H_bal`, or another balancing term.

### Comparative implication

The useful empirical test is not merely

\[
P(\text{polymorphism}) \sim R.
\]

Instead, decompose the association:

\[
R \rightarrow H \rightarrow \text{polymorphism},
\]

\[
R \rightarrow \text{mutation supply/persistence} \rightarrow \text{polymorphism},
\]

and

\[
R \rightarrow \text{sampling effort} \rightarrow \text{recorded polymorphism}.
\]

A positive raw range-size effect that disappears after controlling for heterogeneity
and sampling effort is not evidence that range size itself maintains polymorphism.

## 2. Metapopulation persistence is not the same as local coexistence

The rare-invasion results in `model.py`, `two_patch.py`, `multi_patch.py`, and
`spatiotemporal.py` answer the following question:

> Can each morph increase when rare from the boundary state in the modeled
> population or metapopulation system?

Mutual invasibility implies protected persistence at that modeled scale. It does
**not** imply that both morphs are mixed within every local population.

The same metapopulation-level persistence can be realized as:

1. **geographic mosaic** — different patches are dominated by different morphs;
2. **local coexistence** — both morphs occur within the same local populations;
3. **mixed spatial polymorphism** — some patches are locally mixed while others are
   dominated by one morph.

Therefore

\[
\text{persistence} \neq \text{local coexistence}.
\]

The phase theorem should be interpreted as a persistence theorem. Realized spatial
structure requires the interior equilibrium or long-run distribution of patch
frequencies to be examined separately.

`outcome_structure.py` implements this second descriptive layer.

## 3. Recommended hierarchy of empirical outcomes

Do not collapse all floral-colour variation into a single binary variable when the
literature allows finer coding.

Recommended state hierarchy:

- `monomorphic`
- `geographic_mosaic`
- `mixed_spatial_polymorphism`
- `local_coexistence`

A coarser binary outcome may still be used for the primary species-level analysis,
but the finer outcome should be retained for secondary analyses.

## 4. The theoretical architecture after robustness checks

The framework now separates four logically distinct layers:

### Layer A — opportunity

Origin, mutation supply, range size, number of populations, and opportunities to
encounter heterogeneous selection.

### Layer B — invasion and persistence

Rare-invasion exponents and the phase boundary. In the static network model:

\[
 b + H_{\rm bal} > |d + D_{\rm asym}|.
\]

### Layer C — realized spatial organization

Geographic mosaic versus local coexistence versus mixed spatial structure.

### Layer D — observation

Sampling effort and morph frequency determine whether true variation is recorded.

The observed species label is therefore the end of a filtering process:

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

This decomposition prevents three common confusions:

1. treating large range size as a mechanism of balancing selection;
2. treating species-wide variation as evidence of local coexistence;
3. treating recorded monomorphism as evidence that no alternative morph exists.

## 5. What the current theory does and does not claim

### It does claim

- explicit conditions for mutual invasibility in the minimal, two-patch, network,
  and time-dependent linearized formulations;
- that connectivity can attenuate the contribution of spatially heterogeneous
  selection;
- that environmental variance alone is insufficient to characterize spatial or
  temporal effects;
- that metapopulation persistence and local mixing are separate outcomes.

### It does not claim

- that large range size directly causes protected polymorphism;
- that mutual invasibility guarantees within-population coexistence;
- that all documented species-level colour variation is adaptive;
- that a binary literature label perfectly measures the underlying evolutionary
  state.

These limits are part of the model, not caveats to hide later.
