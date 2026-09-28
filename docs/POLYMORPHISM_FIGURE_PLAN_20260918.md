# Flower-colour polymorphism canonical figure plan — 2026-09-18

This document describes the **implemented main figures** after prospective-confirmation H2 and visual QA. It supersedes the earlier aspirational multi-panel architecture where they differ.

Authoritative claim ledger:

- `docs/POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260918.md`

Authoritative prospective H2 result:

- `results/polymorphism_h2_third_cohort_prospective_white_axis_20260917/result.json`

Canonical generator:

- `scripts/analysis/make_polymorphism_manuscript_figures.py`

Canonical outputs:

- `docs/figures/polymorphism_20260918/`

## Figure 1 — Data provenance and inferential necessity

### Panel A — descriptive D distributions

Show the high-depth **discovery** and **species-disjoint validation** D distributions:

- discovery n = 369;
- validation n = 363;
- dashed cohort medians;
- `D = 1 - sum_k p_k^2`.

Mandatory annotation:

**High-depth cohorts — not a prevalence sample.**

No percentage of species is labelled polymorphic.

### Panel B — inferential sequence

Show a single top-to-bottom inferential sequence. The boxes are **not nested samples**; they explain why each stage exists.

1. **Stage 0 — outcome-blind opportunity frame**
   - iNaturalist metadata: 42,111 metadata-discovered species;
   - 4,730 species capable of >=100 retained photographs;
   - necessity: define candidate species before flower colour is examined.
2. **Stage 1 — discovery + validation**
   - iNaturalist images: 500 discovery × 100 photos;
   - 500 validation × 100 photos;
   - 369 + 363 D-eligible species;
   - necessity: estimate within-species distributions deeply and hold out species for validation.
3. **Target freeze after discovery**
   - fixed white-versus-nonwhite axis, W statistic and structured null;
   - necessity: a data-derived target cannot be called confirmatory in the data that generated it.
4. **Stage 2 — prospective confirmation**
   - 3,230 unused eligible candidates after prior exclusions;
   - 500 species frozen outcome-blind;
   - 499 species × 100 newly sampled photographs;
   - 377 measurement-evaluable; 158 primary-H2 species;
   - necessity: untouched species/photo test of the frozen target.
5. **Stage 3 — post-H2 annotations**
   - image-highlight reacquisition, WorldClim 2.1 climate and V.PhyloMaker2 phylogeny;
   - necessity: bound image-formation effects and simple alternative explanations;
   - dashed styling because these analyses do not create or upgrade the prospective H2 confirmation.

Solid arrows encode the inferential sequence through prospective H2. The post-H2 annotation box is dashed.

Mandatory top-right note:

**Inferential sequence — cohorts are not nested samples.**

---

## Figure 2 — H1 observer-disjoint reproducibility

### Panel A — first-frozen 200-partition result

Show q05–median–q95 of split-half Spearman rho for discovery and validation.

Validation annotations:

- median rho = **0.7891**;
- q05 = **0.7652**;
- q95 = 0.8109;
- primary median floor = 2/3.

The discovery result is displayed for concordance; validation is the decision cohort.

### Panel B — later deterministic stress test

Show deterministic split rho with bootstrap 95% interval for discovery and validation.

Validation:

- rho = **0.7927**;
- 95% CI = **0.7419–0.8324**;
- CCC = **0.8474**;
- strict floor = 0.80.

Annotations must be vertically offset from CI lines and the floor label must not overlap the validation text.

Visual message:

**primary H1 supported; later stricter stress test constrains the claim but does not overwrite the earlier prospective decision.**

---

## Figure 3 — H2 target localization in the original cohorts

This is explicitly a discovery/audit figure. It must never make q_white look prospectively chosen in the original cohorts.

### Panel A — fixed white contrast

Palette order:

`white, yellow, orange, red, pink, magenta, purple, blue, bronze`.

Display

`q_white = normalize([1,-1/8,-1/8,-1/8,-1/8,-1/8,-1/8,-1/8,-1/8])`.

Mandatory annotation:

**Named axis isolated after original broad H2 was opened.**

### Panel B — legacy targeted W against structured null

Show observed W as diamonds, structured-null medians as points, and structured-null 95% intervals as horizontal bars.

Primary 0.10:

- discovery N = 152, W = **0.514625**, p = **0.001**;
- validation N = 129, W = **0.514586**, p = **0.001**.

Strict 0.20:

- discovery N = 75, W = **0.542355**, p = **0.001**;
- validation N = 65, W = **0.510517**, p = **0.008**.

The encoding legend is placed **outside below the right axis** so that it cannot cover the strict-validation point or its annotation.

Projection-removal and non-white-only falsification are reported in Supporting Information rather than forced into this panel.

---

## Figure 4 — Prospective confirmation of the pre-frozen achromatic–chromatic axis

This is the decisive biological figure.

The figure caption must state:

**Species-disjoint prospective test of excess alignment within the same iNaturalist opportunity universe; not an independent-source replication.**

This boundary is intentionally kept out of the plotting field so the decisive null-versus-observed comparison remains legible at journal width.

### Panel A — primary 0.10 tier

Histogram of 999 frozen structured-null W values.

- N = **158**;
- W_obs = **0.5172457461**;
- null median = **0.4571428150**;
- null 95% interval = **0.4358491120–0.4752987776**;
- p = **0.001**.

### Panel B — strict 0.20 tier

- N = **86**;
- W_obs = **0.5329282123**;
- null median = **0.4593196659**;
- null 95% interval = **0.4328679570–0.4867224043**;
- p = **0.001**.

The measurement gate is described in the legend/text rather than as an additional panel:

- 499 species;
- 49,900 terminal rows;
- 25,788 classifiable;
- 377 measurement-evaluable;
- support minimum 250;
- 256/256 partitions;
- 0 replacements;
- pixels not persisted.

---

## Figure 5 — Spatial organization accompanies D while simple explanations fail

### Panel A — D tracks within-species geographic organization

This is the positive biological result and receives the widest panel in the three-panel layout.

Plot the frozen observed partial correlations against their geometry-preserving null intervals:

- discovery, span + clear technical-failure adjusted: partial rho = **0.1266367**, p = **0.007**;
- validation, span + clear technical-failure adjusted: partial rho = **0.0992877**, p = **0.025**;
- validation matched flower-minus-background: partial rho = **0.1162411**, p = **0.010**.

The reporting source is `results/polymorphism_spatial_organization_clue_20260918/result.json`, which copies previously frozen Step-8/Step-9 values without running a new biological analysis.

### Panel B — no detectable broad tree-wide conservation

Validation-cohort Blomberg K:

The three tree-placement scenarios are displayed as **unconnected discrete points**; no line joins S1–S3 because they are alternative frozen scenarios, not an ordered trajectory.

- S1: K = **0.0710190**, p = **0.2716**;
- S2: K = **0.0601476**, p = **0.4134**;
- S3: K = **0.0707577**, p = **0.2674**.

Annotate that 0/3 raw-D scenarios had p < 0.05; label this as a bounded non-support result, not an equivalence test.

### Panel C — discovery span effect collapses in validation

- discovery: rho = **0.1798786**, p = **0.00089996**;
- validation: rho = **-0.0025855**, p = **0.9586021**.

Mandatory boundary:

**The D–spatial association is structural rather than causal; H3a is not an equivalence test; sampled span is not true biological range size.**

---

## Supporting figures

### Fig. S1
Full H1 validation partition diagnostics.

### Fig. S2
H1 discovery concordance and deterministic stress-test details.

### Fig. S3
Broad pre-target H2 geometry: leading concentration, frozen-axis validation transport and discovery–validation axis alignment.

### Fig. S4
Construction-preserving null audit.

### Fig. S5
q_white projection-removal and non-white-only H2 diagnostics.

### Fig. S6
Prospective-confirmation chain of custody: selection, metadata freeze, 256 partitions, support gate and H2_COMPLETE.


### Fig. S7
H3a full sensitivity panel.

### Fig. S8
H3b full sensitivity panel.

---

## Visual QA contract

Before a figure package is marked submission-ready:

1. inspect the actual generated PNG/PDF, not only the plotting code;
2. no annotation, legend or title may obscure a data mark or another label;
3. cohort-flow arrows must encode the real chronology/topology;
4. main-text figure legends must match the implemented panel count;
5. all frozen numerical values must remain unchanged;
6. both PNG and PDF must be regenerated from the same frozen generator;
7. the manifest must retain `scientific_claims_changed = false`;
8. long inferential boundary statements belong in the caption when placing them inside the plotting field reduces legibility;
9. alternative tree-placement scenarios must not be connected by a line that implies an ordered continuum.

The scientific narrative is:

`Fig.1 data provenance + phenotype definition -> Fig.2 measurement validity -> Fig.3 target localization -> Fig.4 prospective confirmation -> Fig.5 positive spatial clue + bounded alternative explanations`.
