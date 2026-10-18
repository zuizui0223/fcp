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

### Panel B — original-cohort targeted W against structured null

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

## Figure 5 — Species-wide flower-colour variation is geographically partitioned

This remains the main ecological figure. **Panel A now displays all five truly measured original-photo cohorts and their pooled mean**, while Panels B–C explicitly retain the original three-cohort technical/environmental comparison (no inferred falsification for the new U100 photos). Chronology is checked against the genuine completed five-cohort original [Release](https://github.com/zuizui0223/fcp/releases/tag/fcp-all4730-original-method-actual-colour-results-20261017) and Git-tracked source receipt `results/fcp_original_447500_five_cohort_verified_reporting_20261018.json`.

### Panel A — five independent taxon-disjoint photo cohorts and pooled local depletion

Show six horizontal rows: five *actual cohorts* and a visually distinct species-equal **pooled mean** row, explicitly labelled not a sixth independent cohort. Display the original 50-km global-minus-local mean mismatch effect, true >=30 local-pair eligible taxon count and 199-original-vertex-permutation upper p:

- discovery: n=**166**, depletion **0.0205293**, p=**0.005**;
- validation: n=**181**, **0.0186730**, p=**0.005**;
- prospective third cohort: n=**204**, **0.0146849**, p=**0.005**;
- post-freeze wave4: n=**190**, **0.0172936**, p=**0.005**;
- newly measured fifth U100: n=**1,020**, **0.0175987**, p=**0.005**;
- **pooled five-cohort**: n=**1,761**, depletion **0.0176149**, p=**0.005**; **not a sixth cohort**.

The pooled line represents the same original fixed 50-km test, selected *post hoc*, applied to 447,500 actual original-photo terminals from 4,475 taxon identities (with 3,297 >=40 classifiable but only 1,761 having enough 50-km pairs). The null fixes every species' original colour-state composition and coordinates. The 25-, 100- and 250-km comparisons are post hoc sensitivities, not replacements for the unchanged fixed 50-km primary comparison.

### Panel B — falsification tests preserve distributed polymorphism

At 50 km plot mean local-depletion effects for three representations/controls **in the three original study cohorts only (discovery, validation, prospective)**. These tests were not transported to all newly measured U100 taxa:

1. **different-observer pairs only**
   - discovery 0.01998, p = 0.005;
   - validation 0.01690, p = 0.005;
   - third 0.01457, p = 0.005.

2. **nonwhite-only three-state analysis**
   - discovery 0.01947, p = 0.005;
   - validation 0.00729, p = 0.025;
   - third 0.00891, p = 0.040.

3. **continuous nine-colour JSD**
   - discovery 0.01518, p = 0.005;
   - validation 0.01366, p = 0.005;
   - third 0.01301, p = 0.005.

Calendar-quarter-stratified and cross-year controls are added to SI once terminal.

### Panel C — geographic distance dominates a smaller BIO5-associated residual

Compare equal-species mean partial-rank effects in the continuous nine-colour analysis:

| Cohort | IBD-like: colour ~ geography \| BIO5 | IBE-like: colour ~ BIO5 \| geography |
|---|---:|---:|
| Discovery | **0.02760**, p=.005 | **0.00925**, p=.005 |
| Validation | **0.02414**, p=.005 | **0.00826**, p=.010 |
| Third | **0.02225**, p=.005 | **0.01019**, p=.005 |

The panel must visually make clear that IBD-like > IBE-like in all three cohorts.

Mandatory interpretation boundary:

**The BIO5 component is phenotypic and technically bounded: validation flower-minus-background and strict same-observer sensitivities do not support a flower-specific, observer-independent environmental effect. Neither component establishes local adaptation.**

H3 phylogenetic and sampled-span panels move to Supporting Information.

---

## Figure 6 — Shared features of photographed FCP across species (new five-cohort synthesis)

**Reporting-only from genuinely measured, SHA256-pinned 447,500 original photos**. Analysis CI [#39091197295](https://github.com/zuizui0223/fcp/actions/runs/39091197295) separately recomputes frozen four-state 50-km nulls; Figure 6 reads only `results/fcp_original_447500_polymorphism_commonality_20261018.json`, not unverified photo IDs or a new image model.

The figure is an **exploratory FCP commonality synthesis**. It does not replace the original prospective H2 axis confirmation in Figure 4 and does not substitute a predeclared 50-km analysis for the original post hoc scale.

### Panel A — species-wide skew in photographed morph-state composition

Across **3,297** taxa with >=40 classifiable images, plot number of observed coarse photo states (1/2/3/4: **268/1411/1080/538** species), annotate D median **0.172** and photographed modal-state share median **90.6%**. The quantity is not a population-wide biological morph-frequency distribution.

### Panel B — local sorting across different modal hue families

Plot mean global-minus-local four-state 50-km pair diversity for species with majority white (807; 0.018552), yellow/orange (521; 0.010326), red/pink (203; 0.025927), or blue/purple (230; 0.023502); each original 199-null upper p=0.005. These are four exploratory strata from the SAME photograph archive, not four independent prospective tests.

### Panel C — amount of D versus spatial allocation

Four D rank quartiles of the 1,761 true 50-km eligible taxa show mean effects **0.002904/0.011920/0.018917/0.036752**. Emphasize that raw depletion grows with potential amount of within-species variation; it does not demonstrate selection.

### Panel D — source-conditional hue-free and technical sensitivities

Compare the original four-state 0.017615 (n=1761), within-species majority-versus-other 0.016485 (n=1761), equal-genus 0.015251 (n=832 genera), different-observer 0.017070 (n=1745), cross-year 0.015845 (n=1643), and nonwhite-only 0.006595 (n=763). All corresponding pooled upper p=0.005, with an **explicit nonwhite exception**: prospective upper p=0.14 and fourth upper p=0.72. The caption must say 50 km is post hoc, the source remains iNaturalist, and neither the photo-state categories nor local neighbourhoods establish genotype/population/fitness mechanisms.

Canonical output stem: `polymorphism_figure6_fcp_commonality`.

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


### Fig. S9 — bounded mechanism evidence

Three-panel reporting-only figure generated from frozen or source-reported secondary evidence.

**Panel A — BIO5 sorting, observer sensitivity and transport**
- prospective primary cohort: median within-species white-minus-nonwhite BIO5 = +0.0690 SD; Holm-adjusted p = 0.0354; conditional OR = 1.073;
- observer-balanced sensitivity: median +0.0541 SD, Wilcoxon p = 0.0750;
- observer-paired sensitivity: median 0.000 SD, Wilcoxon p = 0.485; conditional OR = 0.787;
- discovery transport: median -0.0068 SD, p = 0.743;
- validation transport: median +0.0662 SD, p = 0.0541;
- visual message: the frozen primary within-cohort association is observer-sensitive and does not transport as a common cross-cohort rule.

**Panel B — published *Silene* molecular anchor**
- structured extraction from Casimiro-Soriguer et al. (2016), not raw-read reanalysis;
- plot all source-reported p<0.05 bud-stage pigmented-versus-white expression contrasts among 29 ABP-related loci;
- dark/white: F3h1 49.0× (p = 0.039), C4h2 36.2× (p = 0.013), Myb1a 5.1× (p = 0.009);
- light/white: F3h1 42.2× (p = 0.049), F3′h 4.5× (p = 0.047);
- F3h1 is the only locus significant in both pigmented-versus-white comparisons;
- visual message: one high-frequency PAL system has a quantified regulatory-pathway anchor near F3h1/Myb1a, but this is not cross-species molecular validation.

**Panel C — biochemically anchored PAL/WAL persistence contrast**
- *S. littorea* phenotype definition is supported by HPLC-DAD-MS^n tissue profiles: PAL lacks anthocyanins in petals but retains them in photosynthetic tissues, whereas WAL lacks anthocyanins in both;
- cross-system entries retain Del Valle et al. Table S1 classifications and are not independently tissue-profiled by the current study;
- *Silene littorea* PAL median positive frequency = 15.5% (8–21%);
- *S. littorea* WAL median = 0.21% (0.05–0.86%);
- cross-system PAL median lower bound = 5%;
- cross-system WAL median numeric upper bound = 0.1%, maximum quantified upper bound = 1.4%;
- visual message: flower-restricted pigment loss can reach substantially higher reported natural frequencies, but the source literature is ascertained and heterogeneous.

Hard boundary: **secondary empirical evidence, not a universal causal explanation of H2.**

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
