# Flower-colour polymorphism paper architecture — 2026-09-18

This document supersedes `docs/POLYMORPHISM_PAPER_ARCHITECTURE_20260913.md` where the evidential status changed after completion of the prospective third-cohort H2 test.

Authoritative current claim ledger:

- `docs/POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260918.md`

Working manuscript:

- `docs/POLYMORPHISM_MANUSCRIPT.md`

Current figure plan:

- `docs/POLYMORPHISM_FIGURE_PLAN_20260918.md`

## 1. One-sentence paper

**Species-level flower-colour polymorphism is reproducibly measurable under high-depth observer-disjoint photographic sampling, and a recurrent white-versus-nonwhite colour-space axis identified in the original cohorts was prospectively confirmed in a pre-frozen species-disjoint third cohort from the same iNaturalist opportunity universe.**

The H3 results are retained as bounded alternative-explanation tests: broad tree-wide conservation is not detected in reserve, and the discovery sampled-span association collapses under species-disjoint replication.

## 2. What changed after the 2026-09-13 architecture

The former main limitation was that the named white-versus-nonwhite H2 axis had been isolated only after the original discovery/reserve H2 geometry was opened.

That limitation has now been directly addressed.

A new outcome-blind third cohort was selected and frozen before biological opening. The already fixed:

- white/nonwhite contrast `q_white`;
- statistic `W`;
- 0.10 primary tier;
- 0.20 strict tier;
- construction-preserving structured null;
- support gate;
- one-shot/no-rerun execution contract

were then applied without refitting.

The resulting machine-readable terminal verdict is:

`H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED`

Therefore the paper can now distinguish three evidential stages:

1. **measurement validation** — H1;
2. **target discovery/localization** — legacy H2 discovery/reserve;
3. **untouched prospective confirmation of the already frozen target** — third cohort.

This chronology is the central paper architecture.

## 3. Central conceptual move

Most macroecological trait studies reduce each species to one mean or one category. This paper instead treats **within-species trait diversity itself as a continuous species phenotype** and separates three logically distinct questions:

1. **Measurement:** can a species-level polymorphism phenotype be recovered reproducibly from independent observer sets?
2. **Geometry:** is the within-species variation directionless, or does it repeatedly align with a common colour-space direction?
3. **Explanation:** do broad phylogenetic or geographic-opportunity predictors account for species differences in that phenotype?

The evidence is asymmetric:

- H1: yes, under the first-frozen observer-disjoint validation rule, with a stricter split-invariance caveat;
- H2: yes, for a narrow achromatic–chromatic white-versus-nonwhite axis, now prospectively confirmed;
- H3a: no support for broad tree-wide phylogenetic signal under the frozen reserve test;
- H3b: the discovery sampled-span association does not replicate.

That asymmetry is the result, not a gap to be rescued.

## 4. Working title

Primary:

**A recurrent achromatic–chromatic axis structures within-species flower-colour polymorphism across plant species**

Safer alternative:

**Within-species flower-colour polymorphism is reproducible and repeatedly aligned with an achromatic–chromatic axis**

The primary title is now allowed because the already frozen axis passed the untouched prospective third-cohort test.

Do not use language implying:

- global prevalence;
- independent-source replication;
- pigment loss/gain;
- adaptive mechanism;
- evolutionary direction.

## 5. Abstract logic

### Sentence 1 — problem

Flower colour is usually summarized by a species mean or dominant state, obscuring within-species diversity.

### Sentence 2 — measurement question

Define the continuous four-state polymorphism score `D = 1 - sum p_k^2` and ask whether independent observer sets recover stable species differences.

### Sentence 3 — H1

Reserve repeated observer-disjoint partitions give median rho = 0.789, q05 = 0.765 and median CCC = 0.855; a later strict deterministic split gives rho = 0.793 and misses an 0.80 floor.

### Sentence 4 — H2 chronology

Original discovery/reserve geometry localizes the construction-controlled recurrent component to a white-versus-nonwhite axis, but that named axis was post-audit in those cohorts.

### Sentence 5 — prospective test

A pre-frozen third cohort completes 49,900 terminal rows from 499 species; 377 species pass the support gate.

### Sentence 6 — prospective result

Primary 0.10: 158 species, W = 0.51725, p = 0.001.
Strict 0.20: 86 species, W = 0.53293, p = 0.001.

### Sentence 7 — alternative-explanation filters

Fresh reserve tests show no detectable broad tree-wide conservation of D, while the discovery sampled-span association collapses to essentially zero; broad ancestry and sampled extent are therefore insufficient as simple explanations under the frozen designs.

### Sentence 8 — conclusion

Within-species flower-colour diversity is measurable as a species phenotype, and its strongest confirmed recurrent geometry in these data is achromatic–chromatic rather than a general hue axis.

## 6. Introduction spine

### P1. The missing phenotype

Macroecology commonly collapses within-species variation into a single species state.

The key gap is not only where colours occur, but whether **the amount and geometry of within-species diversity** can themselves be treated as macroecological traits.

### P2. Measurement before explanation

Introduce D as a continuous four-state diversity score.

State H1: independent observer sets should preserve species rankings if D is a measurable phenotype.

### P3. Geometry of polymorphism

Two species can have similar D but differ in which colours define their within-species modes.

State broad H2: variation may be geometrically constrained rather than arbitrary.

Explain that any shared axis must be tested against a null preserving coarse-state composition and the global state-to-palette construction.

### P4. Discovery versus confirmation

Explain chronology explicitly:

- original cohorts reveal and audit the geometry;
- the white/nonwhite target is then frozen;
- a new species-disjoint cohort tests that target prospectively without refitting.

This paragraph is now essential because it turns the former weakness into the strongest design feature.

### P5. Alternative-explanation filters

Introduce H3a and H3b as fresh-data tests of two simple explanations for species differences in D: broad ancestry and sampled geographic opportunity. Preserve that H3a is not an equivalence test and H3b uses sampled span, not true range size.

## 7. Methods order

### 7.1 Global frame and cohort roles

Start with the 42,111-species opportunity frame.

Then distinguish:

- discovery;
- reserve;
- P500 measurement-transport execution;
- third-cohort prospective H2 confirmation.

Do not present all four as equivalent biological replications.

### 7.2 Four biological states and D

Define:

- white;
- yellow/orange;
- red/pink;
- blue/purple.

`mixed_uncertain` is not a fifth biological state.

Define D and the >=40 classifiable gate.

### 7.3 H1 observer-disjoint validation

Describe the first-frozen 200-partition protocol before the later stricter deterministic stress test.

Chronology determines inferential role.

### 7.4 Continuous palette geometry

Define:

- nine-colour normalized palette;
- Hellinger transform;
- deterministic unlabeled two-means;
- species displacement vector;
- primary 0.10 tier;
- strict 0.20 tier.

### 7.5 Construction-preserving structured null

The structured null is main-text material, not a supplement-only technicality.

It preserves species × coarse-state row counts while permuting continuous nine-colour rows within coarse morph across selected species.

### 7.6 White-axis target discovery and freeze

Define

`q_white = normalize([1,-1/8,-1/8,-1/8,-1/8,-1/8,-1/8,-1/8,-1/8])`

and

`W = mean_i (u_i dot q_white)^2`.

State clearly that the legacy cohorts are post-audit for this named axis.

### 7.7 Third-cohort prospective confirmation

Describe:

- 3,230 outcome-blind candidate species;
- deterministic selection of 500;
- fresh metadata;
- 499 authorized species ×100 rows;
- no replacement;
- technical terminal qualification before biological opening;
- 256 terminal partitions;
- support gate before H2;
- one-shot/no-rerun contract;
- durable serialization/read-back validation.

### 7.8 H3a / H3b bounded tests

Put after H2.

They delimit broad explanation; they are not allowed to rescue or redefine the main positive result.

## 8. Results order

### Result 1 — D is reproducibly measurable

Primary reserve repeated partitions:

- median rho = 0.7891;
- q05 = 0.7652;
- median CCC = 0.8548;
- projected reliability = 0.8821.

Stress-test boundary:

- deterministic rho = 0.7927;
- strict floor = 0.80;
- strict test not supported.

Interpretation: reproducible, not near-perfect or split-invariant.

### Result 2 — original cohorts localize the recurrent geometry

Primary 0.10:

- discovery N=152, W=0.514625, p=0.001;
- reserve N=129, W=0.514586, p=0.001.

Strict 0.20:

- discovery N=75, W=0.542355, p=0.001;
- reserve N=65, W=0.510517, p=0.008.

Projection of q_white removes the excess structured-null signal.

Interpretation: the recurrent component is narrow; there is no supported general non-white hue axis.

### Result 3 — third-cohort support gate passes

- 499 species;
- 49,900 terminal rows;
- 25,788 classifiable;
- 377 measurement-evaluable species;
- required minimum = 250;
- 256/256 partitions;
- replacement = 0;
- persisted pixels = false.

This result must appear before the prospective W values.

### Result 4 — frozen white axis is prospectively confirmed

Primary 0.10:

- N=158;
- W=0.5172457461;
- null median=0.4571428150;
- p=0.001.

Strict 0.20:

- N=86;
- W=0.5329282123;
- null median=0.4593196659;
- p=0.001.

Verdict:

`H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED`.

### Result 5 — broad predictors do not replicate

H3a reserve K:

- S1 p=0.2716;
- S2 p=0.4134;
- S3 p=0.2674.

H3b:

- discovery rho=0.1798786, p=0.00089996;
- reserve rho=-0.0025855, p=0.9586021.

## 9. Figure architecture

Use `docs/POLYMORPHISM_FIGURE_PLAN_20260918.md`.

Main sequence:

1. global frame → D phenotype;
2. observer-disjoint H1;
3. H2 target discovery/localization;
4. prospective third-cohort confirmation;
5. bounded H3 alternative-explanation tests.

The key visual separation is Figure 3 versus Figure 4:

**discovery/audit → frozen target → prospective confirmation**

Do not collapse those stages into one figure that makes the named axis look prospectively selected in the legacy cohorts.

## 10. Discussion order

### D1. Establish the phenotype

High-depth citizen-science photography can recover stable species-level differences in within-species flower-colour diversity.

### D2. Establish the recurrent geometry

The strongest recurrent geometry is achromatic–chromatic rather than a shared non-white hue direction.

### D3. Explain why the prospective test matters

The legacy cohorts identified the target; the third cohort tests it without target selection, threshold changes or axis refitting.

This is the manuscript's strongest inferential upgrade.

### D4. Keep mechanism open

Do not infer:

- pigment pathway;
- transition direction;
- pollinator selection;
- climate selection;
- developmental mechanism.

### D5. Use H3 as alternative-explanation filters

Present the reserve results as showing that broad ancestry is not detectably conserved at the tested scale and that the discovery sampled-span effect collapses out of sample. These results make two simple explanations insufficient without claiming zero phylogenetic or geographic effects.

### D6. State the representativeness boundary

The third cohort is species-disjoint but not source-independent.

The 42,111-species frame is not a prevalence sample.

## 11. Main text versus supplement

### Main text

Must retain:

- 42,111 frame / high-depth distinction;
- D definition;
- H1 chronology and result;
- structured null;
- white-axis discovery chronology;
- third-cohort pre-opening design;
- support gate;
- prospective W result;
- concise H3 negatives.

### Supplement/repository

Move detailed receipts to supplement:

- all H1 partition metrics;
- deterministic split implementation;
- full legacy H2 eigengeometry;
- full structured-null diagnostic tables;
- q_white residual/non-white sensitivity;
- selection hashes;
- third-cohort metadata receipt;
- 256 partition receipts;
- P500 terminal postmortem;
- all H3a S1-S3 sensitivity outputs;
- H3b PGLS details.

## 12. Submission-level claim ceiling

The paper can now state:

> Species-level flower-colour polymorphism is reproducibly measurable in high-depth photographic data, and a white-versus-nonwhite colour-space axis first localized in the original cohorts was prospectively confirmed under an unchanged construction-controlled test in a pre-frozen species-disjoint third cohort.

Required modifier:

> The prospective cohort is species-disjoint but drawn from the same iNaturalist opportunity universe and processed by the same measurement system.

## 13. Hard nonclaims

Do not claim:

- global polymorphism prevalence;
- independent-source H2 replication;
- pigment loss/gain;
- evolutionary direction;
- adaptive causation;
- recurrent non-white hue geometry;
- absence of all phylogenetic structure;
- irrelevance of biological range;
- near-perfect H1 reliability;
- a durable biological P500 H2 verdict.

## 14. Immediate manuscript-production priorities

The remaining work is no longer hypothesis discovery.

Priority order:

1. render Figures 1–4 from frozen outputs;
2. add source-audited literature citations to the Introduction and Discussion;
3. build supplementary evidence tables from machine-readable receipts;
4. add a manuscript-result guard that checks all headline numerical claims against frozen JSON/results;
5. only then choose journal-specific length/format.

Do not reopen H2 target selection or search for additional H3 predictors before these production tasks are complete.
