# Flower-colour polymorphism figure plan — 2026-09-18

This figure plan supersedes the provisional figure architecture in `docs/POLYMORPHISM_PAPER_ARCHITECTURE_20260913.md` where it conflicts with the completed prospective third-cohort H2 result.

Authoritative claim ledger:

- `docs/POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260918.md`

Authoritative prospective H2 result:

- `results/polymorphism_h2_third_cohort_prospective_white_axis_20260917/result.json`

## Figure 1 — From a global frame to a measurable polymorphism phenotype

### A. Sampling architecture

Show the hierarchy without implying prevalence:

```
42,111-species global frame
        ↓ high-depth opportunity / frozen cohort rules
original discovery source cohort: 500 × 100 photos
original reserve source cohort:   500 × 100 photos
        ↓ >=40 classifiable
discovery D cohort: 369 species
reserve D cohort:   363 species

third-cohort outcome-blind candidate frame: 3,230 species
        ↓ deterministic frozen selection
500 selected
        ↓ fresh metadata, no replacement
499 species × 100 rows = 49,900 terminal rows
        ↓ >=40 classifiable
377 measurement-evaluable species
```

The visual must explicitly label the 42,111 frame as an opportunity/sampling frame, not a prevalence denominator.

### B. Four biological states and D

Show the four frozen states:

- white;
- yellow/orange;
- red/pink;
- blue/purple.

Define

`D = 1 - sum_k p_k^2`.

Show `mixed_uncertain` outside the four-state simplex and label it "not a biological fifth state".

### C. D distributions

Plot discovery and reserve full-data D distributions for descriptive context only.

Do not annotate a percentage "polymorphic" because the paper does not estimate prevalence.

### D. Observer-disjoint split logic

Diagram observers assigned to A or B with no observer crossing the split.

Purpose: establish why H1 is a measurement-validity test rather than a simple image resampling exercise.

---

## Figure 2 — H1: reproducibility across disjoint observers

### A. Representative split scatter

Reserve D_A vs D_B with 1:1 line.

### B. Distribution across 200 frozen partitions

Show reserve Spearman rho across all 200 observer-disjoint partitions.

Required annotations:

- median rho = **0.7891**;
- q05 = **0.7652**;
- q95 = 0.8109;
- primary median floor = 2/3;
- primary q05 floor = 0.5.

### C. Agreement diagnostics

Show reserve CCC distribution or paired D differences.

Headline:

- median CCC = **0.8548**;
- median Spearman-Brown reliability = **0.8821**.

### D. Later strict stress test

Single clearly separated point/interval:

- rho = **0.7927**;
- bootstrap 95% CI = **0.7419–0.8324**;
- strict floor = 0.80;
- decision = stress-test not supported.

The panel must visually communicate that the later stress test constrains the claim but does not overwrite the chronologically earlier primary H1 result.

---

## Figure 3 — H2: how the recurrent axis was identified

This figure is the discovery/audit figure. It must not visually imply that q_white was prospectively selected in the original cohorts.

### A. Nine-colour continuous representation

Palette order:

`white, yellow, orange, red, pink, magenta, purple, blue, bronze`.

Show Hellinger transformation and deterministic unlabeled two-means.

### B. Species displacement geometry

Illustrate the two-mode displacement vector and unit axis `u_i`.

Show the primary 0.10 and strict 0.20 construction/admissibility gates.

### C. Fixed white contrast

Display

`q_white = normalize([1,-1/8,-1/8,-1/8,-1/8,-1/8,-1/8,-1/8,-1/8])`.

Label:

**target isolated after the original broad geometry was opened**

### D. Legacy targeted W results

Plot observed W against each cohort's structured-null distribution.

Primary 0.10:

- discovery: N = 152, W = **0.514625**, p = **0.001**;
- reserve: N = 129, W = **0.514586**, p = **0.001**.

Strict 0.20:

- discovery: N = 75, W = **0.542355**, p = **0.001**;
- reserve: N = 65, W = **0.510517**, p = **0.008**.

### E. Falsification by white-axis removal

Show that construction-controlled directional support collapses after q_white projection.

If space is limited, move the non-white-only diagnostic to Supplement but retain the projection-removal panel in the main figure.

Purpose: establish the narrow target and its retrospective discovery status before the prospective figure.

---

## Figure 4 — Prospective third-cohort confirmation of the frozen axis

This is the paper's decisive biological figure.

### A. Outcome-blind third-cohort chronology

Use a compact timeline:

```
3,230 outcome-blind candidates
  → frozen deterministic 500-species selection
  → fresh metadata
  → 499 × 100 authorized rows
  → pre-opening technical qualification
  → one-shot biological run
  → 256/256 terminal partitions
  → support gate
  → H2 opened once
```

Include "0 replacements" and "pixels not persisted".

### B. Measurement-support gate

Show:

- terminal rows = **49,900**;
- classifiable rows = **25,788**;
- measurement-evaluable species = **377**;
- required minimum = 250;
- support = **PASS**.

Do not show W in this panel; the chronology is support-before-H2.

### C. Primary prospective result

Structured-null distribution for primary 0.10 tier with observed line:

- N = **158**;
- W_obs = **0.5172457461**;
- null median = **0.4571428150**;
- null 95% interval = **0.4358491120–0.4752987776**;
- p = **0.001**.

### D. Strict prospective sensitivity

Structured-null distribution for 0.20 tier:

- N = **86**;
- W_obs = **0.5329282123**;
- null median = **0.4593196659**;
- null 95% interval = **0.4328679570–0.4867224043**;
- p = **0.001**.

### E. Cross-cohort summary

A compact effect summary for discovery, reserve and third cohort.

Plot W or observed/null-median ratio, but distinguish evidence status:

- discovery/reserve = target localization / post-audit;
- third cohort = untouched prospective confirmation.

Do not pool p-values or call the three cohorts independent data-source replications.

---

## Figure 5 — Broad explanatory tests do not replicate

### A. H3a reserve phylogenetic signal

Three points/intervals for S1-S3 Blomberg K with permutation-null context:

- S1: K = **0.0710190**, p = **0.2716**;
- S2: K = **0.0601476**, p = **0.4134**;
- S3: K = **0.0707577**, p = **0.2674**.

Optional secondary annotation: lambda estimates are near zero and unsupported.

### B. H3b discovery versus reserve sampled-span association

Place effect estimates side by side:

- discovery: rho = **0.1798786**, p = **0.00089996**;
- reserve: rho = **-0.0025855**, p = **0.9586021**.

The graphic should emphasize replication contrast, not "significant vs nonsignificant" alone.

### C. Interpretation boundary

Small text panel:

- phylogenetic test = broad tree-wide signal only;
- span predictor = sampled photographic span, not true geographic range;
- negative H3 results do not negate H1/H2.

---

## Supplementary figures

### Figure S1 — H1 full partition diagnostics

All 200 reserve partitions; paired species counts, rho, CCC, bias and absolute difference.

### Figure S2 — H1 discovery concordance and strict split details

Discovery repeated partitions plus deterministic stress-test diagnostics.

### Figure S3 — Broad pre-target H2 geometry

Leading eigenvalue / axis alignment / reserve transport diagnostics from the original label-free analysis.

### Figure S4 — Structured-null construction audit

Demonstrate which quantities are preserved under the coarse-state-preserving null.

### Figure S5 — Residual and non-white-only H2 diagnostics

Full q_white projection-removal results and low-N non-white-only tests.

### Figure S6 — Third-cohort chain of custody

Selection hashes, metadata freeze, 256 partition receipts, support gate and terminal H2 receipt.

### Figure S7 — P500 terminal postmortem

Measurement support passed, but H2 has no durable verdict because serialization failed after calculation. This figure is provenance/chronology only and must not display reconstructed H2 results.

### Figure S8 — H3a full sensitivity panel

Raw D, unbiased D, opportunity-adjusted residuals, K and lambda across S1-S3.

### Figure S9 — H3b full sensitivity panel

Raw, unbiased, partial-rank and rank-PGLS reserve results.

---

## Main-text visual narrative

The five main figures should read as one argument:

```
Fig. 1  define a species-level polymorphism phenotype
   ↓
Fig. 2  show that D survives observer separation
   ↓
Fig. 3  identify the narrow white/nonwhite geometric target
   ↓
Fig. 4  confirm that already-frozen target prospectively
   ↓
Fig. 5  show that broad phylogeny/span do not explain the result
```

The key transition is Fig. 3 → Fig. 4: **target discovery/audit is visually separated from prospective confirmation**.
