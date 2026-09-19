# disttrait methods-paper architecture — 2026-09-19

Status: active architecture after validated v0.10 promotion.

Package:

`packages/disttrait/`

Current validated version:

`0.10.0`

## 1. Paper-level problem

Large opportunistic datasets often contain many individual observations from
many species, but the scientific question is frequently **within-species trait
organization**.

A pooled analysis can confuse:

- between-species geographic turnover;
- unequal species sampling depth;
- within-species phenotype structure;
- observer/measurement instability;
- model assumptions about whether species share the same response direction.

The methods paper should therefore not be framed as a new correlation
coefficient or diversity index.

The target problem is:

> How can repeated individual observations be converted into validated
> species-level trait distributions and spatial-organization summaries without
> confusing between-species turnover with within-species structure?

## 2. Central methodological contribution

The general architecture is:

[
	ext{individual observations}
ightarrow
	ext{validated species-level distribution}
ightarrow
	ext{within-species geometry / dissimilarity}
ightarrow
	ext{species-conditioned spatial statistic}
ightarrow
	ext{matched null or model-based cross-species inference}.
]

The reusable contribution consists of four separations.

### A. Measurement validity versus biological inference

A species-level distributional phenotype is validated before ecological
interpretation, for example with observer-disjoint reliability.

### B. Within-species structure versus between-species turnover

Spatial association is calculated within species before any cross-species
aggregation.

### C. Estimand choice versus test machinery

The package distinguishes:

- direction-invariant spatial organization;
- common signed responses;
- species-specific signed responses and directional heterogeneity.

These are different biological questions rather than competing estimates of one
universal parameter.

### D. Statistic versus calibration

A statistic may be useful while its asymptotic reference distribution is poor.
v0.8 demonstrates this directly: random-effects slope mean/Q statistics are
informative, but their small-sample analytic calibration can fail and can be
replaced by matched within-species permutation calibration.

## 3. Evidence chain

### Evidence 1 — exact implementation transport from FCP

Compact frozen fixtures establish exact algorithmic equivalence for:

- observer-disjoint species-level diversity reliability;
- deterministic Hellinger two-mode geometry;
- fixed one-vs-rest contrast;
- species-conditioned spatial vertex null;
- matched focal-minus-background null.

This establishes that `disttrait` is a genuine extraction of the validated FCP
inference layer rather than a post-hoc rewrite.

### Evidence 2 — pooled geographic confounding

In the targeted synthetic benchmark, species baseline trait frequencies turn
over geographically while within-species spatial dependence is absent.

Result:

- naive pooled false-positive fraction = **1.00**;
- species-conditioned matched-null false-positive fraction = **0.00** in the
  original 40-world benchmark.

The paper-level point is not that every pooled analysis fails. It is that a
between-species spatial gradient can answer a different question while looking
like within-species organization.

### Evidence 3 — FPR/power surface

Across a 24-cell surface varying effect size, observation imbalance and MCAR
missingness:

- naive pooled inference rejects in every deliberately confounded null cell;
- matched species-conditioned inference has maximum null rejection **0.025**;
- detection increases from **0.225–0.675** at effect 0.8 to **1.00** at the
  strongest effect.

This shows calibration/power behavior beyond a single synthetic example.

### Evidence 4 — aggregation is not the main innovation

Comparator benchmarks show:

- equal-species matched-null maximum null FPR = **0.025**;
- pair-weighted matched-null maximum null FPR = **0.025**;
- species-rho t-test maximum null FPR = **0.05**.

Thus **species-conditioning** is the main protection in this benchmark.
Equal-species weighting is a defensible default under unequal sampling but is
not uniquely optimal.

### Evidence 5 — model-based efficiency when the model is correct

A species fixed-intercept logistic model is substantially more powerful under a
binary-logistic data-generating process that matches its assumptions:

- weak-effect logistic detection = **0.975–1.00**;
- matched-null detection = **0.225–0.675**.

Its worst null cell reaches **0.10**, while the matched-null method remains
lower in that benchmark.

The interpretation is an assumptions/power trade-off, not a winner.

### Evidence 6 — external non-flower empirical transport

The generic continuous-trait layer was applied unchanged to a frozen external
San Francisco street-tree dataset:

- 20 taxon labels;
- 80 trees/taxon;
- 1,600 observations;
- trait = log(DBH);
- within-taxon dissimilarity = absolute log-DBH difference.

Result:

- equal-taxon mean spatial rho = **0.09219**;
- matched-null p = **0.01**.

The purpose is external transport, not causal urban-tree ecology.

### Evidence 6b — a second external transport can be null

A second independent empirical transport uses ShareTrait
`Gammarus insensibilis` metabolic-rate data:

- 375 individuals;
- three Adriatic populations;
- standardized metabolic-rate unit = Joule/day;
- primary trait = log(metabolic rate);
- 999 within-species vertex permutations.

Result:

- log-rate spatial rho = **-0.00683**, matched-null p = **0.68**;
- raw-rate sensitivity rho = **0.01093**, p = **0.198**.

This is non-support for the generic distance–dissimilarity spatial estimand in
this three-population fixture. That null outcome is useful: the external
transport layer is not being populated only with datasets that produce
significance.

### Evidence 7 — direction heterogeneity changes the estimand

When species share the same signed direction, a common-slope model is highly
efficient.

When half the species reverse direction:

- effect 0.4 common-slope detection = **0.05–0.10** in v0.7;
- distance-dissimilarity matched-null detection = **0.475–0.90**.

Thus a shared signed slope and direction-invariant spatial organization are not
interchangeable estimands.

### Evidence 8 — flexible species slopes require calibration

v0.8 adds species-specific signed slopes and random-effects meta-analysis.

The first analytic calibration failed:

- analytic mean maximum null FPR = **0.175**;
- analytic heterogeneity maximum null FPR = **0.30**;
- analytic omnibus maximum null FPR = **0.30**.

After 99 matched within-species trait permutations:

- calibrated combined omnibus maximum null FPR = **0.05**.

With effect 0.4 and 50% reversal:

- common-slope detection = **0.025–0.05**;
- distance-dissimilarity matched-null = **0.40–0.90**;
- calibrated slope-heterogeneity detection = **0.75–1.00**;
- calibrated slope-meta omnibus = **0.525–1.00**.

This is the strongest general methods result after v0.8:

> model flexibility and test calibration are separate problems.

### Evidence 9 — conditional calibration is not protection from joint MNAR selection

v0.10 keeps the latent within-species biological process null and changes only
the observation mechanism.

Across MCAR, trait-only and position-only selection:

- all benchmark cells retain 40/40 evaluable worlds;
- equal-species matched-null rejection remains at or below **0.05**;
- species-level spatial rho remains near zero.

When observation depends jointly on trait state and local position:

- equal-species matched-null rejection rises to **0.80–1.00**;
- species fixed-effect logistic rejection is **1.00**;
- the observed rows acquire mean absolute within-species trait-position
  correlations of **0.21–0.32**.

This should not be described as ordinary type-I error inflation conditional on
the observed sample. The observation process has created real association among
the retained rows.

The methods lesson is:

> matched randomization calibrates the statistic to the observed rows and
> observed sampling geometry; it cannot identify a counterfactual latent
> biological null after an unmeasured trait-by-location observation process has
> changed which rows are present.

Thus observation-process validity and inferential calibration are separate
problems.

### Evidence 10 — nonlinear response shape is another estimand choice

v0.11 generates a centered quadratic within-species response.

When all species share curvature direction at weak effect 0.4:

- correctly specified common quadratic detection = **1.00**;
- distance-dissimilarity matched-null detection = **0.10–0.375**;
- common linear detection = **0–0.025**.

When half the species reverse curvature direction:

- weak-effect common quadratic detection falls to **0.025–0.05**;
- distance-dissimilarity matched-null remains **0.15–0.35**;
- at effect 0.8, matched-null detection is **0.85–1.00** while common
  quadratic detection is **0–0.075**.

This extends the direction-heterogeneity result beyond linear slopes.

The lesson is:

> a correctly specified nonlinear model is highly efficient when species share
> the modeled response shape, while a direction-invariant dissimilarity
> estimand retains spatial-organization information when nonlinear response
> signs differ among species.

## 4. Recommended paper story

### Question 1

**Why is pooling dangerous when the target is within-species trait
organization?**

Answer: between-species turnover can generate overwhelming pooled significance
even when within-species spatial organization is absent.

### Question 2

**What is the minimal robust alternative?**

Answer: define the within-species estimand first, keep species identity during
randomization, then aggregate species-level evidence.

### Question 3

**When should model-based alternatives be preferred?**

Answer: when their response structure and shared-direction assumptions match the
scientific estimand, they can gain substantial power.

### Question 4

**What if response direction varies among species?**

Answer: either use a direction-invariant dissimilarity estimand, or explicitly
model species-specific slopes. Do not interpret cancellation of a common slope
as absence of spatial organization.

### Question 5

**How should flexible model summaries be calibrated?**

Answer: asymptotic calibration should be checked rather than assumed.
Species-conditioned matched randomization can calibrate final statistics to the
observed sampling geometry when the exchangeability null is appropriate.

### Question 6

**What does matched randomization not solve?**

Answer: it conditions on the observations that exist. If observation itself
depends jointly on trait and location, the observed sample can contain spatial
trait structure that neither permutation nor a parametric response model can
attribute uniquely to biology versus selection without additional
observation-process information.

### Question 7

**What if the within-species response is nonlinear?**

Answer: use a nonlinear model when the shared response shape is itself the
scientific target. If species can express equally strong organization with
opposite curvature or other shape reversals, a direction-invariant
dissimilarity estimand answers a different and often more stable question.

## 5. Main claim

A defensible one-sentence claim is:

> **Species-conditioning is the key inferential step for separating
> within-species trait organization from between-species geographic turnover;
> matched randomization provides an assumption-light calibration route
> conditional on the observed sample, while model-based signed-response methods
> offer complementary power and interpretation when their estimands and
> assumptions are appropriate.**

A slightly broader architecture claim is:

> **Repeated individual observations can support comparative inference on
> species-level trait distributions when measurement validity, within-species
> estimands, cross-species aggregation and null calibration are treated as
> separate stages.**

## 6. Hard nonclaims

Do not claim:

- a wholly new statistical family;
- universal superiority over hierarchical models, GAMs or Gaussian processes;
- universal type-I error control;
- that equal species weighting is always optimal;
- that signed slopes and distance-dissimilarity estimate the same parameter;
- robustness to arbitrary MNAR sampling;
- causal interpretation of the San Francisco tree result;
- that every opportunistic dataset supports a stable species-level phenotype.

## 7. Proposed display structure

### Figure 1 — estimand problem

Cartoon/simulation showing:

- between-species turnover;
- zero within-species structure;
- large pooled signal;
- null species-conditioned signal.

Purpose: make the confounding problem immediately visible.

### Figure 2 — calibration, power and observation-process boundary

Panels A–B: effect size × missingness × imbalance.

Compare:

- naive pooled;
- equal-species matched-null;
- pair-weighted matched-null;
- species-rho test.

Panel C: observation-process stress test.

Compare:

- MCAR;
- trait-only selection;
- position-only selection;
- joint trait-by-position selection.

Purpose: separate conditional calibration of the inferential machinery from
validity of the process that generated the observed rows.

### Figure 3 — assumption/estimand trade-offs

Three panels:

- correctly specified binary-logistic model versus matched-null inference;
- shared linear direction versus 50% direction reversal;
- shared quadratic curvature versus 50% curvature reversal.

Purpose: show that the preferred method depends on the target estimand and on
whether species share the signed response form.

### Figure 4 — calibration of flexible species slopes

Show:

- analytic random-effects null rejection inflation;
- permutation-calibrated random-effects null rejection;
- heterogeneity detection under 50% reversal.

This converts the failed preflight into a methodological result.

### Figure 5 — external empirical transport

Two external examples:

**San Francisco trees**

- selected taxa and spatial coverage;
- per-taxon continuous-trait rho;
- positive equal-taxon matched-null result.

**ShareTrait Gammarus**

- three sampled populations and 375 individuals;
- population metabolic-rate distributions;
- null log-rate and raw-rate matched-permutation results.

Purpose: show transport across unrelated systems and show that the same
inference layer can produce either support or non-support.

### Table 1 — method/estimand map

Columns:

- method;
- within-species conditioning;
- direction-sensitive?;
- distributional assumption;
- null/calibration;
- primary estimand;
- benchmark strength/weakness.

## 8. Suggested title family

Primary working title:

**Species-conditioned inference of spatial trait organization from repeated
individual observations**

Alternatives:

- **From opportunistic observations to species-level trait distributions:
  validated species-conditioned spatial inference**
- **Separating within-species trait organization from between-species
  geographic turnover**
- **Estimand-aware comparative inference for spatial trait distributions**

## 9. Remaining gates before external methods-paper submission

The current package is sufficient for drafting the paper.

The highest-value additions are now:

1. a multivariate continuous-trait benchmark;
2. standalone package repository/licence/release metadata.

The MNAR observation-process gap has been stress-tested in v0.10. Arbitrary
MNAR robustness remains a hard nonclaim rather than an unmet benchmark.

A second external empirical dataset is no longer an unmet gate.

These are strengthening steps rather than prerequisites for beginning the
manuscript.
