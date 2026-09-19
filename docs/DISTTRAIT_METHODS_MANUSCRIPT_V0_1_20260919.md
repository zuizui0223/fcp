# Species-conditioned inference of spatial trait organization from repeated individual observations

Working manuscript v0.1 — 2026-09-19

Target article type: methods / quantitative ecology manuscript.

Package implementation: `disttrait 0.12.0`.

## Abstract

Large opportunistic datasets increasingly contain repeated individual-level
phenotypes for many species, but analyses that pool those observations can
confound between-species geographic turnover with within-species trait
organization. We develop a species-conditioned inference architecture that
separates four stages: validation of a species-level trait distribution,
definition of a within-species estimand, cross-species aggregation, and
calibration of the final statistic. The implementation, `disttrait`, uses
matched within-species randomization as an assumption-light calibration route
and supports categorical, scalar continuous and multivariate continuous traits.
In deliberately confounded null simulations, naive pooled analyses rejected in
all tested cells, whereas species-conditioned matched-null inference remained
near the nominal level across observation imbalance and MCAR loss. Correctly
specified response models were more powerful when species shared their response
form, but common signed effects lost power when species differed in response
direction, nonlinear curvature or multivariate orientation. Species-specific
slope models recovered directional heterogeneity, although standard
large-sample meta-analytic calibration was anti-conservative and required
matched permutation calibration. A separate observation-process stress test
showed the boundary of the approach: randomization calibrates inference
conditional on observed rows but cannot recover a latent biological null after
joint trait-by-location selection changes which observations are present.
Finally, the same continuous-trait layer transported to two external non-flower
datasets, producing one positive and one non-supporting result. The central
contribution is therefore not a new correlation coefficient, but an
estimand-aware workflow for comparative inference from repeated individual
observations.

## Introduction

Comparative ecology often compresses a species to a single trait value. That
representation is convenient, but it discards the within-species distribution
that is increasingly visible in community-science photographs, monitoring
programmes, museum records and other repeated individual observations. When the
scientific question concerns intraspecific trait variation, the statistical
problem is not simply how to calculate another species mean. It is how to
construct and validate a species-level distributional phenotype while
preserving the distinction between variation among species and organization
within species.

This distinction becomes especially important for spatial data. Suppose species
differ geographically in both occurrence and baseline phenotype. Pooling
individuals across species can then produce a strong relationship between
geographic distance and trait dissimilarity even if trait and location are
independent within every species. Such a pooled result is not necessarily
numerically incorrect: it answers a between-species turnover question. The
problem arises when it is interpreted as evidence for within-species spatial
organization.

A second difficulty is that there is no single model that is optimal for all
forms of intraspecific organization. A common signed slope is highly informative
when every species changes in the same direction along a comparable spatial
axis. The same statistic can cancel, however, when species are equally strongly
organized but differ in response direction. A nonlinear model can be powerful
when its response shape is correct and shared, but may summarize a different
quantity from a direction-invariant trait-dissimilarity statistic. Multivariate
traits add an analogous problem: species may show similar strengths of
organization while their response vectors point in different directions in
trait space.

A third difficulty is calibration. Flexible model summaries do not guarantee
that their usual asymptotic reference distributions are reliable at the
available species-level sample sizes. Conversely, randomization tests condition
on the data that were actually observed and therefore do not solve arbitrary
outcome-dependent observation processes. These issues are logically separate:
the choice of estimand, the statistic used to summarize it, the calibration of
that statistic, and the validity of the observation process should not be
collapsed into one decision.

Here we develop and benchmark a species-conditioned inference architecture for
repeated individual observations. The implementation, `disttrait`, was
extracted from a flower-colour analysis but is deliberately independent of
flower-colour biology. We ask six linked methodological questions. First, how
severely can pooling confuse between-species turnover with within-species
organization? Second, does species-conditioning remain calibrated under
observation imbalance and random observation loss? Third, when do parametric
response models gain power relative to matched randomization? Fourth, what
happens when response direction or nonlinear shape differs among species?
Fifth, can flexible species-specific slope summaries be calibrated to the
observed sampling geometry? Sixth, what does the method fail to identify when
observation itself depends jointly on trait and location?

The central claim is intentionally narrower than the introduction of a new
statistical family. We test whether repeated observations can support
comparative inference on species-level trait distributions when measurement
validity, within-species estimands, cross-species aggregation, null calibration
and observation-process assumptions are treated as separate inferential stages.

## Materials and Methods

### 1. General inference architecture

The workflow has five conceptual layers:

1. repeated individual observations;
2. validation of the species-level distributional representation;
3. definition of a within-species geometry or dissimilarity;
4. estimation of species-specific spatial organization;
5. cross-species aggregation with matched randomization or an explicit
   model-based alternative.

The package does not require every application to use every layer. For example,
observer-disjoint reliability is relevant when observer identity is available,
whereas an externally measured trait may enter directly at the within-species
spatial stage.

### 2. Species-level distributional traits

For categorical observations with category proportions (p_k), species-level
diversity is summarized by Gini–Simpson diversity,

[
D = 1 - \sum_k p_k^2.
]

The implementation accepts either category counts or proportions.

For multicomponent non-negative trait representations, rows can be normalized
and Hellinger transformed before species-specific geometry is estimated.
`disttrait` includes a deterministic unlabeled two-mode construction and
fixed-contrast alignment statistics because these were required by the source
flower-colour application, but the package does not require those objects for
spatial inference.

### 3. Observer-disjoint reliability

When observer identity is available, the same observer is never split between
halves. Complete observer groups are assigned deterministically to two balanced
sets without reading trait outcomes. The species-level diversity phenotype is
recomputed in each half, and agreement across species is summarized with
split-half Spearman correlation, concordance correlation and Spearman–Brown
diagnostics.

The generic implementation preserves the frozen logic of the source analysis,
including blank-observer exclusion, explicit classifiability rules when
provided, full-row observer balancing and deterministic hash-based tie breaking.

### 4. Species-conditioned spatial organization

For species (i), let (d^{geo}_{jk}) be pairwise geographic distance among
observations and (d^{trait}_{jk}) their pairwise trait dissimilarity. The core
species-level statistic is

[
\rho_i = \mathrm{Spearman}
(d^{geo}_{jk}, d^{trait}_{jk}).
]

Positive values indicate that geographically more separated individuals tend
to be more dissimilar in the chosen trait representation.

For categorical/compositional traits, the default dissimilarity is
Jensen–Shannon divergence. For a scalar continuous trait, the default is
absolute pairwise difference. For a multivariate continuous trait, the package
provides Euclidean pairwise distance. Multivariate dimensions are not
automatically rescaled; applications are responsible for defining a
biologically defensible metric.

### 5. Matched within-species null

The primary randomization keeps species identity, observed coordinates and the
complete set of observed trait rows fixed. Complete trait rows are permuted
among observed positions within species. This breaks the trait–position
association while retaining the observed sampling geometry, trait distribution
and within-row dependence among multivariate dimensions.

Each species therefore has one observed (ho_i) and a matched set of null
values. Cross-species evidence can be aggregated by taking the equal-species
mean in the observed data and in each matched null world. The upper-tail
Monte Carlo probability is calculated with the plus-one rule.

Equal species weighting is treated as a robust default rather than a universal
optimum. Separate benchmarks compare equal weighting with pair-count weighting
and species-level summary tests.

### 6. Model-based comparators

We compare matched-null inference with model-based alternatives under data
generating processes chosen to expose their assumptions.

For binary traits, a species fixed-intercept logistic model estimates one common
non-negative within-species slope and uses a one-sided profile-likelihood ratio
test.

For continuous traits, a species fixed-intercept common signed slope is
implemented by within-species demeaning. A common quadratic model adds a
centered squared-position term.

For direction heterogeneity, one signed slope is estimated for each species.
Those slopes and their sampling variances are summarized by a random-effects
meta-analysis. Because the usual large-sample normal and chi-square calibration
proved anti-conservative in the benchmark, the final mean and heterogeneity
statistics are recalibrated by matched within-species permutation of complete
trait values.

For two-dimensional continuous traits, a common signed multivariate response
vector is estimated after removing species intercepts. Its performance is
compared with Euclidean distance–dissimilarity inference when species response
vectors differ in orientation.

### 7. Pooled geographic-confounding benchmark

We simulate species whose baseline trait frequencies change strongly with
species geographic centre while the trait is spatially independent within each
species. This deliberately creates a setting in which a pooled geographic
analysis answers a between-species turnover question.

The benchmark compares naive pooled inference with the species-conditioned
matched-null statistic. Forty null and forty signal worlds are used in the
initial benchmark.

### 8. Effect, imbalance and MCAR performance surface

A broader 24-cell surface varies:

- within-species effect size: 0, 0.8, 1.6, 2.4;
- species-level observation imbalance: 1× or 4×;
- MCAR observation loss: 0%, 25% or 50%.

Each cell contains 40 replicate worlds and 20 species per world. The spatial
null uses 39 permutations per species. The benchmark records false-positive or
detection fractions and the median species-conditioned statistic.

### 9. Alternative aggregation and model comparison

On the same synthetic worlds we compare:

- naive pooled pairwise Spearman analysis;
- equal-species matched-null aggregation;
- pair-count-weighted matched-null aggregation;
- a one-sided one-sample test of species-specific (ho_i);
- a species fixed-intercept logistic model.

This benchmark separates the effect of species-conditioning from the effect of
one particular aggregation rule.

### 10. Direction-heterogeneity benchmark

For continuous traits, every species receives a local signed linear effect.
The sign is reversed in 0%, 25% or 50% of species. Effect size, reversal
fraction and MCAR loss are crossed in 24 cells.

This benchmark compares a common signed slope with the direction-invariant
distance–dissimilarity estimand. It therefore tests estimand alignment rather
than model ranking.

### 11. Permutation-calibrated species slopes

The species-specific slope benchmark uses the same direction-heterogeneity
worlds. The initial analytic random-effects calibration is retained as a
diagnostic failure. Final inference uses 99 within-species matched
trait-permutation worlds, re-estimating species slopes, their variances, the
random-effects mean, Cochran (Q) and (	au^2) in every null world.

### 12. Observation-process stress test

We hold the latent biological within-species association at zero and vary only
the observation process. Four mechanisms are tested:

- MCAR;
- trait-only selection;
- position-only selection;
- joint trait-by-position selection.

Two selection strengths are used, yielding eight cells with 40 replicate worlds
each. This benchmark is interpreted as a distinction between conditional
calibration on observed rows and identification of the latent biological
process.

### 13. Nonlinear curvature benchmark

A continuous trait follows a centered quadratic local-position effect within
species. Curvature sign is reversed in 0%, 25% or 50% of species. The benchmark
compares naive pooling, matched-null distance–dissimilarity, a common linear
slope and a correctly specified common quadratic curvature.

### 14. Multivariate orientation benchmark

A two-dimensional continuous trait receives a species-specific local response
vector. Response orientations either coincide, span half of the trait-space
circle, or span the full circle. Effect sizes 0, 0.4, 0.8 and 1.2 are crossed
with 0% and 50% MCAR loss, yielding 24 cells with 40 worlds each.

The matched-null analysis permutes complete multivariate rows and uses Euclidean
trait-space distance. The comparator estimates one common signed multivariate
response vector.

### 15. External empirical transport

#### San Francisco street trees

A fixed historical snapshot of the San Francisco street-tree inventory is used
as an external non-flower continuous-trait application. The frozen fixture
contains 20 taxon labels and 80 observations per taxon. The trait is
(log(DBH)), and pairwise trait dissimilarity is absolute difference in
(log(DBH)). Each taxon uses 99 matched vertex permutations.

This dataset is used as a transport validation rather than as a causal analysis
of urban tree size.

#### ShareTrait Gammarus

A second external transport uses ShareTraitDatabase
`Gammarus insensibilis` metabolic-rate data linked to source dataset
TRADAT039. The fixture contains 375 individuals from three Adriatic
populations, all measured in Joule/day. The primary trait is log metabolic rate,
with 999 within-species vertex permutations. Individuals inherit their sampled
population coordinates, so this example represents between-population rather
than fine-scale within-site geography.

### 16. Software and reproducibility

The validated package version is `disttrait 0.12.0`. Compact fixtures test
exact equivalence with the source FCP/RGFCA algorithms. Benchmark receipts,
external empirical fixtures and dedicated GitHub Actions workflows are frozen in
the repository.

The release gate separately builds a source distribution and wheel, validates
both with Twine, installs the wheel in a fresh Python environment and tests the
installed package outside the repository working directory.

## Results

### 1. Pooling can produce overwhelming evidence for the wrong level of organization

In the targeted geographic-confounding benchmark, naive pooled analysis rejected
in every null world, whereas the species-conditioned matched-null test had zero
false positives in the initial 40-world null benchmark.

Across the broader performance surface, naive pooled inference again rejected
in every deliberately confounded null cell. The maximum null rejection fraction
for the equal-species matched-null method was 0.025.

These results show that the main protection is conditioning on species before
cross-species aggregation.

### 2. Equal species weighting is useful but not uniquely optimal

In the comparator surface, maximum null rejection fractions were 0.025 for both
equal-species and pair-weighted matched-null summaries and 0.05 for the
species-rho one-sample test. Under weak signal and 4× observation imbalance,
equal weighting was somewhat more powerful than pair weighting in several
cells, while the simple species-rho test was also competitive.

Thus the general methodological result is species-conditioning rather than the
unique optimality of one aggregation statistic.

### 3. Correct response models gain power when their assumptions are correct

Under the binary-logistic data-generating process, the species fixed-intercept
logistic comparator detected a weak effect in 0.975–1.00 of worlds, compared
with 0.225–0.675 for the equal-species matched-null omnibus. The model therefore
used the correctly specified response form more efficiently.

Its maximum null rejection reached 0.10 in one cell, while its mean across null
cells was 0.05. The matched-null method had lower worst-cell rejection in the
same benchmark.

### 4. Shared signed responses and direction-invariant organization are different estimands

When continuous-trait species shared one signed direction, a common-slope model
had high power. With 50% reversal, weak-effect common-slope detection fell to
0.05–0.10, while matched-null detection remained 0.475–0.90. At effect size
0.8 with 50% reversal, matched-null detection was 1.00 while common-slope
detection was 0.05–0.125.

The common slope did not fail to estimate its target. Its target cancelled
because species did not share one signed response.

### 5. Flexible species slopes need calibration as well as flexibility

The first random-effects slope analysis used standard large-sample
meta-analytic reference distributions. In the small-sample/missingness
benchmark, maximum null rejection rose to 0.175 for the mean component and 0.30
for the heterogeneity and combined tests.

After matched within-species permutation calibration, the combined omnibus had
a maximum null rejection fraction of 0.05. At weak effect with 50% reversal,
calibrated slope-heterogeneity detection was 0.75–1.00 and the combined
slope-meta omnibus detected 0.525–1.00 of worlds.

Model flexibility and test calibration were therefore separate problems.

### 6. Conditional calibration does not recover an unobserved biological null

Under MCAR, trait-only and position-only selection, equal-species matched-null
rejection remained at or below 0.05. Under joint trait-by-position selection,
rejection rose to 0.80–1.00, and the species fixed-effect logistic model rejected
in all worlds.

The joint selection process created a trait–position association among observed
rows. Matched randomization correctly tested exchangeability conditional on
those observed rows but could not identify whether the association was
biological or selection-induced.

### 7. Nonlinear response shape produces the same estimand distinction

Under a shared weak quadratic effect, the correctly specified common quadratic
model detected the signal in all worlds, while matched-null detection was
0.10–0.375 and the common linear model had almost no power.

When half the species reversed curvature sign, weak-effect common quadratic
power fell to 0.025–0.05. At effect size 0.8, matched-null detection was
0.85–1.00 while common quadratic detection was 0–0.075.

### 8. Multivariate orientation heterogeneity preserves trait-space organization while cancelling the common vector

Across effect-zero multivariate cells, maximum null rejection was 0.05 for both
the matched-null and common-vector methods, while naive pooled inference
rejected in every cell.

With effect size 0.4 and shared orientation, matched-null detection was
0.975–1.00 and common-vector detection was 1.00. When species response vectors
spanned the full trait-space circle, matched-null detection remained
0.95–1.00 but common-vector detection fell to 0.075. At effect size 0.8 under
full orientation spread, matched-null detection was 1.00 while common-vector
detection was 0.10–0.25.

The same conclusion therefore extends from signed scalar responses to
multivariate trait-space orientation.

### 9. The inference layer transported to unrelated empirical systems

For the San Francisco street-tree fixture, the equal-taxon mean within-taxon
spatial rho for log(DBH) was 0.09219 with matched-null p = 0.01. The association
between taxon-level DBH spread and spatial rho was not supported.

For ShareTrait `Gammarus insensibilis`, the primary log metabolic-rate
distance–dissimilarity result was non-supporting: rho = -0.00683, p = 0.68. A
raw-rate sensitivity analysis was also non-supporting: rho = 0.01093,
p = 0.198.

The two examples show executable transport with both positive and null outcomes
rather than selecting only significant applications.

## Discussion

### Species-conditioning is the central inferential step

The most consistent result across the benchmark stack is not that one omnibus
statistic dominates every alternative. It is that the level of conditioning
must match the scientific question. When the target is within-species trait
organization, species identity must be preserved before evidence is aggregated
across species.

This principle is visible most starkly in the pooled-confounding simulations.
A strong between-species geographic gradient can produce overwhelming pooled
significance while every species individually satisfies the within-species
null. Equal weighting is useful under sampling imbalance, but comparator
analyses show that other species-conditioned summaries can also be calibrated.
The main inferential protection is therefore species-conditioning.

### Estimand choice should precede method ranking

The benchmark comparisons also argue against a universal winner. A correctly
specified logistic or quadratic model can be substantially more powerful than a
generic dissimilarity statistic. That is a feature, not a contradiction. Those
models use response-form information that the matched-null analysis does not
assume.

The converse is equally important. When response signs, curvature directions or
multivariate orientations differ among species, a common signed parameter can
approach zero even while every species remains strongly organized. A
distance–dissimilarity estimand retains such organization because it asks about
strength rather than shared direction.

The practical implication is that the analysis should begin with the biological
quantity of interest. If the target is a common directional response, a common
slope or response vector is appropriate. If the target is whether phenotypes
become more dissimilar with spatial separation despite species-specific
directions, pairwise dissimilarity answers a different question.

### Randomization calibrates a statistic, not the observation process

Matched randomization has an important but bounded role. It is attractive
because it keeps the realized species identity, coordinates, trait values and
sampling geometry while breaking only the assignment of traits to positions.
This can provide useful finite-sample calibration without specifying a complete
response distribution.

That conditional logic does not reconstruct data that were never observed. The
joint trait-by-location selection benchmark makes this boundary explicit. Once
the observation process induces trait–position association in the retained
sample, both a permutation test and a parametric response model can detect that
association. Additional observation-process information is required to decide
whether it reflects biology, observation or both.

This boundary is particularly relevant to opportunistic data. Observer
restrictions, fixed sampling budgets, matched controls, outcome-blind
acquisition and validation against external data can reduce some forms of
feedback, but they do not prove the absence of every outcome-dependent
observation process.

### Flexible models still require finite-sample calibration checks

The failed analytic random-effects preflight provides a useful caution. Allowing
species-specific slopes solved the structural limitation of a common slope, but
the usual large-sample mean and heterogeneity tests were anti-conservative under
small per-species samples and missingness. Recalibrating those same statistics
by matched within-species permutation restored the combined benchmark behavior.

Thus model flexibility and inferential calibration should be considered
separately. A more flexible model does not automatically imply a better
reference distribution for its test statistic.

### Generalization across trait representations

The package currently supports categorical/compositional traits, scalar
continuous traits and continuous multivariate traits. These representations use
different dissimilarities but share the same species-conditioned logic.

The multivariate benchmark also highlights an application responsibility.
Euclidean distance is meaningful only after the analyst has decided how trait
dimensions should be scaled and combined. The package deliberately does not
hide that decision behind automatic standardization.

### External transport is validation of execution, not a universal ecological pattern

The two external examples serve a deliberately limited role. Street-tree DBH
showed positive within-taxon spatial organization, whereas Gammarus metabolic
rate did not show the same generic distance–dissimilarity pattern. The contrast
is useful because it demonstrates that the package is not being used to select
only supporting biological examples.

Neither result is interpreted causally. Planting history, age, management and
inventory structure can affect urban tree DBH, while metabolic physiology is
known to depend on predictors not represented by the generic three-population
distance–dissimilarity analysis. The transport claim is that the inference layer
executes coherently outside the original flower-colour system.

### Scope and limitations

The framework should not be described as a wholly new statistical family.
Gini–Simpson diversity, rank correlation, Hellinger transformation,
Jensen–Shannon divergence, linear models, random-effects summaries and
permutation principles are established components. The contribution is their
assembly into an explicit comparative architecture with separate stages for
measurement validation, within-species estimation, cross-species aggregation
and calibration.

The framework also does not guarantee universal type-I error control,
robustness to arbitrary MNAR observation, optimal trait-space metrics or
superiority over hierarchical models, GAMs or Gaussian processes. The benchmark
suite instead maps specific conditions under which different estimands and
assumptions become informative.

## Conclusions

Repeated individual observations can support comparative inference on
species-level trait distributions, but only if the level of the scientific
question is kept separate from between-species turnover and from the process
that generated the observations.

The benchmark stack supports three practical rules.

First, condition on species before interpreting spatial trait organization as a
within-species process. Second, define the estimand before choosing the model:
shared signed responses, species-specific directional heterogeneity and
direction-invariant spatial organization are different biological questions.
Third, validate calibration and observation-process assumptions separately.

Matched within-species randomization provides a useful assumption-light route
for calibration conditional on the observed sample. Model-based approaches can
be more efficient when their response forms and shared-direction assumptions
are appropriate. Used together, these approaches provide a transparent route
from repeated individual observations to comparative species-level inference.

## Data and code availability

The working implementation is `disttrait 0.12.0` inside the FCP repository.
Frozen benchmark receipts, exact implementation-equivalence fixtures, external
empirical fixtures and CI workflows are version controlled with the package.

A source distribution and wheel for version 0.12.0 have been built, validated
with Twine and installed successfully in a fresh Python environment. Final
public-release metadata, standalone repository location, software licence and
citation metadata remain to be frozen before external package publication.

## Figures planned

1. **Estimand/confounding problem** — pooled between-species turnover versus
   null within-species organization.
2. **Calibration/power/observation-process boundary** — performance surface plus
   joint-MNAR stress test.
3. **Estimand trade-offs** — binary response model, linear direction,
   nonlinear curvature and multivariate orientation.
4. **Flexible species slopes and calibration** — analytic versus matched
   permutation calibration.
5. **External transport** — San Francisco trees and ShareTrait Gammarus.

## Table planned

**Table 1. Method–estimand map**

Columns:

- method;
- within-species conditioning;
- direction-sensitive;
- trait representation;
- distributional assumptions;
- calibration;
- target estimand;
- demonstrated strength;
- demonstrated limitation.

## Reference work remaining before submission

External literature references are intentionally not finalized in this v0.1
draft. Before journal submission, the introduction and discussion require a
source audit covering at minimum:

- intraspecific trait variation / trait-distribution ecology;
- opportunistic and community-science observation bias;
- permutation/randomization inference;
- spatial trait dissimilarity and distance-based methods;
- hierarchical/random-slope comparative models;
- multivariate trait-distance choices;
- conditional versus observation-process inference.

The external empirical data sources already frozen in the repository must also
be cited in the final reference list, including the San Francisco Public Works
street-tree source and the ShareTrait / Shokri et al. Gammarus records.
