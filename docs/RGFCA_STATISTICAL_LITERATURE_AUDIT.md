# RGFCA statistical interpretation: primary-source audit

Audit date: 7 September 2026. This bounded reading supports interpretation and
method positioning, not a new analysis or a systematic review. It changes no
frozen result, permutation, measurement, admission rule or reserve control.
No reserve outcome file was opened. Source summaries below are deliberately
brief; their application to RGFCA is our interpretation, not a result of those
papers. See also the [image/ecology audit](RGFCA_IMAGE_ECOLOGY_LITERATURE_AUDIT.md).

## 1. Repeated subsampling is established methodology

John M. Drake (2015). *Range bagging: a new method for ecological niche modelling
from presence-only data.* Journal of the Royal Society Interface 12(107):20150086.
[10.1098/rsif.2015.0086](https://doi.org/10.1098/rsif.2015.0086).

Read: indexed [primary full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC4590497/),
Methods 2.1 and Discussion 3.2. Direct PMC access encountered a browser check.
Range bagging subsamples records and environmental dimensions without replacement,
fits marginal niche support and aggregates membership votes. Its target is
environmental niche estimation, not a multispecies colour-boundary field.
This establishes an ensemble/subsampling precedent, not RGFCA's priority or
validity. Repeated map construction alone is not a new-method claim.

## 2. Monte Carlo p-values have a finite reporting resolution

Belinda Phipson and Gordon K. Smyth (2010). *Permutation P-values should never be
zero: calculating exact P-values when permutations are randomly drawn.*
Statistical Applications in Genetics and Molecular Biology 9:Article 39.
[10.2202/1544-6115.1585](https://doi.org/10.2202/1544-6115.1585).

Read: [author manuscript](https://gksmyth.github.io/pubs/PermPValuesPreprint.pdf),
sections 5–7; [published bibliographic record](https://pubmed.ncbi.nlm.nih.gov/21044043/).
The paper treats finite random-permutation sampling and nonzero tail probabilities,
including the conservative plus-one calculation. This supports distinguishing
RGFCA's 999-draw p-value resolution from exhaustive enumeration. A reported
0.001 is the minimum attainable with its fixed formula, not proof that the
exhaustive tail probability equals exactly 0.001. More draws improve Monte Carlo
resolution, not biological sample size. The formula does not establish
exchangeability or remove an observation-process confounder.

## 3. One spatial field and association between two spatial fields are different nulls

Gilles Guillot and François Rousset (2013). *Dismantling the Mantel tests.*
Methods in Ecology and Evolution 4:336–344.
[10.1111/2041-210X.12018](https://doi.org/10.1111/2041-210X.12018).

Read: [publisher full text](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210x.12018),
model definitions and the sections on invalid and valid uses. The paper separates
testing one variable's spatial structure from testing independence of two
autocorrelated variables. Unrestricted permutations can invalidate the latter
by destroying marginal spatial structure; partialling geographic distance is
not a general remedy. It does not pronounce every distance-colour test invalid.
RGFCA's simple spatial-association result cannot establish environmental causation.
Matched coordinate geometry does not mean the null retains colour autocorrelation.
Future environmental association needs its own qualified null; this reading is
not permission to retune existing tests or relabel their outcomes.

## 4. A holdout tests a specified kind of transfer

David R. Roberts, Volker Bahn, Simone Ciuti, Mark S. Boyce, Jane Elith,
Gurutzeta Guillera-Arroita, Severin Hauenstein, José J. Lahoz-Monfort,
Boris Schröder, Wilfried Thuiller, David I. Warton, Brendan A. Wintle,
Florian Hartig, and Carsten F. Dormann (2017). *Cross-validation strategies for
data with temporal, spatial, hierarchical, or phylogenetic structure.*
Ecography 40:913–929.
[10.1111/ecog.02881](https://doi.org/10.1111/ecog.02881).

Read: [publisher full text](https://nsojournals.onlinelibrary.wiley.com/doi/full/10.1111/ecog.02881),
simulation/case-study comparisons and Discussion. Those comparisons show why
random validation splits can underestimate prediction error in structured data.
Blocking should match the intended generalization target and can introduce
extrapolation. This is methodological precedent, not a numerical calibration
of RGFCA. Its taxon-disjoint reserve still shares a platform, regions and a
measurement model; it is not independent of every systematic source of error.
No alternative reserve split is selected after reading outcomes.

## 5. Thinning does not have a universally validated setting

Valerie A. Steen, Morgan W. Tingley, Peter W. C. Paton, and Chris S. Elphick (2021).
*Spatial thinning and class balancing: Key choices lead to variation in the
performance of species distribution models with citizen science data.*
Methods in Ecology and Evolution 12(2):216–226.
[10.1111/2041-210X.13525](https://doi.org/10.1111/2041-210X.13525).

Read: [publisher full text](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13525),
Methods, independent-survey evaluation and Discussion. In models for 102 bird
species, performance effects of thinning/class balancing varied by model,
prevalence and evaluation metric. This supports qualifying sampling decisions
for their target, not claiming one universally optimal thinning distance.
Bird distribution-model results do not validate RGFCA's 0.25–0.5-degree cells,
minimum 40 classifiable images or observer cap of two. Those remain explicit
study choices, not literature-derived guarantees of unbiased flower sampling.

## 6. Bias correction and spatial resampling predate RGFCA

William Fithian, Jane Elith, Trevor Hastie, and David A. Keith (2015).
*Bias correction in species distribution models: pooling survey and collection
data for multiple species.* Methods in Ecology and Evolution 6(4):424–438.
[10.1111/2041-210X.12242](https://doi.org/10.1111/2041-210X.12242).

Read: [publisher full text](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12242),
Summary, observation-process model and spatial-bootstrap discussion. The paper
jointly models survey and collection records under explicit shared-bias
assumptions and uses spatial resampling in its analysis. RGFCA does not have
these complementary survey data or the same likelihood. Its opportunity
denominator corrects sampled edge availability, not the entire process that
determines which plants get photographed. Citing this precedent must not imply
that RGFCA implements or inherits its bias correction.

## Search and use limits

Sources were located by exact title/DOI and targeted queries for range bagging,
random permutation p-values, distance-matrix exchangeability, structured
cross-validation and spatial thinning. Primary publisher text, an author-hosted
manuscript and explicitly identified indexed primary sections were used.
Issue years, rather than earlier online dates, are used where these differ.
No unread reference cited inside a paper was adopted as supporting evidence.

The current manuscript can cite the narrow findings above. This does not close
all submission references: software, data-provider, environmental-layer and
remaining measurement-model citations still need final package-level checks.
Automated citation checks only establish internal consistency with these audited
notes; they cannot verify a paper's truth or certify the study's conclusions.
