# FCP original 42,111 species: dated LCVP phylogeny × local photographed colour discordance × rainfall (2026-10-10)

## Source-tested comparison opportunity

This separate exploratory branch reuses already completed direct backbone local photo comparison preflight (Actions 38008375715). Original 250km local congeneric cohort has 872 source photographed species, 342 genuine LCVP direct tree tips; 500km cohort has 1761 photographed species, 649 genuine LCVP direct tips. Photo-site local comparison groups are fixed at maximum pairwise great-circle distance <=50 or <=100km, using only original species ID, genus, original cell and actual source photo latitude/longitude.

| Historical source cohort | Max photo distance | Direct LCVP species in multispecies groups | Original tipped pair dyads | Genera | Local 3+ groups with varied patristic distances |
|---|---:|---:|---:|---:|---:|
| 872 taxa /342 direct tips | 50km | 167 | 244 | 47 | 22 |
| 872 /342 | 100km | 226 | 477 | 65 | 29 |
| 1761 /649 | 50km | 256 | 318 | 65 | 31 |
| 1761 /649 | 100km | 343 | 593 | 85 | 44 |

Within-group 3+ tipped groups with variable genuine patristic distances contain 89,124,120 and 179 original source taxa, respectively. This makes an exploratory direct-tip regression numerically possible, not a full 1761 species phylogenetic correction.

## Fixed source-exploratory model

Construct all undirected dyads between DIFFERENT original source photo species in each complete-link same-genus × same-original-cell local microgroup. A species can contribute to several dyads within its group, so count distinct species, groups and genera separately.

- Outcome: original coarse four-state flower-photo colour mismatch, not measured genetic polymorphism.
- Baseline: log(1 + source photograph great-circle distance), unsigned absolute-latitude difference, elevation difference, unsigned WorldClim BIO1 annual mean temperature and BIO5 warmest-month maximum difference.
- True phylogenetic feature: log(1 + dated-tree patristic distance) between exact original LCVP species tips. Never graft unsupported tips for this test.
- Moisture: unsigned BIO12 annual precipitation and BIO15 precipitation-seasonality differences.

Compare four models on IDENTICAL original photo dyads and training/test folds: GEO+THERMAL, +true dated phylogeny, +rain without phylogeny, +dated phylogeny+rain. GroupKFold with five folds excludes whole genus from training when predicting source-pair mismatch. All 4 historical source cohort × microdistance combinations are reported. Nested cohorts, overlapping photographs and pairwise records are NOT independent replications.

Estimate source photographed heldout binary log loss and AUC. Perform 999 conditional genus-cluster resamplings of original out-of-group loss differences. A positive patristic increment means the dated reference tree improves these source photo predictions, NOT demonstrated evolutionary causation or heritable pigment variants. A moisture increment after patristic control is exploratory only; genetic population fitness and a proper phylogenetic GLMM are absent.

HOLD if below 100 local direct tips, 20 source genera, 30 original local multispecies groups, 100 photo pair dyads, or missing binary outcome class from training folds. No threshold rescue.

## Mandatory overall inference decision

The original full 872/1761-species phylogenetic correction stays HOLD: the real LCVP direct-tip coverage of 39.2%/36.9% fails the previously imposed 50% threshold. Even a successful direct-tip exploratory model cannot rescue the original local-rainfall selection claim, because truly geographically matched 50/100km source photo-colour nulls were unsupported. All analyses are retrospective and source-exposed.

A species photographed once has not had its genetic colour established; a phylogenetic tree is a reference hypothesis whose intra-genus topology can be uncertain. Genus-bootstrap on shared-species dyads does not capture complete phylogenetic uncertainty, photographer bias, soil/biotic covariation or local adaptation.

Original 1499 high-depth H1/H2 and future 2000+730 independent taxa unchanged.


## Executed actual source result and final ecological decision

**Source-verified [Actions 38008814183](https://github.com/zuizui0223/fcp/actions/runs/38008814183)** passed all synthetic tests, exact source 42,111 / 18,457 denominators and 342 / 649 direct LCVP tips, and preserved the pre-existing full-cohort phylogenetic HOLD. Complete [machine result](../results/fcp_global42111_direct_LCVP_local_photo_environment_20261010/result.json).

One photographed congeneric dyad has a binary source **coarse-colour mismatch** outcome. Positive numbers below denote source original **genus-heldout binary log-loss reductions** from adding a predictor family beyond site great-circle, latitude, elevation and thermal differences, NOT a regression coefficient or adapted fitness.

| Parent photographed locality group | Local maximum pair diameter | Original direct-tip source dyads (mismatch) | Add dated LCVP path after geo+heat (95% genus-block fixed-prediction CI) | Add rainfall after geo+heat+true phylogeny (95% genus-block CI) |
|---|---:|---:|---:|---:|
| 250km | 50km | 244 (147) | −0.013616 [−0.021600,−0.003333] | +0.011687 [−0.025528,+0.055415] |
| 250km | 100km | 477 (295) | −0.007161 [−0.022581,−0.003391] | +0.004917 [−0.006539,+0.028613] |
| 500km | 50km | 318 (181) | −0.011836 [−0.017652,+0.003020] | **+0.033123** [+0.007378,+0.084779] |
| 500km | 100km | 593 (352) | −0.001982 [−0.003735,+0.002773] | +0.006716 [−0.001039,+0.013410] |

Phylogenetic distance did **not** improve heldout mismatch prediction consistently; in 250km source cohorts its bootstrap interval for incremental gain was negative. Rainfall after the dated-tree path showed one positive conditional interval (500km parent, <=50km microgroups), but was **not** robust when matching source parent/subgroup sizes changed; these cohorts are nested and selected original direct-tip taxa comprise only 37–39% of the full source parent populations. The 500km/100km source subset did not confirm the 500km/50km one. Do not choose the apparently favourable one condition as the ecological truth.

### Three distinct noninterchangeable layers of support

1. **True tree comparability:** direct 73,420-tip LCVP branch lengths exist and are non-constant across 3+ original source local phylogenetic groups, but coverage is incomplete.
2. **Exploratory conditional predictive model:** additional source photo rainfall improves mismatch prediction in one nested subgroup; genus-heldout CV with genus-bootstrap on shared-species dyads is not independent genetic or causal confirmation.
3. **Strong spatial assumption check:** the more relevant original broader 872/1761-source-species group-conditioned photographic-colour null fails at within-50/100km microgeographical constraints (all four p=.175–.880). This is not repaired by the newly restricted direct-tip photo-pair model; a different analysis scale and subset cannot retrospectively validate earlier p-values.

### Current paper-facing inference

**Supported as retrospective descriptive biogeography:** broad one-photo flower-colour differences across sampled taxa and regions, weak source climate pattern among congeneric species, observable environmental heterogeneity among real photographed sites, and concrete limitations from fine spatial structure and incomplete direct phylogenetic species coverage.

**Not supported:** any spatially robust, lineage-corrected, genetically identified rainfall effect causing local flower-pigment selection. Neither SoilGrids nor direct source LCVP patristic control provides universal support across matched source cohorts. A small conditional improvement in a selected 318-dyad subset cannot outweigh the incomplete phylogenetic coverage, multiple nested comparisons and spatially matched null results.

**Future identification needed:** independent species-level phylogenetic resolution and source taxon replicates (rather than one source image), genuinely colocated plant populations, flower-colour genotype/pigment and pollinator/reproductive success, plus independent taxon/year replication not already used for these exploratory model choices. The untouched 2,000+730 prospective role is not sampled now.

The complete exploratory numerical evidence is kept in PR #143, which remains a draft. Full source PR #142 stays HOLD on 50%-direct-tip coverage. The original high-depth FCP manuscript and main were not altered.


## Additional source-verified LCVP taxonomic tip-availability selection bias

[GitHub Actions 38009071842](https://github.com/zuizui0223/fcp/actions/runs/38009071842) successfully replayed the immutable original 42111 source photos and fixed 872/1761 taxa, comparing **342/649 original directly-tipped source species** to the remaining **530/1112 original nontipped taxa** without changing original photo class labels or taxon identities. [Source numeric receipt](../results/fcp_global42111_LCVP_direct_tip_selection_bias_20261010/result.json).

- In the 250km original cohort, direct vs nontipped original photo taxa differed in source geographic/climate variables by max |standardized mean difference| **0.195** (annual mean temperature), annual precipitation SMD **−0.122**, and up to **8.22 percentage points** difference in four-colour photographed class fraction (yellow/orange: 28.95% direct vs 37.17% untipped).
- In the 500km original cohort, max climate/geography |SMD| **0.220** (annual mean temperature), annual precipitation SMD **−0.139**, and max original photographic colour fraction difference **1.45 percentage points**. Direct tips had lower original photo-site annual rainfall medians (**658mm vs 756.5mm**), lower annual temperature median (**16.30 vs 16.67 in native reported WorldClim units**), and higher source photo-site elevation medians (**325m vs 270m**).
- The original source selection into available LCVP direct tips was TAXONOMY-defined rather than source photo-colour-outcome-defined, but it remains associated with photographic species distributions and climate. In particular the direct-tip-only source 500km/50km rain interval cannot be generalized to the original 1761 taxa under the preexisting full-sample 50% tip coverage gate.

**Scientific decision retained:** no robust, full-sample, spatially and phylogenetically isolated rainfall effect on heritable local flower-colour evolution can currently be reported. A single nested 318-dyad exploratory positive conditional interval is neither full-population evidence nor independent confirmation. This direct-tip selection diagnostic strengthens the requirement for better species-representative phylogenetic and phenotype sampling, not a rationale for replacing original untipped species with artificial grafts.
