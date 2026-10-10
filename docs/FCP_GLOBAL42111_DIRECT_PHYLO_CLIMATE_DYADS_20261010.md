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
