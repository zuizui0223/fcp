# FCP global 42,111 species: scientific claim gate after fine spatial matched null

**Source audit date: 2026-10-09. Status: SPATIALLY ROBUST LOCAL RAINFALL-COLOUR CLAIM NOT SUPPORTED.**

## Original source and estimands

The immutable photographic global atlas contains 42,111 original named taxa (18,457 original one-photo four-colour classifiable), and separate 85,337 original species×regional-cell photo measurements (39,075 classifiable). A within-genus photocolour study found small rainfall-related predictive improvement among different sampled congeneric species, conditioning on genus category means and coarse geography. The historical within-species one photographed pair per 13,416 taxa (3,182 climate-complete fully classifiable, 793 coarse-colour mismatches) did not show stable source-photo rainfall mismatch prediction. SoilGrids had no robust independent incremental signal when constraining original photographed sites to 250/500km.

A dedicated genus×coarse-cell local photo study restricted all original photo sites in a group to 250km max (872 taxa) or 500km max (1,761 taxa), keeping identical original species, photographed flower-colour labels, latitude/longitude/elevation, temperature and rainfall fields and source photo taxon-hash heldout folds.

## Three spatial-null assumptions and their source-verified results

All figures are POSTHOC, nested photographs and cohorts, NOT independent preregistered tests.

| Complete original photo-group max distance | Group photo label exchange constraint | Original or source site-scoring support | Refitted null p | Interpretation |
|---|---|---:|---:|---|
| 250km | Entire genus×162-cell group (may be 250km apart) | 872 | **0.010** | Positive under broad within-group exchangeability ONLY |
| 500km | Entire genus×162-cell group (may be 500km apart) | 1,761 | **0.015** | Positive under broad exchangeability ONLY |
| 250km | Only within <=50km complete-link original photo microgroups, full original cohort scoring | 872 (314 informatively swappable original photos) | **0.350** | Not robust to closer photo-geography |
| 250km | Within <=100km microgroups, full original scoring | 872 (442 swappable) | **0.175** | Not robust |
| 500km | Within <=50km microgroups, full original scoring | 1,761 (505 swappable) | **0.880** | Not robust |
| 500km | Within <=100km microgroups, full original scoring | 1,761 (746 swappable) | **0.295** | Not robust |
| 250km | Only <=50km **multispecies source-photo neighborhoods scored**, exchange within <=50km | 534 colour-blind eligible original photo species | **0.450** | Not robust even among comparable neighbourhood photos |
| 250km | Only <=100km multispecies neighbourhoods scored, exchange within <=100km | 707 eligible | **0.160** | Not robust |
| 500km | Only <=50km multispecies neighbourhoods scored, exchange within <=50km | 906 eligible | **0.910** | Not robust |
| 500km | Only <=100km multispecies neighbourhoods scored, exchange within <=100km | 1,245 eligible | **0.255** | Not robust |

The within-neighbourhood evaluation subsets were chosen only by original species ID, nominal genus, original photograph coordinates and having **at least two different source species**; the original four-state photo-colour outcome did NOT determine which site photos entered the evaluation. All models still trained on the same full 872/1761 source taxon cohorts and the same source species photo folds. In the 50/100km neighbourhood microgroups, labels were permuted while conserving each source photo microgroup's four-state class counts exactly, then the same predictive models were re-fitted 199 times.

The broad original photo Brier improvements +0.008598 (872 taxa) and +0.005115 (1,761 taxa) are still correct descriptive held-out source statistics. They do NOT have sufficient spatially constrained source-label randomization support to be claimed as rainfall-driven local pigmentation sorting.

## Mandatory claims and nonclaims

**SUPPORTED (descriptive observational):** the sampled global distribution of already classifiable original photographed flower-colour classes shows a small climatic correlation among different congeneric nominal species when geography is controlled coarsely or when train-only genus composition baselines are estimated. Within-genus rain predictors can reduce fixed heldout photo-colour prediction error, but this signal is NOT uniquely distinguishable from fine spatial photo structure in the test shown here.

**NOT SUPPORTED:** robust local precipitation–flower-colour association after constrained 50/100km matched geography. Earlier p=.010/.015 are exchangeability-assumption-dependent and must not be presented standalone as evidence. A score on colour-blind exchangeable-neighbour photos also gives no support. SoilGrids independence and original within-species local colour mismatch models are unsupported or uncertain under relevant geographic holds.

**NOT IDENTIFIABLE FROM ORIGINAL SOURCE PHOTOS:** whether rainfall selected heritable flower-pigment variants, phylogenetic evolution of local congeners, genetic population co-occurrence, species' representative colour, real-world pigment spectra, photographic colour misclassification or pollen/fitness outcomes. A genus name is not an explicit original multi-species phylogenetic tree, and a 50–100km microcluster is not a biological population.

The more local spatial null is partly conservative because it cannot exchange a different photographed colour for all source taxa; this is why the separate photo-geometry-only scoring audit was done. Neither method proves geography is the *sole* cause, but together they remove the basis for a spatially robust positive rainfall-selection claim.

## Source-verified complete actions

- Original 250/500km 199-null photographed source: [original broad null](https://github.com/zuizui0223/fcp/blob/analysis/fcp-global42111-local-congeners-20261009/results/fcp_global42111_local_congener_rain_null_20261009/result.json), source-verified in Actions 37894008208.
- Full original cohort 50/100km exact complete-link microspatial null: [actual result](https://github.com/zuizui0223/fcp/blob/analysis/fcp-global42111-microspatial-label-null-20261009/results/fcp_global42111_microspatial_photo_label_null_20261009/result.json), 6 synthetic tests + 4 × 199 matching-null refits, Actions [37943864327](https://github.com/zuizui0223/fcp/actions/runs/37943864327).
- Colour-blind geographically multispecies photo evaluation support null: [actual result](https://github.com/zuizui0223/fcp/blob/analysis/fcp-global42111-microspatial-label-null-20261009/results/fcp_global42111_microspatial_photo_subset_null_20261009/result.json), 10 synthetic tests + 4 × 199 model refits, Actions [37944581127](https://github.com/zuizui0223/fcp/actions/runs/37944581127).

### Next empirical identification gate

Future biological claims need a source-verified genus-level phylogenetic tree or at least explicit within-genus clade assignment, independently replicated flower-colour identities / photo-level quality, genuinely colocated populations and explicit ecological exposure/pollinator/pigment/fitness data. Those are NOT present in the source numeric photo-colour and local climate regression result. Do not alter original deeper 1,499-taxon analyses or untouched 2,000+730 independent taxon selections.

**Scientific decision:** freeze the present global study as an observational colour-biogeography characterization + a clear spatial-confounding sensitivity result, not as a successful test of rainfall-driven flower-colour adaptation.
