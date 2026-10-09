# FCP original local congeneric moisture hypothesis: photo-neighborhood-only evaluation

## Critical result motivating this follow-up

The originally reported apparent rainfall association among different local congeneric photographed species was based on original 250-km and 500-km maximal source-site groups, with source-conditional permutation p=.010/.015 if flower-colour labels were shuffled across the ENTIRE original genus×equal-area-cell group.

The more defensible [complete-link 50/100km microgeographic photo exchangeability audit](../results/fcp_global42111_microspatial_photo_label_null_20261009/result.json) (source-verified Actions 37943864327) refitted 199 nulls at all four fixed combinations:

| Original max source group diameter | Original photographed species | Micro-site allowed exchange diameter | Original species with >=2 source flower-colour labels in exchange groups | Source conditional permutation p |
|---|---:|---:|---:|---:|
| 250km | 872 | 50km | 314 (36.0%) | **0.350** |
| 250km | 872 | 100km | 442 (50.7%) | **0.175** |
| 500km | 1,761 | 50km | 505 (28.7%) | **0.880** |
| 500km | 1,761 | 100km | 746 (42.4%) | **0.295** |

**The old p=.010/.015 does not survive an ecologically better geographically conditional exchangeability null.** This is a crucial robustness failure and should supersede claims of strongly calibrated local source rainfall-photo-colour associations.

However, in the 50km microgroups, a large fraction of original source photos cannot exchange a different original four-state colour. The new null is partly anchored by unexchangeable original observations, and may have limited power. It is not proof that geography is the sole cause of the photographed-colour pattern.

## Next, outcome-blind scoring of only genuinely comparable sites

To diagnose dilution by unexchangeable site structure, calculate the ORIGINAL heldout Brier gain **for every species** using the EXACT previously frozen models, source taxa and taxon-ID-only fivefold plan, then score it twice:

1. Across ALL frozen original photos in the 872/1,761 taxon sample, reproducing +0.0085980 / +0.0051153 original source total moisture gain.
2. On the fixed subset of photo taxa whose source microgeographic 50/100km subgroup has **at least two distinct original taxa**, without looking at photo-colour values. This does not require >=2 different photo-colour outcomes, which would make the evaluation selection response-dependent.

For each micro-site definition, preserve exact photo metadata, original four-colour labels, source genus×cell frequencies, taxon-ID-defined photo folds and fixed train-only genus-cell model fitting across ALL original photo cohorts. Swap labels only within fixed complete-link <=50 or <=100km source neighbourhoods and repeat 199 whole-cohort training/prediction comparisons. Evaluate the already-frozen near-multispecies photo IDs only, for both observed and permuted models; reject/mark HOLD if there are fewer than 100 original site-matched photo taxa or <10% of cohort in this geometry-only evaluation subset.

This is an additional **POSTHOC estimand** motivated by the failed spatial null. It must never replace the full original source result if its p-value happens to be smaller. Report all 4 predefined source cohort/micro-radius combinations, compare results explicitly to the earlier 0.350/0.175/0.880/0.295 p-values, and do not claim phylogenetic adjustment or local adaptation.

The source-photo geographic cells remain broad, and this method still does not know true intra-genus phylogeny, photographer bias, actual flowering state, genetic colour variants, pollination, pigment chemistry or reproductive success. The older 1,499 deep-photo analysis and untouched prospective 2,000+730 species are not altered.
