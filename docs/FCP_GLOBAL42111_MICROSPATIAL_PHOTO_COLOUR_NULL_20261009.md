# Global FCP source-photographic rain–colour association: finer-spatial exchangeability audit

## Why the earlier 199-permutation result needs checking

The [source-verified local congeneric test](https://github.com/zuizui0223/fcp/pull/140) found rain-related predictive improvement among different photographed species from the same genus and original global equal-area cell, with max source-photo group diameter <=250km (872 taxa; 217 original genus×cell groups) or <=500km (1,761 taxa; 417 groups).

A label-shuffle null within each full genus×equal-area cell preserved local four-state colour frequencies and refit the models, yielding exploratory p=0.010 and p=0.015, respectively, from 199 shuffles. **However the exchangeability assumption may be too broad.** Within-group photographs can still be geographically and phylogenetically structured. Shuffling species across sites separated by hundreds of kilometres could break natural spatial gradients and produce an anti-conservative null, even though latitude/longitude/elevation are in the model.

## New original-site microgeographic exchangeability control

Do not change original 42,111 source-photo species, 18,457 classified, 18,413 climate complete, observer/photo IDs, original source WorldClim environmental values, source group membership, original 5-fold taxon-hash validation, or original climate predictor families. Only modify the NULL permitted label exchanges, not the observed model.

Test all FOUR predeclared combinations:

| Original fully source-photo-diameter cohort | Max distance of each label-exchange microgroup |
|---|---|
| <=250km, 872 original species | 50km |
| <=250km, 872 original species | 100km |
| <=500km, 1,761 original species | 50km |
| <=500km, 1,761 original species | 100km |

Within the preexisting genus×original equal-area-cell source group, sort by ORIGINAL numeric taxon ID and partition all photographed sites using a deterministic **complete-link** greedy clustering rule. A photo may enter a microgroup only if its point great-circle distance from **every** earlier assigned microgroup photo is <=50 or <=100km. This selection NEVER opens photo-colour outcomes; cluster sizes 1 are permitted.

Then report the exact original species denominator, total microgroups, microgroups with multiple species, groups with at least two different four-state photograph colours, and the original species count/fraction whose labels can actually change under a within-microgroup shuffle.

**Feasibility gate fixed in advance:** at least 100 original source species AND at least 10% of the cohort must lie in informative microgroups (>=2 different original photo-colour labels). If not, mark `HOLD_INSUFFICIENT_MICROSPATIAL_EXCHANGEABILITY` and return NO permutation p-value, rather than calling the absence of exchanges a biological null.

For an evaluable threshold, carry out exactly 199 colour-label permutations within each qualifying microgroup, conserving **both** microgroup and parent genus×cell source-photo four-state counts exactly. For every null, retrain the unchanged original within-genus×original-cell photo predictor twice: local-group+geography+BIO1/5 temperature only vs adding BIO12/15 precipitation. Save original observed Brier gain (must match +0.0085980344 at <=250km and +0.0051152608 at <=500km), null mean, spread, quantiles and source-conditional one-sided p=(1+#null>=observed)/200. Do not select favourable microgroup radius after seeing the result, do not change sample size or photographic classifications, and do not treat nested 250km and 500km cohorts as independent replicates.

## Biological limitations

This is a **sensitivity to spatial exchangeability**, not an actual species-level tree or population genetics test. Even if original Brier gain remains rare after more local swaps, closely related species within a genus may share flowering traits and inhabit correlated habitats. No within-genus phylogenetic tree, subgenus clade assignments, alleles, pigment assays or real selection/fitness evidence are present in the source photo table. One original photograph per species is not its true modal, genetically identified flower colour.

A microgeographic HOLD is a coverage/identifiability result, not evidence of no rainfall–flower-colour association. The prior unrestricted 199 shuffles remain valid only under their broader null's assumptions, not as proven phylogenetically/locally robust p-values.

The older 1,499 deeply photographed FCP H1/H2 source and untouched prospective 2,000+730 selected taxa remain unchanged. Result receipt is committed under `results/fcp_global42111_microspatial_photo_label_null_20261009/result.json` only after source and synthetic tests pass.
