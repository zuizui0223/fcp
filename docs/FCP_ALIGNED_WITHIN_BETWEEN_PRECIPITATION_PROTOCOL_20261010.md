# Original photo-aligned within- versus between-species rain signal: protocol (2026-10-10)

## Motivating gap and denominator

Different original studies identified (i) a modest between-congeneric-species precipitation prediction effect among 13,307 original one-photo species and (ii) a modest within-species environmental gain among 20,546 different original photo-site observations from 5,127 taxa. These are **not the same observation set nor the same training information**. Thus the two positive source coefficients cannot simply be called a 'shared driver' even when both include BIO12/15.

## Pre-source-execution question

For exactly the SAME original four-state classified, climate-complete, real-site original taxon×geographic-cell photographs among the existing **85,337 original taxon-cell** and exact **100,543 expanded original photo-site** frames, does adding exactly the SAME **annual precipitation BIO12 plus precipitation seasonality BIO15** improve out-of-original-geographic-cell heldout four-state Brier predictions at both conditional scales?

- **Within-species level:** in geographical cell holdout, use source *same-species* photo colour mean from OTHER training cells; regress deviations of colour and all features from train-only nominal species centroids across training photos, then apply those coefficients to the heldout source photo. This already-established source predictor explicitly uses photos of the same named species, but NEVER in the same heldout geographic cell.
- **Between-species-within-genus level:** in the same geographic fold, use *only different named training species of the source genus* and estimate equal-source-species genus centroids + slopes on species mean colour and environment. For each heldout source photo, the source focal **entire species** is ineligible for this model's training (not just its heldout photo). Require at least 3 other genus species in the between source train. Species-intercept within train can still see the focal species in other cells; the differential training-information sets are an explicit limitation and are not exchangeable.
- **Matched evaluation:** source photographs must be simultaneously evaluable in both models. Five fixed whole-original-geo-cell folds × five source-taxon-hash folds restrict exactly the original heldout photo IDs. Record the common photo count, species count and genera count. If <500 heldout photos or <100 original species, **HOLD** and do not report a strong sign.
- **Matched predictors, outcome, and contrast:** both use the same 4 categorical original photos and original real-site features. Baseline both levels uses the same fixed source photograph absolute latitude, longitude sine/cosine, true site elevation plus **BIO1+BIO5**. Compared alternative includes those plus **BIO12+BIO15** only. No invented colour ordering or white/nonwhite summary.
- **Matched metric:** original photo multiclass Brier score, precipitation increment `Brier(geo+temperature)-Brier(geo+temperature+precipitation)`. Differences in incremental Brier computed on each identical heldout original photo, with separate species-, genus-, and original-geo-cell cluster bootstrap intervals (999 fixed-resample draws) for within gain, between gain, and within-minus-between gain. Report all three and any negative/zero intervals.
- Original classifiable photo-source counts, 85,337 originals, 39,075 classifiable, 38,968 climate complete, and exactly matching original source photo-taxonomic/photo ID triple checked. All source measurement and phenotype labels immutable.

## Strict scientific firewall

1. This is exploratory source matched-photograph *predictive* comparison. A positive within and between Brier increment on the same photographs **does not establish identical causal ecological selection**, shared genetic pigmentation or an individual-level fitness consequence.
2. The ***between*** model sees only source species other than the target species. The ***within*** model sees the target species in other original photo cells. Those different training opportunities and target identification estimands mean Brier effect size comparison is not a direct comparison of common regression coefficients or mechanistic effect strength.
3. The fixed photogeography terms are **not** an actual lowrank or full geodesic spatial random field. The source 50/100km fine-geographic label-shuffle null and true dated LCVP ancestral covariance tests remain separate stronger falsification gates and have not been reversed by this analysis.
4. FCP original 42,111 one-photo species, 18,457 original classifiable one-photo species, and 149,900 high-depth 1,499 species are **separate denominators**, never pooled or redefined as this common-photo sample.
5. Species-centred analysis absorbs named species mean ancestry, not a joint phylogenetic covariance. Between-genus model only bounds broad genus identity, not patristic distances.
6. Original photos lack two-person independent expert floral-organ/colour adjudication, original heritable morph genotypes, actual local fitness and pollinator links. This cannot prove adaptive colour evolution.
7. Main and primary H1/H2, historical original source results and untouched future 2,000+730 taxa remain unchanged. Fail closed if full source photos or exact outcome counts are absent.

**Status at protocol creation:** code, synthetic tests and inference gate completed; real matched original photo result not yet observed.
