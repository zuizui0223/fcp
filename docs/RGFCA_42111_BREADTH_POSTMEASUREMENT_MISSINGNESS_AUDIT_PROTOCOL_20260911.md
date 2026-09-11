# RGFCA 42,111 species-breadth postmeasurement missingness audit protocol

Date: 2026-09-11 JST
Status: frozen while the first blind Step-8F partitions were running, before any complete 42,111-species species-colour join or atlas composition result was available.

## Purpose

The primary Step-8F estimand is the species-equal distribution of one frozen observed coarse flower-colour state per discovered species. Because image acquisition, flower localization and colour admission can fail, the classifiable subset may differ systematically from the 42,111-species denominator. This audit diagnoses that selection without redefining or rescuing the primary result.

## Frozen strata

All strata are defined from metadata frozen before the complete colour join:

1. discovery source: V1 only, V2 only, both;
2. current metadata capacity after observer cap: 0, 1, 2–4, 5–9, 10–19, >=20;
3. previously opened species flag: the earlier 1,000-species measured cohort versus species not in that cohort;
4. breadth-rank blocks: 1–5,000, 5,001–10,000, ..., 40,001–42,111;
5. first-token genus size in the 42,111-species frame: singleton versus non-singleton descriptive genus labels.

## Endpoints

For each frozen stratum report:

- denominator species;
- image-acquisition failure fraction;
- ROI/localization non-admission fraction;
- `mixed_uncertain` fraction;
- four-state classifiable fraction.

The primary selection diagnostic is the maximum absolute difference in classifiable fraction across levels within each stratum family. Report Wilson intervals descriptively. Do not choose or discard strata after viewing colour composition.

## Composition sensitivity

The primary colour composition remains the unweighted composition among classifiable anchors, with the full 42,111 denominator and missingness fraction shown alongside it.

Three prespecified sensitivity surfaces are required:

1. **No-assumption denominator bounds.** For each colour state m, with N=42,111, observed classified count k_m and total unclassified/missing count M, report `[k_m/N, (k_m+M)/N]`. These bounds make no missing-at-random assumption and are not required to sum to one across states.
2. **Inverse-probability weighting.** Fit a logistic model of classifiability using only the frozen metadata strata above. Predicted probabilities are clipped to `[0.01, 0.99]` before inverse weighting. Report the weighted four-state composition, maximum nonzero weight, and effective sample size `(sum w)^2 / sum(w^2)` among classifiable anchors.
3. **Complete-stratum standardization.** Standardize colour composition from metadata cells having at least one classifiable anchor back to their full-cell denominators. Report both the standardized composition and the fraction of all 42,111 species represented by complete standardization cells.

Neither sensitivity may replace the primary estimator. No flower-colour value may enter the classifiability model. The raw classifiable composition, bounds, IPW composition and standardized composition must be reported together rather than selecting the most favourable one.

## Hard boundaries

This audit cannot convert one photograph per species into modal species colour, species polymorphism prevalence, within-species D, C*/S*, or a population-level claim. It cannot replace failed anchors, alter palette/ROI thresholds, exclude a taxonomic/geographic stratum because its result is inconvenient, or reinterpret a missingness correction as evidence of a biological mechanism.
