# FCP 10-variable environmental verification — 2026-10-14

**Status:** post-outcome exploratory sensitivity, not new confirmatory H2 or E1–E3 tests. Source data and 2026-09-25 three-factor signed environmental results are frozen. No population-level allele, fitness, pigment-chemistry or local-adaptation inference.

## Question

Do within-species photographed white-versus-nonwhite differences respond consistently to more than the original three broad environment covariates? Is any putative common association still visible after observer and within-location comparisons?

## Fixed analysis before opening these additional climate results

WorldClim 2.1 long-term 10 arc-minute baseline values are read from the existing SHA256-frozen release `fcp-worldclim-2.1-10m-20260925`. The measured prospective photo data use SHA256 `57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186`; frozen technical join, near-clipping values and high-clip exclusion are reused without threshold fitting.

Ten abiotic variables are grouped *descriptively*, **without claiming 10 independent causal predictors**:

| Family | Covariate | Interpretation |
|---|---|---|
| Temperature | BIO1 | Annual mean temperature |
| Temperature | BIO4 | Temperature seasonality (standard deviation ×100) |
| Temperature | BIO5 | Maximum temperature of warmest month (old E1) |
| Temperature | BIO6 | Minimum temperature of coldest month |
| Temperature | BIO7 | Annual temperature range, BIO5–BIO6 |
| Moisture | BIO12 | Annual precipitation |
| Moisture | BIO14 | Precipitation in driest month (old E2) |
| Moisture | BIO15 | Precipitation seasonality (coefficient of variation) |
| Moisture | BIO17 | Precipitation in driest quarter |
| Radiation | Monthly SRAD mean | Long-term solar radiation (old E3) |

All covariates are standardized within species **before** removing response-blind high-clip photographs, exactly as in the original 2026-09-25 main species-contrast estimator. Minimum 5 white and 5 nonwhite classifiable, technically available photographs are required. In the 2026-09-25 source, 281 prospective species met this gate. The new workflow *must* reproduce the original BIO5, BIO14, and SRAD species-level means, medians and denominators at numeric tolerance `1e-10` before any new outputs are accepted.

## Tests, diagnostics, and inference ceilings

1. **Exploratory cross-species sign screen:** one equal-species signed white-minus-nonwhite mean-standardized contrast for each climate covariate, two-sided Wilcoxon test. Apply **Holm correction across all ten** variables whether or not a result passes. Report species count, effect median, and bootstrap 95% interval for the median (2,000 fixed-seed resamples). A positive screen indicates ecological *association*, not selection.
2. **Observer-balanced sensitivity:** average per species × observer × colour, then compute a between-colour difference with equal observer contribution; require at least 3 distinct observers photographing each state. The species set can differ from the original and cannot retroactively rescue the original E1–E3 gate.
3. **Strict within-observer sensitivity:** include only observers photographing *both* colour groups of the same species. Compute the observer-specific difference before averaging within species, and report the number of identifiable species. This estimates a narrower within-observer effect, not the same broad geographic environmental gradient.
4. **Geographically matched sensitivity:** within a species, identify **1-degree latitude–longitude cells** containing both white and nonwhite records. Compute the within-cell standardized environmental contrast and take an equal-cell average for species with at least 2 mixed cells. At WorldClim 10′ resolution, many photographic cells will carry identical climate or no mixed colours: record this loss of support rather than assigning zero.
5. **Measurement and opportunity audit:** compare the proportion and signs of original/observer/local contrasts, retained photograph and species counts, near-clipping availability, source hashes and climate-factor collinearity. Because colours and photo observers are geographically clustered, no posthoc negative-control success is claimed without independent source/background data.

The original 2026-09-25 test was prospectively frozen for **BIO5, BIO14, SRAD only**; those p-values remain authoritative and separate from this ten-variable screen. All ten outcomes are **exploratory and follow exposure of the photo-derived response**. If they pass the Holm screen but fail observer/local sensitivities, a general abiotic colour-selection claim is **not supported**.

## Real input replay

`fcp_pipeline` is not modified by the new analysis. The read-only engine lives in `scripts/analysis/u100_multiclimate_20261014.py`, and the current claim-guard workflow performs an immutable-data smoke after recovering raster and technical inputs. Original status/threshold rules are not refit.

```bash
python scripts/analysis/u100_multiclimate_20261014.py \
  --measured /trusted/polymorphism_h2_third_cohort_measured_photos_v1.csv \
  --technical /trusted/highlight/technical_table.csv.gz \
  --join-key /trusted/highlight/sealed_join_key.csv \
  --high-ids /trusted/highlight/high_clip_ids.csv \
  --bio-dir /trusted/worldclim/bio \
  --srad-dir /trusted/worldclim/srad \
  --frozen results/polymorphism_white_environment_mechanism_20260925/result.json \
  --outdir /trusted/analysis/u100_multiclimate
```

Output `multi_climate_species_effects.csv`, `multi_climate_summary.json`, `multi_climate_receipt.json`. No photo pixels are acquired, saved or distributed. Old frozen H2, old E1–E3, original 50km null p-values and 4,730-species opportunity frame are unchanged.

## Next biological priorities

Report effect-size heterogeneity and observer/local identification failures rather than selecting whichever BIO covariate is significant. To establish any general *environment × within-species colour* rule, repeat in outcome-blind species/photo cohorts with independent colour/background segmentation, local environmental availability models and phylogenetic comparative adjustment; if environmental effect sizes are genuinely consistent, use interaction experiments, genotype and fitness data to distinguish selection from plasticity, demographic history and measurement artefacts.
