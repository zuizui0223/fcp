# FCP global geography synthesis — tested spatial regularity versus unsupported universal gradient
## 2026-10-08, post-outcome cross-PR decision ledger

### One-sentence ecological result

**Across many photographed plant species, species-wide flower-colour variation is recurrently geographically structured at the neighbourhood scale, including after observed flowering month and year are conditioned on. The same data do not recover a universal signed white-colour or colour-diversity cline along absolute latitude or elevation.**

This contrasts two *distinct* hypotheses, not two contradictory tests:

1. **Spatial allocation of within-species variation:** Do photos of the same species within 50 km have more similar visible flower colours than expected from the species' exact overall sampled colour composition?
2. **Universal macrogeographic direction:** Does one altitude or latitude band consistently show more white flowers or higher within-species floral-colour diversity relative to the same species' full sampled range?

### Primary source-backed evidence

The sampled resource contains 149,900 photographs from 1,499 nominal species, split 500 discovery / 500 validation / 499 third. Under the fixed >=40-classifiable species gates, 369 / 363 / 377 species are evaluable for D. High-depth sampled species are selected by repeat-photographic opportunity, not a global random census.

#### Spatial organization positive (main FCP result)

The first-frozen post-outcome 50-km local-depletion test produced species counts 166/181/204 and mean four-state pairwise colour-discordance depletion **0.020529/0.018673/0.014685**, each permutation P=0.005. All three cohorts remained spatially structured after discarding white observations and restricting different observers, and through pre-existing quarter-stratified and cross-year local-pair tests. These are *photographic neighbourhoods*, not verified mating populations.

#### New exact photo-month conditioning positive (PR #129)

Successful Actions workflow run **37736541255** recovered all three source-SHA-pinned measured tables and passed synthetic tests and all output gates. Source reporting receipt:
results/polymorphism_specieswide_space_season_posthoc_20261008/result.json
Artifact ID 11532736238, digest SHA256 942e25eb9282074cceb1ec9b69f5442f5aae00c24f955bec36a412a0b45a01dc.

| Cohort | Photos of one species within 50 km: observed-minus-month-fixed expected local depletion | 95% species-bootstrap interval | Geographically evaluable species | Month-null identifiable species | P |
|---|---:|---:|---:|---:|---:|
| Discovery | +0.017455 | +0.008467 to +0.027287 | 165 | 140 | 0.005 |
| Validation | +0.015287 | +0.006318 to +0.025659 | 180 | 152 | 0.005 |
| Third | +0.013324 | +0.006078 to +0.021121 | 203 | 174 | 0.005 |

The same-month restricted null shuffles colour labels only **within the same species and calendar month**. It preserves each species' number of photos at each location and the exact month-by-visible-colour composition. Both all local pairs and different-observer local pairs pass positive effect and the prospective-to-this-follow-up >=80-identifiable-species threshold in each species-disjoint cohort.

The stricter **year × month** restricted null also retains positive geographic residual mean +0.012092 / +0.008996 / +0.008761 (all pairs; all P=0.005), although only 92/108/122 species allow null identification in that mode. These are technical calendar restrictions on *photo-observed dates* and may not remove within-month phenology, exposure or repeated-site observation confounding.

Nonexchangeable photo-month species remain in the all-geographically-evaluable equal-species denominator with **zero identified additional effect**, reported separately from the identifiable-only sensitivity. This is an explicit conservative *estimand convention* and not an observation of biological zero.

#### New broad latitude/elevation bands show NO replicated signed global cline (PR #127)

Successful workflow **37739944082** executed 12 synthetic tests and retrieved three immutable photo tables plus a checksum-pinned WorldClim 2.1 10-minute elevation raster. Source receipt on the PR #127 branch:
results/polymorphism_global_geography_bands_posthoc_20261008/result.json
Artifact ID 11533940496, digest SHA256 c3a9a6d9f859fe8fceb566c4b6593b67c8046f8c7bc5d0ca145a0b8676d3fe0c.

This analysis tests latitude bands 0–15°, 15–30°, 30–45°, 45–60°, 60–90° and elevation <250 m, 250–1,000 m, 1,000–2,000 m and >=2,000 m. Species must contribute >=8 photos inside each bin and >=8 elsewhere. Their total photo colour-label composition and exact geography are fixed under 199 within-species permutations, with max-T adjustment over eligible bins within each axis/cohort.

**No band met the cross-three-cohort, >=25-within-species-informative-species, same-sign, max-T-corrected support rule for either white-photo-frequency shift or four-state diversity shift.** Tropical 0–15° and polar/high-latitude 60–90° contrasts fail the 25-species coverage rule in at least one cohort. Some bins show the same exploratory sign but cannot support a universal association: e.g., >=2,000 m D shifts +0.008398/+0.022360/+0.001903, while 15–30° white-photo-fraction shifts +0.013746/+0.018352/+0.002404. One-cohort significant adjusted P values cannot rescue generality.

This is not an equivalence test. A failure to detect a **universal** gradient cannot rule out important local slopes, turnover, nonlinear transitions at untested scales or clade-specific patterns. The geographic altitude layer is ~10-arc-minute resolution and is not field-measured elevation. The white-photo classifier is known to be exposure-coupled.

### Scientific integration

The positive result is **scale-conditional spatial organization** rather than the existence of one global direction of floral whitening or colour diversity at particular elevations/latitudes. The observed spatial clustering could be generated by lineage-specific histories, seed/pollen dispersal, local environments, floral phenology, plasticity or adaptation; it does NOT establish which mechanism dominates, and geographic correlations larger than temperature correlations do NOT prove drift dominates local natural selection.

This is scientifically different from a floristic claim that all higher-altitude flowers are more chromatic, and also different from a genetic statement that segregating flower-colour alleles are maintained within a population. Those stronger claims are outside the photograph data.

The synthetic causal counterexamples in PR #128 and external literature are used only to bound interpretation: Moricandia arvensis can switch its flower colour with season without a segregating colour allele, and Linanthus parryae demonstrates local colour differentiation with experimental divergent fitness without a matching neutral-marker cline. Neither organism defines a universal mechanism for the current 1,499-species photo opportunity universe.

### Suitable central claim and hard nonclaims

**Supported now:** In three species-disjoint resources from one photographic and classification system, flower-colour state diversity is geographically allocated nonrandomly within species. The 50-km clustering persists after preserving observed calendar-month and year×month colour composition, and after observer separation.

**Not supported now:** A global high-altitude/high-absolute-latitude increase in chromatic flower states; an elevation/latitude band with reproducibly higher local colour diversity; global frequency of actual genetically polymorphic flowering species; flower-pigment loss/gain history; direct pollinator/UV fitness mechanisms.

**Scientific scope:** The observed generality is across many separately sampled species *within a common photographic opportunity frame*, not an independent-source world flora census. The 42,111-species metadata discovery frame is an opportunity frame, not the biological prevalence denominator. Repeated mapping of a common species-level estimand is the contribution, while sources and fitness claims are kept separate.

### Actions and manuscript boundary

- Do not overwrite frozen New Phytologist numerical evidence, original H1/H2 interpretation or manuscript decision ledgers.
- PR #129 is the living month-conditioned geography evidence route; PR #127 holds the directional and banded latitude/elevation route.
- Both are **post-outcome exploratory analyses** with real SHA-frozen photo source and successful Actions execution. Results should be distinguished from the separately frozen prospective phenotype-space H2.
- Next geographic program after this result should test spatial *scale and region heterogeneity* with verified same-species photo opportunity, rather than searching for a favourable high-latitude colour rule post hoc.
