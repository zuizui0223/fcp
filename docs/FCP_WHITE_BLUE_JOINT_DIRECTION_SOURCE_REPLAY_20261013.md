# FCP: original-source white-minus-blue/purple rainfall direction — posthoc joint bootstrap (2026-10-13)

**Scientific status:** post-outcome, same iNaturalist source, **draft only**. Parent branch `analysis/fcp-shared-genus-ecoregion-transport-20261011` and frozen New Phytologist H1/H2 remain unchanged. This extends PR #164 with an opt-in derived contrast, not a new phenotype, independent cohort, causal experiment or prospective test.

## Ecological question
Does the **slope of white photographic colour probability relative to blue/purple** along annual rainfall BIO12 differ between the source-photographed **90–180° E** sector and all other longitude sectors, within genera represented by distinct species on both sides?

Define `C = (beta_white,E - beta_blue_purple,E) - (beta_white,other - beta_blue_purple,other)`. Source slope units are **photographed four-state class probability per ORIGINAL global one SD BIO12 (~736 mm)**. `C` is a difference of marginal probability slopes in four mutually exclusive categories, NOT a two-colour conditional odds ratio, pigmentation biochemistry, selection gradient, fitness, or within-population genetic morph frequency.

The source universe contains 100,543 unique originally photographed IDs / 42,111 taxon names, including 45,953 four-state-classifiable photo hues; true photo-site climate + classified hue are available for 45,836 source photos / 25,028 taxa. The dominant source-matched comparison has **394 genera, 9,938 species, 21,176 classified photos**. The source positions, hue labels, genus matches and country/longitude grouping are held unchanged.

## Actual original-source independent local replay
Downloaded original GitHub Actions source ZIP artifact `11587115220` and verified the exact source ZIP SHA256:
`2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0`.

The separate independently written source reader + species aggregation + genus-sufficient-statistics bootstrap recalculated all **16 original-eligible sample × minimum distinct-species × weighting** comparisons using the parent seeds and **999 paired genus draws**. For **each condition**, the original four hue point estimates, east/rest estimates, each marginal 95% percentile interval and the original bootstrap fraction were matched to the sealed previous parent artifact `11690391164` (run `38106584765`), with maximum error below 1e-9.

### Exact jointly resampled contrast (values in **percentage points per BIO12 SD**)

| Source photo condition | Threshold species per sector per genus | Weight | n genera | C point | Joint genus bootstrap 95% CI | Six hue-pair exploratory maxT simultaneous ~95% CI |
|---|---:|---|---:|---:|---|---|
| **All classified** | **3** | **species-equal** | **394** | **+7.628** | **[+3.613,+11.388]** | **[+2.517,+12.739]** |
| **All classified** | **3** | **genus-equal** | **394** | **+5.309** | **[+1.101,+9.618]** | **[-0.156,+10.775]** |
| All classified | 5 | species-equal | 199 | +8.811 | [+3.858,+13.608] | [+2.429,+15.194] |
| All classified | 5 | genus-equal | 199 | +6.071 | [+0.938,+11.482] | [-0.872,+13.015] |
| All classified | 10 | species-equal | 72 | +9.949 | [+3.696,+15.556] | [+2.305,+17.593] |
| All classified | 10 | genus-equal | 72 | +6.156 | [-0.700,+13.014] | [-3.108,+15.420] |
| Common observed rainfall, 356–1,832mm | 3 | species-equal | 275 | +18.665 | [+8.485,+28.731] | [+5.641,+31.688] |
| Common observed rainfall, 356–1,832mm | 3 | genus-equal | 275 | +10.847 | [+0.145,+21.840] | [-3.205,+24.898] |
| Common rainfall | 5 | species-equal | 128 | +24.111 | [+10.406,+37.069] | [+6.867,+41.355] |
| Common rainfall | 5 | genus-equal | 128 | +22.843 | [+7.459,+37.019] | [+3.305,+42.381] |
| Common rainfall | 10 | species-equal | 48 | +24.395 | [+8.081,+42.043] | [+2.745,+46.045] |
| Common rainfall | 10 | genus-equal | 48 | +23.953 | [+5.854,+43.904] | [-0.727,+48.634] |
| **Strict original photo quality** | **3** | **species-equal** | **77** | **+9.107** | **[-1.361,+18.040]** | **[-3.434,+21.647]** |
| **Strict original photo quality** | **3** | **genus-equal** | **77** | **+11.938** | **[-1.846,+23.625]** | **[-4.244,+28.119]** |
| Strict original quality | 5 | species-equal | 31 | +12.919 | [+3.039,+24.146] | [-1.238,+27.077] |
| Strict original quality | 5 | genus-equal | 31 | +17.209 | [+3.697,+28.810] | [-0.492,+34.910] |

Strict-photo minimum 10 taxa per side per genus does not meet the parent **>=30 shared genera** feasibility gate: **HOLD, no result**. Conditions share source images and are NOT independent replications. Common-rainfall restriction and photo-quality restriction change the evaluated taxa/genera.

**Main interpretation:** Exact joint paired-genus bootstrap improves on the previous algebra-only comparison and supports a positive *photo-colour coefficient contrast in the already selected primary data*. The original-frozen species-equal versus genus-equal weighting and the six-contrast multiplicity sensitivity do NOT yield the same support: only species-equal all-photo retains a positive interval under the within-region six-pair maxT correction. On strict-photo primary support, both weighting methods' nominal 95% CIs include zero. Do not promote this to a universal humid-white/purple-loss rule.

## What the extra uncertainty calculations mean
For four photographed hue labels, there are **six unordered hue-pair contrasts**. The source code evaluates all six using the **same matched genus bootstrap draws**, calculates within-pair SD, then applies the 95th percentile of the maximum absolute *point-centered standardized bootstrap deviation* across the six. This is a descriptive **approximate simultaneous confidence envelope for the selected eastern sector**, not a selection-adjusted CI over four alternative sectors, many models, photo filters and potential post-outcome traits. The ordinary joint percentile CIs above also do not adjust for those selections.

A second covariance-free envelope is derivable from the archived hue-marginal 95% percentile intervals: `[white_lower - blue_upper, white_upper - blue_lower]` contains at least 90% of paired bootstrap draws under the empirical percentile interpretation (union bound), and the lower endpoint is individually a conservative one-sided 95% bootstrap quantile bound. That envelope is **not** an exact 95% two-sided joint CI. The paired-source percentile CIs here replace it.

The printed `bootstrap_positive_fraction` is a **resampling frequency conditional on the observed data**, NOT a permutation-null or prospective hypothesis-testing p-value.

## Inferential ceiling / next falsifiers
- Both **90–180 E** and the white–blue hue comparison were emphasized **after seeing source outcomes**, creating post-selection risk. This contrast cannot be described as preregistered, confirmatory, independent ecological replication or FWER-significant among all possible regions/colours.
- The source photo white classifier has measured image-exposure coupling; photo quality and classifiability selection are region-dependent; original 54,590 four-state-unclassifiable photos retain unknown hue. MNAR and original flower-organ identity remain open.
- Equal genus membership is not direct dated phylogeny. Spatial structure, genuine source locality and taxonomic turnover remain potential alternatives. The same-name species east-versus-rest comparison of white frequencies did **not** confirm a common positive rainfall slope.
- No genotypes, colour chemistry, exact pollinator identities, fecundity or relative colour-morph fitness were measured.
- Independent follow-up must pre-fix *both* region and hue direction before opening **previously unexposed taxa or genuinely independent photo/field sources**, and audit flower organs, missingness and shared multidimensional environmental support.

## Executable relationship
The legacy script's return and numerical results are unchanged by default. Opting in with `--posthoc-white-blue` adds a `posthoc_white_blue_region_rain_contrast` object to each eligible pairwise whole-genus bootstrap, without refitting the source models or changing genus resamples. The original four-colour keys remain byte-for-byte unmodified numerically.

Original primary analysis: [PR #164](https://github.com/zuizui0223/fcp/pull/164).
