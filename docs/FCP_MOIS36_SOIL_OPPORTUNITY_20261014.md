# FCP × mois: 36-variable climate/soil opportunity and within-species sensitivity

**Date:** 2026-10-14. **Status:** source-hash-checked, post hoc exploratory; **no confirmed soil-colour effect**. FCP within-species measured photographs are not the same source population as `mois` global between-species photographs. This companion analysis does not amend New Phytologist H1/H2 results or any signed macroecological coefficient in `mois`.

## Exact input lineage

1. **FCP:** immutable prospective 49,900-photo table from commit `7e538e5c51c05a7cc47b2fcf53eea92634c8a863`, `data/derived/polymorphism_h2_third_cohort_measured_photos_v1.csv` SHA256 `57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186` (499 species). Three preexisting technical inputs: `archive/fcp_submission_20260925/highlight/{technical_table.csv.gz,sealed_join_key.csv,high_clip_ids.csv}` and their tracked Git blob hashes. Frozen exposure eligibility: 12,583 photo rows in 281 species with >=5 white and >=5 nonwhite classifiable photographs, exact response-blind high-clip exclusion and near-clip availability.
2. **mois:** original environmental ZIP from historical FCP GitHub Actions artifact `11653183273`, SHA256 `6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca`. Its member `fcp42111_original_photo_all_climate_soil/all_original_100543_photo_sites_expanded_environment.csv.gz` preserves **100,543 unique photo/site IDs, 42,111 original source taxa, site coordinates and 36 real abiotic attributes**. Variable names, family blocks and source/selection limitations reproduce [`mois` methods](https://github.com/zuizui0223/mois/blob/main/scripts/analysis/fcp_abiotic_hue_direction_macro_20261013.py). That source is a sampled **reference set of environmental sites**, not a downloadable wall-to-wall raster for the high-depth FCP observations.

## Estimand and filters

- Compare visible-photo **white versus pooled nonwhite** environmental values **within** each FCP species. Per species, a variable is estimable only with >=5 white and >=5 nonwhite valid site-proxy values, positive within-species SD, and correct frozen classifiability/technical gates. The response is `mean(z_environment | white) − mean(z_environment | nonwhite)`; give species equal weight. Photograph categories are not genetic morph states.
- Join exact photos by **photo ID and matching taxon ID**, with coordinates checked to 100 m. Otherwise choose the **nearest original `mois` environment source site before considering soil completeness**, by 3-D spherical geodesic distance. The source site's original coordinate uncertainty AND the target photograph's positional uncertainty must each be <=1,000 m. Results are separated into **exact identity**, nearest source <=1 km and <=5 km. Proximity/soil coverage is never retrospectively optimized to a colour result; a nearest site missing soil remains missing.
- The 26 climatic/topographic axes are WorldClim elevation, BIO1–BIO19, solar annual mean/monthly CV, wind annual mean/monthly CV, atmospheric vapour pressure annual mean/monthly CV. The **10 soil proxies** are pH, soil organic carbon, total N, clay, source `available_water_proxy`, source CEC, sand, silt, bulk density, and coarse fragments (`soil_cfvo`). `available_water_proxy` is *not* a verified plant-available water-capacity field measurement. Source raw variables are retained, not silently converted into new units.
- For each of three bands, report all 36 variables, a two-sided Wilcoxon across equal-species deltas if at least **30 eligible species**, and **Holm correction over 36 variables including non-estimable tests as p=1**. Also compute all available same-observer and within-one-degree-grid-cell colour contrasts, separately Holm-corrected (1 degree is a *coarse negative control*, not local population identity). Bands, methods and tests are post hoc; no multiple-band familywise selection correction is claimed. Those controls are descriptive and cannot rescue the main null.

## Frozen empirical result (correct photo-position accuracy required)

| Site evidence | photos mapped | named species with >=1 mapped photo | complete 10-soil rows | species with >=5 per colour for soil pH |
|---|---:|---:|---:|---:|
| Same photo ID | 195 | 104 | 122 | **0** |
| Within 1 km | 3,235 | 274 | 1,867 | 28 (below 30-species test gate) |
| Within 5 km | 6,530 | 281 | 3,834 | 94 |

- All **36** fixed environmental-variable contrasts failed Holm p < .05 in both 1 km and 5 km proxy bands; exact-photo contrast has **0** estimable variables under the >=30-species gate. Three original frozen environmental-effect findings in `fcp` are **not** recalculated as confirmations by this nearest-site design.
- **Most extreme 5-km climate association:** BIO8, n=160 species, median white-minus-nonwhite +0.092 SD, nominal p=0.00711, Holm p=0.256. In 1 km, smallest nominal p is BIO3, n=57, nominal p=0.0502; Holm p=1.00.
- **Soil `cfvo` at 5 km:** 94 species, median −0.076 SD, nominal p=0.0300, Holm p=1.00. Soil pH: 5-km n=94, no firm direction. **Within-observer eligible species are no more than one per 5-km axis**, so observer-controlled soil effects are **not estimable**. No same-one-degree-cell factor survives Holm across the 36 variables.
- The earlier proximity analysis lacking an explicit target-photo positional-accuracy cap is **superseded**. Its stronger `soil_cfvo` nominal signal and counts (5 km 6,891 photos; 108 soil-pH species) must not be reported as primary evidence.

## Scientific decision and next defensible escalation

> `NOT_SUPPORTED_UNDER_EXACT_AND_ACCURACY_SCREENED_SITE_PROXY_TESTS`: Neither exact-photo matching nor 1/5-km nearest-environment-site *proxies* justify a recurrent within-species soil direction in the photographed FCP species.

A nearest source record can represent another species, season, slope, land cover or soil parent material. Even a 1-km transfer of SoilGrids values is an **environmental proxy**, not proof of measured environment for the target specimen. Repeated target photographs may map to the same proxy site, reducing independent information. BIO layers are long-term climate, not weather at flowering. Spatial and observer associations persist, and outcomes were opened before designing this analysis.

**The next valid test requires exact-coordinate sampling** of the original SoilGrids rasters for *all FCP high-depth photograph locations*, an independently versioned source raster/scale/checksum and coverage ledger, within-species matched observed versus available habitat, plus observer- and background-aware spatial controls. Retain multiple imputation or accessibility-weighted sensitivity only if specified before looking at new outcomes. A positive historical `mois` between-species white/moisture coefficient does not transfer causally to these within-species photos.

## Reproduce (the original archives must be supplied)

```bash
python -m pip install numpy pandas scipy pytest
python -m pytest tests/test_u100_mois_soil_transfer_20261014.py -q
python scripts/analysis/u100_mois_soil_transfer_20261014.py \
  --source-archive /trusted/original_100543_all36_artifact.zip \
  --prospective /trusted/polymorphism_h2_third_cohort_measured_photos_v1.csv \
  --technical archive/fcp_submission_20260925/highlight/technical_table.csv.gz \
  --join-key archive/fcp_submission_20260925/highlight/sealed_join_key.csv \
  --high-clip archive/fcp_submission_20260925/highlight/high_clip_ids.csv \
  --outdir results/u100_mois36_soil_transfer_20261014
```

GitHub Actions in the FCP claim guard verifies original archive/table hashes, fixed eligibility and the reproduced counts, executes this source-backed diagnostic and uploads full per-species results for audit. The repo keeps the small result table + receipt and does not check in original photographs, near-site source IDs or raw location records.
