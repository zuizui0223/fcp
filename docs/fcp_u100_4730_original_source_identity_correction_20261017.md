# FCP original U100=4,730 photo-source identity correction (2026-10-17)

## Status and scope

This is a **source-identity and analysis-denominator correction**, **not** a new flower-colour statistical result. The original FCP H9 source selection and all 512 already-running original ROI-v4/EfficientSAM workers remain immutable. The original 149,900 measured photographs and previously frozen 50-km, observer-disjoint D and fixed white-axis W implementations are preserved.

The original 4,730 iNaturalist taxon-ID opportunity frame yielded 4,478 source taxa with 100 original photo rows each. These produce **447,800 source-photo occurrences**, but they are **not 447,800 globally unique photo IDs** across all five cohorts, as an earlier release description incorrectly suggested.

## Independently observed collisions

Evidence: [five-cohort real-source check #39006956313](https://github.com/zuizui0223/fcp/actions/runs/39006956313), with exact taxon/photo attribution in [#39007299802](https://github.com/zuizui0223/fcp/actions/runs/39007299802). All old discovery, validation, prospective, wave4 and globally reconciled new-U100 source files were byte-SHA authenticated first.

| Original discovery source | Reused new-U100 source | Reused photo ID | Reused observation ID |
|---|---|---:|---:|
| *Pyrola elliptica*, taxon 82515 | *Melampyrum lineare*, taxon 62367 | 45948761 | — |
| *Phlomoides tuberosa*, taxon 909773 | *Astragalus cicer*, taxon 75633 | 311814882 | — |
| *Rosa multiflora*, taxon 78882 | *Rosa luciae*, taxon 1045161 | 517553374 | 287821814 |

There are three repeated photo IDs and one repeated observation ID; no source cohort shares an iNaturalist taxon ID or binomial scientific name. The collisions are **exclusively discovery versus new-U100**, not wave4 or the other historical measured cohorts.

These observations alone do not establish why iNaturalist associated the same picture with different taxonomic source records. Do not infer the underlying botanical identity or recode the flower morph from the collision.

## Conservative original-method inference denominator

The original analysis requires exactly **100 *independent, authentic* sampled photo IDs per analyzed taxon**, and species are equally weighted. Keeping 99 rows for any of the three affected taxa would violate the fixed-100 rule; inventing or selecting replacement photos after opening classification results would bias the source procedure.

The reproducible implementation in `scripts/analysis/fcp_u100_original_photo_reuse_20261017.py` therefore:

1. SHA/authenticates all original photo terminal and input source files **before** any correction; then checks the precise known discovery-only collision photo/observation IDs against the old frozen records.
2. Excludes **all 100 original ROI-derived photo rows** of each of the three affected *new* taxa **only at pooled statistical assembly**; does **not** remove or edit any original discovery photographs, rerun photos, change ROI-v4, or impute colours.
3. Fails closed on any *additional* overlap and publishes a separate analysis-eligibility ledger. The signed historical H9 4,730-source ledger is never overwritten.

| Original study accounting | Taxa | Source/analysis photo rows |
|---|---:|---:|
| Immutable previous measured 3 cohorts | 1,499 | 149,900 |
| Original new wave4 source | 498 | 49,800 |
| Newly H9 complete U100 **before** collision check | 2,481 | 248,100 |
| New U100 taxa excluded from scientific pooling | **3** | **300** |
| **Original fixed-100 independent analysis eligible** | **4,475** | **447,500** |
| Original 4,730 taxa without 100-photo H9 source | 252 | Not analyzed |
| Additional photo-source identity-excluded taxa | 3 | Not analyzed |
| **Original 4,730 source-frame taxa accounted for** | **4,730** | — |

Distinct IDs in uncorrected five-cohort metadata: **447,797**. The **447,500** output is the more conservative, complete-100-taxa, nonduplicated **analysis-eligible photo denominator**, not an assertion that these photographs have all completed ROI classification. Image completion and the exact original 50-km, H1/D, H2/W results still require both genuine 256-partition worker pipelines to finish.

Source release (unchanged, preserving the 447,800 raw-source occurrences): [U100 source H9 release](https://github.com/zuizui0223/fcp/releases/tag/fcp-original-U100-4730-reconciled-source-20261017). Interpret the older release's word “unique” as **within its original selection ledger**, not as a verified statement of zero duplicates against the old 149,900 measured discovery photos; the full five-cohort audit above supersedes that broader claim.

## Scientific nonclaims

This is quality assurance for observational photo availability, not a measurement of genetic or population colour polymorphism, an independent data provider, or a new ecological effect. No larger spatial-depletion effect, p-value, causal association, change in flower morph classification, or reconstructed 100-photo sample is claimed until the original complete image analysis produces verified statistical outputs.
