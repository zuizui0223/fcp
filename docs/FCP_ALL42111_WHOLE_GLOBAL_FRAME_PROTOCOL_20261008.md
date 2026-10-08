# Global FCP: all 42,111 discovered plant species, no high-depth exclusion
## Whole-world sampling frame and nested inference protocol — 2026-10-08

### Scope decision

**Every one of the 42,111 previously discovered iNaturalist plant species is in scope.** Earlier analyses' 4,730 species are NOT a biological upper limit; they were a restrictive historical photographic opportunity subset with >=100 observer-capped eligible source records. This study explicitly does not restrict its sample to those 4,730 species. The other 37,381 species remain in the frame with their actual observation opportunities, including 789 taxa with no eligible photograph in the historic metadata scan.

The 42,111 species are an **outcome-blind discovery frame of plants with iNaturalist record opportunities**, not all described angiosperms or a probability sample of the world's species. This metadata-only whole-frame step provides a denominator and workplan, **not** 42,111 measured flower-colour phenotypes.

### Original fixed data and numerical audit

Read source files from the immutable pre-photo FCP branch head `7665d1fe38d7146925b8f3450ed3cb5c2736195a`. Verify SHA256 before any analysis:

- 42,111 species identity frame: `data/frozen/global_monte_carlo_species_discovery_combined_v2_species_v1.csv` / `5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc`.
- all-species observer-capped iNaturalist photo-opportunity frame: `data/frozen/global_monte_carlo_capacity_scan_species_audit_v3.csv` / `12d94100d6343597aa2b87670555e801ce0681f4b69e3ad634fa6d80409008e8`.
- two independent original global taxon/geographic-cell indexes: `data/frozen/global_monte_carlo_species_discovery_observation_index_v1.csv.gz` / `b60a8b1b98fcf6745344d6acc394bb6ff905b9d5bded134031402892ff2efa2e`; `data/frozen/global_monte_carlo_species_discovery_v2_observation_index_v1.csv.gz` / `d7f08601ccc4d7dc9ff1943ea3e58ed903b9759277f7b536bc7423bc08189b3a`.
- all original 162 equal-area discovery cells retained. Cells indicate occurrence/opportunity, not confirmed colour-morph frequencies.

The frozen frame yields these observer-cap photographic-support strata:
 
| Records/species | Species available |
|---|---:|
| >=1 | 41,322 |
| >=2 | 33,944 |
| >=5 | 24,612 |
| >=10 | 18,301 |
| >=20 | 12,985 |
| >=30 | 10,330 |
| >=40 | 8,753 |
| >=50 | 7,632 |
| >=60 | 6,770 |
| >=80 | 5,558 |
| >=100 | 4,730 |
| zero | 789 |

These counts are historical *metadata opportunity*, and are NOT the count of photographs successfully retrieved/processed through the flower ROI colour classifier. A species with one historical photo cannot support a within-species polymorphism test, but its occurrence and photo opportunity remain valid components of the complete census.

### Practical full-universe acquisition and biological targets

All 42,111 retain unique taxon IDs and must be included, irrespective of whether an image exists or is classifiable. A photometric atlas should have a separate field **MEASURED vs UNMEASURED**, never classify unobserved as white or monomorphic.

The single growing world atlas will have nested sampling budgets, using the same original fixed photo opportunity:
- **breadth tier 1**: up to 1 candidate photo/species, with 41,322 historically reachable photo-species; record a strictly *observed visible colour*, not its natural frequency.
- **breadth tier 5**: up to 5 candidate photos/species, with 24,612 candidate species reaching 5 (but every photographed 1–4-photo species retained).
- **reliable descriptive tier 20**: up to 20 candidates/species, 12,985 candidate species at full depth. Variation absence in 20 photos is still not proof of genetic monomorphism.
- **within-species ITV tier 40**: up to 40, 8,753 candidate species; must also pass actual classification, observer/location/time and null exchangeability gates.
- **spatial high-depth tier 100**: 4,730 candidate species; original sampled 2,000 species are historical and the remaining 2,730 have prospectively preselected new taxa with two different study roles (2,000 primary and 730 replication). Preserve their selection/analysis independence.

These are **nested** quotas: one photo at stage 1 counts toward the 5, 20, 40 and 100 photo levels, not five independent datasets. Every source photo/observation ID is deduplicated across stages and prior experimental cohorts. Automatic selection of only colourful photos is forbidden, as is excluding poor-quality / unclassifiable taxa from denominator reports.

Fresh iNaturalist observer-cap source capacity and actual photo URLs can have changed since the census; the updated photo draw must be separately source-verified. Large-scale API use should be batched with the historical source's rate guard and complete failed-request reporting, and sequenced to avoid overloading the concurrent high-depth acquisition. The opportunity ledger itself requires no new API calls.

### Inference hierarchy

- **42,111 species global coverage denominator**: report ecological and geographic scope, metadata photo depth, spatial species occurrence across the 162 equal-area cells and missingness, *not* allelic polymorphism or measured colour frequencies yet.
- **>=1 successfully classifiable photos**: possible global colour-map presence/phenotypic occurrence, with observation-conditioned sampling weights and camera/organ classification errors shown.
- **>=5 / >=20 successfully classifiable photos**: progressively stronger descriptions of coarse visible colour diversity, with explicit detection-power bounds for rare states. No white fraction from one photo should count as true fixed-pigment population frequency.
- **>=40 measured photos + geographic/date opportunity**: candidate species-level within-range geographic colour segregation under colour-composition-preserving permutation, requiring <=50km near pair support and independent observers.
- **>=100 selected photography cohort**: higher spatial replication precision. Its 4,730 historical opportunity species are a subset of the 42,111, NOT a substitute for the whole 42,111 universe.
- **independent prospective result**: only the formerly unexposed species (2,000 primary, 730 predeclared replication) can confirm old photo patterns independently; original 2,000 selected species stay explicitly labelled outcome-exposed historical, with two historic selected taxa having failed 100-photo acquisition.

### Global biological question

Do different plant species repeatedly exhibit within-species photograph colour-state variation with geographically structured allocation, and how does measured observational diversity change among geography, sampling depth and lineages when the actual photographic detection opportunities are explicitly modelled?

This is a macroecological flower-colour atlas with a core independent geographical validation, **not** an experimental identification of a universal net selective advantage or metabolic pigment cost. Morph-specific balancing selection needs genotype, biochemistry and reproductive payoffs, which are not present in the raw iNaturalist photo census.

### Stop/claim rules

1. Never write that 42,111 taxa were *flower-colour measured* until pixels have been fetched and the fixed classifier has successfully processed them, with source IDs, denominators, failures, dates and spatial coordinates retained.
2. A zero-photo species is an honest missing-data species, not monomorphic. A one-photo species is observed but cannot be scored for within-species ITV. Failing to see a rare variant among 5 or 20 photos does not prove absence.
3. Any benefit-cost, selfing, genotype, natural-selection and floral-pigment-loss inference must be separately source-verified; photo-white is exposure/lighting coupled.
4. Do not select species based on previously observed H1/H2, white axes, climatic effects or any favourable flower-colour outcome; retain all identities.
5. The 42,111 frame itself is not a census of all described flowering plants. Conclusions generalize to this sampled photographic opportunity frame, with geographic/citizen-science bias stated.

### Original 42,111 photo-ID URL transport correction (before any colour opening)
An initial original-ID URL resolver invocation (GitHub Actions run 37762416621) returned a successful **technical job status** while **421 of 422 API batches failed** and only 11 original photo URLs were verified. The failed route used `/v1/observations/{100-ids}`, which is known to reject many-ID batches with HTTP 422 Too many IDs. This is a **transport failure**, NOT evidence that 42,100 biological observations are missing or unlicensed. A separate, source-preserving technical recovery uses `GET /v1/observations?per_page=200&id=<100 comma-separated IDs>`, with a 100-ID smoke test before the remaining batches. No species, observation/photo IDs, original species names or biological outcomes are changed. Current observation taxonomy is checked against frozen original taxon ID and taxonomic mismatches are flagged, not repurposed as valid focal-species flower colour.

### Source-respecting technical recovery outcome separation

The original `results/fcp_all42111_original_photo_url_resolution_20261008/` is an immutable **unsuccessful transport record**: 421/422 batches failed, leaving only 11 checked photos; it must not be interpreted as biological image absence. A smoke test with the corrected GET query retrieved exactly 100/100 historical observation IDs, but the corrected full run 37770580281 stopped solely because a safety guard refused to overwrite the original failed `result.json`. No further batch requests occurred in that aborted run. Recovery now writes a different directory `results/fcp_all42111_original_photo_url_recovery_v2_20261008/`, retaining old failure provenance and reusing precisely the originally selected 42,111 taxon/photo IDs. A green completion also requires <=21 failed of 422 batches and >=1,000 licence/URL-valid original IDs; these are transport/coverage floor checks, not a flower-colour effect. No photo pixels have been opened.
