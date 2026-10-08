# Global FCP geographic-gradient atlas — latitude and elevation
## 2026-10-08 — three-cohort, species-equal exploratory extension

### Biological scope and question

The target is a **comparative geographic law across sampled plant species**, not an account of one taxon's fitness mechanism. We have 149,900 source photographs of 1,499 sampled plant species in three species-disjoint discovery, validation and prospective-photo resources. Previous within-species signed latitude/elevation contrasts did not transport as one universal direction; this does not exclude nonmonotonic geography, elevation bands, hemisphere-specific sampling differences, or different colour-diversity distributions.

**Primary question:** where, across absolute latitude and approximate elevation, are photographed species locally more white/coloured and more/less four-state-diverse than expected from each *same species'* exact range-wide visible-colour composition and the exact numbers of photos available in those bands?

No climate/UV/altitude adaptation mechanism is inferred by finding a photograph-gradient pattern.

### Input bytes, population and geographical resolution

Use only immutable measured-photo tables at:
- Discovery SHA256 ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4.
- Validation SHA256 0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6.
- Third SHA256 57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186.

The current PR already downloads an exactly ZIP-SHA-pinned WorldClim 2.1 10-arc-minute elevation layer (ZIP SHA256 6a9f4e9a37f289d594ed1567b63be344a4dacf247f0b69756e81df031087e3bd). This is **coarse raster elevation**, not GPS altitude or microclimatic elevation at the actual photographed plant. Latitude is derived directly from source coordinates.

Bands fixed before this route opens results:
- absolute latitude: 0–15°, 15–30°, 30–45°, 45–60°, 60–90°;
- coarse elevation: <250 m, 250–1000 m, 1000–2000 m, >=2000 m;
- north/south latitude-band presence table is DESCRIPTIVE only and explicitly confounded by different species pools.

### Species equal weighting; separating species turnover from within-species geography

The mean white fraction in northern versus southern plant observations can change merely because different species occur in those zones. A fair *intraspecific* FCP comparison therefore requires the **same species** to be observed on both sides of a band distinction.

Use >=40 classifiable records per species globally. To contribute to a given band, a species must have >=8 classifiable photos *inside* that band **and** >=8 classifiable photos outside the band, all along the current geographic axis. Species with >=8 photos in a band but no outside support are retained in a coverage ledger but contribute NO within-species effect. No photos or taxa may be substituted based on colour.

The two within-species responses are:
1. **white shift**: fraction photographed white within band minus fraction photographed white species-wide;
2. **four-state colour-diversity shift**: unbiased within-band pairwise morph-discordance minus the same unbiased species-wide pairwise discordance.

The unbiased finite-sample pairwise discordance is [n^2 - sum_m n_m^2]/[n(n-1)]. This guards against spurious larger diversity solely from larger sample size. Mean band effects give each species one equal weight, irrespective of photographs or genera. Also describe the fraction of band-supported species with >=3 photographed white and >=3 photographed nonwhite records, separately from a **natural-population FCP prevalence** claim. Visible red/pink, blue/purple, yellow/orange are not a common validated biochemical pigment state.

### Exact conditional null and uncertainty

For each species, preserve its exact 4-class photographed flower-colour counts, actual coordinates, raster-derived altitude, bin membership and bin observation count. Under no geography-specific morph assortment, exchange the colour labels among *that same species' sampled photo locations*. Run 199 predetermined seeded permutations; aggregate the band-specific null with exactly the same species set and equal-species weights as the observed estimate.

Report means, null means, species-level 999-bootstrap 95% intervals, within-band n species/n genera/n photos, genus-balanced means, and per-band fraction of observed photographed white+nonwhite groups. Adjust exploratory two-sided band P values using max-T across all adequately sampled bands *within an axis and cohort*, separately for white shift and total D shift. Keep absent bands/status explicit.

### Generality decision rule

At least **25 species with within-species band contrasts** in a cohort/band are required to label that band "coverage-eligible exploratory". Fewer species = **HOLD_SPARSE_WITHIN_SPECIES_BAND_OPPORTUNITY**, not a geographic zero. A geographically directional/midlatitude/mountain band effect should only be presented as a cross-cohort candidate if the signed effect agrees in **all three nonoverlapping species cohorts**, no cohort fails coverage, and max-T p<0.05 in each. Since candidate bands are reviewed after outcomes and all cohorts share the same photo provider/measurement algorithm, even that is exploratory same-source transport rather than untouched experimental confirmation.

Previous first-frozen 50km distributed-polymorphism result (0.02053/0.01867/0.01468, matched p=.005 each) is *different*: it assesses nearby pair colour homogeneity relative to species composition, not macroregional frequency of white or species-local alpha diversity in 15-degree/elevation bands. Do not combine the estimands or imply the long-span band effect identifies mating-population adaptation.

### STOP rules and interpretation

- The 42,111-species opportunity frame and 1,499 high-depth sample cannot be used as an unbiased denominator of global biological white-colour polymorphism incidence.
- Observed white is **known to be coupled to digital photo exposure**; a white image classification is not proof of petal anthocyanin loss.
- Regional composition shifts driven by species turnover must never be renamed within-species shifts. The hemispheric coverage is a map/description, not a test for plant evolution.
- No assertions about pollinators, adaptive thermal protection, UV absorption or selection coefficients are possible from banded photographs.
- This is a self-contained exploratory analysis, not a modification of existing frozen H1/H2/ITV-spatial claims or their source hashes.
- If there is insufficient within-species cross-band geographic support, report this as the central empirical coverage limitation and do NOT hunt additional bins until a positive sign appears.

### Files and CI output

Script: scripts/analysis/run_polymorphism_global_geography_bands_20261008.py.
Tests: tests/test_polymorphism_global_geography_bands_20261008.py.
Reuse the immutable-photo and WorldClim-elevation retrieval from .github/workflows/polymorphism-white-chromatic-clines-20261008.yml.
Output: global_geography_atlas/result.json and species_geographic_band_estimates.csv inside the existing run artifact. These outputs are not frozen biological results until Actions executes, tests pass, shape is audited, and the evidence result is copied into a source-backed Git tracking receipt.
