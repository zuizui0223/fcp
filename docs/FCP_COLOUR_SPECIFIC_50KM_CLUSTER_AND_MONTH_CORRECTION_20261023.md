# Colour-specific geographic clustering of photographed flowers: 50-km source-locked test

**Date:** 23 October 2026. **Status:** completed on original 447,500 photographs with verified photo and release SHA256; [successful authenticated workflow #39887641626](https://github.com/zuizui0223/fcp/actions/runs/39887641626).

**Scientific distinction:** This study tests *photographic spatial autocorrelation* of four coarse flower-colour **labels**, not the existence, discreteness, sympatry or inherited identity of biological colour-morph clusters. Its 50-km radius is a sampling and null-comparison scale, not a directly delimited breeding population or a DBSCAN cluster. All analyses below were designed **post hoc** after earlier geographical colour results.

## Source and conditions

- Exactly **447,500 authenticated source photos**, including **227,107 classified** photos of **4,475 taxa**; all other photos remain unclassified. SHA256 of measured photo archive CSV: `8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6`.
- For each species-specific colour label: **at least 40 classifiable photos** in the original sample, **30 pairs** of geospatial-valid conspecific photos within 50 km, **at least five geo-valid photos of that colour and five of other colours**.
- Observed metric: for each species and each colour separately, proportion of all local (at most 50-km) photo pairs in which *both photos* have that colour label.
- **Null A:** shuffle all four colour labels over fixed original coordinates within each species, preserving species label counts and spatial/photo opportunities (199 permutations).
- **Null B:** additionally preserve within-species **calendar-month** colour composition; labels shuffle only among photos originally observed in the same calendar month (199 permutations). This controls broad seasonal label composition, not precise year-month, true flower stage or locale/camera differences.
- **Observer control:** repeat the test after selecting **one photo per observer per species**, using observer/photo IDs only, *not* the photographed colour, and refit the matched nulls rather than extrapolate effect-size changes.
- Multiplicity: four-colour **Holm correction** within each condition. Four observer/month panels are post hoc sensitivities; the Holm values do not claim an independently preregistered 16-test family.

## Verified results

The value shown is **excess local same-label pair probability over the corresponding matched null**. Values are *absolute probability differences*; for example, +0.01243 means 1.243 percentage points, not 12.43%.

| Photographic label | Species in unthinned comparison | Composition-null excess | Calendar-month-null excess | Species, one per observer | One per observer + month-null excess |
|---|---:|---:|---:|---:|---:|
| White | **787** | **+0.01452** | **+0.01243** | **700** | **+0.01053** |
| Yellow / orange | **403** | **+0.01140** | **+0.00883** | **349** | **+0.00845** |
| Red / pink | **313** | **+0.01431** | **+0.01132** | **282** | **+0.01436** |
| Blue / purple | **377** | **+0.01488** | **+0.01323** | **352** | **+0.01106** |

**Every colour in all four sensitivity panels:** upper-tail Monte Carlo P **0.005** (floor of 199 worlds), within-panel four-colour Holm-adjusted P **0.020**.

**Source replay:** Before testing label-specific clustering, the script independently reproduced the original strict `>=40 classifiable / >=30 near-pair` **1,761-species** 50-km global-minus-local mismatch depletion at **+0.0176149021685411**. The synthetic test confirmed that a colour distinction that occurs **only in different calendar months** cannot produce an excess above a calendar-month-preserving null.

Machine-derived exact panel results: [verified JSON](../results/fcp_colour_specific_50km_month_control_20261023.json). Reproduction code: [analysis script](../scripts/analysis/fcp_colour_specific_50km_month_null_20261023.py) and [GitHub Actions workflow](../.github/workflows/fcp-colour-specific-50km-month-control-20261023.yml).

## What was and was not corrected

**Corrected or controlled:** species-wide colour abundance/rarity for each test, fixed observation vertices and uneven photo density, minimum support for minority colours, multiple colour-category comparisons (Holm per panel), repeat photographs by individual observers, and calendar-month photo-colour mix.

**Not corrected or not identified:** true petal-versus-background identity; illumination/white balance; flower age within the same calendar month; native versus cultivated or introduced individuals; observer camera and photograph-quality differences correlated with geography; shared plant identity; genotype; true morph number; inherited stability; morph fitness, costs or selection. A month-adjusted positive comparison does not remove interannual phenology or regional month-dependent photographic bias.

**Scope of the positive result:** At 50 km, **all four coarse photographic colour labels** independently show a modest, positive species-equal local same-label excess after these stated corrections. The positive effect is **not limited to white-vs-other label comparisons**, although the species denominators differ strongly among the four labels. This supports a **colour-specific broad-scale geographic autocorrelation of photographic labels**, not four discrete geographic populations per plant species.

Earlier [original 4,475-species colour-error audit](FCP_PHOTOGRAPHIC_MULTICOLOUR_ERROR_AUDIT_20261023.md) found **4,189** species with at least two raw coarse labels, but only **2,400** with at least two labels supported by five photos each, **435** with three such labels and **19** with all four. This prevents interpreting raw 2-, 3-, or 4-colour detections as verified genetic morph richness. Moreover, the earlier strict **5-km, same-year/month, distinct-observer** spatial null was unsupported; this new 50-km result **does not overturn that negative test or establish local within-breeding-population colour-morph maintenance**.

## Next decisive ecological test, not yet done

To identify actual *separate patches or clusters* rather than aggregate spatial autocorrelation, one would need a per-species, per-colour spatial connected-component/marked-point-process analysis with a composition-, opportunity-, observer- and month-matched null, followed by blind petal-hue verification on original photographs and confirmation of morph stability/genetic identity for selected plants. No such real genotype-linked cluster conclusion follows from the present data.
