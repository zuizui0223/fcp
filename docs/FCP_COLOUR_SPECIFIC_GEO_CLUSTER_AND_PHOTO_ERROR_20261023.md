# FCP: are individual photographic colour categories geographically clustered?

**23 October 2026 — new source-locked hue-specific hypothesis test; results not yet adjudicated.** This adds a different question to the already completed [original photo label/error audit](FCP_PHOTOGRAPHIC_MULTICOLOUR_ERROR_AUDIT_20261023.md): are *white*, *yellow/orange*, *red/pink*, and *blue/purple* photographs individually overrepresented as local same-colour pairs, once their species' exact four-colour totals are fixed?

## Established empirical starting point (not new results)

- **447,500 source photos, 4,475 species; 227,107 classifiable.** Source SHA256: `8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6`.
- Photo-state richness is highly sensitive to single images: 4,189 species with >=2 labels, **2,400** with >=2 labels supported by **5 or more photos each**; 814 with all four labels, only **19** with >=5 photos *in all four* (and none with >=10 in all four). Under >=40 classifiable photos, **1,592 / 309 / 11** species have >=2/3/4 categories supported by >=5 photos each.
- **The original 50-km mismatched-photo difference is reproducibly above its exact species-fixed-composition/coordinate null** (3,600 species, mean +0.016116, P = 0.005). Its stricter 5-km, year-month, different-observer spatial null **was not supported** (194 species with >=5 pairs, P = 0.300).
- In the [independent observer/source-photo reanalysis](../results/fcp_photo_label_error_audit_20261023/aggregate_results.json), at least two colours with >=5 records and >=40 classifiable source photos plus >=30 local pairs retained **813 species** (P2 depletion +0.028550, P=.005). Taking one colour-blind selected photo per observer retained **727** (+0.026545, P=.005). Requiring >=10 records per colour yielded **462** (+0.037440) and **406** (+0.036319) respectively, both P=.005.
- Small sample photo-colour recurrence is **not** genotype/morph confirmation: the 1-km, same-month, four-observer, two-year screen found **8** candidate species but only **1** passed even the exploratory >=65% primary-palette-fraction screen for all four photographs.

## Additional per-colour geographic statistic and conditional nulls

For species `i` with `n_i` classified original source photographs and category `c` represented by `n_ic`, define the **within-50-km same-colour pair excess**:

[
E_{ic}(50)=
\frac{\#\{(a,b):d_{ab}\le50\,\mathrm{km},y_a=y_b=c\}}
{\#\{(a,b):d_{ab}\le50\,\mathrm{km}\}}
-\frac{n_{ic}(n_{ic}-1)}{n_i(n_i-1)}.
]

The second term is the *exact* species-composition-conditioned expected fraction of category-`c` matching pairs under unrestricted exchangeability, given the original finite photo counts. **Algebraic identity:** `sum_c E_ic = D_species-global − D_species-local`. Four colour groups are **nominal image labels**, not inferred four independent evolutionary morphs, and their aggregate tests are statistically dependent.

**Analysis units:** each source taxon has at most 100 photographs and equal species weight. Main opportunity gate requires >=40 classified photographs, >=30 within-50-km pairs, and at least two colour categories supported by >=5 photographs each. A specific colour enters its *own* inferential mean only if its own count is >=5; the four category null contrasts are studentised and adjusted together using the **maxT** across 199 shared within-species shuffle worlds. Independent observer/photo selection is specified without using hue labels.

**No source photo coordinates, source user IDs, photo IDs or observation IDs are exported.** Original photo/archive/ledger SHA are authenticated before the real test.

## Diagnostics designed to challenge image and season confounding

The [source-verified implementation](../scripts/analysis/fcp_colour_specific_neighbourhood_null_20261023.py) and [CI](../.github/workflows/fcp-colour-specific-neighbourhood-null-20261023.yml) apply separately:

1. Original source classes, >=2 categories with >=5 photos, all available observers (expected **813 eligible taxa**).
2. Same gate after **one colour-blind selected photo per observer** (expected **727 eligible taxa**).
3. More demanding >=2 colours with **>=10 photos** each, with/without one-per-observer selection (expected **462 / 406 eligible taxa**).
4. Retain original sample, but shuffle colours **within each species and calendar month**; this preserves season-specific image-label frequencies while testing additional geographic sorting. A species whose colour labels do not overlap among months is not an exchangeable test of spatial selection under this null. The method is **not** the previous 5-km same-year-month result.
5. Exploratory compositional QC of ROI labels with >=65% share of the original declared colour group and a >=20 percentage-point margin above the largest competing group, with/without observer thinning. Such purity is **not** calibrated confidence in floral morph identity.

The color categories' significance is defined against the relevant category-specific conditioned null, not against an assumption that photographs are error-free. The paired null tests remain **post hoc**, as are the photo-quality thresholds. Multiple analysis tiers represent sensitivity choices and are not independent replicated sources.

## Decision rules and limitations

- **If matched hue-specific nulls remain positive** after observer thinning and calendar-month conditioning, conclude only that within-species *photographed hue labels* have category-specific geographic organisation at this spatial scale. This would challenge single-photograph and some phenology confounding mechanisms but not calibrated geographic camera/flower-part error.
- **If effects disappear under month controls**, report season/source confounding as a plausible explanation; do not treat the original broad geographic association as local adaptation or genotype frequency divergence.
- **If the strict ROI composition filter leaves insufficient eligible species**, record NA and opportunity loss rather than a biological zero or a new assumed local morph count.
- **Do not interpret species with two, three, or four photographed label classes as species with that many inherited morphs.** Distinct petal phenotypes require organ/stage-matched blind image adjudication; inherited morphs require source plants/genetic/common-garden controls. No isolated photographed colour cluster alone establishes local coexistence, pollinator selection, pigment maintenance cost, dispersal barrier, or speciation.

**Provenance and current status:** the new hue-specific analysis runs on an immutable archival release (the raw photographer coordinates and identifying observation metadata remain confined to the CI runner). Small, non-identifying aggregate outputs will be placed in the PR's Actions artifact after source and synthetic checks. Until that genuine-data run succeeds, no hue-specific cluster test has been reported as a biological positive.
