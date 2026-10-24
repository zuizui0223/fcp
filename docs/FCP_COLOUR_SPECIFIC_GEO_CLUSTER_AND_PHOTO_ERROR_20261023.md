# FCP: are individual photographic colour categories geographically clustered?

**23 October 2026 — COMPLETE: authentic hue-specific 50-km cluster audit and error controls verified.** This adds a different question to the already completed [original photo label/error audit](FCP_PHOTOGRAPHIC_MULTICOLOUR_ERROR_AUDIT_20261023.md): are *white*, *yellow/orange*, *red/pink*, and *blue/purple* photographs individually overrepresented as local same-colour pairs, once their species' exact four-colour totals are fixed?

## Established empirical starting point (not new results)

- **447,500 source photos, 4,475 species; 227,107 classifiable.** Source SHA256: `8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6`.
- Photo-state richness is highly sensitive to single images: 4,189 species with >=2 labels, **2,400** with >=2 labels supported by **5 or more photos each**; 814 with all four labels, only **19** with >=5 photos *in all four* (and none with >=10 in all four). Under >=40 classifiable photos, **1,592 / 309 / 11** species have >=2/3/4 categories supported by >=5 photos each.
- **The original 50-km mismatched-photo difference is reproducibly above its exact species-fixed-composition/coordinate null** (3,600 species, mean +0.016116, P = 0.005). Its stricter 5-km, year-month, different-observer spatial null **was not supported** (194 species with >=5 pairs, P = 0.300).
- In the [independent observer/source-photo reanalysis](../results/fcp_photo_label_error_audit_20261023/aggregate_results.json), at least two colours with >=5 records and >=40 classifiable source photos plus >=30 local pairs retained **813 species** (P2 depletion +0.028550, P=.005). Taking one colour-blind selected photo per observer retained **727** (+0.026545, P=.005). Requiring >=10 records per colour yielded **462** (+0.037440) and **406** (+0.036319) respectively, both P=.005.
- Small sample photo-colour recurrence is **not** genotype/morph confirmation: the 1-km, same-month, four-observer, two-year screen found **8** candidate species but only **1** passed even the exploratory >=65% primary-palette-fraction screen for all four photographs.

## Authenticated colour-specific geographic results (real original photos)

The [original-photo CI run #39887066227](https://github.com/zuizui0223/fcp/actions/runs/39887066227) completed successfully after authenticating the exact release, original measured photo CSV and permanent species ledger. The new hue-specific test covered **813 species** with >=40 classified photographs, >=30 50-km pairs and at least two colour classes supported by >=5 photos each. Colour category-specific species counts vary because only sufficiently supported categories enter each category's mean.

| Colour label | Species with >=5 photographed records of that colour | Observed minus matched-null local same-colour pair fraction | Four-colour familywise maxT upper P |
|---|---:|---:|---:|
| White | **734** | **+0.013637** | **0.005** |
| Yellow/orange | **385** | **+0.010747** | **0.005** |
| Red/pink | **290** | **+0.014296** | **0.005** |
| Blue/purple | **359** | **+0.013535** | **0.005** |

The total four-colour global-minus-local mismatch in these **813 species** was **+0.0285495**, exactly reproducing the separate October 23 authenticated error audit's restricted P2 estimator. Each row tests an aggregated photo-colour component **conditional on its own source-photo support**. Because the species panels differ by colour, the four conditional means above **must not be summed** to obtain the overall P2 mismatch. No category is a validated independently heritable genetic morph.

### Does the result survive controls?

| Source constraint | Eligible species for the overall P2 contrast | Category-specific result, maxT-adjusted |
|---|---:|---|
| >=5 photos per category; one photo per observer | **727** | **All four categories supported, P = 0.005 each** |
| >=10 photos per category, original observers | **462** | **All four supported, P = 0.005 each** |
| >=10 photos per category, one per observer | **406** | White, yellow/orange, red/pink P = 0.005; blue/purple **P = 0.020** |
| Month-of-year-preserving label permutation, otherwise same 813 species | **813** | **All four categories supported, P = 0.005 each** |
| >=65% main ROI palette colour, >=20-point separation from next group | **183** | White P = 0.005, yellow/orange P = 0.005, blue/purple P = 0.020; **red/pink not supported, P = 0.105** |
| Same ROI rule plus one photo per observer | **154** | White, yellow/orange and blue/purple P = 0.005; **red/pink not supported, P = 0.065** |

The original source month's photo-colour composition itself explained *some* apparent geography: the 813-species calendar-month null mean mismatch depletion was **+0.004349** instead of the near-zero unrestricted-null mean, but the observed **+0.028550** still exceeded that season-conditioned null (P = 0.005). Importantly, preserving **month of year** does *not* preserve the exact **year–month** or demonstrate local coexisting morphs. The older **5-km same-year/month** distinct-observer matched null remains negative (P = 0.300 for 194 species).

Among the 813 eligible species, **674** had two photo colours with >=5 records, **136** had three, and only **three** had all four. Under the exploratory ROI purity screen (183 eligible species), **170** were supported by two categories, **13** by three and **zero** by four. Thus **strongly populated three/four-colour categories are relatively uncommon among spatially analyzable taxa**, and the model cannot infer three/four genetic morphs from these labels.

**Critical ecological boundary:** these are **species-equal aggregated clustering signals**, not positive tests for each individual species. The fraction of source species showing positive category-specific excess was only approximately **50–59%**, depending on hue and sampling gate; mean effects can be driven disproportionately by stronger species. Source-cohort/observer geography, floral age, within-flower bicolour tissue and region-specific image segmentation remain plausible causal alternatives. The harsher photo-composition screen cannot be equated with a calibrated classification-error rate. The 199-world Monte Carlo floor is **P = 0.005**, not evidence of a more precise tail probability.

**Machine evidence:** [compact verified result receipt](../results/fcp_colour_specific_neighbourhood_null_20261023/aggregate_results.json); [unredacted aggregate statistical output artifact (no exact geo, photo or observer IDs)](https://github.com/zuizui0223/fcp/actions/runs/39887066227/artifacts/12473127210), SHA256 of detailed JSON $\mathtt{b53de01d630d565325757064ba4c913fc45048cabe0c082069a7c4e2d88d5336}$. The original raw individual-photo input is not included.

## Additional per-colour geographic statistic and conditional nulls

For species `i` with `n_i` classified original source photographs and category `c` represented by `n_ic`, define the **within-50-km same-colour pair excess**:

$$
E_{ic}(50)=
\frac{\#\{(a,b):d_{ab}\le 50\,\mathrm{km},\; y_a=y_b=c\}}
{\#\{(a,b):d_{ab}\le 50\,\mathrm{km}\}}
-\frac{n_{ic}(n_{ic}-1)}{n_i(n_i-1)}.
$$

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

**Provenance and current status:** the new hue-specific analysis **passed** against immutable original source photo bytes (source/archive SHA checks, 4 synthetic tests, 4-colour-familywise source-null, and exact reproduction of previously authenticated P2 opportunity-gated means). No exact original coordinates, photo IDs, observation IDs or user IDs are exported. Its positive category-specific **photographic** geographic results are recorded in the preceding table, not promoted into a claim about biological morph clusters.
