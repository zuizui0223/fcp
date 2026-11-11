# FCP: independent botanical descriptions versus the east–other rainfall white–blue/purple direction

2026-10-13. Frozen plan, before opening independent botanical species-level flower-colour labels from NHM raw data for this FCP analysis. This is post-outcome relative to FCP PRs #164/#188–190, not prospectively selected geographic/colour hypothesis.

## Biological question

Can a separate **botanical-description source** (not image photo hues and not photography/classification success) detect the signed difference in rainfall association between **white versus blue/purple** described flower species in the source-defined East 90–180°E sector and other 90° world sectors? Distinguish three levels:

1. Independent measurement of flower colour from botanical literature.
2. Original FCP iNaturalist photograph **locations and long-term climate as environmental sampling opportunities**, reused in this analysis; not independent environment or new occurrence source.
3. True genetic colour morph frequencies, plant fitness, pollen limitation, dry-stress mechanisms: completely unobserved.

## Truly external phenotypic provenance

NHM, *Angiosperm functional traits extracted from botanical descriptions*, DOI **10.5519/p3dm31kc**, downloaded by SHA-pinned GitHub public archive `NaturalHistoryMuseum/plant-traits`, commit **4843c22bb2fc91018d9828dcddebf067c83b1097**, `output/raw/ocr-traits/combined.csv`, 17,279,412 bytes, raw Git blob SHA-1 **e1cfd2d5c518c2061ec91d1a9e8b64fc7ca69389**, licensed **CC BY 4.0**. Source authors: Araujo, Hartley, Brummit. Unlike the parent FCP source, its flower colours derive from machine-parsed botanical descriptions and may be wrong or ambiguous; agreement with independent morphological expert scoring is not assumed.

The original FCP georeferenced 100,543 source photograph records remain exactly SHA256 **2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0** from historical Actions artifact `11587115220`. Combine photo IDs and admit **all georeferenced climate-complete original photo opportunities, even if FCP images could not be colour-classified**. No original FCP photo `morph`, `flower_fraction_*`, classification-success label, or FCP PR190 hue outcome enters NHM colour classification, model building or outcome-based filtering.

## Outcome mapping, frozen before source colour join

- Match **exact botanical binomial species strings only** between FCP photo name and NHM `taxon`. No fuzzy, synonym, authority, infraspecific expansion or FCP-photo-hue-based taxon repair.
- Inspect `petal colour`, `corolla colour` and `flower colour` and require **all observed organ colour fields to identify exactly the same single hue family**. White is only the standalone lexical word `white`. Blue/purple is the standalone word `blue`, `purple` or `violet`. Words denoting other colours (including cream/ivory/lilac/lavender) are intentionally `other` and do not silently enter one of the two states.
- A colour field with no recognized hue (but nonblank text), two recognised hue families, contradictory flower-organ categories, duplicate binomials with conflicting labels, or absent colour fields is **HOLD**. Other unmixed hues are descriptively counted but are **not** in the white–blue binary contrast.
- This is a species-level **described flower colour trait**, not within-species variation or local colour-genotype frequency. Every site record of one taxon carries the same botanical label in this design.

## Frozen estimand

For each of 4 equal-width longitude sectors (`[-180,-90),[-90,0),[0,90),[90,180]`), compare (i) climate at observed FCP sites carrying botanically white species and (ii) climate at sites carrying botanically blue/purple species. Use distinct source species means **within each sector versus all other sectors**, require >=3 independent named species on each side of a common genus, and centre all measured species-level predictors and response by genus *within each side*. Include annual precipitation BIO12 plus BIO1, BIO5, BIO15, elevation, absolute latitude and longitude sine/cosine (8 axes). Standardize predictors using **all eligible FCP photo-site opportunities**, not only colour-classifiable FCP images or independently labelled taxa.

The primary contrast is:

`C=(beta_white,east-beta_blue,east)-(beta_white,rest-beta_blue,rest)`, in species probability units for one **global FCP photo-opportunity SD of annual rain BIO12**. Because this independent description panel only includes white/blue species, the binary probability denominator differs from the four-hue parent FCP outcome; **a matching sign would be convergent support of environmental-direction association, not replication of the same numerical estimand or magnitude**.

Estimate species-equal and genus-equal weighting separately. Paired genus bootstrap uses 999 deterministic multinomial draws, for the same source genera present on both sides. The **four sector-vs-rest tests are one pre-fixed maxT family** if and only if all four pass the support gate; otherwise report the individual conditional estimates only as exploratory and HOLD confirmatory four-sector inference.

## Pre-fixed opportunity/support gates

Primary: >=3 species per side per common genus, >=30 qualified genera, >=100 independently labelled species on EACH side, and >=20 **white and blue/purple named species** on EACH side. If any of these fail, output `HOLD_INDEPENDENT_BOTANICAL_COLOUR_SUPPORT`; do not reduce gates post hoc to extract significance. Additional secondary diagnostics use >=2 species per side and/or taxa with >=2 original FCP photo-site opportunities, all separately labelled posthoc support sensitivities. If missing estimated design rank, or a singular paired bootstrap draw, HOLD. Report unmatched/ambiguous/mixed NHM rows and FCP-NHM exact overlap counts. Never infer biological no-effect when underpowered.

## Scientific decision before outcomes

The FCP parent familywise + genus-equal + strict photo effect did **not** provide stable support. Independent NHM botanically described colours must have a positive, source/weighting-stable east contrast and four-region familywise confidence interval that excludes zero before considering an additional observational replication. Even so, genetic pigmentation adaptation / plant fitness remain untested because sampling sites and environment still come from original FCP iNaturalist distribution. A NULL, negative, or eligibility-HOLD result must be retained unchanged and cannot be replaced by later favourable trait remapping, longitude definitions, post-outcome filtering or using PR #133 2,000 taxa as unexposed outcome (their one-photo colour already appeared in 42,111 atlas).

## Attribution and independence

Araujo, A.C., Hartley, B. & Brummit, N. (2026). *Angiosperm functional traits extracted from botanical descriptions*. NHM Dataset DOI **10.5519/p3dm31kc**, CC BY 4.0; repository `NaturalHistoryMuseum/plant-traits`. The text-derived phenotype is independent of the FCP photo colour measurement algorithm, but both data sources remain taxonomically nonrandom and the FCP environmental opportunity source is reused.
