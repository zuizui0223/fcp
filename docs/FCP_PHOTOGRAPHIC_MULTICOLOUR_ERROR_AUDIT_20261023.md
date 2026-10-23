# FCP colour-label error audit — real 4,475-species source, 23 October 2026

**Current status:** authentic 4,475-species / 114-column joined CSV audited at the species level. A new GitHub Actions workflow (source archive SHA-checked before use) is configured to test **original-photo observer support, same-observation discordances, and 50-km matched nulls after colour-blind observer thinning**. Raw-photo-stage outcomes are **not yet established by this document**.

## Why this is essential

**Four-state photographic label diversity is not verified inherited flower-colour polymorphism.** At least two recorded colours can arise from one wrongly segmented flower, background pixels, organ/age differences, horticultural forms, illumination or true within-taxon biological differences. The number of *observed* colour classes is especially unstable when a colour is seen in only one or two source photos. Neither a very small permutation P nor 100 photos per species is a substitute for independent petal-organ and genotype checks.

Input: `results/fcp_np_all4475_unified_species_ledger_20261021.csv` (source-locked SHA256 `23d3daadc9799c2e7b5867d54ccc1ed98369ba79a81b72dc7f831e5197e63c49`), from photo CSV SHA256 `8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6`.

## 1. Actual colour-class richness is highly sensitive to rare photographs

| Observable count, not genotyped morphs | Species among all 4,475 |
|---|---:|
| Exactly one observed coarse colour | 286 |
| Exactly two observed coarse colours | 1,743 |
| Exactly three observed coarse colours | 1,632 |
| Exactly four observed coarse colours | 814 |
| At least two colours, >=5 photos **per colour** | **2,400** |
| At least three colours, >=5 photos **per colour** | **435** |
| All four colours, >=5 photos **per colour** | **19** |
| At least two colours, >=10 photos **per colour** | 1,223 |
| All four colours, >=10 photos **per colour** | **0** |

The original **4,189 species with >=2 recorded labels** is consequently **not** a defensible estimate of 4,189 species with genetic colour polymorphism. The 814 taxa displaying all four labels likewise do not warrant a four-inherited-morph conclusion: only 19 of these satisfy even five photographed records in each colour class.

With the older eligibility gate **>=40 classifiable photos/species**, 3,297 taxa remain: at least two 1-photo labels in 3,029, at least two 5-photo-supported labels in **1,592**, three 5-photo-supported labels in **309**, and four in **11**. With >=60 classifiable photos, 1,425 taxa remain, of which 458 have two >=5-photo colours and **none** has four >=5-photo colours.

The **minimum number of label substitutions required to make a photographed species monochrome** equals classified-photo count minus the largest original colour-class count. This exceeds five photos for **2,470** species and exceeds ten for **1,489** species, so the entire multicolour pattern cannot always be attributed to *one* incorrectly labelled image per species. However, these are **arithmetic edit distances**, not estimated model error rates or proof that all colours are biological.

## 2. Classification success is strongly coupled to apparent D

The permanent primary source reports **227,107 / 447,500 classifiable photographs** (50.75%): **220,393 are not assigned one of the four labels**. Species with fewer than 40 classifiable photos have mean conditional photo-colour D **0.44401** (1,178 species), versus **0.23751** with >=40 photos (3,297 species); Spearman correlation between the number classifiable and measured D is **−0.58989**.

Even allowing every missing photo *hypothetically* to have one of the same four labels, attainable species-equal all-100-photo D lies between **0.15282 and 0.67063**. These are adversarial mathematical completions, **not** confidence bounds on biological diversity: an unclassifiable record may not contain a valid petal organ at all. Thus **measurement/classifiability bias is a high-severity barrier to species-level colour-richness prevalence and any proposed 3/4-morph comparisons**.

## 3. Preliminary P2 spatial effect among stronger photo-colour-count species

Using the exact 4,475-taxon permanent ledger's *already estimated* per-species global-minus-local 50-km mismatch, **descriptive** filtering gives:

| Source-colour observation filter | Taxa with measured 50-km mismatch | Mean mismatch depletion |
|---|---:|---:|
| At least two colours, >=5 photos per colour | 1,902 | **+0.02307** |
| At least two colours, >=10 photos per colour | 989 | **+0.03540** |
| At least two colours >=5 each and >=40 classifiable photos | 1,274 | **+0.02837** |
| At least two colours >=10 each and >=60 classifiable photos | 198 | **+0.03625** |

**Important:** These means are *observed estimates*, **not new matched-permutation P values**. Filtering on the same photo-colour counts can change effect magnitude, variance, geographic pair opportunity and the species population. This is not evidence for stronger biological selection in well-photographed species.

The original source-conditioned 50-km test across **3,600** species gave mean +0.016116 and matched 199-world P=0.005; the older >=40-classifiable / >=30-local-pairs gate across **1,761** species gave +0.017614902, P=0.005. A stricter same-year/month, <=5-km and distinct-observer comparison **was unsupported** (P=0.300 for 194 species with >=5 eligible pairs; P=0.275 for 14 species with >=30). Therefore geographic photo-colour clustering is observed on the broad scale, but **within-population inherited morph coexistence or true localized differentiation remains unverified**.

## 4. Original-photo analysis being executed

The new audit script `scripts/analysis/audit_fcp_photo_labels_original_20261023.py` and workflow `.github/workflows/fcp-original-photo-label-error-audit-20261023.yml` do the following:

1. **Fail closed on original ZIP and measured-photo SHA**, source table count and exact joining to all 4,475 original taxon IDs.
2. **Count observers and distinct observation events for each colour class**, and separately report how many taxa retain 2/3/4 labels at >=5 photos and >=3 independent observers and observation events *for each label*. Observer-repeat evidence is necessary, not sufficient, for a true inherited morph.
3. **Flag disagreeing colours within a single observation ID**, where two photographs of the same recorded observation receive different coarse labels. This could represent multilobed or bicolour flowers, different photographed organs or image/age errors—not necessarily false labels.
4. **Reproduce the exact original 50-km per-species mismatch and 199-world composition-preserving source null**, then repeat source-conditioned permutations on *fixed-count* robust colour subsets and on **one-photo-per-observer** samples chosen *without inspecting colour*. The latter reduces repeated observer-photo pseudoreplication but still cannot rule out geography-specific camera/flower-stage biases.
5. Upload **only aggregate JSON and a per-taxon observer-support panel with no exact photo IDs, observer IDs or geographic coordinates**. Keep all original measurements and existing results unchanged.

The robust subset P values are **new exploratory tests** and are permitted only because the corresponding photo-colour-count thresholds remain **fixed under each within-species label permutation**. Raw support based on observer counts is reported descriptively and does not automatically receive a P value unless selection is also recalculated under the null. Under photo thinning, selection is observer/photo-ID based and colour-blind prior to the matched null.

## Interpretation rule

- **Supported now:** multiple photographic colour labels exist, but much of the apparent 3/4-class richness has very small label support. Broad-scale same-species photo labels have nonrandom spatial organisation under the original null.
- **Pending actual original-photo audit:** independent-observer support for each colour, within-observation disagreement and repeat-analysed 50-km null conditional on colour-blind observer thinning.
- **Not established:** individual-genotype inherited morph frequencies, true morph counts per wild population, the fraction of geo sorting attributable to biological colour morphs, pigment cost or selection.

The central manuscript should use **“photographic colour-state variation”** where species-level genetic colour polymorphism has not been validated.
