# FCP colour-label error audit — real 4,475-species source, 23 October 2026

**Current status: COMPLETE, original 447,500-photo audit verified.** Both the permanent 4,475-species ledger and original photo CSV were validated by exact SHA256. Independently repeated source-conditional geographic null tests (199 worlds) passed after thresholding minor colour classes and after colour-blind one-photo-per-observer selection. Machine result: [aggregate_results.json](../results/fcp_photo_label_error_audit_20261023/aggregate_results.json); full run [#39879387172](https://github.com/zuizui0223/fcp/actions/runs/39879387172). **These tests do not establish true inherited flower-colour morphs.**

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

## 4. Authenticated 447,500-photo observer audit and new locality nulls: completed

The exact archived original photo bytes passed their SHA checks; the complete CSV contained 447,500 unique source photo IDs and 227,107 coarse-colour classifiable records. Every one of the 4,475 original taxon's colour counts and each originally estimated P2 50-km depletion matched the permanent species ledger.

**Independent observer counts per colour:** original source sampling limited each observer to at most **two source photos per taxon**. Therefore, **each recorded colour supported by five photographs was automatically represented by at least three distinct observers and three separate observation IDs**. This held in all **2,400 / 435 / 19** candidate species with ≥2 / ≥3 / all 4 coarse colours supported by ≥5 photos each. **Crucially, this additional observer-count filter is *mathematically redundant* with the five-photo threshold under the existing acquisition cap.** It is *not* independent proof of true floral morph identity, different genotypes or independent field populations.

**Same-observation disagreement test:** the source contained **zero distinct observation IDs represented by two or more colour-classified photos** (zero eligible multi-photo same-observation comparisons), so **zero observed within-observation label conflicts is NOT evidence of perfect classification**. The comparison is unidentifiable with this source photo selection.

### Matched 50-km spatial null after photographic support and observer thinning

All 199 within-species label-permutation worlds preserve each analysed species' original observed colour counts and fixed geographic sample vertices. Count-based filter thresholds are therefore invariant across the permuted worlds. One-photo-per-observer selection retained a fixed photo per observer using observer/photo ID order **without looking at photo colour**, then recomputed the matched null, not just an original effect mean.

| Real original photographic sample | Species with 50-km comparison | Observed mean (global minus local) | Upper P, 199-world null |
|---|---:|---:|---:|
| All original source photos, baseline reproduced | 3,600 | **+0.016116** | **0.005** |
| One photo per original observer | 3,592 | **+0.015900** | **0.005** |
| ≥40 classifiable, ≥30 local pairs, ≥2 colours × ≥5 photos | 813 | **+0.028550** | **0.005** |
| Same ≥5/photo filter **and one per observer** | 727 | **+0.026545** | **0.005** |
| ≥40 classifiable, ≥30 local pairs, ≥2 colours × ≥10 photos | 462 | **+0.037440** | **0.005** |
| Same ≥10/photo filter **and one per observer** | 406 | **+0.036319** | **0.005** |
| ≥40 classifiable, ≥30 local pairs, ≥3 colours × ≥5 photos | 139 | **+0.031084** | **0.005** |

**Interpretation:** Positive nonrandom **broad-scale spatial structuring of photographed colour labels survives** removing within-observer repeat photographs and restricting analyses to colour classes with several observations. Therefore the 50-km result is unlikely to be *entirely* the product of a single isolated misclassified photo per taxon or two repeat shots by one observer. This is a meaningful **photographic robustness result**, not an empirical estimate of a biological morph-frequency or image-error process.

**Limits:** The thresholds, observer-thinning choice and tests were designed **after** the original results; these are exploratory post hoc sensitivities. P=0.005 is the **minimum obtainable value from 199 worlds**, not evidence of a precisely estimated smaller P. The cohort structure, flower organ/stage composition, camera calibration, conservation of true petal reflectance, intraspecific crop/native confusion, geographically uneven model errors and same-location observer behaviour are **not eliminated**. Colour labels were retained as originally assigned; redoing a photo-level permutation cannot prove any of those labels botanically correct. The existing strict same-year/month ≤5-km null remains unsupported (P 0.300/0.275), so local genetic coexistence is still not established.

**Output and code:** [aggregate JSON](../results/fcp_photo_label_error_audit_20261023/aggregate_results.json) contains authentic photo/archive SHA receipts and the exploratory per-subset observed/null results. [Code](../scripts/analysis/audit_fcp_photo_labels_original_20261023.py) and [workflow](../.github/workflows/fcp-original-photo-label-error-audit-20261023.yml) regenerate it. The Actions upload additionally contains a **taxon-level support panel with no photo IDs, observer IDs or coordinates**. The output is not a substitute for human blind original petal checks.

## Interpretation rule

- **Supported by actual reanalysis:** 4,189 taxa contain ≥2 original photo-colour labels, but higher apparent 3/4-colour richness is vulnerable to scarce minority labels; the 50-km colour-label spatial null survives multiple photographic support tiers and one-photo-per-observer controls.
- **Not identifiable from this source:** same-observation discordance (zero eligible duplicate-observation photo groups), direct independently annotated petal hue, actual classifier false-positive/false-negative rates, and the share of geographic structure attributable to real biological morphs versus spatially systematic imaging/source errors.
- **Not established:** individual-genotype inherited morph frequencies, reliable 3/4 genetic morph count, local sympatric morph maintenance, photo-independent climatic adaptation, pigment costs or differential selection.

The central manuscript should use **“photographic colour-state variation”** where species-level genetic colour polymorphism has not been validated.
