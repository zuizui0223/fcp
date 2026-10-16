# FCP × RGFCA — 552 REAL original-photo ROI-v4 technical robustness audit (2026-10-16)

## Completed empirical work; explicitly NOT the main 1,000-species study

- Original audited pre-colour **1,000-species RGFCA enrolment**: [GitHub Actions #38828417100](https://github.com/zuizui0223/fcp/actions/runs/38828417100), success. Its 200 × 250 species realization schedule has exactly 50 appearances/species.
- **Actual newly downloaded and original ROI-v4 / EfficientSAM-measured photographs** from six enrolled taxa: [Actions #38832829706](https://github.com/zuizui0223/fcp/actions/runs/38832829706), success, artifact **12018829689**; unique original source-photo IDs **552**, no historical colour measurements reused, original weights and source verified.
- Source photo IDs for this TECHNICAL test were acquired from **live iNaturalist observations API**, filtered for `created_d2 <= 2026-09-27`; photo IDs were **NOT reconciled** to the exact, separate **2026-09-27** Open Dataset archive. Only the 1,000-enrolment *taxon IDs and identities* have original dated-source proof. Do not merge these photo observations with the actual main 1,000-taxon dated-source biological dataset.
- Actual classification: **202 of 552** photographs yielded a classifiable four-state original flower ROI (**36.594%**). Others were explicitly retained as ROI/colour ambiguity terminal statuses, not silently replaced.
- The existing **552 real original ROI measurements** were independently reanalysed with **200 deterministic photo re-draws/taxon**, conserving observer-disjoint splits and requiring >=40 classifiable images and >=20 from each split. Depths per observer split: **20, 32, 40, 44, 48**. Audited [Actions #38838940794](https://github.com/zuizui0223/fcp/actions/runs/38838940794) SUCCESS; real output artifact **12020038776**; result CSV and JSON carry the actual six-taxon morph counts and source SHA256.

### Empirical technical data (all source photo outcomes real)

| Species | Original photos | Classifiable | Coarse photo labels | P(label pair), 40+40 | P(label pair), 48+48 |
|---|---:|---:|---|---:|---:|
| *Matricaria chamomilla* | 96 | 64 | white 56 / yellow-orange 8 | **0.830** (200 evaluable) | **1.000** (200 evaluable) |
| *Nassauvia axillaris* | 96 | 43 | white 42 / yellow-orange 1 | 0.000 (4 evaluable) | 0.000 (200 evaluable) |
| *Combretum grandidieri* | 89 | 31 | red-pink 16 / yellow-orange 14 / white 1 | unevaluable | unevaluable |
| *Acourtia nana* | 96 | 31 | white 27 / yellow-orange 4 | unevaluable | unevaluable |
| *Chloris formosana* | 96 | 19 | white 13 / yellow-orange 6 | unevaluable | unevaluable |
| *Grevillea decurrens* | 79 | 14 | white 10 / yellow-orange 4 | unevaluable | unevaluable |

`P(label pair)` here is a **technical conditional fraction among evaluable re-draws from the exact same original technical-pilot photographs**; it is NOT the 1,000-species formal `P_poly`, independent biological validation, or hereditary FCP.

For *M. chamomilla*, detection depended on photo budget: 32+32 yielded 158/200 evaluable re-draws and `P=0.63291139`; 40+40 yielded 200 evaluable and `P=0.83`; 44+44 gave `P=0.975`; 48+48 (the whole 96 raw-photo pool) gave `P=1.0`. **Apparent random-draw robustness is not sufficient evidence of among-individual polymorphism.**

### Discovered false-positive mechanism

For *M. chamomilla*, among its **64 classifiable** original ROI measurements, **58/64 = 90.625%** included both **white and yellow fractions >=0.10 within the SAME PHOTO ROI**. Thus a high proportion of ROI images were intrinsically dual-coloured even when their one selected coarse morph label changed with photographic context.

This is biologically interpretable as a genuine **floral-architecture confounding risk**, because in *M. chamomilla* the white ray florets and yellow disc florets are normal parts of **one capitulum** (Flora of China: https://www.efloras.org/florataxon.aspx?flora_id=2&taxon_id=242331663; Kew: https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:154715-2/general-information). Those standard tissue colours alone do not establish one all-white-morph individual and one all-yellow-morph individual. Actual pixel crop images were discarded by the pilot so these particular photos still require visual organ-level review; **we do not assert that all 58 corresponded to disk/ray organs solely from the palette measurements**.

### Immediate post-pilot protocol correction (implemented)

1. The **image-first species screen** `fcp_pipeline/photo_polymorphism_screen.py` now retains within-ROI white + yellow co-occurrence summaries.
2. The **RGFCA species repeated photo-resampling engine** `fcp_pipeline/rgfca_photo_repeatability.py` keeps both original flower colour fractions (alongside the four coarse labels). When >=40 original classifiable ROI photos support white *and* yellow-orange coarse labels, and >=50% of those SAME ROIs contain white+yellow fractions >=10% each, a formerly `repeatably_photo_polymorphic_phenotype` species is instead placed in **`repeatable_colour_LABELS_WITHIN_ROI_STRUCTURAL_QC_HOLD`**. The raw photo-label repeatability estimator is preserved for inspection, but the species is not silently admitted as FCP.
3. These thresholds were established **after** the separate six-taxon technical pilot, and must be identified as a **post-pilot measurement QC rule established before the main 1,000-taxon original-photo dataset**. They cannot be claimed to have been preregistered before this pilot. The rule deliberately triggers **manual review, not a declaration that the plant lacks polymorphism**.
4. No result from the six source-unreconciled technical taxa is promoted into the main 200×250 replicated species estimator; the original 1,000 species and hash-balanced enrolment are unchanged.
5. Manual review must determine whether distinct floral-display morphs represent different plants/local populations after accounting for flower architecture, developmental stage, repeat photographs of the same plant and observer/photography preferences. Same within-ROI mixed colours are a flag, **not proof of misclassification**.

## Interpretation and operational boundary

The current empirical result is **a calibration of picture-level detectability and structural false positives**, not a claim of six natural polymorphisms or global FCP frequency. The next main biological blocker is genuine original photograph IDs and newly measured ROI-v4 images drawn under the dated 2026-09-27 frame for all 1,000 enrolled taxa. The archive itself is not retained on permanent disk after the earlier count-only hosted-runner job, and no worker-wide 128,000-photo measurement may be declared complete until appropriate storage/provider access and genuine terminal photo receipts exist.

**Scientifically important result already supported:** even **200** image resampling rounds with observer-disjoint halves can reproducibly rediscover differences among **photo labels** that may not correspond to discrete morphs among plants; therefore an organ/developmental validation tier is necessary in addition to sample-size and observer bias controls.
