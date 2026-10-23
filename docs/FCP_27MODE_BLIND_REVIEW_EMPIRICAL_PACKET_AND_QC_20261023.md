# FCP: 27 observed photo-mode taxa, 27 matched controls — source-photo validation packet

**Completed 23 October 2026.** Source file verification, candidate/control assignment, reviewer packet generation and an independent image-technical audit have all passed in GitHub Actions. **This is not a completed human annotation or verified genetic morph observation.**

## Empirical selection and real-photo retrieval

- Frozen source: 447,500 original images of 4,475 plant taxa; 227,107 four-state classifiable photos. Frozen compressed photo source SHA256: `8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6`.
- **27** previously observed reproducible photographic hue-PC1 density-mode candidate taxa versus **27** non-detected mode taxa matched 1:1 on source photo count, 50-km pair opportunity, global nine-hue JSD, PC1 variance, photographed colour category count and cohort. Controls were required to have >=40 original geovalid classified photos; they are *not* natural monomorphic controls.
- 3 photos near each of low/centre/high **continuous hue PC1** per taxon, all selected from different original observers when possible. An additional 2 original **unclassified** source photographs per taxon expose technical/missing-label failures. Total **594** selected authentic photo IDs: **486** classified + **108** unclassified.
- The source original measured CSV did **not** include direct image URLs. Using iNaturalist's [documented open-data photo ID and extension scheme](https://github.com/inaturalist/inaturalist-open-data), verifying the source photographic Creative Commons licence and requiring actual successful image decoding, **589/594** original-source images (99.16%) were recovered. **5** remained unavailable and were left unreplaced.
- Two independent, order-randomized reviewer CSV packs and actual downloaded thumbnail/contact-sheet files were produced. Reviewer sheets hide classification arm, taxon name, original photo ID, original colour, PC1 stratum, geography, observer and date. PI mapping retains source IDs, photo-page links, source licence, class/arm and continuous hue stratum. It must be distributed separately. The GitHub artifact separation is not a cryptographic access control or confidential vault.

Authenticated source acquisition / reviewer artifacts: https://github.com/zuizui0223/fcp/actions/runs/39931967531
- **Reviewer-only image sheets and blank annotations:** artifact ID **12492480632**.
- **PI source-photo mapping (do not give to blind reviewers):** artifact ID **12492660357**.
- **Metadata/status manifest:** artifact ID **12491273672**.

All images retain their original licences, including CC BY, CC BY-NC and other variants. Since attribution and reuse obligations vary, review source photo-page attribution before distributing image collages or publishing them. iNaturalist [usage requirements](https://help.inaturalist.org/en/support/solutions/articles/151000169918) apply.

## Independent technical photograph screen on all 589 recovered thumbnails

The separate [source-only QC run #39932618166](https://github.com/zuizui0223/fcp/actions/runs/39932618166) read 589 **actual recovered thumbnails**, not previous image-classification outputs.

| Criterion | 27 candidate taxa | 27 matched control taxa |
|---|---:|---:|
| Images successfully recovered | **296** | **293** |
| Mean whole-frame luminance | **0.4312** | **0.4382** |
| Mean whole-frame saturation | **0.3686** | **0.3759** |
| Strong bright-clipping >20% of pixels | **0** | **0** |
| Strong dark-clipping >20% | **1** | **0** |
| Very low centre-region luminance variation (SD<0.04) | **0** | **0** |

Source photo licences: 489 CC-BY-NC, 67 CC-BY, 15 CC-BY-NC-SA, 7 CC-BY-SA, 16 CC0 across 594 originally selected photo IDs. The five failed image retrievals were **four original unclassified photos and one central-PC1 classified photo**.

Comparing the *absolute* whole-frame luminance difference between low and high hue-PC1 photographic samples **within the 27 matched pairs** gave candidate-minus-control mean **−0.0194**, paired 499-bootstrap 95% interval **[−0.0461, +0.0053]**. This does not support a clearly larger *global exposure/luminance* contrast among candidate photographs, but it **does not rule out photographic colour artifacts**. Whole-frame brightness is not calibrated petal reflectance, and the analyses are exploratory.

All original thumbnail EXIF fields were stripped before anonymous images were exported. The image QC is nonbiological: no algorithm claimed true flower target organ validity, floral senescence absence, genetic polymorphism or pigment chemistry.

## Pending decisive human review

1. **Two independent masked reviewers** need to record genuine flower/petal visibility, petal vs sepal/bract/background, flower age, perceptual hue, lighting/white balance, cultivation/taxonomy concerns and confidence on each 589 real image; note the 5 missing records without imputation.
2. After both review passes are frozen, combine them to estimate agreement and remeasure how many **of the 27 photographed mode candidates** still have multiple supported *botanically verified* petal colour peaks. Repeat for the 27 matched controls.
3. Only species with botanically valid morph differentiation should proceed to population locality/time and genotype-based differentiation/fitness hypotheses. The previously unsupported strict 5-km same-month independent-observer matched null is not overturned here.

**Claim ceiling:** 27 photographed mode candidates and matched controls are prepared for high-quality biological inspection. The current biological outcome remains **unknown**. Photo review alone still would not establish inherited morph number, spatially stable breeding populations, local adaptation or pigment investment costs.

Machine-readable file: `results/fcp27_botanical_blind_review_verified_20261023.json`. Code: `scripts/analysis/build_fcp27_original_blind_photo_review_20261023.py`; `scripts/analysis/audit_fcp27_retrieved_photo_image_quality_20261023.py`.

