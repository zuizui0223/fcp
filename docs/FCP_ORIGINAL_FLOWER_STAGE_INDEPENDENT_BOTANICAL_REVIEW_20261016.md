# Original iNaturalist flower-stage validation before seasonal white-frequency estimation

Status: 2026-10-16. Scientific provenance: the genuine [2026-09-27 original-source 1,000-taxon census](https://github.com/zuizui0223/fcp/actions/runs/38884223530) and [first genuine original ROI photo unit](https://github.com/zuizui0223/fcp/actions/runs/38901816890). These establish photo IDs and automated ROI classes, **not** that each ROI is a botanically comparable open flower.

## Phenology is missing from the provider archive

iNaturalist's official monthly Open Dataset has four bulk metadata tables: observations, photos, observers and taxa. It includes `observed_on`, but no plant phenology annotation. The public iNaturalist API also supports optional observation annotations: term_id=12 (Plant Phenology), value 13=Flowering, 14=Fruiting, 15=Flower Budding and 21=No Evidence of Flowering. These refer to an **observation**, possibly with several different photographs, and may be absent, incomplete, or different from the pictured structure. They are supplementary evidence; never use them by themselves to assign an open flower to an individual photo.

## Mandatory photographic checks

Every selected original photo, including classifier failures and suspected same-ROI white+yellow cases, receives **independent** stage and ROI-organ review from two different reviewers, blinded to predicted four-state colour label, geolocation/season, outcome model and climate values. The original image URL and taxon ID/name are provided for botanical interpretation.

Predeclared stages: `open_flower`, `buds_only`, `fruit_only`, `no_flower_visible`, `uncertain`. ROI tissues: `petal_tepal`, `other_open_flower_organ`, `spikelet_glume_bract`, `fruit_seed`, `leaf_stem`, `background`, `uncertain`.

For the primary comparable white/chromatic flowering photo outcome, require two independent reviewers to agree that a fully open flower is visible **and** that the original colour-extraction region represents petals/tepals (or a separately predeclared homologous display structure, if analysed as a secondary stratum). Flag spikelet bracts/glumes, fruits, and ambiguous mixed floral organs separately, not as a white/colour alternative. Inter-reviewer disagreement requires a third independent adjudicator. Unknown/absent review = `unreviewed`, **never a biological negative**.

Even when open anthesis is verified, failed/unclassifiable automatic ROI-v4 colour output is not a white or chromatic outcome. Repeated photographs from the same original observation are not independent plants, and separate observation records do not guarantee distinct individual plants.

## Current actual original 64-photo technical evidence

The *Cyperus glomeratus* + *Leptomeria drupacea* original-image unit measured 64 photos, 27 with coarse automatic four-state class, 4 with same-ROI white/yellow organ mixing and 37 unclassifiable images. In *C. glomeratus*, yellow/orange labels typically have large bronze contributions from the conspicuous spikelet/bract region; this does not establish petal pigmentation. *L. drupacea* has tiny pale flowers and differently coloured fruits; white ROIs are extremely size-variable. The archive has 64 original photo URLs and observation UUIDs but **zero independently verified per-photo anthesis/ROI tissue judgements**, hence no botanically established seasonal white rate can be reported from it.

A blank photo-ID/source-verified two-reviewer sheet is available as conversation artifact `first64_open_flower_review_template.csv`; reviewer fields are intentionally blank, and it does not expose the ROI colour/morph label. The tested local validator records all 64 as `unreviewed` until two independent human review records are supplied, retaining explicit uncertainty rather than imputing a no-flowering outcome.

## Observer-disjoint seasonal photo wave

From all 3,430,502 source photo IDs, the real dated source census supplies 268 high-depth matched site-year/taxon comparisons across 39 taxa (minimum 24 independent source observations per season). The independent observer audit showed only one of those candidate site/years supports at least 64 **non-overlapping original observer IDs in each of two seasons**: *Asclepias syriaca* (2025; 0.25-degree cell lat 529 / lon 435), with source autumn/summer observers 80/75 after excluding six reporters active in both seasons. All 128 pre-colour selected photo IDs have different source observation IDs and different observers. Frozen first wave = 32 summer and 32 autumn originals; reserve wave = another 32+32. The original 106,376-image 1,000-taxon randomized repeated-sample cohort is unchanged.

Exact production source selection (zero newly opened pixels): [run #38908308446](https://github.com/zuizui0223/fcp/actions/runs/38908308446); [five passing regression tests](https://github.com/zuizui0223/fcp/actions/runs/38908006010). The first real 64-photo matched-season original image ROI workflow is [run #38909164068](https://github.com/zuizui0223/fcp/actions/runs/38909164068). Do not treat its raw pixel colour labels as verified flowering until per-photo independent botanical QC is complete.

## Inference after open-flower stage pass

Only stage-confirmed, classifiable, observation-independent source photos with site/season/year overlap may enter a descriptive white rate. Report (i) all original photo opportunities (including unreviewed/missing dates), (ii) independently stage-checked images, (iii) genuine open flower ROI photos, (iv) valid colour classes, (v) original independent observations and observers, by taxon×site×season. Any size or direction of seasonal white contrast must be checked for differential flowering visibility, fruiting, image quality and observer coverage. Real observed-month TerraClimate is an additional explanatory layer, not needed to assert a photograph is open flower; do not fill 2026 weather with 2025 or long-term averages.