# A location-blind automated image workflow did not detect flower-colour spatial organization in three iNaturalist species

## Abstract

Large citizen-science image archives could reveal how continuous flower colour varies across geography, but a global atlas is credible only if image measurement is separated from location and if failed measurements remain visible. We hash-froze a fully automated, location-blind workflow before the complete development and spatial outcomes. It used a pinned CLIPSeg model and five fixed text prompts to derive soft flower-candidate weights and continuous CIELAB summaries. Six species, 480 development encounters and 886 photographs were used only for image-measurement gates. Three species passed all gates and entered an independent locked partition of 360 encounters and 717 photographs. All three passed the locked completeness gate, yielding 306 admitted encounters. Coordinates were joined only after those decisions. Within each species, we tested the Spearman association between all pairwise great-circle distances and robust-standardized colour distances using 9,999 whole-vector random-mark permutations, Benjamini–Hochberg correction, a leave-top-observer sensitivity and a flower-minus-background negative control. Spatial organization was not detected in *Erythranthe lewisii* (rho = 0.0077, q = 0.5969), *Hesperis matronalis* (rho = -0.0582, q = 0.9637) or *Orchis mascula* (rho = 0.0200, q = 0.5969). No species passed the combined primary, observer and background-control rule. These non-rejections do not prove spatial randomness. They show that a reproducible zero-shot flower-candidate measurement can complete strict image gates without yielding confirmatory spatial evidence. We therefore stop the planned 20,200-image atlas before bulk pixel opening and require an independent flower-tissue localization benchmark and signal-recovery analysis before scale-out.

**Keywords:** citizen science; colour; fail-closed workflow; iNaturalist; location blindness; negative result; random-mark test; zero-shot segmentation

## 1. Introduction

Georeferenced photographs make it possible to consider flower-colour variation as a continuous spatial field rather than as a small set of literature-derived categories. The attraction is substantial: a species-free public map can display colour directly, while species-conditioned inference can ask whether transition surfaces recur across taxa. The main risk is equally substantial. Camera processing, background vegetation, flower localization and uneven sampling can create geographic structure before any biological signal is measured.

We designed the image-first FCP atlas around ordered gates. Metadata and geometry are frozen before image pixels; image measurement is completed without coordinates; coordinates are released only for species passing fixed completeness and technical gates; and every null permutation remains within species. This paper reports the first independent locked test of that measurement route. It is an admissibility study for the atlas, not the atlas itself.

The measurement model is CLIPSeg, a text-prompted zero-shot segmentation model (Lüddecke & Ecker, 2022). Its output has not been validated as botanical flower tissue in this corpus. We consequently call the output a model-consensus **flower-candidate weight** and analyse continuous CIELAB summaries, not named colours or biological morphs. We asked whether species that pass location-free measurement gates subsequently reject a preregistered species-conditioned random-mark null. The stopping rule was symmetric: a supportive result would justify the next atlas gate; non-detection would be retained and would stop bulk image measurement until the estimator itself was independently benchmarked.

## 2. Materials and methods

### 2.1 Frozen design and information firewall

The complete analysis contract was frozen before the complete development outcome, before the locked partition was opened and before any spatial colour outcome was inspected (`fcp-inaturalist-automated-colour-state-v2`). Human judgement was not used for flower localization, colour extraction, admission or inference.

The development partition contained 80 encounters per species for six species (480 encounters; 886 attached photographs). It was location-free: exact coordinates, observer identity and spatial summaries were unavailable to the measurement and gate code. Only species passing every development criterion were eligible for the locked partition. The locked partition contained 120 new encounters per passing species and was never pooled with development data.

### 2.2 Automated flower-candidate measurement

We pinned `CIDAS/clipseg-rd64-refined` at revision `999e0328d9e10b484360c477313983f9afdd7050`; the `model.safetensors` SHA-256 was `d00ca85d6b859f9d07b7cfb8ef26fe9771cb275b34c9368f2ecf603139307f55`. Model input was 352 × 352 pixels. Positive prompts were `flower`, `petals` and `blossom`; negative prompts were `leaves` and `background`.

At each pixel, a positive prompt weight was its positive logit contrasted with the maximum negative-prompt logit and transformed to [0, 1]. The flower-candidate weight was the mean of the three positive weights. Continuous CIELAB mean, standard deviation and weighted 10th, 50th and 90th quantiles were calculated from this soft weight. Encounter features were componentwise medians across all admitted photographs attached to the encounter.

A photograph was admitted only when it had at least 100 effective weighted pixels; at least two valid positive prompts; maximum prompt-to-prompt mean-Lab distance no greater than 10; original-versus-horizontal-reflection mean-Lab distance no greater than 5; reflection soft intersection-over-union at least 0.50; and finite flower features. Every failure received a machine-readable reason. No manual rescue or post-outcome threshold change was permitted.

### 2.3 Development and locked gates

A species passed development only if at least 70% of its 80 encounters were admitted, at least ten multi-photo encounters were available, within-encounter repeatability was significant in a 9,999-permutation lower-tail test, five-fold ridge prediction of colour from technical image descriptors had cross-validated R-squared no greater than 0.80, and background features were available for at least 70% of admitted encounters.

The locked completeness rule required exactly 120 encounters, at least 84 admitted encounters and background features for at least 70% of admitted encounters. Only after a species passed was its coordinate and observer table joined to the already frozen encounter-colour table. Species failing the rule would retain blank coordinates and become `not_evaluable`.

### 2.4 Species-conditioned spatial test

For each evaluable species, CIELAB encounter means were robust-standardized by median and interquartile range. We calculated great-circle distance and Euclidean standardized-colour distance for every encounter pair. The primary statistic was Spearman's rho between the two distance vectors; positive values indicate that geographically distant encounters tend to be more different in measured colour.

The null reassigned complete three-component colour vectors among coordinates within species and recomputed the statistic 9,999 times. The one-sided probability was `(1 + null >= observed) / 10,000`. Probabilities were corrected across evaluable species with the Benjamini–Hochberg procedure at q = 0.05.

A species was called supported only if the primary rho was positive with q < 0.05, rho remained positive after excluding the single most-contributing observer, and the flower-minus-background rho contrast was positive with q < 0.05 under paired whole-vector permutations. Five equal-pair-count distance bins were retained for descriptive visualization only.

### 2.5 Display–inference separation

The public map and photo bar omit species labels. The statistical input retains species and all permutations are within species. The shareable analysis table contains blind encounter identifiers, continuous features and coordinates for admitted encounters. Observer identifiers are replaced with species-scoped, rank-preserving pseudonyms; this preserves the exact leave-top-observer calculation, including deterministic tie handling, without publishing source identifiers.

## 3. Results

### 3.1 Location-free development gate

Three of six species passed all image-measurement criteria. *Erythranthe lewisii* admitted 73/80 encounters (91.25%), *Hesperis matronalis* 68/80 (85.0%) and *Orchis mascula* 77/80 (96.25%). *Digitalis purpurea* (61.25%), *Hepatica nobilis* (62.5%) and *Protea repens* (60.0%) failed the fixed 70% completeness rule; *H. nobilis* also lacked the required number of multi-photo encounters. These failures were retained as `not_evaluable` and their locked partitions were not opened.

### 3.2 Locked image measurement

The three eligible species contributed 360 locked encounters and 717 photographs. All 717 cache records passed the independent integrity audit, with no missing, partial or unexpected record. Automated measurement admitted 508 photographs and aggregated 306 encounters. Locked encounter admission was 101/120 (84.17%) for *E. lewisii*, 97/120 (80.83%) for *H. matronalis* and 108/120 (90.0%) for *O. mascula*. Background controls were available for every admitted encounter. All three species therefore passed the locked gate before coordinates were opened.

### 3.3 Locked spatial test

No species rejected the random-mark null (Figure 2). *Erythranthe lewisii* had rho = 0.00774, one-sided p = 0.3979 and BH q = 0.59685. *Hesperis matronalis* had rho = -0.05822, p = 0.9637 and q = 0.9637. *Orchis mascula* had rho = 0.01996, p = 0.3364 and q = 0.59685.

Observer sensitivity did not rescue support: leave-top-observer rho was 0.00149, -0.05651 and 0.01206, respectively. Flower-minus-background rho was 0.04356 (q = 0.5106), -0.10022 (q = 0.9664) and -0.08759 (q = 0.9664). Thus all three fixed combined decisions were `spatial_organization_not_detected`.

## 4. Discussion

The workflow achieved its engineering objective but not its confirmatory ecological objective. A pinned zero-shot model, deterministic soft weighting, reflection and prompt stability checks, a coordinate firewall and complete cache audit produced a reproducible locked table. Nevertheless, none of the three independently tested species showed spatial organization under the fixed combined rule.

This is a non-detection, not evidence that flower colour is spatially random. At least three explanations remain compatible with the result: spatial organization may be weak in these species; the sampling extent and encounter count may be insufficient for this distance-based statistic; or the flower-candidate estimator may attenuate real flower-colour differences by including non-flower pixels. The current design cannot identify which explanation is correct because it deliberately avoided post-outcome manual annotation.

The result does not overturn the frozen six-species Chapter 1 analysis, which used different species, a different localization/measurement route and a nearest-neighbour statistic and reported within-species organization. The two studies must not be pooled. Instead, their disagreement identifies the next methodological requirement: an independently annotated flower-tissue benchmark and prespecified signal-recovery simulation that quantify attenuation and false structure for the exact atlas estimator.

The metadata-only 50-species atlas cohort and its 100/250/500-km geometry remain frozen. Its 20,200 candidate images remain unopened. Bulk measurement, spatial-field reconstruction and cross-species shared-boundary concentration are stopped until the estimator benchmark passes. This protects the atlas from converting scalable segmentation into scalable measurement error.

## 5. Claim ceiling

The supported result is limited to the following statement:

> In three development-passing species, model-consensus flower-candidate colour states did not reject a frozen, locked, species-conditioned random-mark null.

The study does not verify flower tissue, botanical morph identity, named colour classes, population frequencies, mechanism, adaptation, universal randomness or absence of local transition boundaries.

## Data and code availability

The repository contains the frozen protocol, model and source hashes, development and locked gate manifests, a privacy-preserving 360-row locked analysis table, the complete species-level results, descriptive distance bins, figure source data, analysis scripts and automated tests. Model weights and copyrighted source photographs are not redistributed. The figure manifest retains source-image hashes, licences and attributions for displayed crops.

## Figure legends

**Figure 1. Species-free display of the locked automated-colour pilot.** Each map point is one of 306 admitted locked encounters and is rendered from its model-consensus flower-candidate mean CIELAB value. The 24 photographs below are selected deterministically at equally spaced positions after ordering admitted encounters by longitude; the first admitted photograph blind ID for each encounter is cropped to the central 90% of CLIPSeg flower-candidate soft-mask mass with 12% display padding. Species labels are omitted from this display but retained in the inferential table. The map and photo bar are descriptive and do not establish a spatial boundary. All displayed photographs are CC0 or CC BY; hashes, licences and attributions are in the figure manifest.

**Figure 2. Ordered image and spatial gates.** (a) Location-free admitted-encounter share for the six development species; the dashed line is the frozen 70% minimum and blue bars denote species passing every development condition. (b) Primary Spearman correlation between pairwise great-circle distance and robust-standardized model-consensus colour distance in the three locked species; labels show Benjamini–Hochberg q values after 9,999 within-species whole-vector permutations. (c) Flower-minus-background rho negative control and its corrected q value. No species met the frozen combined support rule.

## Reference

Lüddecke, T. & Ecker, A. (2022). Image segmentation using text and image prompts. *Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern Recognition*, 7086–7096. <https://openaccess.thecvf.com/content/CVPR2022/html/Luddecke_Image_Segmentation_Using_Text_and_Image_Prompts_CVPR_2022_paper.html>
