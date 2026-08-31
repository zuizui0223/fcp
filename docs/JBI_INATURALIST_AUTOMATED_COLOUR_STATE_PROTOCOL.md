# Automated, label-free iNaturalist colour-state protocol

Status: **development protocol v2 fixed before any spatial colour outcome is opened**  
Date: 2026-08-28  
Response: **within-species spatial organization of reproducible model-consensus flower-candidate colour states**

## 1. Why the response changed

No human morph annotation is used. The photograph stream therefore does not claim botanical morph identity, traditional colour names, population morph frequencies or correspondence to literature C/G labels. It estimates a narrower, fully automated response: a continuous colour feature extracted from an image region that is repeatedly identified as flower-like by a fixed zero-shot model under fixed prompts and perturbations.

The global map may still display every admitted continuous colour state. Inference remains within species and asks whether those marks are more spatially organized than expected after holding observation locations fixed. Species is retained in the null even when it is suppressed in the visualization.

## 2. Frozen model and inputs

- Input partition: the already fixed six-species `development_40` encounter set.
- Input files: licensed images and the location-free reviewer packet only. The extractor does not read coordinates, dates, observer identifiers, environment, literature state or the locked partition.
- Model: `CIDAS/clipseg-rd64-refined`.
- Hugging Face revision: `999e0328d9e10b484360c477313983f9afdd7050`.
- Safetensors SHA-256: `d00ca85d6b859f9d07b7cfb8ef26fe9771cb275b34c9368f2ecf603139307f55`.
- Software: `transformers==4.57.6`, `huggingface_hub==0.36.0`, `safetensors==0.6.2`; the validated host uses CPU `torch==2.11.0`.
- Positive prompts, in fixed order: `flower`, `petals`, `blossom`.
- Negative prompts, in fixed order: `leaves`, `background`.
- Perturbation: original image and horizontal reflection. No crop, colour correction, geographic bin or species-specific prompt is tuned.

CLIPSeg was designed for text- or image-prompted zero-shot segmentation, but its reported benchmarks do not validate flowers in this corpus ([Lüddecke & Ecker 2022](https://openaccess.thecvf.com/content/CVPR2022/html/Luddecke_Image_Segmentation_Using_Text_and_Image_Prompts_CVPR_2022_paper.html)). Its output is therefore called a **flower-candidate weight**, not a verified flower mask. Published citizen-science petal-colour pipelines demonstrate feasibility only after species-specific supervision and exclusions ([Perez-Udell et al. 2023](https://doi.org/10.1002/aps3.11505)); that stronger morph claim is not transferred here.

The upstream CLIPSeg repository states that its source-code license does not cover model weights, while the pinned model card reports Apache-2.0. The weights remain in a private local cache and are not redistributed in the repository or document bundle.

## 3. Flower-candidate weight

For each positive prompt (p), let (z_p(u)) be the CLIPSeg logit at pixel (u), and let (z_N(u)) be the maximum logit across the two negative prompts. Define

\[
w_p(u)=\max\{0,\;2[\operatorname{logit}^{-1}(z_p(u)-z_N(u))-0.5]\}.
\]

The ensemble flower-candidate weight is the mean of the three (w_p(u)). The subtraction makes a pixel contribute zero unless a positive prompt outranks the strongest negative prompt. No binary contour or hand-selected threshold is used for colour extraction.

The original RGB image is resized only to the model's declared 352 × 352 input. Weighted CIELAB means, standard deviations and 0.1/0.5/0.9 quantiles are then computed. These are device-dependent ordinary-sRGB descriptors, not calibrated reflectance, UV colour or pollinator colour space.

## 4. Fixed automatic admission gate

A photograph is `automated_colour_state_admitted` only when all conditions hold:

1. total ensemble effective weight is at least 100 model-input pixels;
2. at least two of the three positive prompts independently have effective weight of at least 100 pixels;
3. the maximum Euclidean distance between positive-prompt weighted mean Lab vectors is at most 10;
4. the Lab distance between original and reflected ensemble means is at most 5;
5. original versus unreflected reflected-weight soft IoU is at least 0.50;
6. all flower-candidate model outputs and extracted flower features are finite.

The background negative control is diagnostic rather than part of flower-photo admission. If its effective weight is zero, its features are stored as explicit missing values, never as zeros or non-standard JSON `NaN`. This rule and the later background-completeness gate are included in the contract hash.

These values are image-measurement stability gates, not perceptual just-noticeable-difference claims. A failure is retained as `automated_colour_state_not_evaluable`; no mask filling, prompt replacement, threshold relaxation, geographic rescue or post-result category merging is allowed.

An encounter is admitted when at least one attached photograph passes. Its feature is the componentwise median across passing photographs. Within-encounter photo-to-photo Lab distance is retained as an outcome-blind repeatability diagnostic.

## 5. Development decisions without human truth

Before coordinates are joined, the six species must pass all of the following:

- at least 70% of the 80 development encounters admit an automated colour state;
- at least 10 encounters contain two or more admitted photographs, and their median within-encounter mean-Lab distance is lower than 9,999 species-preserving photo-to-encounter permutations at one-sided \(p<0.05\), using seed `fcp-inaturalist-automated-colour-state-v2`;
- fixed five-fold cross-validation of a ridge model (`alpha=1`) using whole-image luminance, dark and bright clipping fractions, log pixel count, log aspect ratio and log edge variance must have multivariate out-of-fold \(R^2<0.80\) for the encounter mean-Lab vector. A value at or above 0.80 is treated as near-perfect technical prediction and makes that species `not_evaluable`;
- at least 70% of admitted development encounters retain a finite negative/background candidate vector for the later spatial negative control.

The fold assignment is a deterministic hash of species, encounter blind ID and the same fixed seed. Predictors and outcomes are standardized from each training fold only. A species passes the development gate only if every condition passes; missing or undefined diagnostics fail closed.

The development set may establish feasibility and freeze transformations. It cannot estimate spatial randomness. Coordinates remain sealed until the complete development artifact, code, model revision, prompts, gates, seeds and hashes are frozen.

## 6. Spatial test after the image gate

The complete development gate is evaluated before any colour feature is joined to coordinates. A passing species is then processed once in its already frozen 120-encounter `locked_60` partition, after expanding every reusable photo belonging to each encounter. Development and locked encounters are not pooled for confirmatory inference.

The locked partition independently requires at least 70% (84/120) admitted encounters and a finite background vector for at least 70% of admitted encounters. A failed locked completeness gate is `not_evaluable`; development and locked thresholds are never pooled or averaged.

For a locked species, the encounter mean-Lab vector is componentwise centred by the locked median and divided by the locked interquartile range. A zero or non-finite interquartile range is `not_evaluable`; no alternative scaling is substituted. Great-circle distance is calculated between the fixed public coordinates. Exact duplicate coordinates are retained because they are valid local co-detections, but obscured/private coordinates remain excluded.

The primary statistic is the Spearman correlation between all pairwise great-circle distances and all pairwise Euclidean standardized mean-Lab distances. A positive value represents greater colour-state distance at greater geographic separation. Its one-sided random-mark null uses 9,999 permutations of entire three-channel encounter vectors among the fixed locations within that species, with seed `fcp-inaturalist-automated-colour-state-v2`. Six-species multiplicity is controlled by Benjamini-Hochberg at `q=0.05`; a non-rejection is reported as **no detected spatial organization**, not proof of randomness.

Two predeclared protections are required for an affirmative species result:

1. the primary BH-adjusted result is significant and the sign remains positive after excluding the single most-contributing observer;
2. on encounters with a finite background vector, the flower-minus-background correlation contrast is positive and significant under 9,999 joint paired-vector permutations, also BH-adjusted at `q=0.05`.

Five equal-pair-count geographic-distance bins are reported as a descriptive variogram only; they do not create additional tests or tune a distance threshold. Flowering week, whole-image quality variables, observer share and 50-km cell support are reported as observation-process diagnostics. The frozen development technical-prediction gate and top-observer cap are prerequisites; no post hoc cell, season, device or distance filter may rescue a result.

Cross-species boundary stacking remains sealed until species-conditioned locked tests pass. This six-species development/locked panel can show repeated organization in selected species but cannot establish universality; that requires an independently frozen held-out species panel.

## 7. Claim ceiling

Allowed if all automatic gates pass:

> Reproducible, model-consensus flower-candidate colour states do or do not reject a fixed species-conditioned spatial random-mark null in selected species.

Not allowed without external validation:

- the extracted region is botanically verified flower tissue;
- an extracted state is a discrete flower-colour morph;
- any named colour, pigment, reflectance or pollinator-perceived state is recovered;
- photographed encounters estimate population frequencies or morph absence;
- spatial organization identifies selection, maintenance, barriers, history or environmental adaptation.
