# FCP atlas flower-ROI v3: prospective estimator and validation candidates

Date: 2026-08-31  
Scope: prospective work only. The frozen CLIPSeg/Oxford-17 v2 result (`676/848 = 0.7971698`, below the fixed `0.80` admitted-fraction gate) is not changed, retuned, or rerun.

## Decision

The minimal v3 candidate is the ADE20K flower class from **SegFormer-B0**, pinned as follows:

- model repository: [`nvidia/segformer-b0-finetuned-ade-512-512`](https://huggingface.co/nvidia/segformer-b0-finetuned-ade-512-512/tree/489d5cd81a0b59fab9b7ea758d3548ebe99677da);
- model revision: `489d5cd81a0b59fab9b7ea758d3548ebe99677da`;
- weight file: `model.safetensors`, 15,036,944 bytes;
- weight SHA-256: `6ae39addd01de6b1b8bde2cf677d43a5cd733424b8d186de3f95d1c51fee23f9` (also reported by the pinned Hugging Face blob metadata);
- model configuration SHA-256: `209caa9091e4632f7c8883c11170cd08ad29af68b23c09590aa4a5befb1a2a7f`;
- preprocessor configuration SHA-256: `8039d1d210abaa7117ad78e58cdfd6141a2ec72c03dae891b3cd76737e422c6c`;
- output label: `66 = flower`; `17 = plant` is a negative-control label only and must never rescue or enlarge a failed flower ROI. The pinned [configuration](https://huggingface.co/nvidia/segformer-b0-finetuned-ade-512-512/blob/489d5cd81a0b59fab9b7ea758d3548ebe99677da/config.json) lists both labels.

SegFormer-B0 is the only candidate here with a pretrained, explicit semantic `flower` output and no need for a text prompt, point prompt, automatically selected mask, or flower-versus-leaf heuristic. The model card reports 3.75 million parameters and explains that the encoder was pretrained on ImageNet-1K and the full model fine-tuned on ADE20K; the official paper reports B0 as the smallest SegFormer variant ([model card](https://huggingface.co/nvidia/segformer-b0-finetuned-ade-512-512), [ICCV paper](https://openaccess.thecvf.com/content/ICCV2021/papers/Xie_SegFormer_Simple_and_Efficient_Design_for_Semantic_Segmentation_With_Transformers_ICCV_2021_paper.pdf)). Its small 15 MB safetensors checkpoint makes CPU inference in sharded GitHub Actions plausible, but a timed CPU smoke job must still pass before a full locked benchmark is dispatched. Model size is not runtime evidence.

This is a new estimator family, not a changed threshold for v2. Oxford-17 remains closed.

## Frozen v3 pre/post-processing recommendation

Freeze these operations, dependency versions, and hashes before running any locked image:

1. Apply EXIF orientation, convert once to RGB, and preserve the oriented original dimensions for the final upsampling and colour calculation.
2. Use the pinned SegFormer preprocessor unchanged: resize to `512 x 512` with PIL resample code `2` (bilinear), rescale image bytes, and normalize by mean `(0.485, 0.456, 0.406)` and standard deviation `(0.229, 0.224, 0.225)`. The exact values are in the pinned [`preprocessor_config.json`](https://huggingface.co/nvidia/segformer-b0-finetuned-ade-512-512/blob/489d5cd81a0b59fab9b7ea758d3548ebe99677da/preprocessor_config.json).
3. Run one image at a time in evaluation/inference mode on CPU. Upsample the 150-channel logits to the oriented original height and width using bilinear interpolation with `align_corners=False`, then take `argmax` across labels. This follows the model's semantic-class output and avoids a new probability threshold.
4. Define the flower ROI as every pixel whose hard class is exactly `66`. Retain all connected components. Do not use dilation, erosion, hole filling, largest-component selection, a centre prior, a colour prior, or an image-specific threshold.
5. Preserve the existing operational minimum of 100 flower pixels. An image with fewer than 100 class-66 pixels is not admitted; it is not replaced and receives no fallback crop.
6. Run the identical estimator on a horizontal flip, unflip the resulting mask, and retain the already used stability ceilings prospectively: hard-mask IoU at least `0.50` and CIELAB mean-colour difference no greater than `5.0`. Use the original-orientation mask for the recorded colour; the flipped result is a diagnostic only.
7. Define the background-control mask as hard class `17` (`plant`) with at least 100 pixels. It can support the frozen flower-minus-background control but can never count as flower, make a photo admitted, or be replaced by the whole non-flower complement.
8. Keep `minimum_admitted_fraction = 0.80` for estimator qualification. A result of `0.799...` is a failure, not a rounding pass. The atlas-scale completeness thresholds remain separate and unchanged.

The hard-mask rule deliberately sacrifices possible recall to remove the post-hoc degrees of freedom that promptable mask generators and soft cutoffs would introduce. If any of these rules are changed after development images are viewed, the changed rule is a new estimator version and needs a newly untouched locked set.

## What Oxford Flowers 102 can and cannot validate

The official [Oxford Flowers 102 page](https://www.robots.ox.ac.uk/~vgg/data/flowers/102/) contains 8,189 images, segmentation files, class labels, and fixed split identifiers. The fixed splits contain 1,020 training, 1,020 validation, and 6,149 test images. The [Oxford dataset overview](https://www.robots.ox.ac.uk/~vgg/data/flowers/) states that both the 17- and 102-category collections were gathered from various websites, supplemented by the authors' photographs.

There are two reasons not to call Oxford-102 an independent human-mask benchmark:

1. The official 2008 paper says the Oxford-102 image segmentations were generated by the authors' automatic iterative MRF/shape procedure; they are algorithm outputs, not new manual pixel truth ([Nilsback and Zisserman 2008](https://www.robots.ox.ac.uk/~vgg/publications/2008/Nilsback08/nilsback08.pdf), section 2). A later Oxford paper explicitly notes empty or erroneous Flowers masks and describes them as automatically generated rather than human annotations ([Melaskyriazi et al. 2022](https://www.robots.ox.ac.uk/~vgg/publications/2022/Melaskyriazi22b/melaskyriazi22b.pdf), appendix B.1).
2. The two Oxford collections share a web-image acquisition lineage and contain exact duplicates. An audit of the current official archives (`17flowers.tgz` SHA-256 `fe38a60f8b4a95e657551247d3e7d799a3fafcdbc595be504b12839967823d70`; `102flowers.tgz` SHA-256 `2d01ecc807db462958cfe3d92f57a8c252b4abd240eb955770201e45f783b246`) found five Oxford-102 images identical byte-for-byte to images in the full 1,360-image Oxford-17 archive: IDs `3448`, `3456`, `4657`, `4691`, and `6241`. Two are in the Oxford-102 train/validation split and three are in its test split. Relative to the 848 images actually scored in the frozen v2 trimap gate, two exact overlaps are present; the full-archive and v2-subset counts must not be conflated.

Therefore use Oxford-102 as follows:

- development: official `trnid + valid` only, after dropping exact duplicates against the entire Oxford-17 archive and recording every dropped hash;
- locked algorithm-agreement set: select exactly 20 images per class from official `tstid` by a declared SHA-256 rank seed, after the same whole-archive deduplication; write the 2,040 image IDs and source hashes before inference;
- never call the automatic Oxford masks human ground truth;
- never let an Oxford-only pass authorize atlas pixel opening. It establishes agreement with a historical automated flower-foreground proxy on previously unused images, not biological tissue validity.

The 20-per-class design prevents classes with many images from dominating. No class, image, or mask may be removed after predictions are computed. An empty or malformed official proxy mask remains a recorded `source_mask_not_evaluable`; the exact denominator and per-class counts remain visible.

## Independent field localization gate

The most practical unused independent source is the European Commission Joint Research Centre **Flower Detection** dataset, DOI [`10.2905/JRC.2XJ67GR`](https://data.jrc.ec.europa.eu/dataset/caa582b7-7f45-4748-9223-08e5f145a4a6). Its official [README](https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/DRLL/flower_detection/readme.txt) says:

- 500 grassland-vegetation patch photographs from the 2018 LUCAS grassland module;
- all visible flowers manually annotated with bounding boxes in COCO format;
- an official 400-image training / 100-image test split;
- one object category, `flower`.

The images were collected by the LUCAS survey rather than scraped from the Oxford web sources. The official [copyright notice](https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/DRLL/flower_detection/copyright.txt) applies CC BY 4.0. Freeze its current files and hashes before use.

Use all 400 official training images only for code smoke testing and, if absolutely required, selection of a single v3 component rule before the contract is hashed. Keep all 100 official test images sealed and run them once. Because the truth is boxes rather than pixel masks, report at least:

- exact test-image admission fraction;
- object recall stratified by predeclared box-size bins;
- fraction of predicted flower pixels inside the union of all annotated flower boxes;
- per-image and pooled results, with every zero prediction retained.

The numeric box-hit definition and all acceptance cutoffs must be written from the 400-image development set and scientific measurement requirements before the 100 test images are opened. They must not be chosen from test performance.

A JRC pass demonstrates independent field-photo flower localization. It does **not** by itself validate pixel-accurate petal boundaries because boxes contain background.

## Independent pixel-mask options

Two genuinely independent mask sources were found, but neither is currently a drop-in substitute for the two-gate design above.

- **FlowerMaskDataset v1**, DOI [`10.17632/3pw57gdcj2.1`](https://data.mendeley.com/datasets/3pw57gdcj2/1), reports 3,600 natural-condition images from six flower species, each with a manual LabelMe polygon, and a CC BY 4.0 licence. This is the best prospective target-domain pixel benchmark. However, the repository page currently exposes no enumerated files and says to email for the complete dataset. It cannot enter CI until a complete immutable archive is actually obtained, hash-frozen, licence-checked, and made reproducibly accessible. Do not infer availability from the DOI alone.
- **Senckenberg plant-organ segmentation**, DOI [`10.12761/fj4m-zr97`](https://dataportal.senckenberg.de/dataset/plant-organ-segmentation-images-and-annotations-on-digitized-herbarium-scans), provides CC BY 4.0 flower, fruit, leaf, stem, seed, and root masks. It is valuable for flower-versus-other-organ specificity, but the images are digitized herbarium sheets. Passing it cannot establish performance on natural iNaturalist photographs.

If the complete FlowerMaskDataset archive becomes available **before any of its images are scored**, freeze a species-stratified development/lock split by SHA-256: 400 development and 200 locked images per species. The 1,200 locked masks should then be the primary pixel gate, with JRC test remaining the independent field-localization gate. If that archive remains unavailable, record the independent pixel gate as unavailable; do not silently promote Oxford's automatic masks to manual truth.

## Other model candidates

| Candidate | Exact usable artifact | Licence and practicality | Reason not selected for minimal v3 |
|---|---|---|---|
| U2NETP | [`xuebinqin/U-2-Net@ac7e1c817ecab7c7dff5ce6b1abba61cd213ff29`](https://github.com/xuebinqin/U-2-Net/tree/ac7e1c817ecab7c7dff5ce6b1abba61cd213ff29), `u2netp.pth` advertised as 4.7 MB | Apache-2.0; small. The paper describes general salient-object detection and DUTS training ([paper](https://arxiv.org/abs/2005.09007), [DUTS source](https://saliencydetection.net/duts/)). The checkpoint is served from a mutable Google Drive URL, so it requires an acquired-file SHA. | It segments the salient foreground, not specifically flower tissue, and can include leaves, stems, hands, or insects. |
| EfficientSAM-Ti | [`yformer/EfficientSAM@d525f622e6f640acf5a0fc37c7ca1f243da5bde0`](https://github.com/yformer/EfficientSAM/tree/d525f622e6f640acf5a0fc37c7ca1f243da5bde0), `weights/efficient_sam_vitt.onnx` | Apache-2.0; ONNX CPU inference is supported by the official repository. | It needs a point/box/grid prompt and a mask-selection rule. Those extra choices can become post-hoc tuning unless a complete second-candidate contract is frozen now. |
| MobileSAM | [`ChaoningZhang/MobileSAM@f706ad9c4eb7f219c00d9050e46328518ffb65d2`](https://github.com/ChaoningZhang/MobileSAM/tree/f706ad9c4eb7f219c00d9050e46328518ffb65d2), `weights/mobile_sam.pt` | Apache-2.0; the official README reports about three seconds per image on a Mac i5 and documents ONNX export. | Its automatic generator returns multiple category-free masks; selecting the flower still needs a separate frozen semantic rule. |
| IS-Net / DIS | [`xuebinqin/DIS@b6764e20381f6f42a70f83fa3324181529ed1403`](https://github.com/xuebinqin/DIS/tree/b6764e20381f6f42a70f83fa3324181529ed1403) | Apache-2.0; official general-use checkpoint distributed through Google Drive. | Larger, general foreground extraction; no explicit flower semantics and a less convenient immutable checkpoint channel. |

CLIPSeg is absent from this candidate hierarchy because it already has a frozen v2 result. Changing its threshold, prompts, admission rounding, or Oxford-17 selection would be post-outcome retuning.

## Training-overlap risk

SegFormer declares ImageNet-1K pretraining and ADE20K fine-tuning, not Oxford training. That does not prove image-level independence. Oxford's own page states that its images came from websites, and public documentation does not provide a complete duplicate audit against ImageNet or ADE20K. Record estimator/benchmark training overlap as `not_documented`, not `none`.

The same issue is stronger for SAM-family models trained on large web-scale image collections. A benchmark can be unused by FCP while still having unknown exposure to the upstream model. Exact and perceptual duplicate audits should be run wherever both source archives are available, but an inability to audit upstream training pixels remains a claim limitation rather than evidence of no overlap.

## Fail-closed decision tree

1. Freeze source hashes, the Oxford-102 split IDs, JRC official split, model revision, weight/config hashes, dependency lock, pre/post-processing, per-image admission, aggregate metrics, thresholds, and output schema.
2. Run a small CPU timing and determinism smoke test on development images. A timeout, nondeterministic mask, missing source, or hash mismatch is `not_evaluable`, not a scientific failure.
3. Finish development only on Oxford-102 train+validation and JRC train. Then hash the executable contract. Locked images remain unopened.
4. Run the Oxford-102 balanced locked proxy set and all 100 JRC test images once. If a frozen gate fails, save the full rows and return `stop_roi_v3_failed`. Do not adjust morphology, add `plant`, change the label, drop difficult classes, or round the denominator.
5. Atlas image acquisition can proceed only if the repository's final v3 contract explicitly states whether these proxy/localization gates are sufficient or whether an independent natural-photo pixel-mask gate is mandatory. If the latter is mandatory and FlowerMask is not reproducibly available, atlas pixel opening remains blocked.
6. Any second estimator must already have its complete independent contract and untouched locked set. It cannot be designed after observing where SegFormer fails.

## Claim ceiling

The SegFormer code licence is the [NVIDIA Source Code License for SegFormer](https://github.com/NVlabs/SegFormer/blob/65fa8cfa9b52b6ee7e8897a98705abf8570f9e32/LICENSE), not Apache or MIT. It permits research/evaluation use and restricts the work and derivatives to non-commercial use. The Hugging Face model card identifies the licence as `other` and points to this text. Do not redistribute the weights in the FCP repository, do not describe the resulting software as unrestricted for commercial use, and obtain a separate legal/licensing review before any use beyond non-commercial research/evaluation.

Even if every v3 gate passes, the maximum warranted statement is:

> The fixed SegFormer-B0 class-66 estimator met predeclared agreement and localization criteria on a balanced unused Oxford-102 automatic-segmentation proxy and an independently collected, manually boxed JRC grassland test set.

It would not prove pixel-perfect petal tissue, taxon-uniform error, absence of upstream training overlap, iNaturalist-wide domain transfer, calibrated camera colour, pollinator-perceived colour, or biological validity of any spatial/environmental association. Those limitations stay explicit in the atlas manuscript and bundle.
