# RGFCA measurement and display provider audit

Audit date: **7 September 2026**. Scope: the four providers underlying the
fixed image-measurement and map-display stack, not an exhaustive software
bibliography, a new measurement validation, or legal clearance for a submission.
Research was restricted to primary provider documentation, source, licence
texts, and the original EfficientSAM publication. No reserve images, reserve
measurements, or partial colour outputs were inspected; no models were run.

## Outcome and evidence boundary

The implementation can be described precisely, with separate citations for
photo/data provenance, detection, prompted segmentation, and cartographic
display. Provider licences and generic model benchmarks do **not** establish
flower-ROI validity, calibrated colour, taxon-uniform error, pollinator vision,
or ecological causation. Broad mixed-licence data reuse remains distinct from
the separate, CC0-only discovery photo-bar route.

Local evidence was read from the measurement materializer and its immutable
source commit `9fae6ccdf684a46026f72ba12e98de2c5c54bf2a`, specifically
`docs/supporting/jbi_atlas_roi_estimator_contract_v4.json`,
`fcp_pipeline/flower_roi_v4_runtime.py`, and `requirements-atlas-roi-v4.txt`.
The materializer fetches this frozen runtime and custom detector rather than
silently taking current upstream models. See the
[materializer](../scripts/analysis/materialize_random_photo_first_h9_measurement_infrastructure.sh)
and [immutable source tree](https://github.com/zuizui0223/fcp/tree/9fae6ccdf684a46026f72ba12e98de2c5c54bf2a).

### Exact implementation identities

| Component | Fixed identity and role | What this audit verified |
| --- | --- | --- |
| Detector software | `ultralytics==8.4.112`; YOLO11n **detection**, not YOLO11 instance segmentation | Frozen requirements/runtime and version-specific upstream licence/citation files |
| Upstream detector base | `yolo11n.pt`, Ultralytics assets release `v8.4.0`; 5,613,764 bytes | Expected identity recorded by frozen contract; not newly downloaded or rehashed here |
| Actual detector | `data/atlas/qualification/roi_v4_training/jrc_yolo11n_last_v4.pt`; locally fine-tuned flower detector | Materializer and frozen runtime select the custom checkpoint, not the base checkpoint |
| Segmenter | EfficientSAM-Ti box-prompt inference, repository revision `d525f622e6f640acf5a0fc37c7ca1f243da5bde0`; separate ViT-T encoder/decoder ONNX files | Official pinned tree, README, licence, file sizes and Git blob metadata; not newly downloaded or rehashed here |
| Segmenter runtime | `onnxruntime==1.23.2`; frozen runtime validates expected SHA-256 values before loading | Source inspection, not an inference replay |
| Basemap | GeoPandas `0.14.4` bundled `naturalearth_lowres`, used by the publication figure script | All five official tagged files were independently fetched as bytes and their hashes matched the local figure manifest |

Expected model SHA-256 values are **contract identities**, not new binary
verification results from this audit:

| Artifact | Expected SHA-256 |
| --- | --- |
| Upstream `yolo11n.pt` | `0ebbc80d4a7680d14987a577cd21342b65ecfd94632bd9a8da63ae6417644ee1` |
| Custom `jrc_yolo11n_last_v4.pt` | `f1aaeec4664fe2c178e5cf2bc1f508977bef3e4aa7b40613026cb8ae3de789d5` |
| `efficient_sam_vitt_encoder.onnx` | `84ed466ffcc5c1f8d08409bc34a23bb364ab2c15e402cb12d4335a42be0e0951` |
| `efficient_sam_vitt_decoder.onnx` | `a62f8fa5ea080447c0689418d69e58f1e83e0b7adf9c142e2bd9bcc8045c0b11` |

## 1. iNaturalist: photo rights are not observation-data rights

The platform distinguishes observation data, photographs, and sounds, each
with its own licence. Contributors retain ownership; a default licence is not
proof of the licence on a particular selected object, and licences can be
changed. Thus an observation's licence cannot authorize reuse of its attached
photograph. The official help article was last modified 12 September 2024;
the complete relevant article was read on the audit date.
[iNaturalist, What are licenses?](https://help.inaturalist.org/en/support/solutions/articles/151000175695)

The official Open Data repository explains that photo-specific licences govern
reuse and that photographs remain copyrighted unless dedicated to the public
domain. Its attribution examples are useful operational guidance, not evidence
that every observer is the photographer. Preserve the provided photo-level
credit, source photo/observation URLs and licence declaration; label observer
identity separately when authorship is not established.
[iNaturalist Open Data, README and attribution guidance](https://github.com/inaturalist/inaturalist-open-data)

The current candidate parser retains `photo_license`, `attribution`, observer
ID/login, observation/photo IDs and the large-image URL. Its frozen acquisition
contract permits `cc0`, `cc-by`, `cc-by-sa`, `cc-by-nc`, and `cc-by-nc-sa`.
However, the inspected schema does not retain an observation-data licence or an
explicit versioned photo-licence URL. This is a **schema-level gap**, not a claim
that particular rows lack attribution: data rows and current per-photo rights
were not audited here. See the
[acquisition script](../scripts/acquisition/run_global_monte_carlo_candidate_acquisition_shard.py),
[parser](../fcp_pipeline/random_photo_h9_pool.py), and
[frozen allowed-licence contract](supporting/global_monte_carlo_candidate_acquisition_contract_v1.json).

Official CC 4.0 deeds illustrate why these categories cannot be pooled into an
unqualified single reuse permission: BY requires attribution/licence information
and marking changes; NC adds a noncommercial restriction; SA addresses licensing
of shared adaptations. The inspected abbreviated codes do **not** establish that
every underlying photograph uses version 4.0. Check the actual object-specific
declaration before release rather than assigning a version from this paragraph.
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/),
[CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/),
[CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/),
[CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/)

For citation, iNaturalist gives separate guidance for downloaded datasets,
individual records, and the platform, and prefers a corresponding GBIF download
DOI where applicable. A direct-API dataset must not receive an invented GBIF DOI.
The manuscript should report its actual source, acquisition/export dates and
frozen manifest, while keeping photo credits separate from the platform citation.
The complete guidance, modified 4 December 2023, was read.
[iNaturalist, How should I cite iNaturalist?](https://help.inaturalist.org/en/support/solutions/articles/151000170344-how-should-i-cite-inaturalist-)

### Separate publication illustration route

[RGFCA_PHOTO_BAR.md](RGFCA_PHOTO_BAR.md) specifies a fixed 24-slot, discovery-only
CC0 illustration with current photo-level declarations checked before pixels,
exact original-image hashes and retained palette/count checks, source links and
credits, and retention of derived RGBA crops rather than source-image files.
Its execution/QA receipt is a separate deliverable; this provider audit does not
reopen images or independently assert that receipt's completion. Generic API
credit text such as “no rights reserved” must not be silently replaced with an
observer-as-photographer assertion. The illustration is not representative global
coverage or independent validation, and its CC0 filter does not clear the broader
mixed-licence candidate/discovery data or model stack.

CC0 allows copyright reuse and modification but does not waive every possible
privacy, publicity, patent, trademark or third-party right, provide a warranty,
or imply endorsement. The deed is a summary rather than the legal instrument;
it was read in full on the audit date.
[Creative Commons, CC0 1.0 Universal](https://creativecommons.org/publicdomain/zero/1.0/)

## 2. Ultralytics: cite YOLO11, distinguish base and custom weights

The official YOLO11 documentation distinguishes detection checkpoints such as
`yolo11n.pt` from `-seg` checkpoints. Its citation section identifies the software
as **Glenn Jocher and Jing Qiu (2024), Ultralytics YOLO11**, and states that there
is no formal YOLO11 research paper. The family citation is not the installed
package version. Report both YOLO11n and the actual package `8.4.112`; a generic
COCO model benchmark does not validate a fine-tuned flower detector.
[Ultralytics, YOLO11 documentation, models and citations sections](https://docs.ultralytics.com/models/yolo11)

The exact `v8.4.112` source tag carries an AGPL-3.0 licence. The provider's current
licensing page also states that its default AGPL treatment includes trained and
fine-tuned models, with enterprise licensing as a separate option. This is the
provider's stated licensing position, not a legal determination about every
possible downstream integration. Local fine-tuning alone does not justify
labelling the custom checkpoint unrestricted.
[Version-specific source licence](https://raw.githubusercontent.com/ultralytics/ultralytics/v8.4.112/LICENSE),
[Ultralytics licensing, trained-models guidance](https://www.ultralytics.com/license)

**Citation trap:** the inspected package-tag `CITATION.cff` has a preferred
**YOLO26** paper citation, despite this pipeline using YOLO11n. Copying that
preferred citation would misidentify the model. Use the model-specific official
YOLO11 software citation plus package version and checkpoint provenance; do not
invent a YOLO11 DOI.
[Exact `v8.4.112` citation metadata](https://raw.githubusercontent.com/ultralytics/ultralytics/v8.4.112/CITATION.cff)

The frozen FCP contract names the
[upstream base asset](https://github.com/ultralytics/assets/releases/download/v8.4.0/yolo11n.pt)
separately from the custom JRC-trained, one-flower-class last-epoch checkpoint.
Inference loads the latter and passes detected boxes to EfficientSAM. The
contract records 400 JRC development images and 50 training epochs; those local
provenance facts are not upstream YOLO11 performance claims. JRC training-data
rights and a complete AGPL distribution/integration compliance review were not
part of this four-provider audit. Neither is established by citing the detector
software or its base weights.

Access ceiling: official model/citation and licensing sections, exact-tag
licence identity and complete `CITATION.cff`, plus frozen local source. No new
weights download, pretraining-data audit, enterprise entitlement check or legal
compliance certification was performed.

## 3. EfficientSAM: original method, pinned export, separate validity limit

The original publication is **Yunyang Xiong et al. (2024), EfficientSAM:
Leveraged Masked Image Pretraining for Efficient Segment Anything**, CVPR,
pp. 16111–16121. It introduces an efficient segmentation approach with lightweight
image encoders and masked-image pretraining. The official abstract's general
vision results do not establish flower/petal segmentation accuracy in RGFCA or
validate the combined detector-to-mask estimator.
[Official CVPR publication and bibliography](https://openaccess.thecvf.com/content/CVPR2024/html/Xiong_EfficientSAM_Leveraged_Masked_Image_Pretraining_for_Efficient_Segment_Anything_CVPR_2024_paper.html)

The pinned official README documents the separate encoder/decoder ONNX route
and credits the export contribution; the code and checkpoint files are in the
same pinned repository. RGFCA uses the ViT-T files, not the ViT-S checkpoint or
an unpinned current export. The original paper should be cited alongside the
implementation revision/export identity, rather than treating an ONNX filename
as a new scientific method.
[Pinned README](https://raw.githubusercontent.com/yformer/EfficientSAM/d525f622e6f640acf5a0fc37c7ca1f243da5bde0/README.md),
[Pinned weights inventory](https://api.github.com/repos/yformer/EfficientSAM/contents/weights?ref=d525f622e6f640acf5a0fc37c7ca1f243da5bde0)

The root of this pinned tree supplies Apache-2.0. Its redistribution conditions
include carrying the licence and applicable notices and identifying modified
files; its grant is not a trademark endorsement or warranty. No distinct
weights-only licence was found in the inspected root and weights listing, so
the supported provenance is the provider's repository-level Apache-2.0
declaration, not an independently reconstructed chain of rights for all
pretraining data.
[Pinned Apache-2.0 licence](https://raw.githubusercontent.com/yformer/EfficientSAM/d525f622e6f640acf5a0fc37c7ca1f243da5bde0/LICENSE)

The official GitHub API confirmed encoder size **24,799,761 bytes**, Git blob
`6458f72477ae216a1bd68db41ffa14802c8d54f1`, and decoder size **16,565,728 bytes**,
Git blob `f9310202c916fe5a4ec9a6897edae855caf023f4` at the pinned revision. These
are provider file-metadata checks, not fresh SHA-256 checks of the binaries.
[Encoder](https://github.com/yformer/EfficientSAM/blob/d525f622e6f640acf5a0fc37c7ca1f243da5bde0/weights/efficient_sam_vitt_encoder.onnx),
[Decoder](https://github.com/yformer/EfficientSAM/blob/d525f622e6f640acf5a0fc37c7ca1f243da5bde0/weights/efficient_sam_vitt_decoder.onnx)

The local frozen runtime uses detected boxes as prompts, selects masks by
predicted IoU, thresholds logits, clips masks to boxes and unions them. This is
automated flower-candidate masking, not a botanical guarantee that every retained
pixel is petal tissue. Its measurement-validity ceiling must follow the separate
RGFCA qualification evidence, not the EfficientSAM licence or generic benchmarks.
Species-conditioning is by the observation's taxon label: the frozen image-only
estimator has no focal-taxon input and pools every retained flower instance.
The [qualification audit](RGFCA_ROI_QUALIFICATION_AUDIT.md) records the exact source
and artificial-mask pooling tests. A same-scene flower from another taxon can
contribute; its frequency and effect are unmeasured, and background controls do
not independently validate focal-species attribution.

Access ceiling: pinned README, full licence and provider file metadata; official
CVF-indexed abstract and bibliographic record. A later direct CVF page open was
denied, so this audit does not claim a complete reading of the paper's methods,
supplement or benchmarks. No unverified DOI is supplied.

## 4. Natural Earth: reproducible, display-only cartographic context

Natural Earth declares its raster/vector map data public domain and permits
modification and publication, including commercial use, without requiring
permission or credit. An acknowledgement remains useful scientific provenance.
Its terms also disclaim accuracy/fitness guarantees. These terms do not make
cartographic borders an ecological classification or validate any shared floral
transition.
[Natural Earth, Terms of use](https://www.naturalearthdata.com/about/terms-of-use/)

The exact GeoPandas `v0.14.4` loader identifies the bundled
`naturalearth_lowres` dataset and points to Natural Earth's 1:110-million cultural
vectors. That is a cartographic scale, not a 100/250/500-km inferential resolution.
The original upstream Natural Earth release number was not reconstructed here;
the reproducible identity is the GeoPandas tag plus the actual bundled files.
[Exact tagged loader](https://raw.githubusercontent.com/geopandas/geopandas/v0.14.4/geopandas/datasets/__init__.py),
[Exact tagged dataset folder](https://github.com/geopandas/geopandas/tree/v0.14.4/geopandas/datasets/naturalearth_lowres)

The [publication figure script](../scripts/analysis/make_rgfca_publication_figures.py)
requires GeoPandas `0.14.4`, hashes the five files and draws low-resolution
land/country polygons as gray reference context. These polygons are not
ecoregions, pollinator biogeographic units, inferred floral boundaries, or inputs
to the ecological test. Country polygon edges must not be described as only
coastlines. The [figure workflow](../.github/workflows/rgfca-publication-figures.yml)
pins the display environment.

All five files were fetched from the official GeoPandas tag as bytes in memory
and independently SHA-256 checked on the audit date. Each matched
[the existing figure manifest](supporting/rgfca_publication_figure_manifest_v1.json):

| `naturalearth_lowres` file | Bytes | Independently matched SHA-256 |
| --- | ---: | --- |
| `.cpg` | 10 | `09fc313075748ce8ead962229ed89c919d5a9ff71974ee725a0b48bb87d975a2` |
| `.dbf` | 50,285 | `5cfbcaa21ce5fad798abf2ec65ab0db59538f9bb8a37273ef62b6aa8121487fd` |
| `.prj` | 145 | `a02a27b1d1982c8516d83398e85a3c8b1aef1713c13ef4d84d7bde17430c07c4` |
| `.shp` | 180,744 | `1f689e60b357e1e98702d5d9f774e95e77fc6b324487cadf57eb9317d533ce12` |
| `.shx` | 1,516 | `7933917ebd636eb80822c1bddfa31dd7597ca9e6f89b21e23283c594ed6482c5` |

Access ceiling: full provider terms, exact-tag loader and folder/API metadata,
and byte-level verification of the five tagged cartographic files. No inference
or display regeneration was performed.

## Remaining release questions and claim limits

1. **Mixed-licence records:** a broader reusable data/image package still needs
   object-specific rights and credit provenance, including an explicit decision
   about observation-data licensing and versioned photo-licence references.
   Do not relicense the mixed pool as blanket CC BY or interpret the CC0
   illustration as resolving that pool. Whether a particular aggregate derived
   output carries additional legal obligations was not determined here.
2. **Models and training provenance:** document the distinct base/custom
   detector and pinned segmenter, retain relevant notices, and assess the actual
   proposed distribution/integration against provider terms. This audit neither
   certifies AGPL compliance nor audits upstream pretraining corpus rights.
   The subsequent [JRC training-source audit](RGFCA_TRAINING_SOURCE_AUDIT.md)
   checks the separate dataset-specific CC BY 4.0 notice and citation/split
   provenance; it does not settle the model-distribution question. The
   [qualification audit](RGFCA_ROI_QUALIFICATION_AUDIT.md) separately reconstructs
   saved FCP localization gates and documents their measurement limits.
3. **Citation specificity:** cite actual YOLO11 software, EfficientSAM's original
   CVPR paper plus pinned code/export, the actual iNaturalist acquisition, and
   Natural Earth plus the exact bundled display dataset. Avoid automatic YOLO26
   substitution, an invented GBIF download DOI, or an unverified model DOI.
4. **Scientific interpretation:** attribution, an open licence, byte identity,
   and generic provider accuracy are not flower measurement validation. Keep
   camera-space photo measurements, model qualification, the illustrative
   CC0 photo bar, and ecological inference as distinct evidence layers.

This audit deliberately leaves these questions explicit rather than declaring
the entire manuscript, image collection or software/model bundle legally
cleared. All web sources above were consulted on 7 September 2026; mutable
provider pages describe their observed state on that date.
