# H1–H3 image-model reference audit

Checked 2026-09-15. Scope: original model identity, appropriate citations and
limits of transfer to flower-colour measurement. Only documentation and paper
records were accessed; no images, measured outputs, model weights or checkpoints
were downloaded, and no model was executed. Exact FCP settings and provenance
require the separate frozen-source and artifact audit.

## EfficientSAM: published model, not petal ground truth

Xiong, Y., Varadarajan, B., Wu, L., Xiang, X., Xiao, F., Zhu, C., Dai, X., Wang,
D., Sun, F., Iandola, F., Krishnamoorthi, R., & Chandra, V. (2024).
EfficientSAM: Leveraged masked image pretraining for efficient segment anything.
*Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern
Recognition (CVPR)*, 16111–16121.
[Official CVF publication record](https://openaccess.thecvf.com/content/CVPR2024/html/Xiong_EfficientSAM_Leveraged_Masked_Image_Pretraining_for_Efficient_Segment_Anything_CVPR_2024_paper.html).

The CVF record verifies all authors, title, year, venue and pages. Its abstract
describes lightweight encoders pretrained by reconstructing SAM image-encoder
features, followed by segment-anything fine-tuning on SA-1B. Cite this final
conference publication rather than silently treating the 2023 preprint as the
2024 proceedings article. No proceedings DOI was verified in this check; the
official CVF URL is sufficient and no DOI is invented.

The [author model repository](https://github.com/yformer/EfficientSAM) documents
point- and box-prompt examples, EfficientSAM-S and EfficientSAM-Ti checkpoints,
and separate encoder/decoder ONNX availability. This supports naming a
box-promptable EfficientSAM-Ti component. The current README does not prove the
identity of a historical ONNX export or FCP weight bytes; bind those to the
recorded revision and hashes, not to the moving repository homepage.

Access limit: CVF bibliographic record, abstract and indexed original-paper text
were available, but direct opening of the CVF PDF returned HTTP 403. The entire
paper and supplement were not audited. The repository README was inspected;
weights and example images were not opened.

## YOLO11: cite software, not a nonexistent verified paper

Jocher, G., & Qiu, J. (2024). *Ultralytics YOLO11* [Computer software].
[Ultralytics repository](https://github.com/ultralytics/ultralytics).

The [official YOLO11 documentation](https://docs.ultralytics.com/models/yolo11/)
provides this author/year/title software citation, with template version
`11.0.0`. It explicitly states that a formal YOLO11 research paper has not been
published and its citation DOI is pending. Do not invent a peer-reviewed
YOLO11 paper or DOI. The template version is not proof of the installed FCP
Ultralytics package version.

The same documentation distinguishes detection weights such as `yolo11n.pt`
from segmentation weights such as `yolo11n-seg.pt`. Naming YOLO11n does not alone
establish a segmentation model. Exact FCP flower fine-tuning, training labels,
checkpoint identity and downstream prompt construction must be described from
FCP evidence. Official documentation was accessible; no software was installed.

## Claim boundary for the combined measurement pipeline

Inference from the documented tasks: a detector box followed by a prompted mask
is a localization/segmentation procedure, not an independently validated focal
petal measurement. Generic benchmark performance cannot establish which floral
organs enter FCP masks, flower specificity, or transfer accuracy for the sampled
species. Neither model citation supplies exposure/white-balance calibration,
reflectance measurement, pigment identification or biological morph labels.
Those require target-domain evidence; a failed or incomplete FCP validation must
remain visible rather than being replaced by the original model's benchmark
claims. These citations do not justify altering frozen masks or thresholds.
