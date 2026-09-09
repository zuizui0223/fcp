# RGFCA target-domain measurement-validation feasibility

Source audit: **8 September 2026**. **Bounded feasibility only; no model
measurement, annotation campaign or new ecological inference is authorized.**
The initial audit covered five candidate dataset families through primary provider records,
author repositories and papers. No dataset photographs or masks were downloaded
or viewed, no model was run, and no FCP outcome, model or admission rule was
changed. This is not an exhaustive claim about every public dataset. Later
addenda retain their own read depth: Monarda is a sixth candidate, and the final
section fixes a separate six-document annotation-schema inspection. That step
may transfer embedded image strings but never decodes or displays them.

**Later checkpoint:** the user-supplied Monarda archive was subsequently scored
once under its separate authorization; [its completed region-agreement gate did
not pass](RGFCA_MONARDA_REGION_AGREEMENT_RESULTS.md). Historical no-pixel statements
below describe their respective intake stages, not the current Monarda status.
The [independent-reference follow-up](RGFCA_INDEPENDENT_REGION_REFERENCE_AUDIT_20260908.md)
and [next measurement-target design](RGFCA_MEASUREMENT_TARGET_AND_VALIDATION_V2.md)
retain the new USDA filename split conflict without running another benchmark.

## Decision

**No inspected resource establishes a ready-to-run, independent gold standard
for both focal-taxon petal pixels and target-domain iNaturalist photographs.**
Some resources can support narrower flower-region engineering checks. Passing
those checks would not validate FCP's taxon attribution, calibrated colour,
global transfer or an ecological mechanism.

The immediate recommendation is a **measurement-reference feasibility gate**,
not another search for a significant spatial or environmental result. Preserve
the completed reserve result and its failed flower-specific conjunction as
reported in [the independent replication results](RGFCA_RESERVE_REPLICATION_RESULTS.md).
The original 1,000-taxon / 100,000-photo frame is already outcome-opened; it is
not a fresh ecological validation reserve. A separate benchmark could assess
measurement error, but could not retrospectively make these outcomes unopened.

## 1. What each candidate actually measures

In this table, **verified** means explicitly documented by the inspected
provider or paper. **Limit** is this audit's inference from that annotation
target, not a new benchmark result.

| Candidate | Verified annotation target / access | Fit to the required reference |
| --- | --- | --- |
| Oxford Flowers 17 / 102 | 17: a ground-truth subset with trimap annotation. 102: released segmentations are algorithm-generated, not human masks. | 17 is a possible restricted foreground benchmark; 102 is not independent human pixel truth. Neither supplies the required focal-observation attribution contract. |
| USDA multi-species fruit flowers | Binary flower/non-flower images, with author-described freehand seeds followed by algorithmic refinement. Public-domain dataset record. | Potential orchard-region benchmark; not a fully manual petal reference and not a multi-taxon iNaturalist attribution test. |
| FlowerMaskDataset v1 | Provider advertises manual LabelMe polygons for 3,600 images across six flower labels, CC BY 4.0. Complete-dataset access asks for author contact. | Closest new manual-polygon lead, but ontology, full payload, split independence and focal-taxon truth remain unverified. |
| iNatLoc500 | Checked species-interest bounding boxes in validation/test; all 500 species are animals. | **Exclude:** wrong kingdom, and boxes are not petal masks. |
| Color-Cluster-Kit examples | Two iNaturalist plant species; colour-cluster summaries and classification examples. | Relevant acquisition domain, but released colour outputs are not an independent focal-petal mask reference. |
| Monarda fistulosa author archive | iNaturalist study with manual segmentation development; released segmentation archive is model output, and 500-photo human labels concern lighting/genus presence. | Target-domain lead, but the human polygon payload and its independent membership are not verified. |

### Oxford: distinguish genuine trimaps from generated masks

The Oxford 17 provider supplies class splits and an `imlist.mat` mapping for
the ground-truth subset, explicitly spanning 13 categories. Images originated
from web searches and photography. The original segmentation experiment used
753 images after excluding flower fields, very small/subsampled flowers and
unclear foregrounds. Human-labelled trimaps leave difficult boundaries and
ambiguous regions unlabelled. The paper also groups petals, tepals and sepals
under its modelling terminology. These are not exhaustive, anatomically
petal-only, focal-species instance masks.
[Provider README](https://www.robots.ox.ac.uk/~vgg/data/flowers/17/README.txt),
[Nilsback and Zisserman 2007, sections 2.3 and 3](https://www.robots.ox.ac.uk/~vgg/publications/2007/Nilsback07/nilsback07.pdf).

Oxford 102's README distributes `102segmentations.tgz` separately from class
labels and splits. A primary experimental analysis explicitly identifies its
segmentation references as automated and documents unreliable masks. Agreement
with these outputs is not agreement with independent human truth.
[Oxford 102 README](https://www.robots.ox.ac.uk/~vgg/data/flowers/102/README.txt),
[Melas-Kyriazi et al. 2022, appendix B.1](https://www.robots.ox.ac.uk/~vgg/publications/2022/Melaskyriazi22b/melaskyriazi22b.pdf).

**Admission status:** 17 may be considered for a narrowly scoped, fixed-model
trimap diagnostic only after actual mask inventory, unlabelled-pixel semantics,
rights and overlap checks. Its three classification splits are not proof of
observer/event-independent segmentation validation. The consulted Oxford
pages/READMEs do not establish per-image redistribution rights. The 753 figure
is the paper's experiment, not an audited count of today's archive.

### USDA: useful flower masks, but partly algorithmic reference labels

The provider distributes AppleA, AppleB, Peach and Pear images, flower/non-flower
labels and AppleA training/testing filename lists. Its record labels the
dataset U.S. Public Domain. It records orchard acquisition in 2016–2017.
[USDA Ag Data Commons dataset and file descriptions](https://agdatacommons.nal.usda.gov/articles/dataset/Data_from_Multi-species_fruit_flower_detection_using_a_refined_semantic_segmentation_network/24852636).

The paper distinguishes quicker superpixel-based training labels from evaluation
labels made with human freehand positive/negative traces and RGR propagation.
Thus evaluation masks are **human-guided, algorithm-refined**, not fully hand-drawn
petal boundaries. AppleA supplies training images; AppleB uses a background
panel. The paper reports 30 AppleA evaluation images and 18/24/18 for
AppleB/Peach/Pear. Provider descriptions retain earlier counts alongside later
file additions; archive membership must be reconciled rather than assumed.
[Dias et al. 2018, sections III.A, IV and IV.A; table I](https://arxiv.org/pdf/1809.10080).

**Admission status:** possible method-only stress test after a versioned census
and source-tree/event grouping. Do not treat multiple views or image crops as
independent plants. It has no documented per-instance focal-observation taxon
truth, and success in a few orchard crops would not establish global transfer.

### FlowerMaskDataset: a promising lead with an unresolved access/ontology gate

Version 1, published 23 June 2026, claims 600 images for each of six labels:
Butterfly Pea, Caesalpinia Pulcherrima, Rose, Plumeria, Tecoma Stans and Jatropha
Integerrima. It advertises LabelMe polygons and RGB/mask/foreground/background/JSON
files under CC BY 4.0, but directs users to email for the complete dataset.
The inspected record does not document collection-event groups, an independent
test allocation, per-instance scientific-name adjudication, petal-versus-other-organ
rules or annotator agreement. Some names are not species-level binomials.
[Versioned depositor record, DOI 10.17632/3pw57gdcj2.1](https://data.mendeley.com/datasets/3pw57gdcj2/1).

**Admission status: unresolved, not ready.** A public landing page and licence
label do not prove that the full annotation payload is anonymously retrievable
or that polygons represent focal petals. No contact was sent and no complete
archive was checked. This source cannot presently authorize FCP evaluation.

### iNatLoc500: verified focal-object supervision, but no plants

The author repository documents 12,500 validation and 12,500 test images with
checked boxes and exactly one instance of the species of interest. Crucially,
the same README states that **all classes belong to Animalia**. It also provides
source mappings to iNat17/iNat21. Its repository displays an MIT licence; this
does not independently establish rights to every source photograph.
[Author dataset repository, splits and label hierarchy](https://github.com/visipedia/inat_loc).

**Admission status: exclude.** Its careful focal-object procedure is a useful
design example, not a plant reference to repurpose by relabelling its boxes.

### Color-Cluster-Kit: target-domain examples, not independent mask truth

Perez-Udell et al. describe iNaturalist examples for *Geranium maculatum* and
*Linanthus parryae*. The method uses chosen colour ranges and k-means clusters;
supporting datasets contain the pipeline's summaries and classifications. The
inspected data-availability descriptions do not advertise manual focal-petal
pixel masks. Another algorithm's colour summary cannot serve as ground truth
for this purpose.
[Primary paper, methods and data availability](https://pmc.ncbi.nlm.nih.gov/articles/PMC9934523/).

The author code repository displays GPL-3.0. That software licence does not
relicense iNaturalist photographs or resolve per-photo reuse rights.
[Author code repository](https://github.com/atudell/Color-Cluster-Kit).

**Admission status:** contextual example only for this gate. Its observation
and photo IDs would need comparison with all prior FCP frames before any new
use; independent membership has not been established.

## 2. Independence is still an explicit unresolved check

None of these candidates was compared image-by-image with FCP's JRC
training/qualification data or its opened atlas frames. Provider provenance
differs, but provenance descriptions are **not a cryptographic disjointness
audit**, nor proof of absence from detector/segmenter pretraining. Reuse of an
already-qualified JRC test image would not create an independent target-domain
test. The existing [JRC source audit](RGFCA_TRAINING_SOURCE_AUDIT.md) remains
unchanged.

Before admitting any resource, the reference ledger would need frozen provider
version, exact asset/mask hashes, original source identifiers, all relevant
training/qualification exclusions and duplicate/event grouping. Unknown model
pretraining overlap must remain an explicit limitation, not silently become
`independent=true`. No such membership audit was performed here.

## 3. Proposed next gate — a design recommendation, not execution permission

1. **Specify the measurement target first.** Distinguish visible focal-taxon
   petals/tepals, other floral organs, non-focal flowers, non-floral background
   and indeterminate pixels. Whole-flower or trimap agreement cannot substitute
   for a petal-specific endpoint. If a whole-flower endpoint is preferred, name
   it explicitly in a new protocol without relabelling old results.
2. **Resolve one reference resource before accumulating more atlas photos.**
   Oxford 17 can support only a limited trimap diagnostic; FlowerMaskDataset
   needs complete-payload and annotation-ontology evidence first. Neither
   currently closes the iNaturalist focal-taxon gate. Human adjudication in
   an existing independent reference can be compatible with automated FCP
   measurement; commissioning new annotations requires a separately authorized
   workflow and must not be implied by this note.
3. **Freeze assessment independently of geography and colour outcomes.**
   Predeclare source/event/observer-disjoint allocation, all-frame denominators,
   flower-absent and multi-species scenes, ambiguous-reference handling, pixel
   purity/completeness, focal-taxon attribution error and reference uncertainty.
   Set tolerances from the intended measurement-error budget, not observed
   benchmark performance. A coordinate-blind reference panel and a fixed model
   must precede evaluation.
4. **Keep the two validation questions separate.** Pixel accuracy against a
   valid RGB mask can test extraction error, not reflectance, UV or pollinator
   perception. A later ecological replication must use a separately frozen
   fresh frame; measurement validation alone does not rescue the completed
   matched-background non-support. Model development, if needed, gets a new
   version and a genuinely held-out evaluation set, not retuning on this test.

**Current outcome:** `no_ready_joint_target_domain_reference_in_bounded_audit`.
This is a reference-data feasibility limitation, **not a biological negative**,
not proof that no suitable dataset exists, and not authorization to change the
completed study until a significant result appears.

## Read depth and retained uncertainties

- Oxford provider pages and both READMEs were read completely. Relevant original
  paper sections were read as PDF text; no figure photographs were rendered.
- The USDA provider record and the paper's acquisition/annotation sections and
  dataset table were read. Dataset archives and current filename lists remain
  uninspected; no final image count is claimed.
- FlowerMaskDataset's versioned record and the two author repository READMEs
  were inspected. No full annotation manifests or image rights inventories were
  obtained.
- Color-Cluster-Kit's methods/data-availability text was available through the
  primary-paper search extraction; a direct PMC page open encountered a browser
  challenge. Supplementary files were not downloaded or inspected.
- This note reports annotation semantics and access feasibility, not observed
  performance, actual train/test overlap, calibrated colour or botanical truth.

## Addendum — FlowerMaskDataset public metadata gate, 8 September 2026

This later, metadata-only check updates the **access uncertainty** above; the
original bounded audit is retained. **Public annotation-file entries are now
verified, but no usable/manual focal-petal mask count is established.** The
advertised complete 600-per-label dataset remains unverified.

### Public route and exact inventory observed

The provider's API documentation describes public dataset snapshots and file
listings. Its general API gateway returned an authentication error; no
credentials or access restrictions were bypassed. The dataset website's own
client explicitly uses an anonymous `/public-api` route. That route returned
HTTP 200 for the following metadata requests.
[Provider API documentation](https://data.mendeley.com/api/docs/),
[public website client used in this audit](https://data.mendeley.com/datasets/bundle.js?ec45581598eb1337d2e2).

- [Version-1 public snapshot](https://data.mendeley.com/public-api/datasets/3pw57gdcj2/snapshot/1):
  dataset ID `3pw57gdcj2`, version `1`, DOI `10.17632/3pw57gdcj2.1`,
  `is_confidential=false`, `is_metadata_only=false`. These flags do not certify
  full content or an annotation standard.
- [Version-1 folder listing](https://data.mendeley.com/public-api/datasets/3pw57gdcj2/folders/1):
  319 unique folders, with no dangling parent IDs: one root, six flower-label
  folders, six JSON folders, six masking folders, and 300 numbered masking
  subfolders (50 under each label). Root ID:
  `1b00e3b9-1cb2-43e2-81fe-e19139cfd3e4`.
- Each of the six JSON-folder responses listed **50 distinct completed JSON
  resources**. Together they contain **300 unique file IDs**, all advertised as
  `application/json`, with provider size and SHA-256 fields. Their advertised
  payload size totals **589,916,862 bytes**, ranging from 60,355 to 8,404,729
  bytes per resource. **None of those JSON payloads was downloaded.**

The exact folder responses are linked below. Filenames are identifiers only,
not biological replication units; the Jatropha folder spelling is preserved
verbatim rather than silently taxonomically corrected.

| Provider label | JSON folder ID / primary listing | Listed filenames | Entries |
| --- | --- | --- | ---: |
| Butterfly Pea | [49347441-47f8-4a23-8a74-2b922cedcd78](https://data.mendeley.com/public-api/datasets/3pw57gdcj2/files?folder_id=49347441-47f8-4a23-8a74-2b922cedcd78&version=1) | `1.json`–`50.json` | 50 |
| Caesalpinia Pulcherrima | [9167e445-8515-49b1-bf57-a459dbbf8b3e](https://data.mendeley.com/public-api/datasets/3pw57gdcj2/files?folder_id=9167e445-8515-49b1-bf57-a459dbbf8b3e&version=1) | `171.json`–`220.json` | 50 |
| Jatropha Integerimma | [cf445387-4408-4056-9332-7e6bc16bffe7](https://data.mendeley.com/public-api/datasets/3pw57gdcj2/files?folder_id=cf445387-4408-4056-9332-7e6bc16bffe7&version=1) | `401.json`–`450.json` | 50 |
| Plumeria | [ff41bedc-86e4-46b5-8776-ed7aac29226b](https://data.mendeley.com/public-api/datasets/3pw57gdcj2/files?folder_id=ff41bedc-86e4-46b5-8776-ed7aac29226b&version=1) | `501.json`–`550.json` | 50 |
| Rose | [420f7bab-f24c-404d-8f76-6e60b3e5bae4](https://data.mendeley.com/public-api/datasets/3pw57gdcj2/files?folder_id=420f7bab-f24c-404d-8f76-6e60b3e5bae4&version=1) | `501.json`–`550.json` | 50 |
| Tecoma Stans | [1a19a69d-b974-496a-8452-8474e7cbbd12](https://data.mendeley.com/public-api/datasets/3pw57gdcj2/files?folder_id=1a19a69d-b974-496a-8452-8474e7cbbd12&version=1) | `11.json`–`60.json` | 50 |

**Pagination ceiling:** adding `$start=50&$limit=100` to the Butterfly Pea
website-route query returned the same initial entries instead of a terminal
empty page. The public facade does not demonstrate the gateway documentation's
pagination semantics. Thus these are the **300 entries exposed by the six
observed folder responses**, not a certified exhaustive inventory of every
provider-held file. The matching 50 masking-folder names per label support a
public subset interpretation, but do not prove that additional data do not
exist. No absent tail was inferred and no hundreds-of-folders crawl was done.

### One separable annotation-name reference, not pixel validation

The metadata for Butterfly Pea `1_dataset`, folder
`e677f918-a695-4472-9a03-39e6827d007a`, lists `img.png`, `label.png`,
`label_viz.png` and `label_names.txt`. Only the text file was fetched. The
three image-content URLs remained unopened, so even this example's mask/image
identity and alignment are untested.
[Example folder's file metadata](https://data.mendeley.com/public-api/datasets/3pw57gdcj2/files?folder_id=e677f918-a695-4472-9a03-39e6827d007a&version=1).

The 22-byte `label_names.txt`, resource ID
`c77dc62d-5210-4c29-9269-64a473f32b53`, contains only `_background_` and `Flower`,
with CRLF line endings. Its downloaded SHA-256 exactly matches the provider:
`55fcabe129ddc3d8eadaeae09ee913e24dcf4e03bac50e62bbc757b8b8294bc9`.
[Exact public text resource](https://data.mendeley.com/public-files/datasets/3pw57gdcj2/files/c77dc62d-5210-4c29-9269-64a473f32b53/file_downloaded).

**Inference:** this example supports generic flower/background label names,
not a petal-only category, per-instance taxon identity or a correspondence to
an iNaturalist focal observation. No inference about mask quality or manual
labour is made from these two names. The other 299 corresponding mask/image
file inventories were not fetched; annotation JSON schemas and any embedded
image data remain unopened.

The dataset root and six label-container responses contained no direct files.
No standalone annotation protocol, source-event manifest or train/test split
document was identified in this bounded inspection. This is not a claim that
such documentation is absent from every uninspected location.

### What is still needed, and what this gate permits

The provider landing page explicitly directs readers to author contact for
the complete dataset. Contact is **not needed merely to discover these public
entries**. Full 3,600-record access, however, has not been demonstrated here.
The missing evidence is a versioned complete manifest with retrievable paired
annotation/image IDs, a documented anatomical/focal-taxon annotation ontology,
and source event/plant/observer and development/test membership. An image-free
annotation schema example must exclude embedded image data. The owner could
provide a public release or clarify existing public locations; no email or
other request was sent.

This gate supports preserving the public subset as a **candidate reference**.
It does not admit it as an evaluation set, authorize downloads/model execution,
prove 300 usable independent photographs, or repair FCP's completed
flower-specific non-support. Outcome status:
`public_annotation_entries_verified_joint_reference_not_admitted`.

### Local receipts and exact raw-response hashes

Raw public metadata and the one text file are preserved under
`.artifacts/reference-feasibility-20260908/`. These are response-byte hashes,
not recomputed image or mask hashes. File responses include generated download
expiry fields, so a later response can differ bytewise while naming identical
resources. JSON-content SHA fields remain **provider-advertised and unverified**.

In the table, JSON-folder receipts have the filename
`flowermask-v1-json-files-<folder-id>.json`; other names are given in full.

| Receipt | Bytes | SHA-256 |
| --- | ---: | --- |
| `flowermask-v1-public-folders.json` | 50,700 | `fd30ad1ba248c43eca04c1a53cec31e87617dbd815b5a98b752eb461dec2b892` |
| `flowermask-v1-snapshot.json` | 6,704 | `568bf89c541326877c71f0f5f67041ec8b2d71f9495bdec1b0f93c2b24e41e0c` |
| `49347441-47f8-4a23-8a74-2b922cedcd78` | 37,764 | `63017d75b61de0e7ce350c3d8aebcc2a1f3df2247824e8f51fb8b11c09547a6b` |
| `9167e445-8515-49b1-bf57-a459dbbf8b3e` | 37,793 | `5ddfc501ba45249b35f2b63cd4b01b104349cc3ed6019672e62c28a1c54da7a7` |
| `cf445387-4408-4056-9332-7e6bc16bffe7` | 37,758 | `b503c8ed04f41dbabec5954820b2c62f6be1ac79a4c4f1db172faf1626fa0c32` |
| `ff41bedc-86e4-46b5-8776-ed7aac29226b` | 37,803 | `bc21e1eefb7431021991f67e29cdc5d5fc74e2a180c0ea1ad0ce81b5dd2051b2` |
| `420f7bab-f24c-404d-8f76-6e60b3e5bae4` | 37,729 | `efa3e5d212c4ecbbe1bc5232e508d55160729d06515fc532cfaebb98aab4198f` |
| `1a19a69d-b974-496a-8452-8474e7cbbd12` | 37,654 | `38431ea37dddcad5d5fb86a5be159290d4c67bb57f65dc561f90144443be6f77` |
| `flowermask-v1-root-files.json` | 2 | `4f53cda18c2baa0c0354bb5f9a3ecbe5ed12ab4d8e11ba873c2f11161202b945` |
| `flowermask-v1-butterfly-1-files.json` | 3,004 | `a227790321584118998d390850cca1b66ab4f7ac54276ebc0c0262f566c77485` |
| `flowermask-v1-butterfly-1-label_names.txt` | 22 | `55fcabe129ddc3d8eadaeae09ee913e24dcf4e03bac50e62bbc757b8b8294bc9` |

## Addendum — Oxford access and Monarda annotation provenance

Public HTTP HEAD requests on 8 September returned 200 for Oxford 17's
`datasplits.mat` (7,487 bytes) and `trimaps.tgz` (3,915,787 bytes). No archive
body was obtained: reachability is not a current segmentation census. The
provider README was read through the normal public text endpoint. The original
2007 paper's relevant sections were also read at the
[author-hosted PDF](https://www.robots.ox.ac.uk/~men/papers/nilsback_bmvc07.pdf),
confirming the exclusions and trimap limits above. This does not resolve rights,
focal-taxon attribution or independence.

The [Monarda author repository](https://github.com/pmckenz1/monarda_fistulosa_color)
was inspected at commit `49812a431abce78c871f23aff6bc7978573dd502`.
Its recursive Git listing was not truncated (182 entries); three release tags
had no attached assets. The README describes manual flower-pixel annotation
in Roboflow, but identifies `intermediate_data_files/segmentation_results.zip`
(69,126,884 bytes, Git blob `fa133c7d048312c9426f4f696f8ac30be1044bc8`)
as the trained model's output. That archive was not fetched and cannot be
relabeled human ground truth.

The [validation README](https://github.com/pmckenz1/monarda_fistulosa_color/blob/49812a431abce78c871f23aff6bc7978573dd502/validation/README.md)
documents lighting categories and Monarda-presence labels for 500 photographs,
not anatomical pixel masks. Repository MIT licensing does not establish every
source photograph's licence. The linked
[Roboflow segmentation model](https://universe.roboflow.com/patricks-dashboard/monarda_fistulosa_segmentation/model/1)
returned HTTP 403 to a normal page request; no access restriction was bypassed.
The original human polygons, allocation and source-ID overlap therefore remain
unverified. No CSV outcomes, notebook image outputs, model or photograph was
opened. These archive/preprint method descriptions are not a claim to have read
the inaccessible final published Methods. Monarda's known positive ecological
result is not new independent evidence for FCP.

## Fixed next gate — six annotation-document schemas, not a benchmark

Under the ongoing user-directed implementation, the next bounded action is to
inspect the **structure** of six public FlowerMask annotation documents. This
supersedes the preceding metadata-only acquisition boundary only for these six
documents. It does not admit any benchmark or authorize model evaluation.

1. Freeze the 300 entries from the six observed metadata responses in
   [the selection plan](supporting/rgfca_reference_schema_plan_v1.json).
   Preserve provider names, sizes, SHA-256 and URLs. This is not an exhaustive
   public-file census: the pagination caveat remains. Within each label choose
   the lexicographically smallest provider file ID, without looking at size,
   annotation content, images, colour, geography or FCP performance. Selection
   is deterministic, not a representative random biological sample.
2. The exact plan SHA-256 is
   `48436a74160b4700e2d9b2bbc803aeac7fcfe8a1f09885140a1eade034c5aade`.
   The six advertised sizes total 9,442,182 bytes. The public record's CC BY 4.0
   notice supports inspecting its annotation documents; no photograph
   redistribution clearance or independent test membership is inferred.
3. Commit the script, plan, offline artificial-JSON tests and workflow together.
   The push preflight performs no document acquisition. After it succeeds, a
   manual dispatch on the same commit performs the six document requests once,
   with exact size/hash checks, 45-second request timeouts and a 10 MB/file cap.
   Retain every selected ID on error; do not retry or replace it automatically.
4. JSON documents may contain encoded `imageData`. The process necessarily
   transfers and hashes those opaque bytes, but does **not** decode or display
   images, follow image paths/URLs, load a model, join coordinates or extract
   colour. Retain only structural summaries and terminal outcomes: label/type
   counts, point counts (not vertices), dimensions, field names and hashed image
   path. Raw JSON/image payloads are not written or published as artifacts.
5. A successful structure check establishes only that these six exact documents
   are retrievable and have the reported schema. It cannot establish polygon
   accuracy, manual authorship, petal-only boundaries, focal-taxon identity,
   image/mask alignment, independent events, 300 usable images, target-domain
   transfer or flower-specific ecological replication. These all remain separate
   requirements. Failed schema checks are technical/source outcomes, not ecology.

Implementation: `scripts/analysis/audit_rgfca_reference_schema.py`,
`tests/test_rgfca_reference_schema.py`, and
`.github/workflows/rgfca-reference-schema.yml`. The existing atlas measurements,
all four reserve results and the discovery-background STOP stay unchanged.

### Completed six-document gate — 8 September 2026

The selection and implementation were committed **before acquisition** at
`df0838c15f01eb98eaee6b620c3bb4fd5c9bc78e`.
[Push preflight 34182661045](https://github.com/zuizui0223/fcp/actions/runs/34182661045)
passed 16 artificial-JSON tests and the exact plan hash, with the acquisition
job skipped. All five downloaded preflight artifact files byte-match that Git
commit. Subsequently,
[manual execution 34182715550](https://github.com/zuizui0223/fcp/actions/runs/34182715550)
passed its preflight and completed all six document inspections once.

All six provider SHA-256 and size pairs matched, totalling **9,442,182 bytes**.
The result artifact `10039458598` has provider archive digest
`e0e4587d3c05380769e8af48adff8625e82c6494f8de90c3232da49baa41572e`.
The downloaded and committed [terminal JSON](supporting/rgfca_reference_schema_result_v1.json)
has locally verified SHA-256
`70b2fc123dd2ea5e57dff67de84ce9e77cea52a87124243f687f636e9ef8c23c`.
No raw annotation document or embedded image was published in the artifact.

| Provider group | Recorded polygon label | Polygon count | Populated group IDs |
| --- | --- | ---: | ---: |
| Butterfly Pea | `Flower` | 1 | 0 |
| Caesalpinia Pulcherrima | `flower` | 1 | 0 |
| Jatropha Integerimma | `flower` | 1 | 0 |
| Plumeria | `flower` | 1 | 0 |
| Rose | `rose` | 1 | 0 |
| Tecoma Stans | `Flower` | 1 | 0 |

All six contain non-empty embedded image strings. These bytes were transferred
and hashed as part of the JSON but **not decoded or displayed**. No model,
palette, coordinate join or biological inference ran. Field names `description`
and `mask` are present in shapes, but their values were not reported; neither
their presence nor an empty group ID proves absence of all taxon information.
Observed label case is preserved, not silently normalized into a new ontology.

**Decision:** actual annotation-document access and format are verified for the
six fixed examples. They remain generic flower-region candidates; this does not
admit them as verified focal-petal ground truth. The other 294 registered JSON
bodies and corresponding image/mask pairs remain uninspected. Any next census
must retain these six as already schema-opened examples, verify source and
image correspondence and overlap, and explicitly state a narrower diagnostic
target if it uses whole-flower polygons. No new benchmark score or ecological
claim is available. Later-head CI checks the saved ledger without reacquisition;
the completed six-document acquisition job is closed at later commits.

### Subsequent Monarda archive receipt — 8 September 2026

The user supplied the complete v1 COCO segmentation ZIP, resolving the earlier
access restriction. All 115 ZIP members passed CRC reads; the census contains
110 image/annotation pairs and 788 polygon components (77/22/11 images in the
provider splits). All 110 original photo and observation IDs are mapped using
the author's immutable filename/index convention and metadata. Both full
50,000-photo FCP frames have zero photo-ID and observation-ID intersections.
These are documentary linkage and bounded ID-disjointness results, not model
or perceptual independence, global coverage or segmentation accuracy.

The [dedicated intake audit](RGFCA_MONARDA_REFERENCE_INTAKE.md) records all member
hashes, the unannotated image retained as unknown, export dates that must not
be treated as phenology, historical photo licences (98 noncommercial), source
selection uncertainty and the next prospective flower-region diagnostic gate.
No image pixels were decoded, annotations rasterized, model run, colour
measured or geography joined. This does not reopen the completed FlowerMask
acquisition or alter the reserve's flower-specific non-support.
