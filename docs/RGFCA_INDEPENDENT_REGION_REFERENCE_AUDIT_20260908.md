# Independent flower-region references after the Monarda result

Audit date and external-source access date: **8 September 2026**.
Status: **bounded reference feasibility; no new benchmark admitted or run**.

## Decision first

None of the three resource families inspected here establishes an immediately
usable, independently allocated reference for **focal-taxon petal pixels in
iNaturalist photographs**. This is not a claim that such a reference cannot
exist. The distinction between flower foreground and anatomically identified
petals remains decisive. A smaller, well-defined measurement reference can be
useful, but simply replacing Monarda with another collection does not resolve
that distinction.

| Resource family | Evidence-supported use | Current admission decision |
|---|---|---|
| FlowerMaskDataset v1 | Candidate manual generic-flower polygons; six previously fixed documents have verified schema access | **Insufficient evidence:** complete payload, anatomy, source grouping and independent allocation unresolved |
| Oxford Flowers | 102: algorithm-output comparison, not human pixel truth. 17: a possible restricted human-trimap diagnostic | **102 target/reference mismatch; 17 conditional only:** rights, current inventory and source independence unresolved |
| USDA multi-species fruit flowers | Candidate orchard flower/non-flower diagnostic against human-guided, algorithm-refined labels | **Partly accessible but not admitted:** current public split lists overlap, and paired archives are unverified |

The underlying evidence and its limitations are recorded below. These decisions
are this audit's methodological interpretation, not provider certification.

The completed [Monarda result](RGFCA_MONARDA_REGION_AGREEMENT_RESULTS.md) remains
closed: its 110 images are outcome-opened and its operational agreement gate did
not pass. They are not a new holdout for tuning and rescoring. The
[JRC source audit](RGFCA_TRAINING_SOURCE_AUDIT.md) concerns a different,
box-annotated development source. Neither result is altered by this note.

## 1. FlowerMaskDataset: retain the manual-polygon lead, not a petal claim

The depositor's version-1 record advertises 3,600 manually annotated images,
600 for each of six labels, with LabelMe polygons and RGB/mask/foreground/
background/JSON products. It names Butterfly Pea, Caesalpinia Pulcherrima, Rose,
Plumeria, Tecoma Stans and Jatropha Integerrima. These are provider labels, not
six independently adjudicated species-level identifications. The stated task is
foreground/background separation; the record does not give organ-by-organ
inclusion rules. It displays CC BY 4.0 and requests author contact for the complete
dataset. No contact was made. [Depositor record, DOI 10.17632/3pw57gdcj2.1](https://data.mendeley.com/datasets/3pw57gdcj2/1)
(accessed 2026-09-08).

The earlier FCP metadata audit exposed 300 JSON file entries through six folder
responses, with unresolved pagination semantics. It did **not** establish a
3,600-image public census. The corresponding public routes remain the documented
source of that retained observation; they were not re-crawled in this follow-up.
[Provider folder endpoint](https://data.mendeley.com/public-api/datasets/3pw57gdcj2/folders/1),
[example JSON-folder endpoint](https://data.mendeley.com/public-api/datasets/3pw57gdcj2/files?folder_id=49347441-47f8-4a23-8a74-2b922cedcd78&version=1)
(retained metadata accessed 2026-09-08; exact response hashes in
[the existing feasibility audit](RGFCA_TARGET_DOMAIN_VALIDATION_FEASIBILITY.md)).

The already completed, prospectively selected six-document inspection at
`df0838c15f01eb98eaee6b620c3bb4fd5c9bc78e` is reused, not repeated. Its
[terminal ledger](supporting/rgfca_reference_schema_result_v1.json) verifies
9,442,182 bytes and one polygon per selected document. Labels are `Flower`,
`flower` or `rose`; populated group IDs are absent in these six examples.
Embedded image strings were transferred in that earlier JSON inspection but
never decoded. These observations establish neither exhaustive flower coverage
nor absent taxon information in every other field. [Completed source inspection](https://github.com/zuizui0223/fcp/actions/runs/34182715550)
(retained evidence read 2026-09-08).

**Unresolved gate:** an anatomical annotation protocol; versioned image/annotation
correspondence; per-original-image ownership/attribution evidence; collection
event, plant and observer groups; and development/test membership. The dataset
licence label is recorded as a provider declaration, not silently substituted
for a per-photo provenance ledger. No record inspected establishes absence from
FCP, Monarda, JRC or foundation-model training. The six schema-opened examples
must remain identified as such in any future allocation. This is a candidate
for a specifically named flower-region task, not verified focal-petal truth.

## 2. Oxford family: do not substitute 102 generated segmentations for 17 trimaps

### Oxford 102: reject as independent human pixel truth

The provider supplies separate images, segmentations, class labels and
`setid.mat`. Its README describes web-search/author photographs and records the
merge of two Petunia categories. The author's thesis reports 8,189 images in
102 classes, with 1,020 training, 1,020 validation and 6,149 test images. These
are classification partitions, not verified observer/event partitions.
[Provider page](https://www.robots.ox.ac.uk/~vgg/data/flowers/102/index.html),
[complete README](https://www.robots.ox.ac.uk/~vgg/data/flowers/102/README.txt),
[author thesis, section 3.2](https://www.robots.ox.ac.uk/~men/thesis/menilsback_thesis.pdf)
(accessed 2026-09-08; thesis evidence read through indexed primary text).

The original classification paper explicitly describes automatic segmentation
before feature extraction. Agreement with those supplied segmentations is
agreement with another algorithm, not independent manual annotation accuracy.
[Nilsback and Zisserman, original classification paper, introduction and section 2](https://robots.ox.ac.uk/~men/papers/nilsback_icvgip08.pdf)
(accessed 2026-09-08 through indexed primary text). No current archive was
downloaded, counted or assessed. Public download links alone do not establish
the anatomical quality or completeness of their contents.

### Oxford 17: a narrower manual reference, with intentionally unknown pixels

The provider has 17 classes with 80 images each, but supplies segmentation
ground truth only for a subset in 13 categories, mapped by `imlist.mat`. Its
three classification splits are distinct from the segmentation experiment's
allocation. [Provider page](https://www.robots.ox.ac.uk/~vgg/data/flowers/17/index.html),
[complete README](https://www.robots.ox.ac.uk/~vgg/data/flowers/17/README.txt)
(accessed 2026-09-08).

The original experiment retained 753 images after excluding fields of flowers,
very small/subsampled flowers and unclear foregrounds. Its trimaps deliberately
leave boundaries and ambiguous regions unlabelled; only labelled pixels enter
the overlap score. The modelling terminology includes petals, tepals and sepals
under one simplified name. Therefore neither all-pixel accuracy nor petal-only
anatomy follows from this reference. [Original 2007 paper, sections 2.3 and 3](https://www.robots.ox.ac.uk/~men/papers/nilsback_bmvc07.pdf)
(accessed 2026-09-08 through indexed primary text).

The expanded paper explicitly calls the 753 trimaps manual and describes
260 training images and 493 test images for its segmentation experiment.
These published counts are **not** a newly verified inventory of today's
archive or a mapping to present FCP data. [Author-hosted expanded paper,
sections 3.1-3.3](https://www.robots.ox.ac.uk/~vgg/publications/2009/Nilsback09a/nilsback09a.pdf)
(accessed 2026-09-08 through indexed primary text).

**Access and independence limits:** both READMEs were retrieved completely by
ordinary HTTPS with HTTP 200, despite the browsing service failing their
redirects. No image/mask body was requested. The previous audit's Oxford 17
archive HEAD success remains only historical reachability evidence. The
inspected provider pages/READMEs do not supply a per-photo redistribution licence
or observer/event ledger. A third-party software licence must not fill that gap.
Oxford's own duplicate-detection report also lists `Identicals = 4` for its
8,189-image Flowers collection; this audit did not open the illustrated matches
or infer their partitions. [Oxford FII report](https://www.robots.ox.ac.uk/~vgg/software/fii/)
(accessed 2026-09-08). Neither common class labels nor distinct archive names
prove disjointness from Monarda, JRC, FCP or model pretraining.

**Permissible future scope, if admitted:** agreement on the reference-labelled
foreground/background pixels, with a separate accounting of unknown pixels.
Calling unknown pixels background, eroding the reference after seeing results,
or treating this selected easy-foreground subset as representative iNaturalist
petal validation would change the question.

## 3. USDA fruit flowers: concrete source access, but a split conflict

The provider describes four image sets covering apple, peach and pear, with
binary flower/non-flower masks, orchard acquisition in 2016-2017, and a
dataset-specific **U.S. Public Domain** label. This dataset declaration is
separate from any algorithm/code licence. [USDA provider record](https://agdatacommons.nal.usda.gov/articles/dataset/Data_from_Multi-species_fruit_flower_detection_using_a_refined_semantic_segmentation_network/24852636)
(accessed 2026-09-08 via primary indexed text and its live API).

The paper distinguishes superpixel-based training labels from evaluation labels
made using human freehand positive/negative traces followed by RGR propagation.
The latter are human-guided, **algorithm-refined** masks, not fully manual
anatomical boundaries. It reports 100 AppleA training and 30 evaluation images,
and 18 AppleB, 24 Peach and 18 Pear images; AppleA originally contains 147 images.
AppleA selection is described as random images from one acquisition collection,
not plant/event-blocked allocation. AppleB uses a background panel. These facts
support only an orchard flower-region estimand, not focal-observation attribution
or global transfer. [Dias et al., sections III.A, IV, IV.A and table I](https://arxiv.org/html/1809.10080)
(accessed 2026-09-08). No provider performance value is transferred to FCP.

### New metadata-only verification in this follow-up

The live public API returned article `24852636`, version `1`, published
`2024-02-08T21:37:06Z`, modified `2025-11-21T23:22:13Z`, and **11 file entries**:
four source-image ZIPs, five label ZIPs and two split-text files. Current names
include suffixes such as `AppleA_Labels_1.zip` and `val_0.txt`, whereas the prose
retains earlier names and file-addition notes. Provider-advertised checksums are
not verified archive-body checksums. [Version-1 article API](https://api.figshare.com/v2/articles/24852636)
(accessed 2026-09-08; response explicitly reported version 1).

Only the two small split-text bodies were retrieved. Their exact-byte MD5s
match the corresponding provider entries; nonempty lines were counted without
renaming or deduplicating the files. **`IMG_0339.JPG` occurs in both lists.**
There are 129 distinct filenames in the union, not 130. This was checked again
with case-sensitive matching; no image content was opened.

| Public split resource | Bytes | Nonempty / unique names | Exact MD5 | Independently computed SHA-256 |
|---|---:|---:|---|---|
| [train.txt, file 44358653](https://ndownloader.figshare.com/files/44358653) | 1,400 | 100 / 100 | `38c13cfcc80d6e70442fd9c3729d190e` | `cadbea55cf2232b94b44429f7a3f4f72c4b14cfd2170c0d2f844e1ca870f1abf` |
| [val_0.txt, file 44358659](https://ndownloader.figshare.com/files/44358659) | 390 | 30 / 30 | `f13899724a34ad6a0a12cb74eaee83e2` | `959181b42150e82b1479aff5c7221970aa6c1c1c5b43dfad1fe7bf11e62d1643` |

Both text resources were accessed with HTTP 200 on 2026-09-08. Their contents
were inspected in memory, not added to the repository. This exact filename
overlap is a source-list finding, **not** proof of which bytes were used in the
published training or proof of model leakage. It prevents treating the current
published lists as an already verified disjoint split. Do not silently delete
the overlapping entry, invent a corrected split, or infer independent events
from the other distinct filenames.

The subsequent main-task verification retained both exact text receipts, with
the same SHA-256/MD5 identities, and implemented an offline reproducibility
check. See [the retained metadata audit and next measurement design](RGFCA_MEASUREMENT_TARGET_AND_VALIDATION_V2.md).
The original split lists are unchanged; this follow-up still acquired no images
or masks and did not admit a benchmark.

HEAD-only checks of the advertised [AppleA image archive](https://ndownloader.figshare.com/files/44358617)
and [AppleA evaluation-label archive](https://ndownloader.figshare.com/files/44358644)
returned HTTP 403 on 2026-09-08. No GET body was attempted for either archive.
Thus metadata and split-text access are demonstrated, while full paired-archive
access and membership remain **unverified**; the HEAD response is not proof
that a normal GET would also fail. No access restriction was bypassed.

**Admission consequence:** resolve the list discrepancy and exact image/mask
census, recover source-tree/event grouping where available, and audit overlap
with all FCP development/qualification/opened frames before calling any subset
independent. Different orchard provenance from JRC/Monarda is suggestive source
separation, not an image-level or foundation-training disjointness audit.

## 4. Next gate supported by this evidence

These are design recommendations, not acquisition or scoring authorization:

1. Specify whether the new endpoint is focal petals/tepals or generic visible
   flower foreground. Preserve a separate category for uncertain or non-focal
   tissue. A generic reference may support the latter endpoint but cannot
   adjudicate the former by changing its name.
2. Resolve the reference's ontology, rights scope, full paired manifest and
   source-group allocation **before** any model development or scoring. Keep
   missing/ambiguous/duplicate records in an explicit ledger rather than a
   successful-download subset. Preserve the USDA split conflict as reported.
3. Separate prospective development from a genuinely independent validation
   source/group allocation. Record prior use of Monarda, JRC, all opened FCP
   frames and the six FlowerMask schema examples. Unknown pretraining overlap
   remains unknown; do not encode it as `independent=true`.
4. Freeze the error budget and full-denominator analysis for the actual target.
   Do not reuse the Monarda operational floors as externally validated petal
   standards or choose a new benchmark by whichever supplies a passing score.

No examined resource closes both anatomical and target-domain validation in
this bounded audit. Measurement validity, calibrated colour and ecological
replication are separate questions; a future narrow engineering pass would not
retroactively validate the completed atlas or explain its non-support.

## Read depth and actions deliberately not taken

- Read the existing target-domain, JRC-source and completed Monarda-result
  documents; reused the saved FlowerMask terminal ledger without reacquisition.
- Read current provider pages, both Oxford READMEs, USDA article metadata and
  two split-text bodies. Inspected primary-paper text/search excerpts only;
  direct Oxford PDF opens were unreliable, so relevant sections are identified
  with that read-depth limitation rather than a claimed full-paper audit.
- No dataset photograph or mask was downloaded, decoded, displayed or
  rasterized in this follow-up. No embedded FlowerMask image string was fetched
  again. No model was loaded or run; no geography or image-derived outcome was
  joined. Download-link HEAD probes transferred no archive bodies.
- No new annotation campaign, training, login, author contact, paid access,
  rights change, benchmark execution, model selection or ecological test was
  initiated. Only this audit note was created.
