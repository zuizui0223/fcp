# RGFCA JRC training-source audit

Audit date: **7 September 2026**. This primary-source follow-up addresses the
JRC Flower Detection dataset's rights, citation identity and original sampling
units. It does not certify legal clearance or reproduce FCP qualification.
No training, test or reserve images were fetched; no scientific outcome files
were inspected and no training, evaluation, relabelling or retuning was done.

## Outcome

The dataset-specific provider notice explicitly applies **CC BY 4.0** to
copyright and/or sui generis rights in this dataset, requiring appropriate
credit and identification of changes. This is stronger, more specific evidence
than the catalogue's generic European Commission reuse label. It is not a
licence for unrelated detector/segmenter software or a warranty of flower-pixel
accuracy. The complete notice was read directly over HTTPS.
[Dataset copyright notice](https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/DRLL/flower_detection/copyright.txt)

Exact frozen FCP implementation and qualification outcomes are a separate
audit; this note does not infer them from provider evidence.

## 1. Citation identities and conflicting date fields

### Dataset

**European Commission, Joint Research Centre (2026). Flower Detection [Dataset].**
[DOI: 10.2905/JRC.2XJ67GR](https://doi.org/10.2905/JRC.2XJ67GR)

The catalogue's current suggested citation uses 2026. Its other fields report
issued **2022-01-01**, created **26 October 2022**, and modified **8 July 2025**.
The catalogue describes 500 grassland-patch images with one flower category,
distributed as 400 training and 100 test images.
[Official JRC dataset record](https://data.jrc.ec.europa.eu/dataset/caa582b7-7f45-4748-9223-08e5f145a4a6)

This is not merely an assumed dynamic web-page year: the DOI registration API
also returned `publicationYear: 2026` and an `Issued` date of `2022-01-01`, with
European Commission, Joint Research Centre as creator/publisher. Its empty
`rightsList` is incomplete registration metadata, not evidence overturning the
explicit dataset notice.
[DataCite DOI registration metadata](https://api.datacite.org/dois/10.2905/JRC.2XJ67GR)

The provider README is dated **10 October 2022** and attributes the source images
to the 2018 LUCAS grassland module coordinated by Eurostat. It describes manual
CVAT bounding-box annotations and COCO-format data. Its reference to
Elvekjaer et al. (2022) predates the final publication and is not the correct
year for the final journal article.
[Complete provider README](https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/DRLL/flower_detection/readme.txt)

**Citation decision:** use the currently registered 2026 dataset citation with
an explicit note that the source README/issued date are 2022 and access was on
7 September 2026. Do not silently change an immutable earlier citation or
interpret 2026 as evidence of new image collection, a revised dataset, or a
new download. Version/content equivalence requires the frozen acquisition
manifest, which was not re-audited here.

### Final research article

**Elvekjaer et al. (2024), final research article.**
[DOI: 10.1002/2688-8319.12324](https://doi.org/10.1002/2688-8319.12324)

The publisher reports first publication on **23 May 2024**; the JRC repository
displays **21 October 2024**. These dates do not justify citing the final paper
as 2022.
[Publisher record](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1002/2688-8319.12324),
[JRC137937 publication record](https://publications.jrc.ec.europa.eu/repository/handle/JRC137937)

## 2. Original study: points before slices

Methods §2.1.2 partitions 250 expert-surveyed European LUCAS points before
slicing; same-point slices stay together. Flower-present top-down images were
selected; the visible-petal category excludes grasses. Results §3.1 reports one/two/three retained
slices for 6%/88%/6% of original images.
[Methods §2.1.2; Results §3.1](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1002/2688-8319.12324)

| Original provider partition | Original points | Final slices |
| --- | ---: | ---: |
| Training | 150 | 300 |
| Validation | 50 | 100 |
| Test | 50 | 100 |

Combined development represents **200 points / 400 slices**; test:
**50 points / 100 slices**. Metrics pool box detections. Section 2.2.1 loosely
says test-set; §3.2/Tables 1–2 specify validation selection and independent
testing. These are provider Faster R-CNN results, not FCP performance.
[Methods §§2.1.2–2.2.1; Results §3.2 and Tables 1–2](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1002/2688-8319.12324)

### Interpretation for the FCP audit, not a new result

- Flower-present source-image selection does not establish image-level
  specificity on flower-absent scenes.
- A slice count is not a count of independent original field visits. Point
  disjointness is also not, by itself, proof that nearby points are spatially
  independent or that a geographically distant region was held out.
- The public dataset's 400-image training label and the article's original
  training/validation distinction should both be recorded. How FCP uses that
  development pool must be established from its frozen code and manifests,
  not assumed to reproduce the provider's tuning procedure.
- Do not transfer provider Faster R-CNN or Pl@ntNet success rates, thresholds
  or uncertainty to a different estimator. This note intentionally reproduces
  no provider performance estimate and no FCP qualification statistic.

## 3. Rights scope and attribution

The local release's source ledger should retain the provider, dataset title,
DOI/source URL, explicit licence link, accessed notice and a description of any
changes actually made to redistributed material. The retrieved notice bears
European Union copyright for 1995–2026; its current year range is not an image
collection date.
[JRC dataset-specific notice](https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/DRLL/flower_detection/copyright.txt)

CC BY 4.0 permits sharing and adaptation within the licensed rights, subject to
its conditions. When sharing covered material, retain supplied attribution and
notices, link the licence/material where practicable, and mark modifications.
The licence addresses database rights but does not grant every privacy,
publicity, patent or trademark right or imply endorsement; it disclaims
warranties. Whether a particular trained model or derived output constitutes
adapted material is not resolved by this audit.
[CC BY 4.0 legal code, §§2–5](https://creativecommons.org/licenses/by/4.0/legalcode.en)

The licence is attached to the identified JRC Flower Detection dataset. It does
not automatically relicense all LUCAS assets, third-party code, model weights,
iNaturalist photographs, or the entire FCP submission. Keep this training-data
notice separate from the software/model and photo rights documented in
[RGFCA_MEASUREMENT_PROVIDER_AUDIT.md](RGFCA_MEASUREMENT_PROVIDER_AUDIT.md).

## 4. Measurement and generalization ceiling

The source is an object-detection dataset whose advertised reference annotations
are bounding boxes, not pixel-level petal masks or calibrated floral spectra.
The inspected provider description/README supplies no reference standard for
reflectance, UV, or pollinator-perceived colour.
[Dataset description](https://data.jrc.ec.europa.eu/dataset/caa582b7-7f45-4748-9223-08e5f145a4a6),
[README](https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/DRLL/flower_detection/readme.txt)

Accordingly, the following are **audit inferences**, not new measurements:

- Reference-box containment cannot establish petal-mask intersection-over-union
  or pixel purity: a correct enclosing rectangle can include non-floral pixels.
- Box-annotation performance and stable RGB summaries do not establish
  calibrated flower reflectance or human/pollinator visual equivalence.
- A European grassland development resource does not establish global,
  taxon-uniform accuracy across opportunistic photos, camera conditions,
  phenological stages or growth forms. The target-domain error ceiling requires
  separate target-relevant evidence.
- A reuse licence supports lawful reuse within its terms, not scientific
  accuracy, botanical completeness, ecological mechanism or official endorsement.

These limits do not invalidate a separately qualified frozen estimator; they
define what this provider audit cannot add to that estimator's evidence.

## Read-depth and unresolved evidence

| Source | Access/read depth on 7 September 2026 | Remaining limit |
| --- | --- | --- |
| JRC dataset catalogue | Complete rendered record and date/citation fields | No image/annotation downloads or fresh identity check |
| `copyright.txt` and `readme.txt` | Complete text fetched directly by HTTPS | Current notice, not proof of every historical served byte |
| DataCite API | Creator/title/publisher/year/date/rights/URL fields | Registration metadata does not reconcile differing date semantics |
| Final Wiley article | Full Methods §§2.1–2.3, Results §§3.1–3.4, relevant tables, Discussion §§4.1–4.3, data availability and bibliography read through web extraction | No original figure images fetched; no independent reconstruction of source data |
| Supporting information | Publisher advertises `eso312324-sup-0001-Supinfo.pdf`, with parameter/performance-table descriptions | Link retrieval failed; PDF contents **not read** |
| JRC137937 repository | Public abstract, authors, DOI and repository date | No separately downloadable final file used |
| CC BY 4.0 legal code | Relevant scope, attribution, database-rights and warranty sections read | Not legal advice or a submission-wide clearance |

The publicly readable Methods URL is the publisher's
[full-text HTML endpoint](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1002/2688-8319.12324).
Web extraction succeeded; an ordinary direct HTTP request encountered a
Cloudflare 403 challenge. The supplement remains unread, and neither its
contents nor an image-by-image point mapping are silently inferred from the
main article. Exact frozen FCP training/holdout membership, qualification counts,
and any point-cluster uncertainty analysis remain the parent audit's remit.
