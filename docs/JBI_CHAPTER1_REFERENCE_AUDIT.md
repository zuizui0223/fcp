# JBI Chapter 1 — targeted reference audit

Audit date: 2026-09-03  
Scope: contextual and methodological support only. No reference added through this audit may alter the frozen sample, phenotype representation, null model, thresholds, spatial support, branch order or scientific decision.

## Why this audit is needed

The current Chapter-1 draft already has a coherent empirical result, but its bibliography is unusually sparse for a Journal of Biogeography Research Article. The goal is not to inflate citation count. The draft needs stronger support in three specific places:

1. geographic intraspecific flower-colour variation can reflect multiple adaptive and historical processes;
2. community photographs can provide scalable colour phenotypes, but photographic variation and sampling bias require explicit measurement caution;
3. random-labelling inference is a standard way to ask whether marks are spatially organized while holding locations fixed.

The four additions below are sufficient for those gaps. None is used to justify a post-result analysis.

## Priority additions

### 1. Koski & Galloway 2020 — geographic flower-colour variation

**Citation**  
Koski, M.H. & Galloway, L.F. (2020). Geographic variation in floral color and reflectance correlates with temperature and colonization history. *Frontiers in Plant Science*, **11**, 991. https://doi.org/10.3389/fpls.2020.00991

**What it supports**

A range-wide study of *Campanula americana* separated geographic patterns in petal reflectance and pollinator-perceived colour and related them to climate and post-glacial history. It is useful because it explicitly shows that geographically organized floral colour can reflect both contemporary abiotic conditions and historical population structure rather than one universal mechanism.

**Recommended placement**

- Introduction paragraph beginning `Flower colour is a tractable system...` after the sentence listing abiotic stress, gene flow, drift and demographic history.
- Optionally the next paragraph on broad-scale geographic variation.

**Claim ceiling**

Use only as precedent for multiple possible drivers of geographic flower-colour structure. Do not imply that the present Stage-A pattern is caused by temperature or colonization history.

### 2. Luong et al. 2023 — colour extraction from photographs

**Citation**  
Luong, Y., Gasca-Herrera, A., Misiewicz, T.M. & Carter, B.E. (2023). A pipeline for the rapid collection of color data from photographs. *Applications in Plant Sciences*, **11**, e11546. https://doi.org/10.1002/aps3.11546

**What it supports**

A directly relevant plant-image workflow for extracting colour data from photographs, including community-science imagery, while discussing practical variation in image conditions and sampling. This provides method context for treating photographs as scalable but operational rather than spectrophotometric colour measurements.

**Recommended placement**

- Introduction paragraph on community photographs, beside the existing photographic-record precedents.
- Discussion subsection `Community photographs as comparative trait data`, especially the caution that image-derived visible colour is not calibrated spectroscopy.

**Claim ceiling**

Use as methodological precedent, not as validation of the exact Florence/ROI pipeline used here.

### 3. McKenzie, Church & Hopkins 2026 — high-throughput iNaturalist flower-colour phenotyping

**Citation**  
McKenzie, P.F., Church, S.H. & Hopkins, R. (2026). High-throughput iNaturalist image analysis reveals flower color divergence in *Monarda fistulosa*. *The American Naturalist*, **208**, 101–109. https://doi.org/10.1086/739413

**What it supports**

A recent high-throughput study that processed more than 40,000 community-science photographs to identify flowers, extract colour and test geographic colour divergence within a widespread plant. It is the closest published precedent for the scale-and-image logic of the present programme.

**Recommended placement**

- Introduction paragraph on community photographs and geographic colour variation.
- Discussion subsection `Community photographs as comparative trait data` as evidence that large community-image datasets can support geographic flower-colour phenotyping when the measurement path is explicitly validated.

**Claim ceiling**

Use as an external methodological precedent only. Its successful taxon-specific pipeline does not imply that ROI-v4 should transfer across the 200-species terminal atlas; the terminal measurement-completeness NE remains unchanged.

### 4. Law et al. 2009 — random labelling and spatial marks

**Citation**  
Law, R., Illian, J., Burslem, D.F.R.P., Gratzer, G., Gunatilleke, C.V.S. & Gunatilleke, I.A.U.N. (2009). Ecological information from spatial patterns of plants: insights from point process theory. *Journal of Ecology*, **97**, 616–628. https://doi.org/10.1111/j.1365-2745.2009.01510.x

**What it supports**

The point-process distinction between spatial locations and marks, including random-labelling nulls that retain the observed point pattern while randomizing properties over locations. This is the appropriate conceptual citation for the Stage-A logic of fixing species-specific geography while permuting complete colour vectors over those locations.

**Recommended placement**

- Methods subsection `Random-labelling null and prospective gate`, immediately after the description of holding locations and graph edges fixed while permuting complete colour vectors within species.
- Optionally Introduction paragraph explaining why unconditional global shuffling would test species composition rather than colour organization.

**Claim ceiling**

Cite for the null-model concept, not as evidence that the particular nearest-neighbour statistic or equal-species aggregation is uniquely optimal.

## Existing references — retain

The following six references already serve distinct purposes and should remain unless a final reference-by-reference audit finds a bibliographic error:

- Dalrymple et al. (2020): macroecological flower-colour context;
- Delph & Kelly (2014): balancing selection and maintenance of within-species variation;
- Farquhar et al. (2023): crowdsourced photographs and geographic colour polymorphism;
- Jansen et al. (2025): citizen observations and geographic colour polymorphism;
- Narbona et al. (2018): flower-colour polymorphism, maintenance and evolutionary interpretation;
- Trunschke et al. (2021): pollinator-mediated selection and broader flower-colour evolution.

## Do not add merely because they concern a focal species

Species-specific studies of one of the six focal taxa may be useful background, but they should not be added simply because the taxon matches. In particular, post-result citation of a focal-species climatic, pollinator or historical mechanism can make the manuscript look as though mechanism hypotheses were chosen after seeing the Stage-A species-level results. Such references should be included only when they support a genuinely necessary general statement or a clearly labelled species-specific biological note.

## Proposed citation map

| Manuscript location | Existing support | Add |
|---|---|---|
| Introduction: mechanisms of geographic flower-colour variation | Narbona 2018; Trunschke 2021 | Koski & Galloway 2020 |
| Introduction: community photographs / geographic variation | Dalrymple 2020; Farquhar 2023; Jansen 2025 | Luong et al. 2023; McKenzie et al. 2026 |
| Introduction: preserve species geography in the null | conceptual text only | Law et al. 2009 |
| Methods: Stage-A random-labelling null | no general method citation | Law et al. 2009 |
| Discussion: operational image-derived colour | current limitations text | Luong et al. 2023; McKenzie et al. 2026 |
| Discussion: multiple possible drivers of spatial organization | current general text | Koski & Galloway 2020 |

## Integration rule

These citations may be integrated into prose without changing any numerical or inferential statement. After integration, the manuscript guard must still preserve the frozen Stage-A and Stage-B tokens and claim boundaries. The reference audit is complete when:

1. the four priority references are cited only where their evidence is relevant;
2. full bibliographic entries are added consistently;
3. no new mechanistic inference is introduced;
4. manuscript consistency CI remains green.
