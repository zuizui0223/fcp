# Observer-disjoint reproducibility and white–nonwhite geometry of photo-derived flower-colour diversity

Working manuscript — incomplete, not submission-ready. Numerical evidence is
frozen; references, figures and measurement qualification remain outstanding.
The companion is [Supplement](FCP_H123_SUPPLEMENT.md). This is the active H1–H3
draft; older RGFCA and literature-comparison manuscripts are historical records.

## Abstract draft

Repeated photographs may characterize within-species colour diversity, but
sampling reproducibility must be distinguished from biological measurement
accuracy. We evaluated a four-state diversity score in species-disjoint,
high-depth discovery and reserve cohorts, examined recurrent directions of
within-species palette variation, and tested bounded phylogenetic and geographic
associations. In the 363-species reserve, a frozen observer-disjoint analysis
recovered species rankings with Spearman rho 0.811 (bootstrap 95% interval
0.765–0.847). Existing-cohort targeted analyses supported a white–nonwhite
direction beyond a coarse-state-preserving null, including the strict reserve
analysis (p=0.008). This particular direction was identified after inspection
of broader geometry, so it is not an untouched prospective confirmation.
Broad phylogenetic signal and the discovery association with sampled geographic
span were not supported in reserve. These findings support reproducibility and
constrained geometry of the current image-derived measurements; unresolved
localization and highlight confounding prevent treating them as validated
biological morph frequencies or identified ecological mechanisms.

## Introduction

The manuscript separates three questions: whether species-level photographic
colour diversity is reproducible across observers (H1), whether within-species
palette variation shares a recurrent direction (H2), and whether broad
phylogenetic structure or sampled geographic span explains the measured diversity
(H3).

Citizen-science photographs can provide useful visible-colour information, but
their accuracy depends on the measurement target and photographic conditions.
Laitly et al. (2021) compared photographic and controlled measurements and found
that aggregation could improve agreement; this does not establish a universal
sampling threshold for within-species diversity. Luong et al. (2023) demonstrated
landscape-scale flower-colour analysis using manually selected Erysimum petal
pixels, with a separate colour-correction assessment. These precedents establish
feasibility, not the accuracy of our automated multi-species estimator.

The distinction matters because calibrated reflectance measurement requires
additional acquisition information, including linearized image measurements and
reference standards (Troscianko & Stevens, 2015). A reproducible descriptor from
ordinary photographs need not be calibrated reflectance or a pollinator's colour
signal. We therefore evaluate sampling reproducibility and observed palette
geometry separately from their possible biological interpretation.

## Methods — evidence-backed outline

The global frame contains 42,111 species, but the high-depth discovery (369
species) and reserve (363 species) are eligibility-selected validation cohorts,
not a representative global sample. They must not be used to estimate worldwide
polymorphism prevalence. The 34-species literature comparison and six-species
spatial study are not analyses or supporting evidence in this manuscript.

H1 uses D = 1 - sum(p_k squared) for white, yellow/orange, red/pink and
blue/purple image states; mixed/uncertain is excluded from those four states.
The direct analysis uses observer-disjoint halves, a frozen minimum of 20
classifiable observations per half, and reserve decision criteria of rho >=0.80,
bootstrap lower bound >0.70 and permutation p<0.001. Repeated-partition and
later strict-split results are reported separately, not substituted post hoc.

H2 uses the frozen nine-colour palette, unlabeled two-mode construction and
sign-invariant displacement axes. The fixed targeted white contrast and W
statistic are evaluated against a null preserving coarse-state composition and
global coarse-state-to-palette structure. Primary and strict mode thresholds
remain 10% and 20%. Full operational definitions and frozen inputs must be
cross-checked against their original protocols before submission.

H3a uses the frozen S1–S3 phylogenetic placement scenarios. H3b tests sampled
photographic span, not true biological range size. The reserve is the replication
cohort. Failed replications are retained without searching for replacement
predictors.

## Results

### H1: repeatability, not proof of biological accuracy

The direct observer-disjoint reserve analysis includes all 363 eligible species:
rho=0.8109164415, bootstrap interval [0.7648383249, 0.8474976907],
permutation p=0.0000499975, and Lin CCC=0.8582944314. All frozen decision
conditions pass. The finite-sample-corrected D gives rho=0.8116141836.

Separately, 200 observer-disjoint partitions yielded median reserve rho=0.789103
and fifth percentile=0.765165. A later strict deterministic split yielded
rho=0.792728 and failed its 0.80 criterion. These results constrain the claim:
repeatability is substantial, but not perfect or invariant to every partition.
The historical approximately 0.971 value is not evidence for current D.

### H2: a narrow, retrospectively targeted pattern

The targeted white–nonwhite test gives W=0.514625 in discovery (152 species)
and W=0.514586 in reserve (129 species), with structured-null p=0.001 in
both. Under the strict threshold, discovery (75 species) gives W=0.542355,
p=0.001; reserve (65 species) gives W=0.510517, p=0.008. After removal of
the white contrast, the residual structured-null analyses do not support a
general recurrent nonwhite hue direction. This is an existing-cohort targeted
result, not prospective confirmation of the specific white axis.

### H3: bounded explanatory tests do not replicate

Reserve phylogenetic signal is unsupported across S1–S3 (raw Blomberg K
permutation p=0.272, 0.413 and 0.267). The discovery sampled-span association
(rho=0.179879, p=0.000900) does not replicate: reserve rho=-0.002586,
p=0.958602. These tests do not establish that evolutionary history or true
geographic range size is biologically irrelevant.

## Discussion — limits governing interpretation

Observer-disjoint repeatability establishes a property of the measurement under
the sampled design. It does not establish accurate focal-petal localization,
correct image-level biological classification or population morph frequencies.
The present estimator pools retained flower regions; co-photographed nonfocal
flowers can contribute. The completed Monarda region-agreement gate failed and
cannot be repaired by relabeling, selecting a successful subset or retuning on
the same images. White–nonwhite geometry could also be affected by digital
highlight failure; the proposed measurement-control test has no result yet.

External photographic studies do not remove these limitations. In particular,
Luong et al. (2023) manually avoided unsuitable petal pixels and studied a
restricted hue range; that validation cannot certify an automated white–nonwhite
contrast. Laitly et al. (2021) likewise does not establish that more photographs
eliminate all colour-measurement error. Colour-space conversion alone supplies
neither reflectance calibration nor missing ultraviolet information
(Troscianko & Stevens, 2015).

No pigment pathway, adaptive direction, pollinator mechanism or climatic cause
is identified. A new prospective white-axis test would strengthen inference
only if its measurement and chronology gates are satisfied. P500 is currently
metadata-only and must not appear as a completed analysis in this manuscript.

## Outstanding before submission

Complete primary-source citations; verify all numerical source mappings and
eligibility denominators; generate and inspect dedicated H1–H3 figures; complete
the companion supplement; resolve measurement qualification and P500 execution
gates, or explicitly retain their unresolved status and reassess the resulting
claim scope. No submission readiness or acceptance is asserted.

## References — verified initial set

Laitly, A., Callaghan, C. T., Delhey, K., & Cornwell, W. K. (2021). Is color data
from citizen science photographs reliable for biodiversity research? Ecology
and Evolution, 11(9), 4071–4083. https://doi.org/10.1002/ece3.7307

Luong, Y., Gasca-Herrera, A., Misiewicz, T. M., & Carter, B. E. (2023). A pipeline
for the rapid collection of color data from photographs. Applications in Plant
Sciences, 11(5), e11546. https://doi.org/10.1002/aps3.11546

Troscianko, J., & Stevens, M. (2015). Image calibration and analysis toolbox – a
free software suite for objectively measuring reflectance, colour and pattern.
Methods in Ecology and Evolution, 6(11), 1320–1331.
https://doi.org/10.1111/2041-210X.12439

Access details and non-transferability limits are recorded in the
[measurement-reference audit](FCP_H123_MEASUREMENT_REFERENCES.md). Statistical,
phylogenetic, image-model and wider ecological references remain to be audited.
