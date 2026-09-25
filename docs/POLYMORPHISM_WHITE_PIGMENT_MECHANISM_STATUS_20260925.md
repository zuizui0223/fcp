# White/pigmented proximal-mechanism audit status — 2026-09-25

## Third-cohort targeted search

The fixed OpenAlex audit completed successfully in workflow run **36106047401** (artifact **10850239977**).

- target species: 281
- fixed OpenAlex requests: 94
- failed requests: 0
- candidate works: 186
- species with any candidate: 64
- automated candidate A species: 2
- automated candidate B species: 1
- species with temperature-linked candidate text: 10

Automated tiers were navigation aids only.

## Natural-FCP expansion

The follow-up literature-derived natural-polymorphism panel completed in workflow run **36106592537** (artifact **10852105007**, digest `sha256:79077be039562e1e418a19a0037e61d5753617936e7f3de47a540d46eb8acab5`).

- natural-FCP species: **111**
- high-confidence literature-supported species: **97**
- candidate works: **252**
- species with candidates: **51**
- high-confidence species with candidates: **44**
- automated candidate A species: **7**
- automated candidate B species: **3**

The A/B automation was then manually adjudicated at the species level. Nine unique species entered manual review because *Vicia faba* appeared in both A and B.

## Strict natural molecular Tier A

Three species pass the strict boundary of (i) natural population / naturally occurring colour variation and (ii) direct proximal molecular evidence.

### 1. Ipomoea purpurea

Naturally occurring albino and pale/ivory colour variants implicate more than one lesion in the anthocyanin pathway:

- **CHS loss-of-function** can produce albino flowers;
- **bHLH regulatory disruption** reduces expression of downstream anthocyanin structural genes.

Interpretation: different mutations can converge on reduced floral pigmentation through the same biosynthetic network.

### 2. Parrya nudicaulis

A widespread natural purple-white polymorphism across Alaska has direct expression evidence implicating **petal-specific CHS downregulation / cis-regulatory control** near the threshold of the anthocyanin biosynthetic pathway.

Interpretation: white flowers can arise through regulatory reduction of pathway flux rather than obligatory coding loss of a pathway enzyme.

### 3. Pleione limprichtii

Natural Huanglong populations contain rose-purple, pink and white individuals. Metabolomic and transcriptomic analyses implicate:

- **PlFLS** in white formation;
- **PlANS / PlUFGT** in pigmented-flower intensity;
- a candidate **MYB-bHLH-WD40 regulatory complex** in colour variation.

Interpretation: the white endpoint again maps to altered anthocyanin/flavonoid pathway allocation and regulation, but through a different specific molecular route.

## Other manually adjudicated candidates

- **Boechera stricta**: strong natural ecological evidence. Drought can induce pigmentation in white-flowered lineages and colour has a genetic component, but the audited study does not identify a proximal causal pigment gene. Keep as ecological mechanism evidence, not molecular Tier A.
- **Silene littorea**: natural floral anthocyanin/flavonoid variation is directly measured, but no single causal white-state gene is resolved. Biochemical Tier B only.
- **Abronia fragrans**: natural white/pink polymorphism and ecological selection evidence, but the pigment system is **betalain**, not anthocyanin. This is an important bound against claiming a universal anthocyanin mechanism.
- **Digitalis purpurea**: genome/pigment-pathway resource in an ornamental context; strict natural white/pigmented segregation is not established by the audited evidence.
- **Medicago sativa**: white-purple transcriptomic/metabolomic comparison is mechanistically interesting, but strict wild natural-population provenance is unresolved for this cultivated species.
- **Vicia faba**: TT8/TTG1 genetics are direct but come from crop/breeding systems; excluded from strict natural-FCP evidence.

Machine-readable adjudication is stored in `results/polymorphism_natural_fcp_pigment_audit_20260925/manual_adjudication.csv`.

## Current synthesis

The strongest proximal-mechanism statement is now:

> **Natural white/pigmented flower-colour polymorphisms repeatedly map to reduced or rerouted pigment-pathway output, with multiple distinct regulatory or structural routes converging on the achromatic endpoint.**

For anthocyanic taxa, the verified examples converge on the anthocyanin/flavonoid network, but **not on one universal gene**. The *Abronia* betalain case prevents promotion of “anthocyanin loss” to a universal angiosperm white-flower mechanism.

This molecular convergence is compatible with the New Phytologist paper's recurrent white/non-white phenotype-space axis, but it remains a separate mechanism line. It does not establish evolutionary transition direction: the independent OpenTree endpoint test did **not** support recurrent nonwhite→white asymmetry.

## Ecological discrimination remains separate

The current ecological tests do not identify one universal selective cause:

- long-term BIO5 shows a broad-scale white-state sorting association;
- local/same-observer diagnostics weaken a universal local heat-selection interpretation;
- the prospective seasonal-heat full gate fails narrowly at the same-cell test;
- stable bee-community predictors do not explain the broad BIO5 association.

Thus the current best two-level model is:

1. **proximal accessibility:** multiple molecular routes can repeatedly alter pigment-pathway output and reach a white/achromatic phenotype;
2. **ecological sorting:** temperature and other local factors may alter where white states occur or persist, but no universal pollinator or environmental selective cause is confirmed.

The frozen New Phytologist manuscript remains unchanged.
