# Flower-colour polymorphism paper architecture

Frozen working architecture: 2026-09-13 JST

Evidence authority: `docs/POLYMORPHISM_MAINLINE_EVIDENCE_LEDGER_20260913.md`

## One-sentence paper

**Flower-colour polymorphism is a reproducibly measurable species attribute in high-information photographic data, and its within-species variation is geometrically constrained primarily along an achromatic-chromatic white-versus-non-white axis, while broad phylogenetic signal and sampled geographic span fail fresh replication.**

This is the mainline. Do not turn the manuscript back into a global shared-boundary paper or a predictor-hunting paper.

---

## Central conceptual move

Most macroecological flower-colour analyses reduce species to one representative colour or ask where colours occur. This paper instead treats **within-species colour diversity itself as a continuous phenotype** and separates three questions that are often conflated:

1. **Can polymorphism be measured reproducibly at species level?**
2. **When a species is polymorphic, is variation geometrically arbitrary or concentrated along particular colour-space directions?**
3. **Do broad lineage or geographic-opportunity predictors explain that species-level variation?**

The empirical answer from the current frozen evidence is asymmetric:

- yes for measurement reproducibility;
- yes, narrowly, for white-versus-non-white geometric concentration;
- no for the tested broad phylogenetic and sampled-span predictors.

That asymmetry is the paper's result rather than a problem to be rescued.

---

## Working title set

Primary working title:

**Flower-colour polymorphism is reproducible across observers and constrained along an achromatic-chromatic axis**

Alternative, more conceptual:

**Within-species flower-colour diversity is measurable but not geometrically arbitrary**

Alternative after prospective white-axis confirmation:

**A recurrent achromatic-chromatic axis structures flower-colour polymorphism across plant species**

Do not use the third title until the frozen white-axis target is confirmed prospectively beyond the current discovery/reserve cohorts.

---

## Abstract logic

### Background sentence

Flower colour is usually analysed as a species mean or categorical state, leaving unresolved whether within-species colour polymorphism can be measured reproducibly at macroecological scale and whether polymorphism occupies recurrent directions in colour space.

### Approach sentence

Use a 42,111-species global sampling frame and species-disjoint high-depth validation cohorts to define a continuous four-state polymorphism score, test observer-disjoint reproducibility, and characterize label-free within-species colour geometry under structured construction-preserving null models.

### Result 1

In the fresh reserve cohort, observer-disjoint split-half D has median Spearman rho 0.789 (5th percentile 0.765), corresponding to median Spearman-Brown reliability 0.882; discovery is independently consistent.

### Result 2

Label-free colour differences are strongly axis-concentrated and survive a coarse-state-preserving structured null in discovery and reserve. Audit localizes the recurrent component to white versus non-white; removing that axis eliminates structured-null support for residual hue concentration.

### Result 3 / boundary

Broad phylogenetic signal in D is unsupported across three frozen tree-placement scenarios, and a discovery association with sampled photographic span fails species-disjoint reserve replication.

### Conclusion sentence

Within-species flower-colour diversity can therefore be treated as a reproducible species phenotype in high-information photographic data, but its strongest replicated structure lies in the geometry of variation rather than in the broad predictors tested here.

Mandatory abstract qualification until prospective expansion: call the white-axis result "targeted/replicated in existing validation cohorts" or equivalent; do not imply the specific axis was prospectively specified before all current-cohort outcomes were inspected.

---

## Introduction — four-paragraph spine

### P1. The missing phenotype

Start from the conceptual gap: macroecology has rich maps of species distributions and trait means, but intraspecific polymorphism is often collapsed away. Flower colour is an unusually tractable test case because repeated image observations can, in principle, recover species-level distributions rather than a single modal state.

End with the problem: a species-level polymorphism score is only biologically useful if it is not dominated by observer identity or sampling idiosyncrasy.

### P2. Measurement before explanation

Introduce `D = 1 - sum p_k^2` over four frozen biological colour states as a continuous diversity phenotype. Explain why continuous D is preferable to an arbitrary polymorphic/monomorphic threshold.

State H1: independent observer sets should recover similar species rankings in D if the phenotype is measurable.

### P3. Geometry, not only amount

Two species with the same D can vary along biologically very different colour directions. Therefore the next question is whether polymorphism is directionless in continuous colour space or repeatedly organized along a lower-dimensional axis.

State H2 initially in broad, label-free form. Explain that the decisive null must preserve coarse-state composition and global state-to-palette structure, otherwise a geometric signal can be generated mechanically by the colour classification system itself.

Do not introduce white-versus-non-white as if it were the untouched preregistered H2. In the Introduction it can be framed as the axis identified by the frozen audit/targeted decomposition, with its retrospective status made explicit in Methods or Results.

### P4. Explanatory scope and falsification

Broad phylogenetic conservatism or geographic opportunity could organize species differences in polymorphism. Treat these as bounded explanatory tests, not required pillars of the paper.

State that the study uses species-disjoint reserve tests and preserves negative replications. This sets up the logic that a measurable and geometrically structured phenotype need not be captured by the broad predictors tested.

---

## Methods order

### 1. Global sampling frame and high-depth validation cohorts

Report the 42,111-species frame first, then explain why hypothesis-specific repeated-observation requirements lead to smaller high-depth cohorts. Clearly distinguish frame, opportunity, eligibility, discovery, and reserve.

Do not open Methods with "369 species were studied" because that recreates the false impression that 369 is the intended global sample.

### 2. Image-level biological states and continuous D

Define the four biological states and the treatment of `mixed_uncertain`. Define D and the >=40 classifiable gate. Keep this concise and auditable.

### 3. H1 observer-disjoint reliability

Describe observer-level splitting before giving any H2 geometry. Emphasize that balancing uses observer identity and observation counts only, never morph labels or D.

Report the frozen reserve decision thresholds.

### 4. H2 continuous palette and label-free two-mode construction

Define the nine-colour normalized palette, Hellinger transform, deterministic unlabeled two-means, primary/strict smaller-mode fractions, delta vectors, and axial concentration statistic.

### 5. Structured null and white-axis audit

This section is essential, not Supplement-only.

Explain three layers:

- isotropic reference;
- coarse-state-preserving structured null;
- white-axis projection / targeted fixed statistic.

Explicitly label the white-axis targeted test as post-audit for the current cohorts.

### 6. H3 bounded explanatory tests

Put phylogenetic signal and sampled-span replication after the positive measurement/geometry methods. These are scope tests, not the organizing method of the paper.

### 7. Reproducibility and frozen decisions

Briefly summarize species-disjoint cohort roles, input hashes, frozen protocols, and no-rescue rules. Detailed receipts go to Supplement/repository.

---

## Results order

### Result 1. Polymorphism is reproducibly measurable across independent observers

Lead with H1, not the distribution of D.

Core numbers:
- reserve median paired N = 329;
- reserve median split-half rho = 0.789;
- reserve rho q05 = 0.765;
- projected full-estimate reliability = 0.882;
- median CCC = 0.855;
- discovery median rho = 0.812.

Interpretation sentence: species differ consistently in D even when no observer contributes to both estimates.

Then show the continuous distribution of full-data D as descriptive context, not as prevalence.

### Result 2. Within-species colour variation is strongly axis-concentrated

Show broad label-free geometry first:
- discovery lambda1 = 0.541;
- reserve lambda1 = 0.535;
- reserve transport = 0.524;
- independent-axis dot = 0.986.

Then structured-null result:
- primary discovery p=.001;
- primary reserve transport p=.001;
- strict discovery p=.001;
- strict reserve transport p=.001.

This establishes that the concentration exceeds the frozen construction-preserving expectation.

### Result 3. The recurrent component localizes to white versus non-white

Show that 82-85% of primary H2 species have white in the top-two coarse pair.

Then show the crucial falsification:
- after white-axis removal, discovery and reserve structured-null support disappears (primary p=1.0/1.0; strict p=1.0/.962);
- non-white-only discovery concentration is not supported by the structured null.

Then give the fixed targeted white-axis statistic:
- primary discovery W=.515, p=.001;
- primary reserve W=.515, p=.001;
- strict discovery W=.542, p=.001;
- strict reserve W=.511, p=.008.

Interpretation: the supported recurrent geometry is narrow; it is not a generic recurrent hue axis.

### Result 4. Broad explanatory signals do not replicate

Report H3a and H3b together as inferential boundaries.

H3a:
- reserve Blomberg K permutation p=.272/.413/.267 across S1-S3;
- adjusted tests likewise unsupported.

H3b:
- discovery span-D rho=.180, p=.0009;
- reserve rho=-.0026, p=.959;
- adjusted reserve and S1-S3 rank-PGLS also null.

Do not present this section as a disappointing hunt for predictors. Its role is to show that the strong measurement/geometry conclusions do not automatically imply a broad phylogenetic or sampled-span explanation.

---

## Figure architecture

### Figure 1 — From global frame to a measurable polymorphism phenotype

Panels:
A. 42,111-species sampling-frame schematic and hypothesis-specific eligibility funnel.
B. Four biological colour states and definition of D.
C. Distribution of D in discovery and reserve high-depth cohorts.
D. Observer-disjoint split schematic showing that observers, not images, are the unit kept separate.

Purpose: establish what is measured and prevent 369/363 from being misread as a global random sample.

### Figure 2 — H1: observer-disjoint reproducibility

Panels:
A. Representative reserve split: D_A versus D_B with 1:1 line.
B. Distribution of rho over 200 reserve partitions with frozen thresholds marked.
C. Discovery versus reserve distributions of rho / Spearman-Brown reliability.
D. Agreement diagnostic: distribution of |D_A-D_B| or CCC.

Headline number belongs in the panel title/caption: reserve median rho = 0.789, q05 = 0.765.

### Figure 3 — H2: recurrent colour geometry and structured null

Panels:
A. Conceptual nine-palette / Hellinger / unlabeled two-mode delta construction.
B. Leading-axis concentration in discovery and reserve.
C. Observed lambda1 / frozen discovery-axis reserve projection against coarse-state-preserving structured-null distributions.
D. Discovery-reserve axis alignment.

Purpose: establish recurrent geometry before naming its biological direction.

### Figure 4 — H2 localization to the achromatic-chromatic axis

Panels:
A. Loadings / canonical white-versus-non-white contrast.
B. Targeted W against structured-null distribution in discovery and reserve, primary and strict gates.
C. Residual concentration after projecting out q_white, showing collapse relative to structured null.
D. Small non-white-only diagnostic, clearly labeled low-N/sensitivity.

The visual logic should be: **signal exists -> signal is white-axis -> remove white-axis and signal disappears**.

### Figure 5 or Supplement — explanatory boundaries

Preferred main-text treatment if journal space permits:
A. H3a reserve K across S1-S3 with permutation null/interval.
B. discovery versus reserve sampled-span-D effect sizes.

If narrative becomes crowded, move full H3 diagnostics to Supplement but retain one main-text panel or explicit Results paragraph showing that these hypotheses were prospectively tested and failed.

---

## Discussion order

### D1. Establish the phenotype

The first contribution is methodological/conceptual rather than a prevalence claim: repeated photographic observations can recover stable between-species differences in within-species flower-colour diversity even across completely disjoint observers.

### D2. The geometry is constrained, but narrowly

The second contribution is biological/descriptive: polymorphism is not directionless. The recurrent component is predominantly achromatic-chromatic, not a general shared direction among chromatic hues.

This is stronger than simply saying white morphs are common because the targeted statistic exceeds a null that preserves coarse-state composition and the global state-to-palette relationship.

### D3. What the result does not yet identify

The data do not establish whether white-versus-non-white transitions reflect pigment loss/gain, developmental switching, pollinator selection, abiotic selection, or repeated evolutionary transitions in either direction. These are mechanistic alternatives for future work, not interpretations to choose among here.

### D4. Measurable structure without a broad predictor

Use H3a/H3b negatives constructively: a stable phenotype and a recurrent geometric axis can exist without detectable broad phylogenetic conservatism or a replicating relationship with sampled photographic span under the current design.

Avoid the stronger statement that phylogeny/range do not matter biologically.

### D5. Prospective confirmation and scaling

End with the strongest next test: apply the now-frozen q_white, W statistic, mode construction, and structured null prospectively to a larger high-depth subset drawn from the 42,111-species frame. That test distinguishes a discovery-audit pattern from a general macroevolutionary regularity.

---

## Main-text versus Supplement boundary

### Main text

- 42,111 frame / validation-cohort distinction;
- D definition;
- observer-disjoint H1 design and result;
- H2 mode geometry;
- coarse-state-preserving structured null;
- white-axis localization and targeted result;
- concise H3a/H3b negative boundary.

### Supplement / repository

- all 200 H1 partition metrics;
- exact deterministic hashing/splitting details;
- H1 >=15-per-half sensitivity;
- full H2 primary/strict numeric receipts;
- isotropic-null diagnostics;
- non-white-only low-N diagnostics;
- tree-placement preflight and manifests;
- all H3a S1-S3 sensitivity outputs;
- H3b rank-PGLS details;
- workflow/run/artifact/hash provenance.

---

## Submission gate

### Paper is already coherent if frozen now

A defensible current paper exists with:

1. direct independent-observer validation of D;
2. replicated broad H2 geometry;
3. a structured null that rejects construction sufficiency;
4. localization to a white-versus-non-white axis;
5. explicit negative explanatory tests.

The main weakness is not H1 and not statistical robustness. It is the **post-audit status of the specific white-axis target**.

### Highest-value upgrade before ambitious submission

Prospectively evaluate the already-frozen white-axis statistic in a larger high-depth expansion from the 42,111-species frame. Do not change q_white, W, mode construction, primary/strict thresholds, or the structured null after selecting that expansion.

If that prospective test passes, the manuscript can elevate Claim 2 from targeted replicated evidence to prospective confirmation. If it fails, preserve the failure and keep the current cohorts as discovery/validation evidence rather than re-tuning the axis.

### Explicit stop rule

Do **not** spend the next analysis cycle searching for new H3 predictors. H3a and H3b have already served their purpose by bounding explanation. The next unit of information with the highest value is prospective H2 confirmation.