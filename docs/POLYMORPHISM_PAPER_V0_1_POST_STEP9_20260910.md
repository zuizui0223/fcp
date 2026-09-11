# FCP polymorphism paper v0.1 — post Step 9

Date: 2026-09-10 JST
Branch: `feat/polymorphism-paper-v0-1-post-step9`
Evidence base frozen through Step 9.

## Working title

**Flower-colour polymorphism is geographically organized within species without a shared global boundary**

Alternative shorter title:

**Species-specific geography of flower-colour polymorphism**

## One-sentence claim

Across two species-disjoint fixed-photo global frames, species with greater four-state flower-colour polymorphism show stronger within-species geographic organization of flower colour, even though neither a common boundary location nor a privileged colour direction is shared across species.

## What this paper is not

This is not a paper claiming:

- a universal global flower-colour boundary;
- a universal climate threshold;
- adaptation or local selection;
- formal phylogenetic signal;
- a recurrent privileged colour axis;
- angiosperm-wide prevalence of polymorphism.

The paper is about the **amount of within-species flower-colour polymorphism and its geographic organization**.

---

# Abstract — v0.1

Flower-colour evolution is often studied by asking which colours occur in which environments or whether particular transitions recur across lineages. Less is known about a different property: how much discrete flower-colour variation a species carries, and whether that variation is geographically organized within its range. We quantified four-state flower-colour polymorphism in fixed photo-derived global species frames and tested its association with within-species geographic colour structure. In a 369-species discovery frame, flower-colour diversity varied continuously among species, with a second colour morph contributing at least 10% of admitted observations in 46.6% of species. More polymorphic species showed stronger within-species geographic colour organization, but this organization did not correspond to a shared cross-species boundary or privileged recurrent colour direction. The association recurred in a species-disjoint 363-species reserve photo tranche and remained positive under observer and calendar controls, matched-background subtraction, sampled-geographic-span adjustment, and explicit separation of technical ROI/flip failures from ambiguous colour observations. Under a sensitivity model that assigned ambiguous-palette observations to the four existing morph classes, the spatial association persisted at both the exact minimum- and maximum-diversity completion endpoints. Genus-level taxonomic clustering of polymorphism also recurred across tranches but was more sensitive to ambiguity completion. These results suggest that flower-colour polymorphism is repeatably organized at the species level even when the geographic location and direction of colour differentiation are not shared among species.

---

# Introduction — argument structure

## Paragraph 1 — move the question from colour state to variation itself

Most comparative studies of floral colour ask which colour state evolves, how often a transition occurs, or which environmental or pollinator regime is associated with a particular colour. Those questions treat colour largely as a species-level state or transition. Yet many species contain substantial intraspecific colour variation. The amount of such variation may itself be an ecologically and evolutionarily meaningful species property.

**Target last sentence:**

> The unresolved question is therefore not only why one colour replaces another, but why some species carry substantially more flower-colour polymorphism than others.

## Paragraph 2 — separate amount, direction, and geography

Three distinct features of colour variation are often conflated:

1. the **amount** of within-species colour diversity;
2. the **direction** in colour space along which that diversity occurs;
3. the **geographic organization** of colour states within a species.

A species can be highly polymorphic without varying along the same colour axis as another species. Likewise, two species can each show strong geographic differentiation while their transitions occur in different places and directions. A common global boundary is therefore a stronger hypothesis than species-specific spatial organization.

**Target last sentence:**

> Distinguishing these levels makes it possible to ask whether geographic organization is repeatable at the species level even when boundary geography is not repeatable across species.

## Paragraph 3 — observational challenge and design

Large citizen-science image collections make global within-species comparison possible but create three major problems:

- unequal sampling opportunity;
- background/image-processing artefacts;
- ambiguous observations that cannot be safely assigned to discrete morphs.

The design therefore needs fixed photo budgets, species-disjoint replication, spatial nulls that retain each species' sampling geometry, background controls, and an explicit accounting of classification failure and ambiguity.

## Paragraph 4 — study questions and predictions

Frame three questions only:

**Q1. Is flower-colour polymorphism measurable as a continuous species-level property in a large fixed-photo frame?**

**Q2. Are more polymorphic species more strongly geographically organized internally?**

**Q3. Is any species-level tendency taxonomically structured, and how robust is that pattern to measurement ambiguity?**

Do not introduce climate, latitude, pollination, life form, recurrent M1-M2 direction, or common global boundaries as coequal hypotheses. Those enter as controls, rejected alternatives, or secondary diagnostics.

---

# Methods — paper-facing structure

## 1. Fixed-photo global species frames

### Discovery

- 369 species retained for the polymorphism analysis;
- each candidate species entered the measurement campaign with a fixed 100-photo raw denominator;
- four admitted morph classes: white, yellow-orange, red-pink, blue-purple;
- species required at least 40 classifiable observations for the primary species-level D frame.

### Reserve

- independent species-disjoint photo tranche;
- same fixed 100-photo raw denominator;
- 363 species met the corresponding classifiable threshold;
- reserve outcomes were used for replication rather than model/threshold search.

State clearly that species inclusion is conditional on the source/capacity design and that these frames are not random samples of all angiosperms.

## 2. Species-level polymorphism metric

For admitted four-state observations:

`D = 1 - sum(p_k^2)`

where `p_k` is the proportion of classifiable observations in morph k.

Retain finite-sample sensitivity:

`D_unbiased = n/(n-1) * D`.

The second-largest morph fraction is used only for interpretable threshold summaries (>=10%, >=20%). Continuous D is preferred for the main gradient analyses.

## 3. Within-species geographic colour organization

Use the existing within-species spatial statistic based on the association between geographic separation and flower-colour dissimilarity.

Important paper-facing language:

- each species has its own observed spatial rho;
- null randomizations shuffle colour structure within the species while retaining that species' sampling geometry;
- global tests aggregate species equally;
- a positive species-level rho means colour difference tends to increase with geographic separation within that species.

Avoid implying that all species share the same spatial boundary.

## 4. Discovery subset and continuous gradient

Report the original whole-frame spatial result first.

Then explain that the discovery re-aggregation asking whether more polymorphic species had stronger spatial structure was post-outcome exploratory:

- second morph >=10%;
- second morph >=20% sensitivity;
- continuous association between D and species spatial rho.

This transparency matters because the independent reserve tranche supplies the stronger replication evidence.

## 5. Species-disjoint reserve replication

The reserve uses the already-fixed polymorphism thresholds and spatial aggregation rules.

Report four reserve spatial responses:

- primary;
- observer-pair exclusion;
- calendar-quarter stratification;
- matched flower-minus-background differential.

The matched-background differential is important because it tests whether the spatial structure attributed to flower colour exceeds spatial structure in the image background.

## 6. Sampled geographic opportunity

Use sampled span only as an opportunity/control variable, not as true biological range size.

Discovery primary span:

`log1p(maximum_span_km_after_observer_cap)`.

Reserve span:

maximum pairwise great-circle distance among all 100 fixed reserve photo coordinates, log1p-transformed before morph/classifiability selection.

The reserve D-span correlation is essentially zero, which argues against making sampled span a general ecological predictor.

## 7. Measurement-process decomposition

Every species begins from 100 acquired photos. Partition every raw photo into exactly one of:

1. classified four-state morph;
2. clear technical ROI/flip failure;
3. ambiguous palette composition;
4. no biological palette mass.

Only ROI/flip gate failure is treated as clearly technical.

Do not describe `n_classifiable` as sampling effort. It is post-acquisition classification yield and is strongly coupled to ambiguous-palette observations.

## 8. Ambiguous-palette endpoint sensitivity

For a species with four-state counts `c_k` and A ambiguous-palette rows, consider the sensitivity model that each ambiguous row is actually an unobserved member of one of the same four morph classes.

Define exact species-level endpoints:

- `D_min4`: assign all A rows to the currently most abundant morph;
- `D_max4`: allocate A rows to equalize the four counts as much as possible.

These are exact minimum and maximum four-state D under that completion assumption.

Stress that downstream tests using all species at D_min4 or all at D_max4 are **uniform endpoint stress tests**, not exhaustive worst-case optimization over every species-specific latent allocation.

## 9. Genus-level taxonomic clustering

Use the equal-genus pair-weight statistic already frozen in Step 4:

- only genera represented by >=2 species contribute;
- every repeated genus receives equal total weight;
- lower within-genus D difference indicates clustering;
- species labels are permuted to form the null.

Call this taxonomic clustering, not phylogenetic signal.

## 10. Multiplicity and claim hierarchy

The prospectively defined Step-4 species-attribute family used Holm correction.

For later post-Step-4 robustness/replication analyses, do not pretend they were part of the original family. Report their inferential role explicitly:

- discovery exploratory;
- reserve replication;
- measurement/opportunity robustness;
- ambiguity endpoint sensitivity.

---

# Results — evidence order

## Result 1 — polymorphism varies continuously among species

Open with the descriptive species-level signal, not with a map.

Discovery:

- n=369 species;
- D min=0;
- D max=0.707645;
- second morph >=10% in 46.61%;
- second morph >=20% in 26.56%.

Interpretation:

> Flower-colour polymorphism is not naturally reduced to a rare binary exception in the admitted frame; it varies continuously among species.

Do not call 46.61% a global prevalence estimate.

## Result 2 — polymorphism does not follow one recurrent colour direction

Keep this short and negative.

Frozen M1-M2 directionality test:

- rho(D, recurrent M1-M2 variance share) = -0.1647;
- primary test unsupported;
- random orientation test unsupported.

The earlier positive SD(M1)/SD(M2) patterns track total colour spread rather than a privileged recurrent direction.

Interpretation:

> Greater polymorphism means more colour diversity, not expansion along a common global colour axis.

This sets up the key contrast with spatial organization.

## Result 3 — more polymorphic species are more geographically organized

### Discovery

- whole-frame equal-species spatial mean rho = 0.027021, p=0.001;
- >=10% second-morph subset: n=172, mean rho=0.034565, p=0.001;
- <10% complement: n=197, mean rho=0.020435;
- difference=0.014130, p=0.014;
- continuous rho(D, species spatial rho)=0.089213, p=0.034.

State explicitly that the subset/continuous comparison is exploratory in discovery.

### Reserve replication

- >=10% subset: n=151, primary mean rho=0.032299, p=0.001;
- observer-pair exclusion p=0.001;
- calendar-quarter p=0.001;
- matched-background differential mean rho=0.013622, p=0.005;
- >=10% minus <10% difference=0.011672, p=0.041;
- continuous rho(D, primary spatial rho)=0.101601, p=0.025.

Interpretation:

> The association between polymorphism and within-species spatial organization recurs in a species-disjoint photo tranche.

## Result 4 — sampled span and clear technical failure do not explain the spatial relationship

### Sampled-span control

Geometry-preserving Step 6b:

- discovery partial rho(D, spatial | span, n_classifiable)=0.170198, p=0.001;
- reserve primary=0.078186, p=0.080;
- reserve flower-minus-background=0.087180, p=0.037.

Then explain why `n_classifiable` is not a clean effort covariate.

### Pure technical-failure control

After replacing broad classification yield with the clearly operational ROI/flip failure rate:

- discovery primary partial rho=0.126637, p=0.007;
- reserve primary=0.099288, p=0.025;
- reserve flower-minus-background=0.116241, p=0.010.

This is the cleaner measurement-process result to emphasize.

## Result 5 — the spatial result survives ambiguity endpoint completion

Ambiguity interval width:

- discovery median=0.0742; q95=0.2533;
- reserve median=0.0732; q95=0.2790.

Observed D remains highly rank-correlated with both endpoints:

- discovery: rho(D,D_min4)=0.9963; rho(D,D_max4)=0.9730;
- reserve: rho(D,D_min4)=0.9969; rho(D,D_max4)=0.9610.

Reserve primary spatial:

- D_min4 partial rho=0.097078, p=0.029;
- D_max4 partial rho=0.125286, p=0.008.

Reserve matched flower-minus-background:

- D_min4 partial rho=0.116299, p=0.009;
- D_max4 partial rho=0.132785, p=0.006.

This is the strongest robustness result in the paper.

Interpretation:

> The spatial association is not dependent on whether ambiguous-palette observations are completed in the diversity-minimizing or diversity-maximizing direction under the four-state sensitivity model.

Do not say arbitrary ambiguity is fully identified.

## Result 6 — genus-level clustering is real but secondary and ambiguity-sensitive

Discovery Step 4:

- repeated genera=61;
- species in repeated genera=169;
- gain=0.202963;
- raw p=0.002300;
- Holm p=0.011499.

Species-disjoint reserve:

- repeated genera=54;
- species=146;
- raw gain=0.153654;
- p=0.021699.

Cross-tranche genus means among 23 shared repeated genera:

- rho=0.530632, p=0.009900;
- D_unbiased rho=0.521739, p=0.010950.

Sampled span alone does not remove the reserve pattern:

- gain=0.141642, p=0.022999.

Sampled span + clear technical failure:

- gain=0.121211, p=0.044848.

But ambiguity endpoints differ:

- D_min4: gain=0.110212, p=0.061647;
- D_max4: gain=0.141546, p=0.023799.

Meanwhile shared-genus mean concordance remains positive under both endpoints:

- D_min4 rho=0.559289, p=0.005900;
- D_max4 rho=0.571146, p=0.005200.

Interpretation:

> Polymorphism is not completely taxonomically exchangeable, but the strength of genus-level clustering depends more on how ambiguous observations are completed than does the spatial result.

Keep this out of the title and abstract conclusion sentence unless space permits a qualified secondary clause.

## Result 7 — no common geography is required

Use the earlier shared-boundary analyses as a contrast rather than a failed centerpiece.

Key message:

- within-species spatial structure is repeatedly recoverable;
- species can differ in where and in what colour direction transitions occur;
- global shared boundary/common-direction claims are unsupported or not identifiable under the current geometry.

Target result sentence:

> Repeatable species-level spatial organization therefore coexists with non-repeatable cross-species boundary geography.

---

# Discussion — paragraph plan

## Discussion 1 — main answer

Lead with the spatial result:

> The strongest result is not that flower-colour transitions recur in the same place, but that species carrying more colour polymorphism are more geographically organized internally.

Stress discovery + species-disjoint replication + ambiguity robustness.

## Discussion 2 — organization without shared boundary

Explain why this is biologically plausible rather than contradictory.

Different species can respond to:

- different environmental gradients;
- different pollinator assemblages;
- different historical dispersal barriers;
- different demographic histories;
- different genetic architectures.

These mechanisms could all produce within-species geographic structure without producing one global boundary map.

Do not choose among them with the present data.

## Discussion 3 — amount is not direction

Use the failed recurrent-direction result constructively.

The amount of polymorphism is a scalar species property, whereas the direction of variation in colour space can differ among species. This distinction explains why D can predict spatial organization even when M1-M2 recurrence fails.

## Discussion 4 — measurement uncertainty as part of the inference

This is a methodological strength of the paper.

Explain:

- fixed 100-photo denominator;
- technical failure separated from ambiguity;
- background matched control;
- exact species-level D_min4/D_max4 endpoints;
- reserve replication.

Do not claim the ambiguity problem is solved; say the main spatial conclusion is insensitive to two mathematically extreme uniform completions under the four-state sensitivity model.

## Discussion 5 — taxonomic structure, cautiously

Genus clustering suggests that polymorphism propensity is not fully exchangeable across taxonomy, but D_min4 weakens the reserve genus test beyond p=0.05. This makes genus structure a hypothesis-generating result for future phylogenetically explicit work rather than a finished claim of phylogenetic conservatism.

Natural next study:

- formal dated phylogeny;
- measurement model jointly treating ambiguous colour observations;
- phylogenetic comparative model of latent species-level polymorphism.

Do not add this analysis to the current paper unless a reviewer requires it and the necessary phylogenetic source can be frozen prospectively.

## Discussion 6 — connection to Camellia without merging papers

One paragraph only.

The Camellia work addresses how similar visible flower-colour outcomes can be generated repeatedly through partially reusable molecular modules. The present global analysis addresses a later level: how much colour variation is carried within species and how that variation is organized geographically. Together they motivate separating **generation, persistence, and spatial organization** as distinct levels of flower-colour evolution.

Do not imply that the Camellia mechanism causes the global FCP pattern.

## Discussion 7 — limitations

Keep four limitations explicit:

1. photo-derived frames are not globally random species samples;
2. D is a four-state representation, not continuous spectral colour;
3. ambiguous-palette observations remain biologically unresolved;
4. observational spatial structure does not identify selection or adaptation.

Then point out that the core spatial result survives the implemented controls despite these limitations.

---

# Figure plan

## Figure 1 — What is the estimand?

**Goal:** make the paper understandable in 20 seconds.

Panels:

A. workflow: 100 fixed photos/species -> four admitted morphs + missingness classes -> species D;
B. distribution of D across 369 discovery species;
C. examples/schematic of low-D versus high-D species;
D. conceptual distinction: amount of polymorphism vs colour direction vs spatial organization.

Do not start with a world map.

## Figure 2 — Flagship discovery + reserve spatial result

**Goal:** the main claim.

Panels:

A. discovery D vs species spatial rho with nonparametric trend;
B. reserve D vs primary spatial rho;
C. >=10% versus <10% subset difference in discovery and reserve with null intervals;
D. reserve primary / observer / quarter / matched-background effect summary.

Use equal-species visual weighting.

## Figure 3 — Robustness to opportunity and measurement ambiguity

**Goal:** answer the reviewer before they ask.

Panels:

A. measurement-status partition: classifiable / ROI-flip failure / ambiguous / no-biological-palette;
B. adjusted partial-rho estimates: raw, span, span+technical failure;
C. D_min4 / observed D / D_max4 interval schematic;
D. reserve primary and matched-background partial rhos for observed, D_min4, D_max4 with spatial-null reference.

This may be the most persuasive figure after Figure 2.

## Figure 4 — Secondary taxonomic structure and the missing shared boundary

Possible main or supplementary figure depending journal.

Panels:

A. discovery and reserve genus clustering gains;
B. 23 shared repeated-genus mean D concordance;
C. D_min4/D_max4 sensitivity showing genus is less robust than spatial result;
D. conceptual contrast: species-specific geographic organization versus unsupported universal boundary.

If the paper needs only three main figures, move A-C to Supplement and keep only a concise conceptual panel in the Discussion figure or remove Figure 4 entirely.

---

# Supplement architecture

## Supplement S1 — sampling and eligibility ledger

- species discovery/capacity criteria;
- 100-photo fixed denominator;
- classifiable threshold;
- discovery/reserve disjointness.

## Supplement S2 — four-state colour classification QA

- ROI/flip gate;
- mixed/ambiguous status definitions;
- palette nuisance controls;
- background extraction.

## Supplement S3 — rejected recurrent-direction analysis

- M1/M2 loadings;
- total-variance control;
- random orientation test;
- M3 auxiliary result.

## Supplement S4 — shared-boundary analyses

- original RGFCA/sharedness tests;
- qualification failures;
- localization limitation.

## Supplement S5 — species-attribute family

- span;
- latitude;
- family;
- pollination gate;
- life-form gate;
- Holm family.

## Supplement S6 — genus robustness

- raw / D_unbiased;
- span-only;
- n_classifiable-only;
- technical-failure-only decomposition;
- ambiguity endpoints.

## Supplement S7 — ambiguity completion mathematics

- proof/argument for D_min4;
- convex integer equalization for D_max4;
- species interval summaries;
- explicit statement that uniform endpoints do not exhaust arbitrary species-specific assignments.

---

# Evidence hierarchy for wording

## Strongest wording allowed

**Replicated and robustness-supported:**

> Greater flower-colour polymorphism is associated with stronger within-species geographic organization in two species-disjoint fixed-photo frames.

**Flower-specific robustness-supported:**

> In the reserve frame, this association remains when matched background spatial structure is subtracted and after sampled-span, technical-failure, and ambiguity-endpoint sensitivity analyses.

## Moderate wording allowed

**Secondary taxonomic result:**

> Species-level polymorphism shows genus-level taxonomic clustering in discovery and reserve frames, although the reserve clustering strength is sensitive to the ambiguity-completion endpoint.

## Wording not allowed

Do not write:

- `polymorphism is adaptive`;
- `range size causes polymorphism`;
- `phylogenetic signal is significant`;
- `pollinators drive global polymorphism`;
- `flower-colour boundaries occur in the same environments worldwide`;
- `ambiguous rows are intermediate morphs`;
- `the four-state D interval fully identifies latent polymorphism`.

---

# Current stopping rule

The empirical spine is sufficient for v0.1 manuscript assembly.

Do not add another broad predictor search before writing the manuscript. New analysis is justified only if it addresses one of these named remaining vulnerabilities:

1. a formal reviewer request for phylogenetic control;
2. an identified mathematical flaw in the ambiguity endpoint construction;
3. failure to reproduce the frozen discovery/reserve spatial statistics from the paper figure pipeline;
4. a mismatch between manuscript numbers and frozen result receipts.

The next technical task should therefore be **figure-data assembly and manuscript-number synchronization**, not hypothesis expansion.
