# Systematic search protocol for natural floral-colour ITV and FCP

**Protocol version:** 1.0.0  
**Evidence-synthesis type:** systematic map with a nested comparative species dataset  
**Reporting framework:** ROSES systematic-map protocol/report checklist and flow diagram; PRISMA 2020 flow terminology is used where compatible.

## 1. Review objective

Identify, deduplicate and classify published evidence for naturally occurring intraspecific variation in floral display colour across angiosperms, including:

1. discrete flower-colour polymorphism within natural populations;
2. continuous within-population flower-colour variation;
3. geographically structured flower-colour variation among natural populations;
4. systems showing evidence at both within- and among-population scales.

The search is intended to support a reproducible global comparative dataset. It is not designed to estimate the prevalence of colour variation among all angiosperms unless a separate probability-based flora sampling design is added.

## 2. Terminology and conceptual hierarchy

### 2.1 Intraspecific floral-colour variation

**Intraspecific floral-colour variation (IFCV)** is the umbrella category. It includes stable differences in floral display colour among conspecific individuals, genotypes or natural populations. Variation may be discrete or continuous and may occur within populations, among populations or at both scales.

### 2.2 Intraspecific trait variation

**Intraspecific trait variation (ITV)** is the broader ecological concept covering variation in any trait within a species. In this review, ITV is restricted to floral display colour. A publication does not need to use the term `ITV` to be eligible.

### 2.3 Strict flower-colour polymorphism

**Flower-colour polymorphism (FCP)** is coded strictly as the coexistence of at least two discrete floral-display colour morphs in the same natural population. A rare aberrant individual is not automatically treated as a stable polymorphism. The source must establish coexistence or report morph frequencies from at least one population.

### 2.4 Geographic floral-colour variation

**Geographically structured variation** means that floral colour or colour-morph frequencies differ among natural populations, regions or positions along a geographic or environmental cline. Within-population coexistence may be absent, unreported or present in only some populations.

## 3. Evidence units

The retrieval unit is a publication record. The biological evidence unit is a species × study combination. Multiple publications on one species are retained and linked because they may provide different evidence for natural status, variation form and spatial scale.

Reviews are retained as discovery seeds and citation sources but are not treated as primary biological observations unless they report original data that meet the eligibility criteria.

## 4. Information sources

The primary reproducible search uses:

- OpenAlex title, abstract and available full-text search;
- Crossref bibliographic metadata search;
- Europe PMC title/abstract and life-science literature search;
- one-generation backward and forward citation chasing in OpenAlex from prespecified seed reviews.

Optional supplementary searches may use Web of Science, Scopus, BIOSIS, CAB Abstracts, Google Scholar or regional databases when institutional access permits. Results from non-API interfaces must be exported with the exact search date, query, result count and export format. Google Scholar results are supplementary because ranked results and total counts are not fully reproducible.

No publication-date or language restriction is applied. English, Spanish, Portuguese, French, German, Italian, Japanese and Chinese query blocks are included. Additional languages may be added as versioned amendments.

## 5. Search strategy

The complete versioned query list is stored in `literature/itv_fcp_search_config.json`. Queries combine four concept blocks:

1. **Floral display organ:** flower, floral, petal, corolla, tepal, labellum or showy floral bract;
2. **Colour phenotype:** color/colour, pigmentation, pigment, anthocyanin, carotenoid, hue, chroma or reflectance;
3. **Intraspecific variation:** polymorphism, morph, dimorphism, intraspecific variation, conspecific variation, continuous variation, geographic variation, cline or population differentiation;
4. **Spatial or natural-population context:** within population, among populations, morph frequency, wild, natural, field, population or locality.

Search retrieval is deliberately broader than final eligibility. Terms such as `natural population` are not mandatory at retrieval because relevant abstracts often omit sampling context. Cultivated, induced and ontogenetic terms are used for prioritisation, not irreversible automatic exclusion.

Special query blocks target commonly missed literature:

- white-flowered, alba, albiflora, acyanic and anthocyanin-loss morphs;
- continuous hue, chroma and spectral-reflectance variation;
- orchid labellum colour variation;
- showy bract colour variation;
- geographic and interpopulation variation that is not labelled polymorphism;
- selection-mechanism terminology associated with known FCP systems.

## 6. Seed set and citation chasing

Seed reviews are specified before the search in the configuration file. They cover strict FCP, continuous colour variation, orchids, Mediterranean systems and classic flower-colour evolutionary studies. For each seed:

1. resolve the DOI to an OpenAlex work;
2. retrieve its reference list;
3. retrieve all citing works, subject to a documented technical maximum;
4. deduplicate citation-derived records with database-search records.

If a citation query reaches the technical maximum, the query is flagged as truncated and must be partitioned or rerun before the search is declared final.

## 7. Deduplication

Records are merged first by normalized DOI. Records lacking a DOI are merged by normalized title and publication year. The preferred retained metadata record is the one with the longest abstract, followed by citation count. All contributing databases, query IDs and source record IDs are preserved.

## 8. Automated prioritisation

Automated screening produces a review priority only. It does not make final inclusion or biological classification decisions.

Priority signals include:

- floral-display colour terminology;
- explicit variation, morph, polymorphism or population differentiation terminology;
- natural, field or population context;
- within- or among-population wording;
- candidate angiosperm binomials matched to the repository census.

Potential exclusions are flagged for cultivated-only material, induced mutation, transgenic or tissue-culture studies, ontogenetic colour change and colour of non-display organs. Records containing both inclusion and exclusion signals remain in human review.

## 9. Human screening stages

### Stage 1: title and abstract

Include for full-text assessment when the publication may describe variation in floral display colour among conspecific individuals or populations. Do not require proof of natural status at this stage if the abstract is incomplete.

### Stage 2: full text

A study is eligible for the natural-ITV evidence map when all conditions are met:

1. the focal entity is an angiosperm species or an unambiguously identifiable infraspecific taxon;
2. floral display colour varies among conspecific individuals, genotypes or populations;
3. the variation is documented in wild, natural or naturalized populations, or in a common-garden/experimental study using wild-origin material for which natural variation is independently documented;
4. colours are compared at equivalent developmental stages;
5. the phenotype concerns petals, corolla, tepals, labellum or a showy floral bract.

### Main exclusions

Exclude from the main natural-ITV dataset when evidence is limited to:

- cultivars, breeding lines, commercial varieties or horticultural collections with no evidence from natural populations;
- induced mutation, transgenic manipulation, gene editing or tissue culture;
- colour changes within the same flower during development, pollination or senescence;
- environmentally induced colour changes in one genotype without evidence of naturally occurring among-individual or among-population variation;
- interspecific colour differences only;
- fruit, leaf, stem or seed colour;
- pollen, anther or stigma colour only.

Non-display floral-organ colour can be retained in a secondary evidence category but must not be silently combined with petal/corolla colour.

## 10. Biological coding

Each eligible species × study record receives the following fields.

### Natural status

- `wild_field`
- `naturalized_field`
- `wild_origin_common_garden`
- `cultivated_with_independent_wild_evidence`
- `cultivated_only_exclude`
- `induced_exclude`
- `unclear`

### Variation form

- `discrete_morphs`
- `continuous`
- `mixed_discrete_continuous`
- `rare_form_only`
- `unspecified`

### Spatial scale

- `within_population`
- `among_population`
- `mixed_within_among`
- `unclear`

Strict FCP requires `discrete_morphs` and direct evidence of `within_population` coexistence. Continuous within-population ITV is included in the ITV map but is not labelled strict FCP.

### Display organ

- `petal_corolla_tepal`
- `labellum`
- `showy_bract`
- `multiple_display_organs`
- `non_display_floral_secondary`

### Evidence fields

Record the verbatim supporting passage within copyright limits, page or section, population/locality identifiers, morph-frequency information, DOI or stable source ID, and reviewer notes.

## 11. Reviewer calibration and agreement

Before full screening, two reviewers independently screen the same calibration set of at least 50 records spanning all automated priorities. Eligibility rules are revised before the main screen if disagreements reveal ambiguous criteria.

All provisionally included records and all records used to assign within/among/mixed spatial scale should receive duplicate independent review. Report raw agreement and Cohen's kappa separately for:

- title/abstract eligibility;
- full-text natural-ITV eligibility;
- strict FCP status;
- spatial-scale classification.

Disagreements are resolved by documented consensus or third-reviewer adjudication. Automated labels remain separate from adjudicated labels.

## 12. Search completeness diagnostics

The final report should include:

- ROSES/PRISMA-style flow counts by database and search stage;
- query-level reported and retrieved counts;
- numbers of truncated queries;
- duplicate counts;
- contribution of each database and citation chasing to unique included studies;
- recovery of all prespecified benchmark studies and species;
- saturation across successive search blocks and citation-chasing rounds;
- reasons for full-text exclusion;
- species accumulation curves as searches are added.

Failure to recover a benchmark paper triggers query revision before biological analysis.

## 13. Versioning and reproducibility

Every search run stores:

- exact configuration SHA-256;
- UTC start and completion times;
- database and API versions where exposed;
- raw records;
- query log;
- deduplicated records;
- untouched human-screening fields;
- output file hashes;
- code commit SHA.

Search amendments require a new protocol version and a change log. The biological dataset is frozen only after human adjudication; model results must not be consulted when deciding study eligibility or spatial classification.

## 14. Interpretation boundary

The resulting sample represents documented natural cases recoverable through the declared literature-search and screening process. It is not a random sample of angiosperms and cannot estimate the global prevalence of floral-colour ITV or FCP without a separate denominator design.

## Execution log

- 2026-08-02: triggered the first full sharded retrieval after the pull-request smoke test passed.
- 2026-08-02: retriggered the branch push explicitly to start the 15-query, seven-seed full retrieval; verify the push-event run separately from the PR smoke run.
