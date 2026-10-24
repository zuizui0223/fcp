# FCP: do photographic flower colours form biological clusters or measurement artefacts?

**Status: 23 October 2026 — source-locked 4,475-species ledger audit AND the stricter original-447,500-photo palette/GMM and matched-geography quality analysis are complete, with successful CI and a frozen result receipt.** Original source measurements, photographic 50-km matched-null P2, and previous local co-occurrence nulls are immutable.

## 1. Three distinct meanings of cluster — do not confuse them

| Question | What the actual present source establishes | What it does not establish |
|---|---|---|
| **Spatial photographic cluster:** conspecific photographs of the same coarse colour occur nearer each other than an equal-colour-composition, fixed-location label shuffle | Primary P2 at ≤50 km: 3,600 species with any eligible pair; species-equal global-minus-local mismatch difference **+0.0161159**, 199 matched worlds upper P **0.005**. Original strict n ≥40 classifiable and ≥30 local pairs gives 1,761 species, difference **+0.0176149**, P **0.005** | Genetic population clusters, locally exclusive inherited morphs, selection, photo camera/flower-age bias removed |
| **Observed coarse photo-colour labels:** how many of the four broad categories occur at least once in each species' photographic observations | All 4,475 source taxa with original 100 photos each and 227,107 classifiable photos; see exact counts below | Independent stable biological morph count. The nine-palette Hellinger *two-means* W statistic uses an **assumed/fitted 2-group partition**; this **cannot test whether the true number of biological morphs is two** |
| **Independent multicolour modes:** whether measured nine-palette photographs actually show stable separated visual clusters rather than a continuum, image/flower-age changes or part-specific bicolour | **Not yet established** by P1/P2 or the one-row-per-species source ledger. Needs original per-photo nine-palette source, alternative 1–4-component fits, independent-observer replication and botanical organ/age validation | Genetic variants and evolutionarily maintained within-population polymorphism, even after a positive photographic mixture test |

## 2. Exact species-ledger quality and photo-colour richness counts

These counts are based on the **Git-tracked full 4,475 × 114 original-source species CSV**; it has original measured-photo SHA256 `8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6` and species-ledger CSV SHA256 `23d3daadc9799c2e7b5867d54ccc1ed98369ba79a81b72dc7f831e5197e63c49`. They were recomputed from the four actual `miss_coarse_photo_*` counters in each species, not inferred from labels in a paper.

| Coarse photographed colours per source species | ≥1 classified photo in each colour | ≥40 classified photos in the species **and ≥5 classified photos in every detected colour** | Original strict 50-km species (≥40 classified, ≥30 local pairs) | Strict 50-km species after ≥5-in-every-detected-colour gate |
|---|---:|---:|---:|---:|
| 1 | **286** | **268** | 159 | 159 |
| 2 | **1,743** | **488** | 763 | 253 |
| 3 | **1,632** | **134** | 588 | 64 |
| 4 | **814** | **11** | 251 | **3** |
| **Total** | **4,475** | **901** | **1,761** | **479** |

The 901-species total in the third column includes **268 one-photo-colour species**. Requiring **at least two** colours with ≥5 photos each (regardless of other singleton labels), plus ≥40 classifiable photos, retains **1,592** species. Three colours each ≥5 photos: **309** species; all four each ≥5: **11** species. At ≥10 classified photos per colour, at least two colours: **921** species; three: **60**; all four: **zero**.

The ≤50-km strict eligibility among those with ≥2 colours each ≥5 is **813** species; their **observed** species-equal global-minus-local mismatch difference is **+0.02855**, with 61.9% of species positive. This is a **descriptive post hoc subset statistic**, **NOT a matched-permutation-significant result** until the null is recalculated on exactly this subset. It is also distinct from the “all four observed colours with ≥5 photos in each group” stratum (only **3** strict-geography taxa, too small to infer geographic clustering of four stable morphs).

Critically: simply seeing four labels in 814 taxa **does not** show 814 four-morph polymorphic species. Even a five-photo minority support rule has **no botanical validation**; independent observers per *colour*, flower part, anthesis, horticultural/source identity and colour quality remain unconfirmed by a species aggregate alone. Conversely, failed support **does not prove biological monomorphism** because roughly half of the original photo terminal records had no classifiable coarse colour label.

## 3. Already completed strict local-source falsification

Existing [verified original source analysis](../results/fcp_fivecohort_recurrent_local_morph_candidates_20261018.json) required same calendar month in two years, four different observation IDs and observers, within <=1 km complete-link diameter: **8 species** with repeatedly contrasting coarse photographic colour labels. Only **one** retained a four-photo minimum coarse assigned-colour-group fraction ≥65%. With a <=5 km complete-link radius, **17** candidates appeared but **three** met the same exploratory purity condition. In the original eight 1-km candidates, **none** had all four pictures with ≤20% other coarse hue group fractions. A colour-class switch can simply be a two-coloured corolla, flower age or background/ROI contamination.

The distinct 5-km, same-year/month, different-observer *matched label permutation* test was **not supported** (194 eligible taxa with ≥5 pairs, P **0.300**; 14 taxa with ≥30, P **0.275**). This is a separate, more local/source-constrained analysis and must not be overwritten by the broader 50-km positive P2.

## 4. New verification runs — separate ecological and visual estimands

### Stage A: full 4,475-species, source-verified ledger QA

Reproducible script: `scripts/analysis/fcp_observed_colour_cluster_qa_20261023.py`, regression test `tests/test_fcp_observed_colour_cluster_qa_20261023.py`, CI `.github/workflows/fcp-observed-colour-cluster-qa-20261023.yml`. Directly checks the original source ledger SHA256, the photo-label sums and all reported photo-richness thresholds, without accessing private photo ID/location data.

### Stage B: original photo-palette modes and geographic-cluster robustness

Exploratory script: `scripts/analysis/fcp_photo_mode_geo_quality_20261023.py` with genuine source download / SHA verification in `.github/workflows/fcp-original-photo-mode-and-purity-geography-20261023.yml`.

1. Read exact frozen 447,500 original-photo measurements; retain the 227,107 legitimately classifiable observations. Calculate actual coarse-assigned fraction from **nine normalised ROI flower-palette fractions**, excluding “bronze” from the four coarse display hues.
2. For taxa with ≥40 original classifiable photos, apply **post hoc diagnostic** flower-palette group purity ≥0.65 (0.75 sensitivity); require ≥24 retained pure photos, at least two colour labels with ≥5 purity-passing photos and **≥3 distinct observers per such label**.
3. **Visual-mode pilot, not a natural-population census:** Hellinger-transformed nine-palette vectors → PCA2 → Gaussian mixture K = 1–4, selecting model by BIC. Up to 80 hash-ordered taxa per photo-label richness stratum (2/3/4) are evaluated. A photograph-mode candidate must have ΔBIC ≥10 against K=1, every selected cluster ≥5 photos, and fixed-K, independent-observer-half prediction ARI ≥0.75. These arbitrary exploratory gates distinguish reproducible photographic modes from some noise; **even successful candidates may be bicoloured flowers or developmental stages**.
4. **Spatial-sorting purity sensitivity on all eligible taxa:** keep only quality-passing original photos and ≥30 local 50-km pairs also available for **distinct-observer** contrasts. Apply 199 within-species original-label, fixed-location matched null permutations both overall and on observer-disjoint photo pairs; report each supported label-richness stratum with its own null and **do not update original frozen P2**. A significant result can still reflect spatial flower-stage/camera/observer structure, not adaptive morph clusters.

The source-image tarball and original measured CSV were checksum-verified. The [successful original-photo workflow #39880111954](https://github.com/zuizui0223/fcp/actions/runs/39880111954) produced [original result JSON artifact #12469982240](https://github.com/zuizui0223/fcp/actions/runs/39880111954/artifacts/12469982240) and [frozen compact machine receipt](../results/fcp_original_photo_palette_mode_geo_quality_audit_20261023.json). Neither output contains photograph/observation/observer IDs or exact coordinates. **These are post hoc photo-measurement results, not genetic morphs, evolutionary clustering or independently validated flower-stage colour.**

## 5. Actual original-photo cluster-mode and quality-constrained spatial results

**Source:** exactly 447,500 original photo measurement rows, SHA256 `8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6`. Exactly 227,107 four-state classifiable photos. The released tar SHA and source csv SHA both passed the unmodified source gates. The script first measured how much of each assigned coarse colour covered the *normalised nine-palette segmented ROI*, then imposed source-blind fixed diagnostic cutoffs and distinct observer support per label.

| Quality/photo support | ≥2 photo colours, each ≥5 photos and ≥3 observers | ≥3 colours | ≥4 colours |
|---|---:|---:|---:|
| ROI assigned-colour-group share ≥65%, ≥24 quality photos from species initially having ≥40 classifiable | **895 taxa** | **81** | **0** |
| Same, stricter 75% fraction | **445 taxa** | **26** | **0** |

For the **65%** support screen, **306 species** also had ≥30 within-50-km pairs from purity-passing photos with ≥30 pairs contributed by **different observers**. Their mean species-equal global-minus-local photo-label mismatch was **+0.046149**, with an exploratory **199-world fixed-photo-geography, fixed-photo-label-count matched upper P = 0.005**. On distinct-observer pairs, mean **+0.045571**, P **0.005**. Among those 306 species, **279** had exactly two supported coarse hue categories and **27** had exactly three; **zero** supported all four. Separate subgroup conditional 199-world tests gave P = 0.005 for the two-label and three-label strata but **were selected post hoc and have no familywise multiplicity correction**. The stricter 75% cohort's *spatial* null was **not** run: the 75% sensitivity here concerns the *number of supported photographic colour labels* only.

These results strengthen the conclusion that broadly **geographic photo-colour clustering is not entirely driven by singleton minority photo labels or same-observer local pairs**. They **do not** exclude different cameras at different sites, local flower-organ identity, flower-age/phenology, horticultural records or environmentally induced pigment variation. The negative strict ≤5-km same-year/month matched-null result remains intact and limits inference about true local morph clusters.

**Photographic mode structure (exploratory balanced pilot only):** among **160 purposely selected** source species that already showed multiple purity- and observer-supported colour categories, Hellinger/PCA2 Gaussian-mixture BIC selected 2 components in **53** species, 3 in **87**, and 4 in **20**; **95** met the additional fixed-K independent-observer-half ARI ≥0.75, minimum cluster size five, ΔBIC ≥10 criteria. This is a *photo-palette distribution* diagnostic, with sample selection dependent on the existence of at least two observed hue labels. Therefore the lack of best-K=1 selections is not evidence that all flowering plants are multimodal, and even a candidate with best K=4 could simply contain differently lit petals, a multicoloured flower, different floral organs or flower ages. This is **not an independently validated number of biological morphs** and is not a population-wide prevalence estimate.

**Interpretation:** There is robust *photographic geographic structure* at the source-wide 50-km scale and under one strong purity/observer screen, plus exploratory evidence that some photographs form repeatable numerical colour-space modes. The actual number of *biologically discrete, heritable flower-colour morphs* across 4,475 species remains **unresolved** and must not be equated to one, two, three or four by the photo classifier or the original forced Hellinger two-means model.

## 6. What would count as convincing biological clusters?

The decisive next gate is **blinded, photo-by-photo botanical review of the same flower organ and developmental stage**, followed by ≥2 seasons of individually marked field plants, inherited morph frequencies, and independent neutral-genetic versus colour-morph fitness and environment contrasts. A positive visual GMM, 50-km clustering or D-richness relationship is never sufficient by itself.
