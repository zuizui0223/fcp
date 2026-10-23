# FCP / RGFCA: genuine 12-species pixel calibration, repeatability denominators, and photo-reuse identifiability

**Date:** 2026-10-16 JST  
**Status:** two real-photo technical waves complete; **NOT** the main 1,000-species dated-source biological analysis.

## Result now supported by actual original ROI-v4 / EfficientSAM pixels

The RGFCA outcome-blind species roster contains 1,000 sampled taxa selected from the independently verified 2026-09-27 iNaturalist photo-opportunity frame of 42,111 taxa and 93,209,351 eligible SOURCE PHOTO IDs. The 1,000-species species-sampling ledger uses R=200 × 250, exactly 50 appearances per species. It is source metadata only, not an old literature-defined polymorphism list.

We ran **two genuinely nonoverlapping six-taxon technical photo-ID pilots**, each selected independently of its flower-colour results from fixed source-only ranks 1–2 and 3–4 within each of the three source-photo opportunity strata. Every photo was newly evaluated using the same byte-pinned original ROI-v4 + EfficientSAM method, and technical acquisition failures were not silently replaced.

| Source photo stratum | Technical taxa | Real ROI photos | Four-state classifiable ROIs | Classifiable fraction |
|:--|--:|--:|--:|--:|
| 128–255 | 4 | 332 | 110 | 0.3313 |
| 256–1,023 | 4 | 318 | 93 | 0.2925 |
| ≥1,024 | 4 | 384 | 202 | 0.5260 |
| **All** | **12** | **1,034** | **405** | **0.3917** |

The source-only species-equal mean of individual photo-classifiability fractions is **0.3756**, as distinct from the photo-weighted 405/1,034=0.3917. The difference in tier rates is descriptive, with **four species per tier**, not an ecological effect or a demonstrated universal relationship between source photo availability and flower-colour variation.

**Two original pixel sources / exact independent analyses**:
- [Six-taxon real wave 1, GitHub #38832829706](https://github.com/zuizui0223/fcp/actions/runs/38832829706), artifact 12018829689: 552 actual original photos, 202 classifiable.
- [Disjoint six-taxon real wave 2, GitHub #38840579008](https://github.com/zuizui0223/fcp/actions/runs/38840579008), artifact 12022221110: 482 actual original photos, 203 classifiable.
- [Wave 1 taxon/observer/geography/within-ROI audit #38841451431](https://github.com/zuizui0223/fcp/actions/runs/38841451431), artifact 12020598724.
- [Wave 2 re-audit and 12-taxon source-photo comparison #38844397987](https://github.com/zuizui0223/fcp/actions/runs/38844397987), artifact 12022132603; all source identifiers and receipt hashes rechecked.

### Technical evaluability is not the same as 40 good photos

Five of the twelve taxa had **≥40 classifiable** images across both photographer-separated source groups. But only **three** reached ≥30 evaluable Monte Carlo repeats at 48+48 raw-photo depth, because the other two lacked ≥20 classifiable images *per fold*. Those three were:

- *Matricaria chamomilla*, 64/96 classifiable, technically evaluable; two coarse photo labels repeatedly detected at depth 96.
- *Nassauvia axillaris*, 43/96 classifiable, technically evaluable, but no repeated two-label support (42 white / 1 yellow-orange).
- *Thelesperma filifolium*, 67/96 classifiable, technically evaluable, but no repeated two-label support (66 yellow-orange / 1 white).

Thus in the combined 12-taxon technical pilot, **1** taxon had a repeated two-coarse-photo-label signal. That one taxon was **held** for floral-organ-structure QC, rather than admitted as an FCP species. All other nine lacked adequate repeated photo support or showed no repeated two-label signal. **None of the twelve was certified genetically polymorphic, biologically monomorphic, or as carrying natural-population sympatric colour morphs.**

*Matricaria* provides the critical negative-control-like technical observation: **58/64** classifiable original photos contained **both white and yellow components >=10% in the SAME flower ROI**, consistent with the normal ray/disk organ architecture of a composite inflorescence. At 40+40 raw photos, the two-photo-label empirical repeatability was 0.83, and at 48+48 was 1.0. Those are **photo-label** repeatabilities, not separate individual morphs. The technical images were discarded after pixel measurement, so a human floral-organ check of original image crops would still be necessary to certify the causal source of colour labels; the intraflower palette evidence alone does not prove every photographed organ.

### Critical resampling overlap problem

A Monte Carlo photo re-draw only perturbs the sample meaningfully when the pool exceeds the draw. Sampling k images without replacement from an available pool of N images gives expected intersection **k²/N** photos between two independent draws, and expected share of the k photos reused in both draws is **k/N**. This is a property of the *same fixed photo pool*, not biological independence.

| Independent source-photo pool | Photos per draw | Expected shared photo fraction |
|---|---:|---:|
| Main proposed 64+64 = 128 | 20+20 = 40 | 31.25% |
| Main proposed 64+64 = 128 | 40+40 = 80 | 62.5% |
| Main proposed 64+64 = 128 | 60+60 = 120 | **93.75%** |
| Main proposed 64+64 = 128 | 64+64 = 128 | **100%** |
| Real six-taxon pilot 48+48 = 96 | 48+48 = 96 | **100%** |

These numbers assume one distinct photo per selected observation and full balanced source pools. Where photographer/source-availability constraints reduce effective pools, actual overlap may be even higher. An apparently convincing **P=1 over 200 draws at full-pool depth is mathematically unsurprising and supplies no evidence about unmeasured flower individuals**.

The implementation has now been amended to **compute actual photo-ID overlap across draws for each species and all depth sensitivities** by counting how many repeated selections contained each original photo ID. It also writes how many distinct source-ID sets were actually drawn. This is more trustworthy than pretending 200 fixed-pool re-draws are independent.

### Three denominators that must be shown together

For taxon i, let R_i denote the 50 scheduled balanced species inclusions, E_i the number of those repeated photo draws with sufficient classifiable ROIs and independent-observer/observation coverage, and D_i the number of evaluable draws showing independently replicated *photo-label* colour states.

1. **P_eval(i) = E_i / R_i.** Probability the fixed source-photo draw is technically informative at the frozen ROI criteria.
2. **P_photo-colour|eval(i) = D_i / E_i.** Conditional repeated two-photo-label detection, **undefined** if E_i=0.
3. **P_joint(i) = D_i / R_i = P_eval × P_photo-colour|eval.** Repeated original-photo study draw is both technically informative and photo-label positive.

**The main intended FCP positive call additionally needs** within-flower organ QC, independent source observations/observers, spatial/population-aware verification and photo-ID-overlap-aware sensitivity to draw depth. P_joint=0 when E=0 **does not establish a monomorphic plant**. Neither P_joint nor P_photo-colour|eval is a genetic polymorphism probability.

## Exact remaining scientific dependency

All 1,034 measured pixels in the two six-taxon technical waves came from the *live iNaturalist observation API* with a created-before 2026-09-27 filter. Taxon identities were drawn from the dated 2026-09-27 original source roster, **but exact sampled PHOTO IDs were NOT reconciled to the original immutable 36-GB monthly 2026-09-27 archive**. Consequently neither pilot's colour labels are admissible as the full dated-source, 1,000-taxon formal RGFCA P_poly result. The 36-GB archive and 93-million-photo-ID SQLite are not saved in the GitHub-hosted ephemeral source-counting scratch. Further *live-API technical pilots are not a substitute* for that immutable provenance unless a separately declared prospective source cohort and a matching sampling frame are adopted.

**Decision:** move the 1,000-species real-photo screen only when source photo-ID admission and a systematic photo-ROI acquisition protocol are available; use measured P_eval, P_photo-colour|eval, P_joint, empirical photo-ID overlap and floral-organ structural review to decide whether repeated photo colour states represent reliable evidence. Do not expand image counts blindly, treat poor-quality species as colour-negative, or claim that a repeated shared global RGFCA hotspot was supported (historical G1 p=0.070).


## Audited actual repeat-photo-ID sets, and modified main RGFCA criterion

The earlier table gives expected overlaps under complete uniform source-photo sampling. The **new independently replayed real original-photo source-ID audit**, not a theoretical calculation, is [GitHub Actions #38845056079](https://github.com/zuizui0223/fcp/actions/runs/38845056079), successful; artifact \`12023026437\`. It re-read both original 552 + 482 ROI photo-ID manifests, retained original observer folds, and replayed 200 salted source-ID draws per technical taxon without downloading a single new image.

**Actual mean fraction of source photo IDs shared by pairs of re-draws**, averaged equally over 12 technical taxa:

| Photos selected per observer split | Actual mean photo-ID reuse |
|---|---:|
| 20 | 0.4900 |
| 32 | 0.7343 |
| 40 | 0.8774 |
| 44 | 0.9415 |
| 48 | **1.0000** |

At **48 photos per split**, **all 12 taxa** had exactly **one distinct photo-ID set in all 200 repeated Monte Carlo realizations**. This includes the apparently striking Matricaria repeatability=1.0; it did **not** constitute 200 independent original-photo samples. Some taxa with low observer availability already had full-pool repetition at shallower depths. Real source observer constraints explain why mean overlap was higher than the complete-balanced 48-per-fold theoretical values at earlier depths.

### Main 1,000-species RGFCA rule corrected AFTER the two technical pilot waves, BEFORE main original-photo measurements

The 1,000 original species and 200 × 250 exact-balanced *species* inclusion schedule are preserved unchanged. But with **64 source-photo IDs per observer group**, the old primary **60+60 photo** re-draw shares ~93.75% of photo IDs, in expectation, between two draws. It is a high-depth photo-**detectability** readout rather than meaningful independent photo-sampling robustness.

The strengthened, *pre-main-image-measurement, post-technical-pilot* method retains **60+60** as a high-power detection screen but adds **48+48** source-photo re-draws for mandatory robustness qualification:

- Minimum **30 evaluable** lower-depth (48+48) species realizations.
- At least **80%** of those lower-depth evaluable realizations reproduce the two photo labels in both observer-disjoint halves.
- Actual measured pairwise source photo-ID overlap in the 48+48 draws must be **≤0.80**, and there must be **≥30 distinct selected photo-ID sets**. In a full 64+64 pool, the theoretical expected overlap fraction is 48/64=0.75; in a smaller effective pool, the recorded overlap determines whether a repeated signal was actually perturbed.
- Simultaneous **P(photo labels | evaluable) ≥0.80** in the 60+60 high-power draws, and adequate evaluable 60+60 support.
- Photo-label positives dominated by white/yellow in the SAME single flower ROI remain on the **structural-organ manual review hold**.

If the lower-depth protocol cannot be evaluated because there are too few classifiable flower images, the species is **\`repeatability unresolved\`**, not biologically monomorphic. The above fractions are repeated *photo-label* stability, not heritable genetic polymorphism probabilities. The 0.80 overlap cutoff and dual-depth confirmation are **posttechnical calibration choices**: do not retroactively claim they were fixed before the first six-taxon pilot. A subsequent independent cohort would be needed to establish the reliability of these gates.

The implementation reports three separate denominator-protected quantities **P_eval**, **P_photo-colour|eval**, and **P_joint** for the primary 60+60 source-photo draw, along with actual photo-ID overlap and independent photo-set counts for depths **40/80/96/120/128**. A photo ID is never treated as a new biological replicate simply because it reappears in another RGFCA realization.

**Production remains gated on truly source-admitted photograph IDs** and newly measured original ROI-v4 outcomes for the prospective 1,000 source-sampled species. No original 2026-09-27 dated-photo 1,000-species formal \`P_poly\` has been estimated from these 12 technical taxa.
