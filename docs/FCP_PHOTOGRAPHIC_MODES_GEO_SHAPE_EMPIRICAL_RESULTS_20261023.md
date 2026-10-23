# FCP: photographic modes and geographic structure — original-data findings

**Completed and authenticated: 23 October 2026.** This is a new exploratory analysis of original one-provider photographs, not genetic flower-colour morph confirmation. Its source population is 447,500 photographs, 4,475 species, including 227,107 four-classifiable photographs; the strict analyzable geographical/phenotypic sample is **1,761 species**. Exact historical P2 strict mean mismatch +0.0176149021685411 was replayed.

Primary authenticated execution: https://github.com/zuizui0223/fcp/actions/runs/39923312111
Independent joining and audit of two immutable prior original-source result artifacts: https://github.com/zuizui0223/fcp/actions/runs/39924481803
Machine result and provenance: `results/fcp_photo_modes_geo_shape_verified_20261023.json`, `results/fcp_photo_modes_geo_shape_independent_audit_20261023.json`.

## Result A: 4 hard labels are NOT counts of distinct phenotypic density modes

On the dominant Hellinger nine-hue photograph phenotype axis, Gaussian mixture K=1–4 BIC fits were compared against both a unimodal normal and a skew-normal, then screened for actual density peaks, a minimum photo/observer support and repeated observer subsets.

| Evidence status | Taxa (among eligible 1,761) | Interpretation |
|---|---:|---|
| Recurrent photograph density multimodality candidates | **27** | Evidence of at least two recurring photographic peaks on first axis |
| Unstable photographic multimodality candidates | **43** | Density peaks do not replicate reliably on observer subsets |
| No robust multimodality identified by this procedure | **1,691** | Includes truly unimodal, continuous, low statistical power, and PC2/PC3-unresolved cases |

Among the 27 recurring candidates, 23 were fitted by a K=2 mixture and four by K=3, **but all 27 fitted densities displayed only two visible peaks** under the conservative peak definition. *Three fitted components do not imply three natural colour types.*

Species-level observed 2/3/4 coarse photo categories are **not** evidence for corresponding inherited morph K. The nine-palette and camera masks cannot establish chemically calibrated petal reflectance or individual stability.

## Result B: photographed geographical organization does not generally extrapolate as a broad linear cline

Of the 1,761 species, 1,566 allowed observer-disjoint validation and only **449** had sufficient geographic spread, spatially separate blocks, and holdout support for three-region evaluation. Within these 449, 8 satisfied an exploratory transferable broad-cline threshold, one had observer-predictive local nonlinear excess without an established transferable linear trend, and 440 did not meet either threshold. The remaining 1,312 were spatially underidentified under the three-block rule.

Across 449 taxa, median spatial-block transfer improvement of the broad linear model over the month-only baseline was **−0.12395**, with **10.02% positive**. Across the 1,566 observer-eligible taxa, the median observer-holdout gain was **−0.02092**, and the median local nonlinear increment was **−0.00420**. Negative out-of-block predictive skill does not contradict positive, within-50km nine-hue photographic resemblance: the objectives and extrapolation conditions differ.

These signatures are NOT a census of real geographical clines or mosaics. A residual nonlinear pattern can be a curved gradient, a patch, a camera-location artifact or phenology. Sparse long-distance data and restrictive block feasibility make failure to identify spatial shape especially uninformative.

## Result C: photographic modes and geographic sorting are different measurements, but may be related

The 27 recurrent photographic density-mode candidate taxa had original four-colour 50km photo mismatch depletion mean **+0.08375**, versus **+0.01657** across 1,691 without detected photographic multimodality. This is an unadjusted comparison and shares source photographs.

The separate continuous nine-hue **calendar-month-matched** geography estimator had a difference **+0.04947** (mode candidates minus non-detected taxa) before covariate correction. After controlling species-wide nine-hue JSD diversity, classified photo count, and near-pair opportunity, the difference was **+0.03263**, **genus-bootstrap 95% interval [−0.00038,+0.06684]** (499 replicates over 821 genera). Additional adjustment for variance captured by first hue PC and raw category richness gave **+0.02841 [−0.00426,+0.06276]**. Forty-three unstable candidates were excluded from this contrast.

**Decision:** an intriguing **positive effect estimate**, but not a statistically well-separated positive independent rule after continuous diversity and photo-opportunity correction. Do NOT state 'multimodal species are proven to segregate more', 'modes are genetically distinct', or 'mosaic explains the remaining taxa'.

## Immediate botanical validation target

Use the 27 recurrent photographic mode candidates and a size-/diversity-/geography-matched set without detected modes for blinded original-image review. Within each candidate taxon, sample photos from both peaks and intermediate densities. Verify flower organ, developmental stage, native/wild status, taxonomic identity, camera lighting and observer independence before biological discrete-morph terminology is introduced.

The independent audit preserves a candidate list (species names and diagnostic summary, **no sensitive precise coordinates or observer IDs**) as the GitHub Actions artifact in run #39924481803. It also reports nine named species with exceptional predictive spatial signatures, including `Delphinium hesperium` and `Xyris caroliniana` with both a recurrent photographic density mode and a broad transferable trend, and `Lomatium martindalei` with recurrent photographic density mode plus local nonlinear predictive signal. They are candidates for inspection, not field-validated conclusions.

### Limits and claim ceiling

1. Entire inference is post hoc and one provider; 220,393 original photos have no eligible four-state label and remain missing.
2. The most stringent earlier <=5km same-year/month different-observer colour null did **not** support added local structure (P=0.300/0.275). It is preserved.
3. One-dimensional mixture densities can miss complex hue configurations; fitted mixture components need not be true peaks.
4. Mode counts, broad trend counts and nonlinearity counts are **not** population-level frequencies of genetic FCP. Neither selection nor pigment costs/benefits were measured.
5. No common macroecological environmental predictor or genealogical mechanism was established by these models. Next causal step requires validated independent flower-colour morphs and matched geography/flowering stages, plus neutral genetic structure where possible.

