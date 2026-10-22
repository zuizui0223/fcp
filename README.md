# FCP — within-species flower-colour variation

**Verified completed five-cohort status (22 October 2026):** FCP has genuinely measured and analysed **447,500 original photographs from 4,475 species**, with **227,107 four-state classifiable photographs**. The unchanged original 50-km species-conditional test had **1,761 eligible species**, mean range-wide-minus-local mismatch **+0.017614902**, P = **0.005** (199 permutations). This supersedes the *current-status wording* of the old 17 October pending ledger; the historical snapshot remains intact. See [current measured-data status and paper-specific evidence routes](docs/FCP_CURRENT_MEASURED_RESULTS_AND_NP_DRAFT_STATUS_20261022.md) and the [permanent finished-data release](https://github.com/zuizui0223/fcp/releases/tag/fcp-all4730-original-method-actual-colour-results-20261017). The **4,475-species New Phytologist Full Paper** directly tests four **post hoc photographic predictions (P1–P4)** in [the active manuscript](docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md); its prior stacked development is retained as [source PR #252](https://github.com/zuizui0223/fcp/pull/252), and the clean integration is documented in [PR #254](https://github.com/zuizui0223/fcp/pull/254). The paper **has not been submitted**.

**Repository boundary — 14 October 2026:** This **fcp** repository is exclusively the current **within-species flower-colour variation** paper: image measurement validity, D, repeated achromatic–chromatic phenotype-space axis W and geographic distribution within sampled plant species. The **global between-species flower-colour × rainfall/vapour/soil study**, including its original macro scripts, 1,280 direction rows and 576 matched spatial coefficients, has been migrated to [**zuizui0223/mois**](https://github.com/zuizui0223/mois). Historical Git branches and FCP photo-frame provenance remain accessible; the superseded integrated one-paper PR #196 should not be merged into FCP main. [Split scope and all-three-cohort within-species axis note](docs/FCP_WITHIN_SPECIES_SCOPE_AFTER_MOIS_20261014.md).

**Seasonal photo-white × observation-month climate analyses now live in [`zuizui0223/mois`](https://github.com/zuizui0223/mois).** The October 2026 source-only date audit, observed-season matched locality white-rate QC, month-specific TerraClimate matching, and the source-only seasonally balanced photo-selection design are documented in MOIS and passed its independent CI (see [MOIS seasonal plan](https://github.com/zuizui0223/mois/blob/main/docs/MOIS_SEASONAL_WHITE_MOISTURE_RESEARCH_20261016.md), [merged transfer PR #2](https://github.com/zuizui0223/mois/pull/2)). **FCP remains the upstream owner of original source-photo ID provenance, immutable iNaturalist photo/observer identity, blind ROI-v4/EfficientSAM image measurements, and the separate within-species D/W/distributed-polymorphism paper.** Historical FCP season code/Actions are retained for SHA provenance and must not be interpreted as a competing active MOIS seasonal manuscript. As of the transfer no empirical seasonal climate–white association or fully flowering-verified 3.43M images have been established.

The current [New Phytologist Full Paper](docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md) uses all **4,475 measured species** as the common source registry (unavailable estimands remain **NA**). The older prospective tests and historical preextension manuscript remain frozen separately, including [the canonical historical manuscript](docs/POLYMORPHISM_MANUSCRIPT.md). This is a working paper, **not a submitted or accepted article**.

## Repository split: fcp = within-species variation; mois = global interspecific moisture

**This repository (`fcp`) is the active research home for within-species photographic flower-colour variation.** The historical paper tested observer-disjoint reliability of species-wide diversity `D`, geographically structured local differences and an earlier fixed white–chromatic direction (including an historical 158-species prospective test). The **current hypothesis paper** uses one **4,475-species** source registry to test four observable predictions about phenotype direction, geography, climate and their cross-species relationship. The historical 158-species outcome is *not* the current main denominator. The 42,111-species discovery frame describes sampling opportunity, **not** genetically verified natural-population polymorphism prevalence.

**The separate [`mois` repository](https://github.com/zuizui0223/mois) is the active research home for *between-species* global floral colour × rainfall/atmospheric-moisture analysis.** Its source-linked manuscript, 1,280-row environmental coefficient ledger, 576-row matched spatial-support audit, original analysis code, uncertainty/missingness sensitivity records and global moisture figure were copied from their historical FCP branches. The transfer preserves original Git blobs and scientific negative results; it does not rerun the source-photo models. See [the FCP↔mois split receipt](docs/FCP_MOIS_RESEARCH_SCOPE_SPLIT_20261014.md).

**Seasonal original-photo programme (16 October 2026):** [`mois`](https://github.com/zuizui0223/mois) is also the active home for observation-month white-versus-chromatic flower composition, real monthly temperature/radiation/water conditions, observer/site/year matching and independent open-flower-organ review. [MOIS PR #4](https://github.com/zuizui0223/mois/pull/4) was merged after source and botanical tests passed; MOIS independently reproduced the signed 39-taxon/268-cell-year *photo-opportunity* rank from FCP date metadata and generated an **unreviewed** original-photo anthesis verification sheet. Those outputs are **not** botanical white-flower rates or ecological climate effects. FCP remains authoritative for original licensed photo/observation acquisition, blind ROI-v4/EfficientSAM classification and its own within-species W/D polymorphism manuscript. Historical FCP source and season-side receipts are retained for provenance, not two active seasonal manuscripts.

**Publication boundary:** the integrated global-plus-ITV Draft [FCP PR #196](https://github.com/zuizui0223/fcp/pull/196) is historical/superseded, not the active FCP manuscript; the [active 4,475-species New Phytologist hypothesis manuscript](docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md) is the working paper, while the [canonical historical manuscript](docs/POLYMORPHISM_MANUSCRIPT.md) and [stacked source PR #252](https://github.com/zuizui0223/fcp/pull/252) preserve earlier model/analysis chronology. Prior macro development PRs [#193](https://github.com/zuizui0223/fcp/pull/193) and [#194](https://github.com/zuizui0223/fcp/pull/194) remain source-provenance records and are **not** a second active FCP paper. Neither repository has direct pigment chemistry, genotype-frequency or morph-specific reproductive-fitness evidence.

The current Full Paper directly tests four **observational** predictions:

1. **P1 — recurring colour direction:** white–chromatic W exceeds a matched coarse-state-preserving null across all 4,475 species (observed 0.4789; null median 0.4320; exploratory P = 0.005).
2. **P2 — geographic allocation:** conspecific photographs within 50 km have lower colour mismatch than a fixed-species-composition null predicts (3,600 estimable species; mean +0.0161; P = 0.005).
3. **P3 — climate association:** the photographed white/nonwhite contrast correlates weakly with warmest-month maximum temperature (4,112 estimable species; univariate mean +0.0219 log odds per °C), but varies by cohort and fails strong spatial prediction.
4. **P4 — coupled outcome:** the white–chromatic W direction does not strongly predict how much a species geographically sorts its photo colours or the sign/strength of its temperature association (rank CIs include zero).

All tests on the complete source roster are **post hoc**; they concern photographed phenotypes, not direct genetic polymorphism, locally stabilizing selection, or a universally adaptive white flower. All five species cohorts use the same iNaturalist source and ROI procedure: the expansion is **not an independent-source replication**. The historical 158-species prospective test remains an archived confirmation of an earlier specifically frozen target, **not** a criterion for admitting species into the present P1–P4 study.

## Repository scope

This repository is the paper-facing FCP analysis surface. The generalized trait/interaction **spatiotemporal turnover** programme now lives in [`zuizui0223/turnover`](https://github.com/zuizui0223/turnover) and is intentionally excluded from the active FCP manuscript, workflows, results, tests, and development documents here. Earlier pre-split development remains recoverable from Git history.

## Start here

1. **Manuscript:** [`docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md`](docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md)
2. **4,475-species data and hypothesis status:** [`docs/FCP_CURRENT_MEASURED_RESULTS_AND_NP_DRAFT_STATUS_20261022.md`](docs/FCP_CURRENT_MEASURED_RESULTS_AND_NP_DRAFT_STATUS_20261022.md)
3. **Expanded Methods/Results:** [`docs/supporting/FCP_NP_ALL4475_EXTENDED_METHODS_RESULTS_20261021.md`](docs/supporting/FCP_NP_ALL4475_EXTENDED_METHODS_RESULTS_20261021.md)
4. **Data figures:** [`docs/figures/fcp_np_all4475_white_axis_five_cohorts_20261021.svg`](docs/figures/fcp_np_all4475_white_axis_five_cohorts_20261021.svg), [`docs/figures/fcp_np_all4475_figure2_geo_50km_20261021.svg`](docs/figures/fcp_np_all4475_figure2_geo_50km_20261021.svg) and five additional paper/supplemental figures.
5. **Reproducibility contract:** [`CURRENT_PAPER_REPRODUCIBILITY.md`](CURRENT_PAPER_REPRODUCIBILITY.md)
6. **Frozen claim ledger:** [`docs/POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260918.md`](docs/POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260918.md)
   Historical source and claim-to-input map: [`docs/POLYMORPHISM_DATA_LINEAGE_MAP_20260925.md`](docs/POLYMORPHISM_DATA_LINEAGE_MAP_20260925.md).
7. **Figure plan:** [`docs/POLYMORPHISM_FIGURE_PLAN_20260918.md`](docs/POLYMORPHISM_FIGURE_PLAN_20260918.md)
8. **Secondary mechanism evidence hierarchy:** [`docs/POLYMORPHISM_SECONDARY_MECHANISM_EVIDENCE_LEDGER_20260928.md`](docs/POLYMORPHISM_SECONDARY_MECHANISM_EVIDENCE_LEDGER_20260928.md)
9. **Frozen provenance release receipt:** [`archive/fcp_submission_20260925/NP_PROVENANCE_RELEASE_RECEIPT.md`](archive/fcp_submission_20260925/NP_PROVENANCE_RELEASE_RECEIPT.md)

## Current evidence surface

**Integrated positive ecology (source-locked historical/current synthesis):** The older two species-disjoint high-depth cohorts showed that greater photographed species-wide D correlated with stronger geographic colour organization (discovery adjusted rho **+0.1266**, P **0.007**; validation **+0.0993**, P **0.025**). These are distinct from the 4,475-species P2 local pair-null comparisons and were not refitted in all taxa. The first-frozen observer-disjoint H1 reliability test passed, although a later stricter split failed; historically tested H3a broad phylogenetic signal and H3b photographed sampled-span replication were unsupported in their respective smaller validation samples. See [the H1–H3 / P1–P4 ecological synthesis and source receipts](docs/FCP_NP_H1_H3_P1_P4_ECOLOGICAL_SYNTHESIS_20261022.md). None of these results estimates natural genetic-morph fitness or independent-provider replication.

From the **447,500 actually measured photographs**, **227,107** were classifiable into four coarse flower-colour labels. The evidence supports a recurrent photographed achromatic–chromatic colour direction (P1) and nonrandom geographic clustering of photographed conspecific colour states (P2). A small positive mean BIO5 white-label association (P3) is heterogeneous by species/cohort; the size of recurrent W alignment weakly predicts either spatial sorting or signed temperature response (P4). These are observed **phenotypic** statistics, not phylogenetic independence, genotype-verified morph maintenance, pigment production costs, fitness or causes of geographic clustering.

The source frame of **42,111 metadata species** is an opportunity frame, not the prevalence denominator. Species with no eligible locality pair or white/coloured comparison remain in the 4,475-species source registry with **NA** in that outcome. Photo-label missingness is material (220,393 original photo rows unclassified), and multiple controls including the stricter 5-km same-month null were unsupported, as described in the [full measured-data status and evidence routes](docs/FCP_CURRENT_MEASURED_RESULTS_AND_NP_DRAFT_STATUS_20261022.md).

**Historical method continuity:** The first four cohorts total **1,997 species / 199,700 photos**, and their original 50-km result is retained without numerical change. The completed fifth-cohort expansion yields **1,761 original strict-gate 50-km species**, mean +0.017614902 (199 permutation P=0.005). The earlier [17 October precompletion snapshot](docs/FCP_U100_REAL_199700_EMPIRICAL_EVIDENCE_AND_447500_PENDING_20261017.md) is intentionally kept intact as historical evidence and no longer describes the current computation.

## Reproduce the paper

The supported local entry point is:

```bash
python -m pip install -e .
python -m pip install pytest
python -m pytest \
  tests/test_polymorphism_manuscript_claims.py \
  tests/test_polymorphism_new_phytologist_submission.py \
  tests/test_fcp_np_all4475_submission.py \
  tests/test_fcp_latest_source_status_20261022.py \
  tests/test_make_polymorphism_manuscript_figures.py \
  tests/test_repository_layout.py -q
```

For the isolated full-source integration checks, see `.github/workflows/fcp-np-main-integration-guard-20261022.yml`, which reruns the historic/current manuscript safeguards, source-scope test and synthetic estimator/figure self-tests without rereading or replacing photo outcomes.

GitHub Actions additionally replay the frozen primary and secondary numerical analyses. Exact environments, immutable input hashes, release assets, and replay tolerances are documented in [`CURRENT_PAPER_REPRODUCIBILITY.md`](CURRENT_PAPER_REPRODUCIBILITY.md).

## Repository layout

| Path | Role |
|---|---|
| `docs/POLYMORPHISM_*` | Current manuscript, protocols, claims, lineage and supporting information |
| `results/polymorphism_*` | Frozen/current machine-readable study results |
| `scripts/analysis/` | Study analysis and figure code |
| `tests/test_polymorphism_*` | Claim, figure and submission regression guards |
| `archive/fcp_submission_20260925/` | Immutable provenance inputs, checksums and release receipts |
| `.github/workflows/` | Current manuscript validation and reproducibility workflows |
| `archive/workflows/` | Historical and superseded paper workflows retained for provenance but inactive in Actions |

## Reproducibility boundary

Repository organization must not alter frozen numerical results, biological decision rules, immutable source commits, or checksums. When prose and machine-readable artifacts differ, the authority order in [`CURRENT_PAPER_REPRODUCIBILITY.md`](CURRENT_PAPER_REPRODUCIBILITY.md) applies.
