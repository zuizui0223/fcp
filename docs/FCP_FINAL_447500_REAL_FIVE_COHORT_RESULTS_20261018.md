# FCP original 447,500-photo five-cohort analysis: completed empirical evidence

**Evidence reviewed:** 18 October 2026 JST. **Status: genuinely measured, analysed and permanently released**, not a source-only preflight. FCP is **within-species** photographed floral-colour structure; MOIS is the separate global interspecific/seasonal moisture programme.

## Immutable evidence and chronology

- Complete original ROI-v4/EfficientSAM [run #38997169400](https://github.com/zuizui0223/fcp/actions/runs/38997169400): 256/256 actual new-U100 image partitions succeeded; the original measurement/terminal inference also succeeded. Source-authenticated complete 248,100 new-U100 raw photo terminals were kept through initial ROI to allow exact pre-specified cross-cohort reuse exclusions.
- Five-source [join and unchanged original 50-km/H1/W run #39049598596](https://github.com/zuizui0223/fcp/actions/runs/39049598596): all five named jobs completed **success**, including source join, exact 199-label-permutation spatial test, observer-disjoint H1 and structured-null W; the source run is not a PR preflight.
- [Release workflow #39050465077](https://github.com/zuizui0223/fcp/actions/runs/39050465077) independently validated 447,500 actual photo terminals, the original three frozen source hashes and cohort spatial effects, immutable W/H1 execution, and successfully published and round-trip verified the complete [permanent five-cohort real-results release](https://github.com/zuizui0223/fcp/releases/tag/fcp-all4730-original-method-actual-colour-results-20261017). **Do not mistake later green workflow_run jobs with skipped statistical steps for the primary science run.**
- Archive \`fcp_all4730_original_method_actual_measured_data_and_inference.tar.gz\`: SHA256 **\`99f7ee68c3fd009e8ab60fd0cfcdf1fb407e3a8de07edd38341b27889310bc25\`** (108,115,413 bytes). Release receipt SHA256: \`b0b152fcfaa1e0015a6e5524192265bb95c1663c10b5da3a4e885cb7be49747c\`. Source merged measured CSV SHA256: \`8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6\`.
- The release also contains the final auto-authenticated Results packets in JSON (SHA256 \`e4d5244fb7e9de0046d44cf01154573d12f2ea3418fa538022a8dd97b1d60dc6\`) and Markdown (SHA256 \`d129d44770977b7506a3151f47c745936b80737251a00f163434301b64b86515\`). The Git-tracked reporting summary is \`results/fcp_original_447500_five_cohort_verified_reporting_20261018.json\`, compiled from the actual successful science job logs and immutable release metadata, not a recomputation or new estimator.

## Actual denominators, not opportunity-only rows

Of **4,730** candidate U100 taxa, 252 had fewer than 100 authentic eligible original H9 photographs. Another **3 newly sampled taxa** had authentic reused prior original photo/observation IDs and were conservatively removed *whole*, preserving all 1,499 original discovery/validation/prospective taxa and every original measurement unchanged. Thus **4,475 taxa × 100 real original photo terminals = 447,500 rows**.

| Cohort | Distinct 100-photo species | Actual photo terminals | Classifiable four-state rows |
|---|---:|---:|---:|
| Discovery (frozen) | 500 | 50,000 | 25,377 |
| Validation (frozen) | 500 | 50,000 | 24,885 |
| Prospective confirmation (frozen) | 499 | 49,900 | 25,788 |
| Fourth wave (post-freeze, measured) | 498 | 49,800 | 24,585 |
| Fifth wave new U100 | 2,478 | 247,800 | 126,472 |
| **All five** | **4,475** | **447,500** | **227,107** |

The other **220,393** original photo terminals remain nonclassifiable; they are not biological white/coloured states. **3,297 species** passed the original >=40 classifiable-photo D/H1 support threshold. The separate, more restrictive **50-km test includes only 1,761 taxa**, which also passed the >=30 local-photograph-pair gate. Neither 4,475 nor 3,297 may be substituted for its null-test denominator.

## Primary unchanged geographic comparison: species-global minus within-50-km local pair dissimilarity

The effect for each species is the full frozen four-state global mismatching-pair fraction minus the corresponding proportion among photo pairs separated by at most **50 km**. All species receive equal weight. The 199 within-species **complete-vertex label permutations** preserve each species' observed colour-state counts and coordinates; the minimum one-sided plus-one upper p is 0.005. **Fifty kilometres was chosen post hoc**, not a preregistered evolutionary scale.

| Actual image cohort | 50-km eligible species | Mean depletion | Upper permutation p |
|---|---:|---:|---:|
| Discovery | 166 | 0.020529255 | 0.005 |
| Validation | 181 | 0.018672972 | 0.005 |
| Prospective | 204 | 0.014684920 | 0.005 |
| Fourth wave | 190 | 0.017293596 | 0.005 |
| **Fifth new U100** | **1,020** | **0.017598698** | **0.005** |
| **Five-cohort pooled** | **1,761** | **0.017614902** | **0.005** |

In all five cohorts, the local observed photo-colour composition was more homogeneous than expected by conditionally shuffling each species' unchanged observed states among its original sampled coordinates. The new fifth-cohort effect was **0.01760**, close to the old four-cohort pooled **0.01764**. All three original cohort statistics were reproduced byte-for-byte at the numerical inference level; this is a later, same-provider species-disjoint **fixed-method expansion**, not an independent-source replication and not a newly untouched prospective test.

For the final 1,761 species, the mean global mismatching-pair proportion was **0.2243045205** and local mean was **0.2066896183**. The original 199-world pooled null for mean depletion had mean **0.0000445943** and central 95% range **[−0.0020032432, 0.0024735998]**. Approximately **55.88%** (984/1,761) of eligible species had positive individual depletion; most species need not have strong individual effects. The species-level D-to-local-depletion association is positive but small (rho = **0.1479056**, matched upper p = **0.005**).

All four additional originally computed radii (25, 50, 100, 250 km) were kept unchanged; for the newly measured fifth cohort all four mean depletions were positive with upper p = 0.005 (25-km n=710, 100-km n=1361, 250-km n=1772). Those are **post hoc scale sensitivities**, not a search for a new best radius.

## Original observer-disjoint H1/D and frozen white-versus-chromatic W

The original 200 observer-disjoint H1 repeated splits, applied without estimator changes to the full real 447,500-photo table, yielded median cross-split species-D Spearman **rho = 0.8035357424** (5th–95th percentile **0.7931469330–0.8094055921**), median **3,044 paired species**, median Lin CCC **0.8553147636**, and median Spearman–Brown **0.8910671671**. This is a reproducible photographed species-level *D*, not a population morph-frequency census.

The original prespecified zero-sum white-versus-eight-other-colour axis W, with its **999 construction-preserving structured null assignments**, gave:
- primary second-mode threshold 0.10: **1,223 vector-eligible species**, observed **W = 0.5187948657**, null median **0.4553650050**, upper **p = 0.001**;
- strict 0.20: **658 species**, observed **W = 0.5254168292**, null median **0.4612009436**, upper **p = 0.001**.

These fixed-W tests use *squared, sign-invariant axial alignment*; they do **not** identify an evolutionary increase in white flowers, pigmentation costs, molecular anthocyanin gain/loss or adaptive fitness. The original prospective third-cohort H2 confirmatory result (W = 0.5172457461, n=158, p=0.001) remains a distinct, earlier untouched experiment; pooled five-cohort W is not a substitute confirmation.

## What is and is not new

**Genuinely extended evidence:** the historical, unchanged methods and fixed target now run on 447,500 actual source-photo terminals, including 2,478 additionally eligible independent taxa, 1,020 of which supply true local pairs at the 50-km scale. The separate newly measured fifth cohort reproduces the **positive same-magnitude geographic depletion**.

**Unchanged limitations:** community-science photograph presence/detection, image/classification coupling (notably white/near-clipping), species-range sampling and spatiotemporal/observer structure constrain field-population inference. Selection, gene flow, founder history, plasticity and pollinator responses are *competing hypotheses*, not mechanisms identified by this null. Earlier nonwhite-only, same-observer, calendar-quarter and continuous-colour falsifications were established for the initial three cohorts; **they are not claimed to have been rerun for all 2,478 additional taxa**. The earlier four-cohort median/top-tail robustness is likewise a distinct 199,700-photo post hoc analysis and cannot silently be promoted to a full five-cohort test.

**Superseded status:** \`docs/FCP_U100_REAL_199700_EMPIRICAL_EVIDENCE_AND_447500_PENDING_20261017.md\` is a historically accurate precompletion snapshot, **not** the current study status. Preserve that document for provenance while using this document and the final actual source release for current claims.
