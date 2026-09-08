# Completed independent RGFCA reserve replication

8 September 2026. **Weak photo association replicated; flower-specific robust
replication did not pass. Not submission-ready.** This reports the already
completed fixed tests. It neither changes their design nor opens a new analysis.

## 1. Full cohort and measurement gate

The reserve contains the complementary 500 species and 50,000 photographs from
the original 1,000-species candidate frame. It is species-disjoint from discovery
and has no observation/photo-ID overlap with discovery or the three other prior
sets listed in the metadata freeze. It shares iNaturalist, measurement models,
regions and potentially observers; disjoint IDs do not remove systematic error.

All 256 terminal partitions (512 CSV/receipt files) and all 50,000 unique photo
and observation IDs were verified before authorization of inference. Of 24,885
classifiable photographs, **20,903 from 363 species** passed the fixed
at-least-40-photo-per-species rule, exceeding the required 250 species. All
non-evaluable terminal records remain: 22,681 ROI/flip-gate failures, 2,415
ambiguous palettes and 19 without biological-palette mass. These are measurement
outcomes, not ecological frequencies.

Discovery (500 species / 50,000 photographs; 369 eligible species / 21,424 photos)
and reserve thus total **1,000 species / 100,000 measured photographs**. They are
not pooled retroactively for inference. All existing figures remain discovery-only.

## 2. Every prespecified result

All 20 inference shards completed the same four tests over all 363 admitted
species, with 999 stored randomizations per test. All null distributions were
nondegenerate and evaluable.

| Fixed test | Equal-species mean rho | Upper-tail p | Positive-support gate |
|---|---:|---:|---|
| Primary distance-colour association | 0.0254826061 | 0.001 | Passed |
| Observer-pair exclusion | 0.0252180230 | 0.001 | Passed |
| Calendar-quarter-stratified randomization | 0.0254826061 | 0.001 | Passed |
| Matched flower-minus-background differential | 0.0044772519 | 0.087 | Did not pass |

The primary result is a replication of the discovery's weak positive association
(rho = 0.0270213, p = 0.001), not a claim of a large ecological effect. The
primary mean's conditional species-bootstrap 95% interval is
**[0.0170055402, 0.0343321275]**, using the already-fixed 4,999 draws and seed
202609070903. This is **not spatially or phylogenetically independent uncertainty
for all plants**. Shard-order species are retained when reconstructing that
specific bootstrap; changing row order would change the fixed draws.

The quarter control uses the same observed pairs and coefficient as the primary
test, while changing the permitted randomization. It is a calendar proxy, not
measured phenology. The background statistic is
`Spearman(geographic distance, flower12 JSD - background12 JSD)`, averaged equally
over species; it is **not subtraction of two Spearman coefficients**. Flower and
background counts were taken symmetrically from the first image decode under
the prospective reserve protocol.

The frozen decision is:

- `directional_photo_association_replicated = true`;
- **`flower_specific_robust_replication = false`**, because all four gates were
  required and the background differential did not pass;
- `cause_or_shared_boundary_inferred = false`.

The background p = 0.087 is **non-support**, not proof of no effect, equivalence
to zero or background causation. No extra permutations, new seeds, altered
thresholds, replacement species or selected successful subset were used. A
passing observer/calendar check does not eliminate every observational bias.

Species-conditioning refers to taxon-labelled photographs. The detector pools
retained flower regions without assigning every region to the observation's
focal species. Neither this replication nor a background contrast validates
taxon-to-mask attribution. The [ROI audit](RGFCA_ROI_QUALIFICATION_AUDIT.md)
documents that unresolved measurement limitation.

## 3. Distinguish two different background outcomes

The older postoutcome discovery reacquisition was **not evaluable**: 85 of
21,424 records failed exact reproduction, and no adjusted statistic was computed.
The reserve's prospectively specified, first-decode differential was **evaluable
but not supported**. It is a separate design, not a repair or replacement of
the discovery recovery. Both outcomes remain in the
[manuscript](RGFCA_MANUSCRIPT.md) and
[discovery recovery audit](RGFCA_BACKGROUND_RECOVERY_COMPLETION.md).

## 4. Immutable execution and audit trail

| Stage | Immutable evidence |
|---|---|
| Design and metadata freeze | `9fd4ae98632c74caa9e66dd24e7390ec221efb4c` |
| Inference implementation preflight, 27 tests | [34091922722](https://github.com/zuizui0223/fcp/actions/runs/34091922722), head `708b305e314087047a39d27956a76c9327585497` |
| Complete blind measurement | [34091091640](https://github.com/zuizui0223/fcp/actions/runs/34091091640), execution `1f80af7f5db81d61f28ae2818c130ef058a9f850`, result `9f5abe7b45fcdc20ba83adf75d1a8d4a640f622c` |
| Complete-census inference authorization | `391caecaa1a6d3c3e2407f5f19d0bad7057922e0`, binding 13 tested code/design files and three completed measurement files |
| All four completed inference tests | [34178957447](https://github.com/zuizui0223/fcp/actions/runs/34178957447), result `ef00a78a2eda7bc79f7e6b88f719a8a696dc8b43` |

The final inference artifact is `10038248014`,
`rgfca-reserve-replication-full-result-v1`, archive SHA-256
`874d71daec34407adf5809f67d3ee6bb08867ef27ff0e2437bf2cb4f18dc2a8b`.
The [fixed result JSON](supporting/rgfca_reserve_replication_result_v1.json)
has exact source-byte SHA-256
`1cb9de6feb8bda3bf05f3fbe51abd8eeebb41774a9cd8253034d618b78dcb3f3`.
It records both output CSV hashes and all 20 shard receipts, including every
species allocation and shard table/null-array hash.

The [measurement census audit](supporting/rgfca_reserve_complete_census_audit_v1.json)
reassembles all terminal files through the frozen completeness checks. The
[inference artifact audit](supporting/rgfca_reserve_inference_artifact_audit_v1.json)
verifies every shard and the full 363 x 4 x 999 saved null tensor, reconstructs
the reported summaries and fixed bootstrap, and reconciles 2,541 direct SciPy
checks (maximum difference 1.1102230246251563e-16). It uses absolute tolerance
2e-12 only for numerical reconstruction; file hashes are exact, never relaxed.
It obtains no new measurements or inferential permutations.

Reconstruction scripts are
`scripts/analysis/audit_rgfca_reserve_complete_census.py` and
`scripts/analysis/audit_rgfca_reserve_inference_result.py`. Run from the repository
root with `PYTHONPATH=.` after downloading the literal runs' artifacts. Their
module docstrings specify the local artifact layout and exact-source dependency
requirements. Publication tests independently check the committed full result,
all 363 species rows, all 999 global null rows per test and the fixed bootstrap;
they do not require image acquisition or repeat the inference.

## 5. Next scientific gate

The present data establish a reproducible, weak photograph-level spatial
association and an unresolved flower-specific interpretation. They do not meet
the active goal of an independently validated ecological signal. Preserve this
result as a bounded research output. Before a stronger ecological claim, close
the target-domain measurement/attribution gap and specify a genuinely fresh
validation route. These 1,000 species' image outcomes are now opened; they
cannot be relabelled as unused validation. Environmental or pollinator questions
require a new estimand, documented source layers, appropriate spatial nulls and
prospective validation, not repeated testing until this failed gate passes.
