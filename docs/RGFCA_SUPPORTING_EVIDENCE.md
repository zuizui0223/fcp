# RGFCA supporting evidence and outstanding submission gates

Discovery audit, 7 September; independent reserve update, 8 September 2026.
**Share with caveats; not submission-ready.** This is the supporting index for the active
[manuscript](RGFCA_MANUSCRIPT.md), not a new analysis or a replacement for any
original execution contract. Independent reserve inference is complete, but
the required flower-specific robust-replication gate did not pass. No reserve outcome
was opened for the original discovery audit; the separately identified reserve
result below was subsequently opened after complete-census authorization.

## S1. Populations, estimands and completed evidence

The following files are fixed at commit
`730612e4af6c64178b1c079d23bb6dc802e1dac0`; their committed Git blob identities are
checked by `tests/test_rgfca_supporting_evidence.py`. The check uses exact Git
object bytes, not a claim of byte identity between platform-specific working
copies. Figure input SHA-256 values remain separately enforced by the figure
builder. No new permutations are run by these publication checks.

- [Discovery measurement](supporting/global_monte_carlo_measurement_result_v1.json):
  500 species, 50,000 terminal photos; 25,377 classifiable and 24,623 uncertain.
  The at-least-40 rule admits 369 species; 21,424 classifiable photos enter
  inference. Classifiable images from excluded species are not added to that
  denominator. The candidate frame has 1,000 species, not all flowering plants.
- [Within-species omnibus](supporting/global_rgfca_within_species_spatial_omnibus_result_v1.json):
  species-equal all-pairs rho = 0.0270213, p = 0.001. This is postoutcome discovery.
  The eight original execution tests, complete 20-shard reconstruction and 12
  direct calculation checks belong to its original receipt; current publication
  tests do not claim to repeat that computation. Null quantiles are not an effect
  confidence interval, and nine BH-detectable species do not estimate prevalence.
- [Primary repeated field](supporting/global_rgfca_g1_result_v1.json):
  p = 0.070 from 69 exceedances in 999 null draws; not supported. The 200 balanced
  map realizations are repeated subsamples, not additional biological replicates.
- [Prespecified robustness](supporting/global_rgfca_prespecified_robustness_result_v1.json):
  coarse 1,000-km p = 0.043, fine 250-km p = 0.006 and minimum-ten-species-support
  p = 0.260 are descriptive sensitivities. Odd/even p = 0.115/0.038. The full
  strong-stability claim remains not evaluable in that receipt. None replaces
  the primary p = 0.070 or establishes species-disjoint transfer.
- [Species-disjoint commonness](supporting/global_rgfca_species_disjoint_commonness_result_v1.json):
  median of five fold correlations = -0.0880319, p = 0.856; only two folds have
  positive coefficients. Shared boundary architecture is not supported.

## S2. Environmental effects, interactions and heterogeneity

[Five-block mean-effect panel](supporting/global_rgfca_expanded_environmental_panel_final_result_v1.json):
all five predeclared blocks are retained. Their raw/Holm p-values are thermal
0.010/0.050, atmospheric energy/dryness 0.026/0.104, edaphic 0.148/0.444,
terrain 0.150/0.444 and water balance 0.338/0.444. The frozen strict p < 0.05
gate admits none. Soil-layer coverage is only 58.45% of edge occurrences in the
reported coverage audit; all five blocks must not be described as having
identical data support. These postoutcome associations do not identify causes.

[Ten-interaction family](supporting/global_rgfca_environmental_interaction_family_verification_v1.json):
all ten pairs and 369 species per pair were checked. No pair passes the fixed
Holm correction. The strongest nominal pair, water balance × terrain, has raw
p = 0.020 but adjusted p = 0.200. Nominal results do not authorize selective
variable-level decomposition or rescue a parent decision.

The separate **heterogeneity inference** is stored on the parent branch, not
as a local file on this branch. The former status entry incorrectly implied a
local path. The completed [run 34039210812](https://github.com/zuizui0223/fcp/actions/runs/34039210812)
and its three-file artifact were directly checked against immutable commit
`f2f9c58e1d857a5d4b5b35a78adab8e8324b25eb`:

- [Result JSON](https://github.com/zuizui0223/fcp/blob/f2f9c58e1d857a5d4b5b35a78adab8e8324b25eb/docs/supporting/global_rgfca_environmental_species_heterogeneity_inference_result_v1.json),
  SHA-256 `c9b8e29c6a9b8f8588417031c9394309be973a1c17aca97c677ce89f91defd36`.
- [All 999 null rows](https://github.com/zuizui0223/fcp/blob/f2f9c58e1d857a5d4b5b35a78adab8e8324b25eb/data/derived/global_rgfca_environmental_species_heterogeneity_inference_null_v1.csv),
  SHA-256 `810bb23a89f651fa6c184ccaa60ee5352f27f6429db134ec729c4448413600fd`.
- [All 369 species rows](https://github.com/zuizui0223/fcp/blob/f2f9c58e1d857a5d4b5b35a78adab8e8324b25eb/data/derived/global_rgfca_environmental_species_heterogeneity_inference_species_v1.csv),
  SHA-256 `ecb95ed7e92145df44957cf8877b602e9d648278bf5042d59e158083f620b138`.

Artifact ID: `9992993814`; GitHub-reported archive digest:
`1412543c383531b82cbdd8c3f34208f5e85b8860af5fae22292f22daa1a28d6f`.
Downloaded file bytes match all three committed Git objects. This is a source
reconciliation, not a rerun. Five main-effect variance tests, ten interaction
variance tests and two coordinated-syndrome statistics all remain unsupported
after their respective corrections. Publication tests reconstruct all 17
reported tail probabilities and the three Holm families from saved null values,
and all 15 observed variances from saved species coefficients. They do not
rerun environmental extraction or the original shared-null permutations.

## S3. Failed qualification is not a biological negative

- [Real-climate synthetic qualification](supporting/hypervolume_real_climate_synthetic_qualification_result_v1.json):
  moderate full-sharing recovery was 0%, 3.2% and 6.8% in the three tested blocks.
  The qualification gate failed. Those are synthetic operating characteristics,
  not estimates of the frequency of real ecological boundaries.
- [Sharedness-specific predictive qualification](supporting/global_rgfca_sharedness_specific_predictive_result_v1.json):
  moderate full-sharing recovery was 0.8%; qualification failed. This does not
  open a new empirical sharedness claim or demonstrate absence in nature.

These results and the links in [current status](RGFCA_RESEARCH_STATUS.md) preserve
the tested alternatives. This index is not an assertion that the entire historic
repository is one prospectively specified statistical family. The six-species,
34-species and earlier small photo-first studies remain separate legacy lanes.

## S4. Current validation assessment

Numerical reconstruction and source identity support sharing the bounded
discovery result with its caveats. Figure source/denominator and double-render
checks are documented in [the figure audit](RGFCA_PUBLICATION_FIGURES.md).
The [image/ecology](RGFCA_IMAGE_ECOLOGY_LITERATURE_AUDIT.md) and
[statistical](RGFCA_STATISTICAL_LITERATURE_AUDIT.md) readings limit interpretation;
neither software checks nor precedent establishes this signal's ecological cause.
The [real photo-bar release](RGFCA_PHOTO_BAR.md) now provides 24 verified CC0 ROI
crops with source links, exact reproduction checks and a separate figure receipt.
This closes the missing-photo illustration item, not the measurement-validity
or independent-replication gates. Four [core software references](RGFCA_SCIENTIFIC_SOFTWARE_AUDIT.md)
have also been checked against project-owned citation records.
The [four-provider audit](RGFCA_MEASUREMENT_PROVIDER_AUDIT.md) adds exact model
and basemap identities, separates photo from observation-data licences, and
records the remaining model/training-data and mixed-licence reuse questions.
The separate [JRC source audit](RGFCA_TRAINING_SOURCE_AUDIT.md) now confirms the
dataset-specific CC BY 4.0 notice and the source article's point-before-slice
partition. It preserves the dataset's conflicting 2026 citation/2022 issued
dates rather than treating them as new image collection. The
[ROI qualification reconstruction](RGFCA_ROI_QUALIFICATION_AUDIT.md) checks all
100 saved test rows and all eight original gates without model execution.
Box precision/recall are 0.730447/0.795563; 85/100 images pass admission.
The large-object denominator is only four boxes, and mask containment is not
petal segmentation accuracy. These audits do not close target-domain measurement
validation or the remaining model-distribution/mixed-licence questions.

### Completed independent reserve evidence

The [complete result report](RGFCA_RESERVE_REPLICATION_RESULTS.md) and
[artifact audit](supporting/rgfca_reserve_inference_artifact_audit_v1.json) bind
result `ef00a78a2eda7bc79f7e6b88f719a8a696dc8b43`, inference run
[34178957447](https://github.com/zuizui0223/fcp/actions/runs/34178957447), and all
20 shard receipts. The full reserve has 500 species / 50,000 terminal records;
363 species / 20,903 photos meet the fixed eligibility rules. Primary rho is
0.0254826 (p = 0.001), observer rho 0.0252180 (p = 0.001), quarter rho 0.0254826
(p = 0.001), and matched-background differential rho 0.0044773 (p = 0.087).
All four tests were evaluable; the flower-specific robustness conjunction failed.
This is not proof of no effect or of background causation. The conditional
species-bootstrap primary interval [0.0170055, 0.0343321] is not spatially or
phylogenetically independent uncertainty for all plants. Publication checks
reconstruct all four results from the exact committed tables and stored nulls;
no new permutation test is performed.

The material submission blockers are scientific, not cosmetic:

1. The reserve measurement, full census and all fixed inference tests are now
   complete. Retain the weak replicated photo association and failed
   flower-specific gate together; do not rerun opened cohorts until that gate
   passes. Stronger ecology requires target-domain measurement validation and
   a genuinely fresh, prospectively specified validation route.
2. Preserve the now-completed discovery background recovery as not evaluable:
   21,339 exact rows and 85 reproduction failures, no background-adjusted test.
   The [complete artifact audit](RGFCA_BACKGROUND_RECOVERY_COMPLETION.md) closes
   the execution/accounting task, not flower-specific validation. Stronger
   interpretation is not established by the now-completed reserve controls and
   cannot be certified by the discovery p-value or a successful recovery subset.
3. Complete remaining software-release, data-provider, model and environmental-layer references,
   and audit the final manuscript/SI/reproducibility package together.

The independent photo-association replication is now documented, but no
flower-specific robust replication, shared global boundary, pollination
mechanism, completed submission bundle or submission is claimed.
