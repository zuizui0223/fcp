# Polymorphism H2 third-cohort result and manuscript-claim freeze — 2026-09-17

## 1. Status

This document freezes the durable biological result and the manuscript claim boundary for the third-cohort prospective H2 test.

The authoritative machine-readable biological result is:

- `results/polymorphism_h2_third_cohort_prospective_white_axis_20260917/result.json`

The authoritative measurement/support receipt is:

- `results/polymorphism_h2_third_cohort_prospective_measurement_20260917/result.json`

If prose, chat text, notes, or later summaries disagree with those committed machine-readable files, the committed machine-readable files control.

Final frozen verdict:

- `H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED`
- stage: `H2_COMPLETE`
- status: `untouched_prospective_test_of_previously_frozen_axis`

## 2. What was tested

The test was the previously frozen white-versus-non-white direction

`q_white = normalize([1, -1/8, -1/8, -1/8, -1/8, -1/8, -1/8, -1/8, -1/8])`

with the previously frozen statistic

`W = mean_i (u_i dot q_white)^2`.

The analysis used the frozen construction-controlled structured null: within the third cohort and frozen coarse morph, normalized nine-colour rows were permuted across already selected H2 species while preserving species × coarse-morph row counts, followed by the unchanged label-free Hellinger two-means fit.

No alternative axis, threshold, null construction, or test statistic was selected after biological opening.

## 3. Prospective execution identity

- branch: `analysis/h2-third-cohort-preopening-20260916`
- biological execution head: `f7582767da237aeb344dfc072452c6daf92e662b`
- workflow run: `35177668182`
- terminal reassembly/test job: `105194354979`
- immutable biological result commit: `7e538e5c51c05a7cc47b2fcf53eea92634c8a863`
- result commit message: `Complete third-cohort prospective biological H2 test`
- metadata-freeze commit recorded at execution start: `ca69930986e39e3ea9b2d2f97ab247045ba9db0b`
- validated measurement-source commit recorded at execution start: `9fae6ccdf684a46026f72ba12e98de2c5c54bf2a`

The run was executed under the frozen one-shot/no-rerun contract. The terminal job explicitly required all 256 terminal receipts before opening the metadata-colour join, ran only the frozen prospective q_white/W test after the support gate, independently read back the durable terminal result, verified the terminal archive, uploaded the final artifact, and committed the immutable result.

## 4. Measurement and support gate

Frozen terminal measurement receipt:

- selected/authorized biological cohort entering the run: 499 species
- target rows per species: 100
- terminal rows: 49,900
- unique measurement IDs: 49,900
- duplicate measurement IDs: 0
- terminal partition receipts: 256 / 256
- classifiable rows: 25,788
- nonclassifiable rows: 24,112
- minimum classifiable rows per species for measurement evaluability: 40
- measurement-evaluable species: 377
- minimum required measurement-evaluable species: 250
- support decision: `PASS`
- replacement species: 0
- replacement rows: 0
- persisted image pixels: `false`

The support gate therefore passed before the H2 statistic was opened.

## 5. Frozen H2 result

### Primary admissibility tier — 0.10

- vector species: 158
- observed W: `0.5172457461053418`
- structured-null repetitions: 999
- structured-null mean: `0.45664389702916275`
- structured-null median: `0.45714281500506854`
- structured-null 2.5–97.5% interval: `[0.4358491119845972, 0.4752987775665408]`
- observed − null median: `0.06010293110027326`
- observed / null median: `1.1314751739007576`
- upper-tail structured-null p: `0.001`
- decision: `PASS`

### Strict sensitivity tier — 0.20

- vector species: 86
- observed W: `0.5329282123135909`
- structured-null repetitions: 999
- structured-null mean: `0.45967400889662496`
- structured-null median: `0.45931966586726264`
- structured-null 2.5–97.5% interval: `[0.43286795703924286, 0.48672240425829794]`
- observed − null median: `0.07360854644632825`
- observed / null median: `1.1602555952124205`
- upper-tail structured-null p: `0.001`
- decision: `PASS`

Both the primary test and the pre-specified stricter sensitivity tier support the frozen white-versus-non-white axis.

## 6. Chain of custody

### Outcome-blind candidate and selection freezes

- third-cohort outcome-blind candidate frame: 3,230 species
- candidate CSV SHA256: `7fc0337074b55a91c0b50773a8d5c3ef82e877074c8b3cacc21f9af58e0ba77e`
- frozen selection size: 500 species
- selection salt: `FCP_H2_THIRD_COHORT_20260916_V1`
- selected CSV SHA256: `4ad1191f39068e0fb2229f84361b1004803e24793190566d21e5c474aef2002a`
- selection-freeze commit: `2768b2dd0f8baf4ed1185a1b64c1060768b36c00`
- selected-species manifest: `results/polymorphism_h2_third_cohort_selection_20260916/selected_species_manifest.tsv`
- manifest SHA256: `16db6a2fc265f0ab20e7dae4e6f9058b907f5c19af3ddab5b6077797dd1def59`

The selection itself used only the frozen deterministic hash ordering and did not use biological outcomes, palette vectors, polymorphism geometry, D, H2 results, geography, climate, phylogeny, family, genus, or other outcome-bearing variables.

### Pre-opening terminal qualification

- qualifying head: `7fcabcda61c90aaac3ad84e7f53ad5cd6c8c5663`
- synthetic qualification run: `35065814180`
- synthetic qualification job: `104695722267`
- qualification tests: 17 / 17 PASS
- synthetic qualification artifact ID: `10434700102`
- synthetic qualification artifact digest: `sha256:2727d23a516cbf162b60606ba36932bd4fc1a4942ddc4242a2a03d18a6692e74`
- synthetic qualification evidence-freeze commit: `f8756913c82952382fbe05951a92f6174aaaac14`

This qualification was technical-only and did not constitute biological evidence.

### Final prospective terminal artifact

- artifact ID: `10496492307`
- artifact name: `h2-third-cohort-prospective-result-20260917`
- artifact digest: `sha256:319c040aaffc30d3cd97dcbcc217409e75010ec0f67bcc86c6bb82d8275779c7`
- artifact size: 11,526,750 bytes
- workflow run: `35177668182`
- execution head: `f7582767da237aeb344dfc072452c6daf92e662b`

## 7. Manuscript claim freeze

### Allowed main claim

> In a pre-frozen species-disjoint third cohort drawn from the same iNaturalist opportunity universe, the previously specified white-versus-non-white polymorphism direction was prospectively supported by the unchanged location-blind measurement pipeline and construction-controlled structured null. The result was supported at both the primary 0.10 admissibility tier (158 species, W = 0.51725, p = 0.001) and the pre-specified strict 0.20 sensitivity tier (86 species, W = 0.53293, p = 0.001).

### Allowed interpretation

This is a prospective species-disjoint transport/confirmation of the previously frozen white-versus-non-white direction **within the same iNaturalist source/opportunity universe**. It is stronger than a post hoc reanalysis of the legacy cohorts because species selection, the tested direction, the statistic, the support gate, both thresholds, the structured null, and the one-shot execution contract were fixed before biological opening.

### Required boundary language

The result is **not** an independent-source replication. It does not by itself establish:

- global flower-colour polymorphism prevalence;
- pigment chemistry or a pigment-loss mechanism;
- evolutionary direction of white/non-white transitions;
- pollinator, climate, or other adaptive causation;
- a recurrent non-white hue axis.

## 8. Relation to P500

The earlier P500 prospective execution remains **not durably evaluable** because its biological H2 calculation reached an in-memory terminal object but failed during post-calculation serialization before a durable terminal result was written.

The successful third-cohort test does not retroactively certify P500 as confirmed or non-confirmed, and no P500 replay can be re-labelled as untouched prospective confirmation under the frozen one-shot contract.

## 9. Post-freeze rule

Any analysis that changes the axis, construction gate, admissibility thresholds, statistic, null, cohort definition, or inferential decision rule after this result freeze must be labelled post-confirmatory/exploratory unless it is separately prospectively frozen before opening new biological data.

This document does not replace the machine-readable result files. It freezes their manuscript interpretation and claim boundary.