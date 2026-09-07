# RGFCA research status

Updated 2026-09-07 (Japan). RGFCA is the active FCP mainline. The six-species
Chapter 1 and 34-species literature study are legacy and are not input to the
current ecological-signal search.

## Current evidence

| Stage | Evidence | Interpretation |
|---|---|---|
| Measurement | 500 species; 50,000 terminal records; 25,377 classifiable records | Automated public-photo colour measurements |
| Eligible RGFCA frame | 369 species; 21,424 photos; at least 40 classifiable photos per species | Conditional inference frame, not a census of all plants |
| G1 repeated field | Primary p = 0.070; fine-scale sensitivity p = 0.006 | Primary not supported; scale-sensitive candidate structure |
| Species-disjoint commonness | p = 0.856; median fold correlation = -0.0880 | No supported transfer of boundary geography across held-out species |
| Five environmental process blocks | None passes the fixed five-block Holm gate; thermal partial rho = 0.00836, adjusted p = 0.050 | Weak exploratory thermal candidate, no confirmed mechanism |
| Environmental effect heterogeneity | No corrected main-effect variance, interaction variance or syndrome support | Species/context-specific environmental mechanisms are hypotheses |
| Real-climate synthetic qualification | Failed; moderate full-sharing recovery 0%, 3.2%, 6.8% across three blocks | Insufficient recovery under the tested design; no biological result |
| Sharedness-specific predictive qualification | Failed; moderate full-sharing recovery 0.8% | Method development, not evidence of absence in nature |

Result files are `docs/supporting/global_monte_carlo_measurement_result_v1.json`,
`global_rgfca_g1_result_v1.json`, `global_rgfca_prespecified_robustness_result_v1.json`,
`global_rgfca_species_disjoint_commonness_result_v1.json`,
`global_rgfca_expanded_environmental_panel_final_result_v1.json`,
`global_rgfca_environmental_species_heterogeneity_inference_result_v1.json`,
`hypervolume_real_climate_synthetic_qualification_result_v1.json`, and
`global_rgfca_sharedness_specific_predictive_result_v1.json` in that same directory.

## In progress: species-level spatial information

Complete the frozen 369-species within-species spatial omnibus under
`global_rgfca_within_species_spatial_omnibus_contract_v1.json`. This tests an
equal-species mean association between geographic distance and colour
dissimilarity. It is exploratory because earlier descriptive G3 outcomes are
known; its 999 permutation outcomes have not yet been computed successfully.

Run `34085861343` stopped at `ModuleNotFoundError: No module named 'fcp_pipeline'`
in all species shards, before numerical execution. Technical recovery exposes
the checkout through `PYTHONPATH` in all workflow jobs and smoke-tests both real
entry points. The data, runner, fixed seeds, statistics and thresholds are unchanged.

## Route to an ecological result

1. Finish that test and inspect effect size, species information and measurement
   limitations before choosing further analyses.
2. Develop bounded, registered exploratory questions using the existing atlas
   and qualification work. Local patchiness, monotone distance effects, shared
   boundaries and environmental mechanisms are distinct targets.
3. Register promising candidates for independent validation. Every existing
   outcome set is discovery data for hypotheses motivated by it; re-splitting
   previously explored species does not create an untouched validation set.
4. Produce species-free maps, species-conditioned inference, ecological
   interpretation, manuscript and reproducibility artifacts. Existing legacy
   positives cannot stand in for a new RGFCA result.

The active goal is in `RGFCA_RESEARCH_GOAL.md`. Technical and synthetic successes
alone do not complete it. All completed negative and non-evaluable outcomes stay
in the research record.
