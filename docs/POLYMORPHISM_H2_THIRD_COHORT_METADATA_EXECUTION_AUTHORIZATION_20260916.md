# H2 third-cohort fresh-metadata execution authorization

Date authorized: 2026-09-16 JST
Status: AUTHORIZED FOR ONE OUTCOME-BLIND METADATA DRAW ONLY

Implementation parent HEAD: `191c05760c97990388a4eb21ece18069b2c8a2b4`

This authorization permits exactly one GitHub Actions attempt of `.github/workflows/h2-third-cohort-fresh-metadata.yml` on the frozen 500-species third cohort. The authorization commit must differ from the implementation parent by this file only. Any rerun attempt is forbidden and must fail before network requests.

The preceding authorization run `35107243536` failed during the pre-request dependency/test step because `pytest` was not installed. All three frozen fingerprints passed before that failure; the metadata acquisition step was skipped, so zero third-cohort metadata requests were made and no biological or colour outcome was opened. This authorization supersedes that unconsumed technical preflight authorization after the dependency-only workflow repair.

Frozen inputs and fingerprints:

- selected species manifest: `results/polymorphism_h2_third_cohort_selection_20260916/selected_species_manifest.tsv`
- selected species manifest SHA256: `16db6a2fc265f0ab20e7dae4e6f9058b907f5c19af3ddab5b6077797dd1def59`
- canonical selected CSV SHA256: `4ad1191f39068e0fb2229f84361b1004803e24793190566d21e5c474aef2002a`
- fresh-metadata protocol SHA256: `3cb1ad1126c460b2b9dd8dd837662023ca252e941dbb84fa3d05e94f5cfc9c59`
- acquisition executor SHA256: `5bbc725aa677b1f61b10cdd0eacb371f45dd8ab6748405edcfd76f468076fb54`
- prior observation/photo exclusion union: exactly 178,462 / 178,462 IDs
- target rows per species: 100
- observer cap: 2
- minimum full-target species for metadata PASS: 300
- request-error ceiling: 0.05
- request retries: 0
- no species replacement
- no target relaxation

The draw may read only opportunity/observation metadata required by the frozen sampler. It must not open or compute image pixels, flower-colour outcomes, morph, palette, D, H2 vectors, W, structured-null outcomes, H3 predictors, or recovered P500 H2 outcomes.

A metadata PASS does not authorize biological pixel opening. After the resulting exact species/row denominator is durably committed, a separate one-shot biological execution authorization is required.

A metadata transport/capacity failure is not a biological negative and must not be repaired by rerun, species replacement, or threshold retuning.
