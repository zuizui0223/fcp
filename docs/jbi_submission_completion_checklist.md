# Journal of Biogeography submission completion checklist

This checklist contains only **current submission gates** for the frozen 34-species paper. Historical PR numbers, S1–S19 numbering, expired workflow artifacts and superseded sensitivity analyses are intentionally excluded.

## Analysis package

- [x] Canonical manuscript is `docs/jbi_manuscript.md`.
- [x] Canonical downstream input is `data/frozen/frozen_34species_five_metric_dataset.csv`.
- [x] Dataset checksum is fixed in `data/frozen/freeze_manifest.json` and enforced by CI.
- [x] Frozen scope is 34 species, 25 families, 20 within-population / 14 geographically structured cases.
- [x] Five climatic-niche metrics are analysed symmetrically with `among ~ metric_z + effort_z`.
- [x] Family-clustered inference, 9,999 row-order-invariant label permutations, Holm context and leave-one-family-out refits are reproduced by the canonical workflow.
- [x] Collinearity diagnostics, Open Tree sensitivity and dated V.PhyloMaker2 S1–S3 sensitivity are reproduced.
- [x] CR2/Satterthwaite and design-based power/precision diagnostics are reproduced.
- [x] `tests/test_manuscript_consistency.py` prevents stale result values and removed exploratory analyses from re-entering the canonical manuscript.
- [x] Supporting-data roles are indexed in `docs/jbi_supporting_information_index.md` without historical S-number accretion.

## Classification review

- [x] Frozen source-traceable classification manifest retained.
- [x] Classification rule audit and review protocol retained.
- [x] Blinded review sheet and separate rule key prepared.
- [ ] Complete independent blinded review of the 34 classifications.
- [ ] If reviewer labels disagree with the frozen rule label, document adjudication explicitly before changing any analytical input.
- [ ] Update manuscript wording only after the completed review record exists. Until then, retain **source-traceable, rule-derived classifications**.

## Occurrence-data citation

- [x] GBIF taxonomic-resolution audit retained.
- [x] Preparatory GBIF DOI / Derived Dataset bundle retained for provenance and external registration.
- [ ] Obtain a citable GBIF occurrence identifier or Derived Dataset DOI for the occurrence data underlying the frozen climatic summaries.
- [ ] Insert the final GBIF citation into the manuscript Data Accessibility Statement.

The current manuscript analysis does **not** depend at runtime on a GitHub Actions occurrence artifact. The climatic summaries are committed in the durable frozen 34-species downstream dataset.

## Repository archive

- [ ] Freeze the final accepted submission commit after all manuscript edits are complete.
- [ ] Archive that exact repository state at a permanent DOI (for example Zenodo or another approved archive).
- [ ] Insert the permanent repository DOI into the Data Accessibility Statement.
- [ ] Verify that the archived release contains the same frozen dataset checksum used by CI.

## Author-controlled submission fields

These fields require author confirmation and must not be inferred from repository history:

- [ ] final author list and order;
- [ ] corresponding author;
- [ ] affiliations and current contact metadata as required by the journal;
- [ ] ORCID identifiers;
- [ ] CRediT roles;
- [ ] funding and grant numbers;
- [ ] conflict-of-interest declaration;
- [ ] acknowledgements;
- [ ] any journal-required biosketch or author statement;
- [ ] final cover-letter declarations.

Templates can remain separate from the analytical pipeline because these fields do not affect the scientific results.

## Final scientific checks

Before changing the PR from draft to ready for review:

- [ ] verify the latest `Frozen 34-species paper pipeline` run is green on the final manuscript commit;
- [ ] verify Table 2 values match the generated five-metric model output;
- [ ] verify Table 4 values match the generated Open Tree and dated-tree outputs;
- [ ] verify the finite-sample text matches the generated CR2 and power/precision outputs;
- [ ] confirm no manuscript claim treats occupied climatic breadth as physiological tolerance or morph-specific climatic adaptation;
- [ ] confirm moisture breadth is described as the strongest observed contrast, not a uniquely established causal mechanism;
- [ ] confirm no unreviewed expanded-set, matched-control, fragmentation or environmental-turnover result appears in the active submission narrative.

## Remaining external blockers

At present, the only blockers that cannot be completed from repository code alone are:

1. completed human classification review / adjudication;
2. external GBIF data citation registration;
3. permanent repository archival DOI;
4. author-controlled authorship and declaration fields.

Everything else in the analytical paper package should be required to pass the canonical repository CI.