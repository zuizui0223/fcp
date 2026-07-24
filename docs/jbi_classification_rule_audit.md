# Classification-rule audit

## Scope

This audit records technical changes to the deterministic text-screening rules. It is separate from the biological classification correction log because the change described here did not alter any frozen baseline label or model input.

## 24 July 2026 — plural `polymorphic populations`

### Trigger

The blinded-review packet compared each frozen baseline label with a fresh rule application to the resolved queue text. Thirty-three species matched immediately. For *Solanum carolinense*, the frozen label was `within_population`, whereas the queue text recomputed as `unclear`.

The classification source and queue best source were identical (`10.5091/plecevo.128527`). The retained title was:

> A field investigation into potential reproductive and resistance advantages of purple-petalled plants in polymorphic populations of horsenettle (*Solanum carolinense*)

The rule recognized singular `polymorphic population` but not plural `polymorphic populations` because of the terminating word boundary.

### Change

The within-population expressions in the following classifiers were aligned to use `polymorphic populations?` and the corresponding plural forms of `within population` and `same population`:

- `analysis_evidence_spatial_scale.py`
- `analysis_enriched_spatial_scale.py`

`analysis_evidence_spatial_scale_enriched.py` already used the plural-capable expression.

### Effect

- Frozen baseline label before the code correction: `within_population`.
- Recomputed queue label after the code correction: `within_population`.
- Number of frozen baseline classifications changed: **0**.
- Number of species added or removed from the 34-species baseline: **0**.
- Climatic values and model results changed: **no**.

This is a code-parity and reproducibility correction, not a post hoc biological reclassification. Human verification is still required before describing the labels as independently human-reviewed.

## Guard

`validate_jbi_classification_transparency.py` and the S18/S19 review-packet builder verify that:

- the manuscript calls the labels source-traceable and rule-derived;
- the resolved queue remains `unreviewed` unless new review metadata is deliberately added;
- the current rule can reproduce all frozen labels from the retained queue text;
- the blinded human-review sheet hides the frozen rule label;
- any future mismatch is surfaced rather than silently accepted.
