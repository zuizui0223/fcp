# Hesperis conditional colour-diversity hierarchy

Retrospective published-summary analysis. This integrates the unchanged analysis, independent verifier, tests and input from the conversation package `FCP_Hesperis_群落間差は単色群落だけか_20261012.zip` (48,641 bytes, SHA256 `997bafee8278cce17aa9dbbffc614161ae79a7abfc4be035599c55837b92aab5`). Only new repository integration and CI are added. Parent PR173 remains unchanged.

From repository root:

```sh
python -m pytest -q packages/fcp_hesperis_hierarchy_20261012/tests packages/fcp_hesperis_entropy_20261012/tests
python packages/fcp_hesperis_hierarchy_20261012/scripts/hierarchy.py
python packages/fcp_hesperis_hierarchy_20261012/scripts/verify_independent.py
```

## Interpretation

Shannon between-stand entropy is decomposed into differences among polymorphic stands, differences among monomorphic stands, and differences between their pooled colour compositions. All stands have equal weight. No stand-level rows are reconstructed.

For 2021 and 2023, conservative rounding enclosures place 66.675–73.821% and 55.883–63.347% respectively of the total between-stand entropy among the polymorphic stands themselves. These are not sampling confidence intervals. The analysis does not distinguish biological differentiation from finite-sampling noise or measurement error and does not identify selection, drift, urban effects, adaptation or a temporal trend.

The 2022 published all-stand and polymorphic-only white frequencies cannot simultaneously satisfy ordinary nearest-rounding bounds. Preserve `HOLD_INCONSISTENT_PUBLISHED_SUMMARIES`. Sensitivity with wider ±0.001 discrepancies is explicitly separate, not a repair to the source and not confirmatory evidence.

## Exact input source

Maunder et al. (2026), *Annals of Botany* 137:1290–1304, DOI `10.1093/aob/mcag035`. Table 2 supplies all-stand and polymorphic-only mean frequencies and monomorphic counts; Results Neutral predictions 1 & 2 supplies mean Shannon entropy. Original page: https://academic.oup.com/aob/article/137/5/1290/8488381 . Source values were rechecked against the publisher text during this integration. The raw stand-year counts, sampling denominators, and repeated-site IDs remain unrecovered.

## Verification scope

Before integration: 22 tests and 16,774 independent Fraction/60-digit Decimal checks passed after fresh execution. The workflow verifies that integrated code and input hashes match the delivered package, reruns this analysis and the parent summary analysis, and compares output hashes. Complete-repository pytest is attempted in a separate diagnostic step after installing the declared core dependencies; its exit status and full log are preserved, not converted into a claim that the whole suite passed. CI does not attempt raw-data recovery and is not independent biological replication.
