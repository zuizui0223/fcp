# Hesperis observation-subsampling bound

A source-exposed, conditional finite-census sampling audit. Start with docs/REPORT_JA.md and docs/PROOF.md. No raw stand counts were recovered or reconstructed.

Python 3.11+; analysis and independent verifier use only the standard library. Tests require pytest.

```
python sampling_bound.py
python verify_independent.py
python -m pytest -q
```

The copied parent_entropy.py is unchanged from PR173; it only supplies entropy and the prior rounding decomposition. The published input CSV is checksum-locked. Results are supplementary; no claim of selection, rejection of drift, corrected effect sizes, independent field replication, or a temporal trend. See docs/DESIGN.md for all declared conditions.
