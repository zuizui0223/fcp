# Hesperis colour-classification information sensitivity (2026-10-12)

Standalone observational, source-locked diagnostic. Does not recover raw stand-year counts; does not change any 2026-10-12 parent PR.

Run from this directory:

```
python -m pytest -q tests
python src/information_bound.py --parent data/parent_entropy_result.json --outdir results
python src/verify_independent.py --parent data/parent_entropy_result.json --result results/result.json --outdir results
```

`docs/PROOF.md` gives the proof, `docs/REPORT_JA.md` the bounded ecological interpretation.
`results/annual_sensitivity.csv` has all scenario values, and `results/result.json` all metadata/nonclaims.

The source DOI is `10.1093/aob/mcag035`; the parent frozen result file SHA256 is `980b50bfe99f69b31b57d273a2995c814e280cf552a52a1eefaf17e828b4faf1`.
