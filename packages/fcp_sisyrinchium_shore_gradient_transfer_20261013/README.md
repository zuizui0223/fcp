# FCP source-locked shoreline-gradient transport in *Sisyrinchium*

**Status:** retrospective source-exposed exploratory reanalysis of Takahashi et al. 2016, DOI `10.1371/journal.pone.0164381`. New numbers concern **predicting withheld field lines after borrowing slope estimates from other *named populations***. They are not selection, fitness, climate, genetic FST or new independent biological data.

1. Original Dryad data DOI `10.5061/dryad.j8q34`, official public mirror Zenodo record 5003124. The included `source/Sisyrinchium_data.zip` is original, hash-verified: 127,502 bytes, MD5 `008ec196aa5acb769dbfd36e835e32ba`, SHA256 `4859b803af66db1f8f18796fe08ce2581363258444de715f1e9f2f988c5` (CC0 original deposit). Internal `frequency.txt` SHA256 `a12cdf46483a6d8b78bd00110c196a25f3c8849967ae8a2c37d6aee3d29f516e`.
2. Parent PR185 immutable result and 185 predicted quadrat records: `source/parent_result.json` and `source/parent_predictions.csv`. Six latest named sites, 32 heldout field lines, 168 scored quadrats, 6,692 counted M plants. Violet/white M counts are distinct from other white L/S species.
3. Protocol frozen on GitHub before new transport outcomes: `docs/PROTOCOL.md`; GitHub `docs/FCP_SISYRINCHIUM_SHORE_GRADIENT_TRANSFER_PROTOCOL_20261013.md` commit `636d073b288a8fd048bf53c19e5ccf32c0e9678d`. Initial source preview had already shown slope signs; this is **not externally preregistered**.
4. Primary comparison is a source-site equal median of five other populations' ridge-logit coefficients **per original recorded distance unit**, adjusting only the target site's training-line intercept. Fixed site-equal binomial log loss. Secondary prespecified arithmetic mean of five. **Posthoc only:** distance-standardization alternative `posthoc_standardized.py`, separate JSON and ledger, not used to overwrite the primary claim.
5. `src/verify_independent.py` re-reads the original CR-delimited source and independently refits all 6×5 source/target models using Newton/Hessian without importing primary modules, and verifies every prediction and both spatial aggregations.

## Reproduce

```bash
python -m pip install -r requirements-lock.txt
python src/transport.py --source source/Sisyrinchium_data.zip --parent-result source/parent_result.json --parent-predictions source/parent_predictions.csv --outdir results
python src/posthoc_standardized.py --source source/Sisyrinchium_data.zip --primary-result results/result.json --outdir results
python src/verify_independent.py --source source/Sisyrinchium_data.zip --result results/result.json --posthoc results/posthoc_distance_scale.json --outdir results
python -m pytest -q tests/test_transport.py
python src/make_figures.py --root .
```

Frozen baseline and previous PR185 results unchanged. Main H1/H2 protected; Draft branch only. Complete ecological caveats and numeric results in `docs/REPORT_JA.md`.
