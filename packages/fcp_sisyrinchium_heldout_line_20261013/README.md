# *Sisyrinchium* held-out field-line prediction: distance vs co-flowering white species

Explore an ecology question in the original Takahashi et al. 2016 field frequency data **without pseudoreplicating quadrats or conflating M-white with L/S white**. This package is posthoc and source-exposed. Read `docs/PROTOCOL.md` (frozen on GitHub before this analysis) and `docs/REPORT_JA.md` (all outcomes and limitations).

## Original author source

Takahashi Y, Takakura K-I, Kawata M (2016), *PLOS ONE* 11:e0164381, DOI https://doi.org/10.1371/journal.pone.0164381; Dryad original DOI https://doi.org/10.5061/dryad.j8q34; legitimate same-DOI public Zenodo mirror https://zenodo.org/records/5003124 (CC0 original data). Exact original ZIP 127,502 bytes MD5 `008ec196aa5acb769dbfd36e835e32ba`. Data from the original author are never recreated from summary percentages. The original source is in `source/Sisyrinchium_data.zip` in the conversation package; on GitHub the workflow downloads it via the pinned parent source-recovery code.

## Reproduction

```bash
python -m pip install -r requirements-lock.txt
python src/predictive.py --source source/Sisyrinchium_data.zip --outdir results
python src/verify_independent.py --root .
python -m pytest -q tests/test_predictor.py
python src/regularization_sensitivity.py --root .  # POSTHOC diagnostic
python src/plot.py --root .
```

The primary site-equal held-out-line log-loss comparison has **6 named sites, 32 left-out lines, 6,692 target-species counted individuals**. DIST outperforms BASE (0.6098 vs 0.7201 logloss). Adding concurrent other-white-species local relative abundance to DIST does not stably improve prediction (BOTH 0.6264; ridge strength sensitivity includes near-zero/positive). Neither of these predictive associations identifies ecological causality. Source observations are nested in the same original study; this is not new wild sampling or extrapolation to unseen ecosystems.
