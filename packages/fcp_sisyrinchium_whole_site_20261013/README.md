# FCP: Water-edge-distance versus congener prevalence — whole-site transfer audit

This is an **out-of-site stress test**, not an independent new field study or a mechanism test. It uses exactly the previous PR185 primary 168 quadrats / 6 named sites and compares the old *leave-one-line-out within the same site* validation with the new *leave-one-entire-site-out* validation. Site name and all its years/lines must be excluded from training in the latter.

**Main result:** The heldout-whole-site model with water-edge distance alone did not consistently improve over a source-site-mean baseline (equal-site log loss 0.75676 vs baseline 0.71303), in contrast to the strong within-site line holdout (0.60977 vs 0.72011). A contemporaneous observed congener-abundance covariate was the lowest-loss single predictor of the four fixed candidates on the selected six sites (0.67258), but it is NOT evidence of a causal ecological interaction. Failure of water distance is especially sensitive to the heldout Shiga site (37.88% of its M individuals outside the range of training-site distances); balancing total training weights across other source sites changes the sign of the mean distance improvement. Thus source-conditioned regional transfer remains **unestablished**.

The source is Takahashi, Takakura & Kawata (2016), PLOS ONE DOI `10.1371/journal.pone.0164381`, original Dryad DOI `10.5061/dryad.j8q34`, Zenodo mirror 5003124. Original-source 512 field-quadrat archive and parent PR185 frozen result are retained **unaltered**. See `PROTOCOL.md` (frozen GitHub commit before new scores), `docs/REPORT_JA.md`, `results/result.json` and source QA.

To run from the unpacked conversation package:

```bash
python -m pip install -r requirements-lock.txt
python -m pytest -q tests
python src/transfer.py --archive parent/source/Sisyrinchium_data.zip --parent parent/results/result.json --parent-code parent/src/predictive.py --outdir results
python src/verify_independent.py --archive parent/source/Sisyrinchium_data.zip --parent parent/results/result.json --result results/result.json --out results/independent_verification.json
python src/make_figures.py --root .
```

The package reuses the **exact original-source checksum gates** and the same fitted likelihood from PR185 while withholding all field lines of the evaluated source site. It supplies a separate original-source reader and independent Newton optimizer. The full FCP repository pytest suite is not equivalent to the 12 focused tests. All outputs are *retrospective, source-exposed exploratory*, not externally preregistered.
