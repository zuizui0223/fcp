# FCP wild-quadrat colour allocation: *Sisyrinchium* and *Boechera*

Original public dataset mirrors preserving original Dryad DOI and author-source MD5:

- Sisyrinchium Takahashi et al. (2016): `10.5061/dryad.j8q34` Zenodo record 5003124, `source/Sisyrinchium_data.zip` (CC0).
- Boechera Vaidya et al. (2018): `10.5061/dryad.q0032` Zenodo record 4993384, original individual + plot TXT and README (CC0).

**Read `docs/REPORT_JA.md` first.** This is a retrospective source-exposed exploratory reanalysis, not previously unreported biological data or proof of selection. Keep its actual field units, nested repeated timepoints, site-dependent Boechera result, and zero-target vs white distinctions.

Run from this directory:

```bash
python -m pip install -r requirements.txt
python -m pytest -q tests/
python src/field_partition.py --source-dir source --outdir results --reps 9999
python src/verify_independent.py --root .
python src/posthoc_site_fragility.py --source-dir source --outdir results --reps 9999
python src/make_figures.py
```

Frozen pre-outcome protocol committed separately at `docs/FCP_TWO_SPECIES_FIELD_MORPH_SPATIAL_PARTITION_PROTOCOL_20261013.md` on the source branch, SHA `1e7a01a876218e9d51cf071872aeb4bf032ada14`. Original source data author credits/licenses should travel with any further distribution; original Zenodo mirrors are open CC0.

Outputs include original source-authenticated quadrat ledger, per-siteyear outcomes with HOLDs, exact fixed-opportunity null replicates, independent arithmetic check, explicitly posthoc leave-one-siteyear sensitivity, and reproducible plots. Do **not** report a single positive Boechera p without its Judd Falls leave-one-site-out failure.
