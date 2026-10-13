# Sisyrinchium M-colour × other white-flowered species, within original field lines

This is a **retrospective reanalysis** of Takahashi, Takakura & Kawata (2016), *PLOS ONE*, DOI `10.1371/journal.pone.0164381`, archived in Dryad DOI `10.5061/dryad.j8q34` and publicly mirrored in Zenodo record `5003124`. The original authors already reported the population/quadrat association of other white species and M-species purple flowers, as well as collinearity with shoreline distance. We do **not** claim this positive association is a new ecological phenomenon, or that its cause has been proven.

The new, tighter question: do the original counts show (a) an association **inside each recorded line of quadrats**, when M-purple/M-white totals are fixed within every line, and (b) residual flower-colour aggregation *inside* those same lines? Source-exposed protocol committed at `docs/FCP_SISYRINCHIUM_WITHIN_LINE_COGENERIC_EXPOSURE_PROTOCOL_20261013.md` **before new endpoints were calculated** (but after previous source/results were seen). Main/H1/H2 and old PR outputs unchanged.

## Original source

* `source/Sisyrinchium_data.zip`: original data.zip, 127502 bytes, MD5 `008ec196aa5acb769dbfd36e835e32ba`. Public source `https://zenodo.org/api/records/5003124/files/data.zip/content`.
* `source/frequency.txt`: original ZIP member `data/frequency.txt`, SHA256 `a12cdf46483a6d8b78bd00110c196a25f3c8849967ae8a2c37d6aee3d29f516e`. 512 field quadrats, 16 siteyears at 13 named populations. M-species violet `v=5081`, white `w=8723`; *other species* L `l=4463`, S `s=29`—**do not pool with white M**. 126 quadrats had no target M plants.
* `humidity` missing in 391/512 and semantically mixed; no imputed humidity effects. `line` is source's recorded line label; groups may contain a **single quadrat** and cannot then support within-line inference.

## Main findings (original site/year sources only)

* **Purple minus white M morph difference in L/S individual share WITHIN same recorded line**: +0.12554956 (12.55 percentage points), 8 eligible site-years / 7 named sites, all eight positive, one-sided 9999 fixed-inventory conditional shuffles **p=0.0001**, Holm family-of-two p=0.0002. Data source original 16 site-years, not eight independently evolved populations. The comparator is not a genetic or fitness measure. Latest-year-only (six named sites) +0.09112136, p=0.0001.
* **Within-line quadrat clumping**: 11 eligible site-years/10 named sites; mean original site-level B=0.44947, line-level between=0.19375, exact random-label expected within-line=0.03801. Excess after fixing line morph totals = **0.21770757**, i.e. +21.77 percentage points; p=0.0001, Holm family-of-two p=0.0001. All 11 siteyears positive; each leave-one-out still positive. Latest-year-only eight sites +0.25072, p=0.0001. This is a phenotype segregation benchmark on the source quadrat scale, **not** genetic F_ST or causal selection.
* **Water-distance sensitivity**: within-line L/S share residual after fitting n-weighted intercept+linear distance has average difference **+0.00995631** (+1.00 pp), nominal one-sided p=0.0378. **Posthoc leave Banpaku Y13 out → -0.00058275, p=0.5548**, so not robust. Comparing exact same line subset before and after projection, unadjusted +0.13304 and adjusted +0.00996 (92.52% numerical attenuation), but this is **not a causal mediation fraction**. Only five of eight siteyears have positive residual difference.

## Run from the directory containing this README

```
python -m pip install -r requirements-lock.txt
python src/recover_original_sources.py --outdir source
python -m pytest -q tests
python src/within_line.py --source source/frequency.txt --outdir results --reps 9999
python src/verify_independent.py --source source/frequency.txt --result results/result.json --outdir results
python src/posthoc_matched_distance.py --source source/frequency.txt --result results/result.json --outdir results
python src/make_figures.py --result results/result.json --out docs
```

Source recovery is strict: check original ZIP MD5, byte size, ZIP CRC, internal member SHA256. No AppleDouble metadata parsed. Outputs are versioned and source-frozen; CI reruns all outputs and independent checks before making claims. For scientific details and limitations see `docs/REPORT_JA.md` and `docs/PROOF.md`.
