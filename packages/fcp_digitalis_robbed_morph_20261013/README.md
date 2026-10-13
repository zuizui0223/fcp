# Independent Digitalis field-flag analysis (2026-10-13)

Source: Lozada Gobilard et al. (2026), *AJB*, DOI 10.1002/ajb2.70186. Figshare original raw data DOI 10.6084/m9.figshare.29095451.v2, file 54632993. CC BY 4.0, original author attribution required.

This is an exploratory analysis of recorded `Robbed_flower` yes/no flags in selected 2024 flowers clustered by recorded plant/patch (`COUNTRY`, `Year`, `GPS-code`, `Plant_ID`) in five original study locations. It is not a census of population morph frequencies and cannot estimate fitness or pollinator choice. Manuscript H1/H2 and previous FCP outputs unchanged.

Run:
```bash
python -m pip install numpy==2.3.5 pandas==2.2.3 pytest==9.0.2
python -m pytest -q tests/test_robbery.py
python src/robbed_morph.py --source source/digitalis_field_individuals.csv --outdir results
python src/verify_independent.py --source source/digitalis_field_individuals.csv --outdir results
```

The program aborts unless original source bytes have length **345652**, MD5 **a1200e0197b282ee126700ae234eb95e** and SHA256 **723c7a03b6ffa010e6093996d8203b6770ab1ee1bc6ddc21ec9c3d5b47d3b8c9**. They decode with CP1252, not necessarily UTF-8. If the raw source file was not packaged, download precisely `https://ndownloader.figshare.com/files/54632993` and verify hashes before use. The user package contains original authenticated Figshare source with attribution and logs.

`docs/REPORT_JA.md` is the interpretation. `results/plant_patch_ledger.csv` is a key-and-QC ledger, not true census plants. `results/site_primary_contrasts.csv` is the per-site estimate. `results/result.json` contains the 9999 permutations, 1999 bootstrap samples and 2024 date-controlled sensitivity; raw replicate draws are not retained. `results/independent_verification.json` checks original source via separate CSV/Fraction implementation.
