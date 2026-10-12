# FCP Linanthus: out-of-year geographic frequency transfer

Source-exposed ecological field-count reanalysis (2026-10-13) based solely on 28 published Appendix rows of Schemske & Bierzychudek (2007), doi:10.1111/j.1558-5646.2007.00219.x. Four fixed large sections within **one** ravine site, seven years. Not genetic F_ST, local 2m subplot polymorphism, local adaptation test, independently replicated species, or a prospective prediction.

## Run

```bash
python -m pytest -q tests/test_spacetime.py
python scripts/field_spacetime.py --source source/linanthus_counts.csv --outdir results
python scripts/verify_independent.py --source source/linanthus_counts.csv --result results/result.json --out results/independent_verification.json
python scripts/make_figures.py        # optional matplotlib; the core analysis is Python stdlib
```

Expected input SHA256 `fa408c7ad6e521a8410252d903b024a057b1ee4f7e27406927e77fcc22bf4ac7`.

## Read first

- `docs/PROTOCOL.md` — eligibility, estimands, leakage guards, inference boundaries
- `docs/REPORT_JA.md` — biological interpretation including negative/untestable cases
- `results/result.json` — full fixed numerical record
- `results/heldout_section_model_predictions.csv` — 60 individual heldout estimates
- `results/yearly_side_summary.csv` — actual white fractions and within-site side contrast
- `docs/public_population_dataset_discovery.csv` — independently published sources and **real** raw-byte-recovery status
- `results/independent_verification.json` — exact Fraction reimplementation independent of the main calculator

The only included raw observational file is the manually transcribed **published aggregate Appendix table**, not new plant records. Percentages are calculated from integer counts; its inconsistent 1992 south-west printed white percentage remains flagged. The source does not provide person or plot coordinates or individual plant IDs.
