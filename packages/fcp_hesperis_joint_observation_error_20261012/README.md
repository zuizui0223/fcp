# FCP Hesperis joint observation-error bound (2026-10-12)

This reproduces a **conditional null tail probability ceiling** for the simultaneous presence of finite colour-independent within-stand subsampling and arbitrary stand-dependent flower-colour classification errors. It does **not** estimate actual classification errors and is not a novel biological observation.

Primary Japanese report: `docs/REPORT_JA.md`. Full proof, including the *joint event* versus *unconditional null* distinction: `docs/PROOF.md`.

Reproduce with Python >=3.11; the analysis and independent verification require only standard library; focused tests require pytest:

```bash
python src/joint_bound.py --parent data/parent_entropy_result.json --outdir results
python src/verify_independent.py --root .
python -m pytest -q tests
```

Fixed input: PR173 annual Shannon result JSON with original SHA256 `980b50bfe99f69b31b57d273a2995c814e280cf552a52a1eefaf17e828b4faf1`. Earlier sibling results PR174 (sampling-only) and PR176 (error-only) are not multiplied or deemed independent. The zero-error scenario matches PR174; the full grid has 120 rows. Source-exposed retrospective exploration, no raw stand-year data, no measured error rate, no fitness or population genomic inference.
