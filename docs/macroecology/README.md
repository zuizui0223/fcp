# FCP macroecology companion paper v0.1 (2026-10-14)

**Status:** Manuscript drafting and source-claim reconciliation completed. **Not submitted and not claimed as prospectively confirmed.** This is a separate, source-overlapping companion to the current FCP within-species *New Phytologist* Full Paper; it is not a replacement of that paper's H1/H2.

## Open these first

1. `MANUSCRIPT_V0_1.md` — 5,364-word English complete research-article draft (Abstract–Conclusions–References–Tables/Figure captions).
2. `FCP_MACROECOLOGY_MANUSCRIPT_WORKING_V01.docx` — editable coauthor Word version, 16 pages, rendered and checked; attached at conversation level, **not necessarily part of Git PR source**.
3. `MAIN_TABLES.md` — real Table 1/2 and supplementary botanical measurement table, with non-comparable denominators kept separate.
4. `FIGURE_TABLE_PLAN.md` — main figures 1–4 and supplementary limits. `figures/` contains **preview** data-driven Figures 2 and 3 (SVG + PNG; not yet the final manuscript figures).
5. `EVIDENCE_AND_CLAIM_LEDGER.md` — published-vs-observational-vs-not-supported claim statuses and links to the exact PRs/results.
6. `PAPER_STRATEGY_AND_GATES.md` — Japanese publication strategy, duplication disclosure, human metadata and required future evidence.
7. `SOURCE_LINKS.json` and `NUMERICAL_RECONCILIATION.json` — controlling source map and local numeric-source check.

## Reproducibility / document tests

```
python -m pytest -q tests/test_companion_paper_claims.py
python verify_paper_numbers.py \
  --macro-ledger <PR193-source-verified-1280-row-CSV> \
  --matched-ledger <PR194-source-verified-576-row-CSV> \
  --photo-source <original-photo-archive-ZIP> \
  --abiotic-source <all36-abiotic-archive-ZIP> \
  --output NUMERICAL_RECONCILIATION.json
python make_prototype_figures.py \
  --macro-ledger <PR193-source-verified-1280-row-CSV> \
  --matched-ledger <PR194-source-verified-576-row-CSV> \
  --outdir figures
```

The first test currently passes **9/9**. The source-number audit passes **30/30** rounded manuscript claims; 1,280+576 source rows and both ZIP SHA-256s match the historical FCP replay artifacts. These are computational consistency tests, not new observations. Parent datasets/statistical results were already replayed from source under GitHub Actions #38502660130.

**Do not merge the retrospective companion into the existing active NP manuscript or main until authors approve the scientific claim boundaries and data overlap/disclosure.** Figure prototypes are not publication-final; Figure 1/4 remain to be assembled from the stated documented sources. Proper manuscript word limit, full literature and permissions review, independent genotype/fitness data and coauthor approval remain outstanding.
