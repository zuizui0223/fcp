# Figure QC — FCP polymorphism paper v0.1

Date: 2026-09-11 JST

## Status

**PASS — final 168-mm visual inspection completed.** All four main figures render as both PNG and vector PDF from the synchronized `paper_numbers.json` and frozen figure-data tables. A 168-mm, 300-dpi preview set was generated after the final renderer changes and every main figure was inspected at that target width for text legibility, clipping, panel-title collisions, tick-label collisions, legend placement and annotation overlap.

The final renderer uses lowercase parenthesized panel labels `(a)`, `(b)`, … and a 9.5-pt source typography floor, providing margin over an 8-pt final-size requirement. No data values, statistical results, axis limits or inferential decisions were changed during visual QC.

## Final visual inspection

- **Figure 1 — PASS.** Panel (a) uses a vertical fixed-photo → measurement gate → four colour states → species-diversity flow, avoiding the crowding seen in the earlier horizontal layout. Panels (b–d) are legible at 168 mm, with the discovery-frame prevalence caveat visible and no panel-title collision.
- **Figure 2 — PASS.** Discovery and species-disjoint reserve remain visually distinct. Panel (c) uses abbreviated `Disc.` / `Res.` tranche labels to prevent tick-label contact. Panel (d) nuisance-control labels and P-value annotations remain readable without clipping.
- **Figure 3 — PASS.** The measurement-status partition, ambiguity ECDF and completion interval are legible. Panel (c) retains the three D variants but removes redundant pointwise P-value text; exact P values remain in Results and the frozen figure-data table. This resolves the only remaining annotation-density problem without altering the statistics.
- **Figure 4 — PASS.** The secondary genus result remains visually qualified by the unsupported Dmin4 endpoint. The six deliberately labelled genera in panel (b) use fixed point offsets; the previous `Salvia`/`Ipomoea` overlap is removed and the cross-tranche concordance remains readable.

## Claim-to-figure contract

- Figure 1: species-level polymorphism estimand and descriptive distribution only.
- Figure 2: replicated association between polymorphism and within-species geographic organization.
- Figure 3: sampled-opportunity, technical-failure, matched-background and ambiguity-endpoint robustness of the spatial association.
- Figure 4: qualified genus-level taxonomic clustering; not formal phylogenetic signal.

No main figure is permitted to imply adaptation, a universal climate threshold, a shared global boundary, or that ambiguous observations are known hidden morphs.

## QC provenance

Earlier renderer audits identified small text, long panel titles and overlapping annotations. Those issues were resolved through layout-only passes, followed by JBI-oriented typography/panel-label normalization and two 168-mm artifact inspections. The final 168-mm preview set was produced by the dedicated final visual-QC workflow after the Figure 2(c) tranche-label and Figure 4(b) genus-label placement fixes. The science bundle was revalidated against the regenerated figure hashes before the final figures were committed.
