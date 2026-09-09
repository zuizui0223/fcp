# RGFCA surviving flower-side candidate leave-one-genus audit — result

Date: 2026-09-09 JST
Protocol freeze: `19168f1a7bc8fc01b410108855ba522952d58fb9`.

The sole candidate that survived both species-disjoint and genus-disjoint partitions, `(-23,11,M2+)`, was audited by deleting each represented genus in turn within every supported recurrence realization.

Across all 200 realizations:

- the full-cell equal-species flower-only M2 score was positive in **95%** of realizations;
- **87%** of realizations remained positive after **every** eligible single-genus deletion;
- median support was **19 species from 18 genera**;
- the median minimum leave-one-genus-out score was `+0.01704`;
- its 2.5–97.5% range was `-0.01487` to `+0.03910`;
- median maximum absolute influence of any single genus was `0.00954` (97.5% `0.02071`).

The predeclared single-genus robustness gate (>=0.70 of supported realizations preserving the positive sign under every deletion) therefore **passed**. No single genus can readily account for the recurrent positive M2 signal, although some realizations remain sensitive and the lower-tail minimum LOO score crosses zero.

All five output CSVs were reproduced byte-for-byte in a complete rerun; upstream recurrence support matched to `4.44e-16`.

## Claim boundary

This is still a post-discovery influence audit. The result supports within-reserve taxonomic robustness only. It is not an independent-data confirmation, a phylogenetic analysis, a pollinator mechanism, or proof of a biological transition boundary.
