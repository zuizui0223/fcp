# Two repositories, two non-overlapping scientific questions — 2026-10-14

**Status: active research-scope split; documentation-only.** No FCP biological photos, colour classifier, H1/H2 protocol, species gate, old local-depletion null, confidence interval or original result has been changed by this split. All original source history remains available in the originating Git commits and historic PRs.

## fcp: within-species photographic flower-colour variation

**Main research question:** Across plant species with sufficiently repeated colour observations, do within-species photo-colour distributions exhibit recurrent displacement directions in colour space, and is the sampled colour variation geographically structured beyond the same species' overall colour composition?

| Active within-species dataset | Frozen result / precise limit |
|---|---|
| High-depth source, 500 discovery + 500 validation + 499 prospective species | 149,900 photographs /1,499 named species; sampled photographs, not verified genetic morphs |
| D reliability | First-frozen observer-disjoint median Spearman rho 0.789; fresh-photo same-source rho 0.968 (136 overlapping species), not external spectral/biological measurement |
| Prospective white–chromatic direction | 158 eligible species, W=0.517 vs source structured-null median 0.457, p=.001; 439 total across original 152+129 and prospective 158 are **descriptive**, not 439 prospective confirmations |
| Local geographic colour structure | Species-wide minus local 50-km colour-pair diversity 0.02053/0.01867/0.01468 across three disjoint cohorts, p=.005 each; photo-site geography, not census of breeding populations |
| Exposure/artifact guard | Photographic white is associated with near-clipping OR 1.444; exclusion retained 142/158 =89.9%, under frozen 90% retention gate: status INDETERMINATE |
| Environment in high-depth sample | Small posthoc BIO5/IBE-like signals are vulnerable to observer/flower-background controls; **no universal natural-selection mechanism inferred** |

**Active manuscript:** `docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md` on FCP main; `docs/POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260918.md` and the FCP data-lineage map retain authority for H1/H2 and all source results. The manuscript may describe global photographic sampling **as a frame**, but does not borrow the global interspecific rain/pigment results as an FCP causal mechanism. *Flower-colour polymorphism* requires verified heritable floral morphs in appropriate natural populations; photographic source colour diversity by itself does not meet that genetic/fitness standard.

## mois: global *between-species* photographed floral colour × moisture

**Main research question:** Which directions of classifiable photographed floral hue correspond to long-term rainfall, vapour and joint environment **among** plant taxa worldwide? How much of the observation remains under genus, geography, strict photo quality and alternate phenotype measurement?

**Active manuscript:** `zuizui0223/mois/docs/MOIS_GLOBAL_MOISTURE_MANUSCRIPT_20261014.md`.

| Migrated source result | Original FCP authority, now copied to mois |
|---|---|
| Global photo opportunity and classifiable photo colours | 100,543 photo IDs /42,111 named taxa, 45,953 classified four photo hues, 54,590 unknown (not biologically imputed) |
| Global climate26 congeneric signed photo-white/BIO12 | 20,694 taxa/1,934 genera, +2.255 pp/source SD with within-scope maxT95 [+0.319,+4.190]; vapour +3.952 [+1.199,+6.704] |
| Strict image / local geographic cross-species limits | BIO12 strict-photo +1.229 [−2.046,+4.504], genus×same coarse cell +0.808 [−3.556,+5.172]; interval includes zero |
| Fixed support geography contrast | Matched same species East B−C +3.178 pp, 24-contrast maxT95 [−0.258,+6.613], not familywise supported |
| Externally described flower hue | Botanical labels distinct from photo classifier; same FCP site climates; primary East comparison +3.56 pp, interval [−15.34,+22.60] |
| Historical source chain | FCP PR #193 (global direction), PR #194 (identical-species spatial decomposition); source result CSVs and scripts **byte-preserving** in mois |

**Source integrity** is checked by the standalone mois workflow `.github/workflows/mois-source-audit.yml` on every change, with original SHA-1 Git blob matching and numerical receipt invariants. **Important limitation:** large historical original photo and environment ZIP archives are not physically copied into mois, so the CI is an original-source-receipt audit, not a new complete photo/environment analysis replay. Source photo archive SHA256 `2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0`; abiotic environment SHA256 `6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca`. The original full source-run link is preserved at https://github.com/zuizui0223/fcp/actions/runs/38502660130.

## Historical PR and duplicate-output handling

- **FCP main remains the same original within-species frozen paper.** Old D/H2/space-analysis source receipts and figure assets are not moved, altered, or deleted.
- **FCP PR #196**, which attempted a single combined macroenvironment + within-species Full Paper, is a historical editorial attempt **superseded by this split**, not an active New Phytologist submission candidate.
- **FCP macro PRs #193/#194 and related historical branches** remain valuable source-provenance identifiers. Their original outputs are mirrored into mois and should not be merged into FCP as a newly active interspecific paper. Historical work is not deleted to avoid breaking source citations and immutable replay chains.
- Overlap in the original metadata universe, photographic image IDs and colour-classification pipeline **does not give independent biological replication between papers**. Cross-citation can acknowledge the shared photography and contrasting question, but no morph-fitness or genotype equivalence follows.
- The NEXT-2000 metadata-only prospectively selected candidate species and 730 reserve taxa are **not** analysed new H2 observations. Original independent-opening gates remain frozen.
