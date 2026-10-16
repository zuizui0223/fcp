# FCP × mois: separate inferential mainlines from related photographic source data

Date: 2026-10-16. **Active FCP paper: within-named-species photographed flower-colour variation.**

The October 14 repository split is retained. This document and `study_scope.json` add an **executable boundary**, not a new numerical reanalysis or an alternative manuscript. The original frozen H1/H2 and post hoc results are not changed.

## Two frozen populations, not one common denominator

| | fcp (this repository) | mois (other repository) |
|---|---|---|
| Question | How is photographed colour variation **within a named species** arranged through geography and phenotype space? | How do photographed species-level white and other hue frequencies vary **among species and genera** with global rainfall, vapour, and other environment? |
| Primary fixed image support | **149,900** photographs, 1,499 species, discovery 500 / validation 500 / prospective 499 | **100,543** unique originals, 42,111 named taxa, 45,953 hue-classifiable photos |
| Primary unit | Conspecific observations or colour pairs conditioned on species, plus species-specific W | Named-species photographic hue composition with genus and geographic conditioning |
| Main response | `D`, local-versus-specieswide colour-pair depletion, construction-controlled `W` | Four photographic colour proportions, particularly white-versus-rain/vapour regression directions |
| Primary test | Species-composition-preserving spatial null and prospectively frozen H2 null | Within-genus slopes with whole-genus bootstrap/maxT and fixed-support spatial/photo-quality tests |

Both pipelines originate from the same iNaturalist source and measurement family. A checksum-verified original-photo-ID audit (FCP workflow run **38911659260**, artifact **12049893936**) identified **3,411 exact matching photo IDs**, representing **2.28%** of FCP's 149,900 and **3.39%** of MOIS's 100,543 original photos. Cohort-specific overlaps were discovery **1,248/50,000**, validation **1,326/50,000**, prospective **837/49,900**; FCP's three frozen cohorts themselves have **zero shared photo IDs**. Of 1,693 common photos classifiable in both frozen outputs, 1,693 had the same coarse four-hue label (same source system, **not** independent validation). One matched photo has discordant recorded taxon IDs; MOIS original two-ledger union contains one duplicate photo row with differing incomplete classification, both retained as provenance warnings. Neither paper is independent-source replication of the other. A shared upstream data source is fine: the biological estimand and evaluation population, not repository or photo duplication, define the analyses.

## FCP primary claims and exclusions

1. Observer-disjoint and fresh-photo within-system reproducibility of sampled species-wide diversity `D`; no global genetic-polymorphism prevalence estimate.
2. Local geographic partitioning: at the fixed post hoc 50-km scale, conspecific colour pairs are less mixed than a null preserving each species' exact sampled colour composition (166, 181, and 204 species across three cohorts; respective mean depletion ~0.02053, 0.01867, 0.01468). This does not identify dispersal, selection, genotype, or local breeding-population composition.
3. Prospective confirmation of an achromatic/chromatic displacement **in photo phenotype space**, `W=0.517246` vs frozen coarse-morph-preserving null median `0.457143`, 158 eligible species, p=.001. The 439 species across three cohorts are descriptive, **not** 439 prospective confirmations.
4. BIO5 and IBE-like follow-ups remain bounded by observer/background/photographic controls and must not be presented as adaptive temperature mechanisms.
5. **Excluded:** importing the global congeneric white–BIO12 coefficient from mois as evidence for a within-species physiological effect, fitness balance, universal local coexistence, or pigment chemistry.

The older JBI literature-derived within-population `C` / geographically structured `S` evidence is retained as archived comparative provenance in this repository. It is not equivalent to photographed high-depth `D`, and absence of evidence for `C` is not evidence of monomorphism. The current New Phytologist manuscript remains `docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md`.

## Next FCP-specific analysis, independently of mois

For each well-sampled species, make within-species geographic `D_{s,g}`, expected photo-colour pair diversity and colour-vector turnover conditional on its own range-wide colour counts. Stratify or adjust for observers, photo colour classifiability/near-clipping, observed flowering dates, spatial sampling opportunity and site distance. Test nulls that preserve species membership and, where appropriate, observer/site/season sampling blocks. Do not use a between-genus coefficient as a within-species effect.

The **season dimension here** is the temporal heterogeneity of **conspecific morph/photocolour arrangements**, not global monthly white-photo abundance. `C` at actual breeding-population level will require field/individual repeat verification; photo proximity alone is not verified morph coexistence.

## Engineering gate

`python scripts/check_study_scope.py` reads the immutable H2 and 50-km receipt files, verifies this inferential level, the fixed original counts and absence of the migrated mois primary result/script paths. CI executes it separately from the repository layout regression tests. It **does not** refit source-photo models, recompute p-values, or resolve the unknown cross-paper photo-ID overlap. The script deliberately guards primary result locations, not historical archived source references. It also requires the frozen cross-paper original-photo overlap receipt at `results/fcp_mois_original_photo_id_overlap_20261016.json`. Audit source code and SHA-locked original-photo workflow live at `scripts/audit_fcp_mois_photo_id_overlap.py` and `.github/workflows/fcp-mois-original-photo-id-overlap.yml`.

See the complement in `zuizui0223/mois/docs/STUDY_SCOPE_FCP_MOIS_20261016.md`.
