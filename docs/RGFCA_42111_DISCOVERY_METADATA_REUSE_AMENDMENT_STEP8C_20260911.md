# RGFCA Step 8C amendment — reuse of pre-outcome V1/V2 discovery metadata

Date: 2026-09-11 JST
Status: frozen before any new Step-8C image pixel is opened.

This amendment applies only to the new 42,111-species Step-8C breadth/depth program. It does not modify, rerun, or reinterpret the earlier 1,000-species RGFCA discovery/reserve experiments.

## Rationale

The V1 and V2 global discovery observation indexes were frozen as metadata-only products before their image pixels or flower-colour outcomes were opened. They therefore provide an already outcome-blind pool of exact observation and photo IDs for all 42,111 discovered species.

For Step 8C, these pre-outcome IDs may be reused as the first candidate source instead of discarding them and forcing a new API draw solely for procedural novelty. This reduces API load while preserving outcome-blind selection.

## Fixed reuse rule

1. Combine V1 and V2 observation indexes and deduplicate exact observation/photo pairs.
2. For each species, rank candidate pairs by SHA256(`20260911|inat_taxon_id|observation_id|photo_id`).
3. Take the first `combined_target_depth` candidates up to availability.
4. For a fresh-validation surface, first remove all observation/photo IDs recorded in prior opened/exclusion ledgers, then apply the same hash ranking.
5. Any remaining target shortfall may be filled only by a separately frozen metadata query executed before new colour outcomes are opened.
6. No colour, image content, climate, C*/S*, D, significance, family or literature evidence may affect candidate ordering, target depth, or replacement.

## Breadth anchor

Every one of the 42,111 discovered species has at least one V1/V2 observation/photo pair. Therefore the descriptive breadth surface may retain all 42,111 species with one pre-outcome discovery anchor, including species whose later independent capacity scan returned zero currently eligible records.

A capacity-scan zero is interpreted as failure of that later query to recover an eligible candidate, not proof that the species had no eligible discovery photograph.

## Independence labels

- IDs already opened by an earlier colour experiment remain `previously_opened` and cannot constitute fresh validation.
- Unopened V1/V2 IDs remain prospectively usable for Step 8C because their identity was fixed without colour.
- Descriptive pooled summaries must preserve opened/fresh status explicitly.

## Failure handling

Failure to resolve or download a frozen candidate is terminal for that candidate. It is not replaced after seeing image or colour outcomes. Any technical recovery must be specified before inspection of the affected colour outcome and must preserve the original species denominator.
