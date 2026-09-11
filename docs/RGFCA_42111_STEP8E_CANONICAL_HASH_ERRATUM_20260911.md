# RGFCA Step 8E erratum — canonical content hashing for regenerated gzip products

Date: 2026-09-11 JST
Status: frozen after a pre-pixel firewall failure and before any Step-8F image pixel or flower-colour outcome was opened.

## Failure

The first Step-8F firewall attempt failed closed with `Step 8E denominator hash drift`. The Step-8E scientific denominator and transport result had not changed: 42,111 species remained in the denominator, 42,110 exact anchors remained resolved, and one anchor remained unresolved.

The drift occurred because Step 8E wrote `breadth_denominator_42111.csv.gz` and `breadth_source_manifest_resolved.csv.gz` and then used SHA-256 of the compressed file bytes as the authorization identity. Rewriting an otherwise identical gzip stream can change gzip-container bytes such as the timestamp header. Thus raw gzip-byte identity is not a stable identity for a regenerated canonical table.

No image byte or flower-colour outcome was opened before this diagnosis.

## Frozen repair

The scientific tables, row-selection rules, 42,111-species denominator, 42,110 resolved anchors, one unresolved anchor, 0.95 transport floor, and all downstream measurement rules remain unchanged.

Step 8E now records two identities for each gzip table:

1. a **canonical content SHA-256**, computed from the complete dataframe after deterministic sorting by `breadth_rank`, serialized as uncompressed UTF-8 CSV with `index=False` and `lineterminator='\n'`;
2. the raw gzip-file SHA-256, retained only as a diagnostic receipt.

Only the canonical content SHA-256 may bind Step-8F authorization. Raw gzip-byte SHA-256 may not be used as a scientific-lineage gate for a regenerated table.

## Outcome firewall

This is a technical identity repair only. It uses no image pixels, flower-colour state, ecological effect, taxonomic outcome, climate variable, literature label, or downstream inference. It cannot add, remove, replace, or reorder species based on colour results.
