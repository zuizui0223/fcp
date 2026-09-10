# RGFCA Step 8D — batch resolution of frozen discovery anchors

Date: 2026-09-11 JST
Status: frozen before any Step-8D image pixel is downloaded.

## Purpose

Recover current iNaturalist photo metadata and large-image URLs for the exact V1/V2 observation/photo IDs already selected outcome-blind in Step 8C. This is metadata resolution only; returned image URLs must not be opened by the resolver.

## API strategy

iNaturalist recommends fetching multiple observation records by concatenated observation IDs. Use batches of at most 200 IDs per request against `https://api.inaturalist.org/v1/observations/<comma-separated-ids>` with a project-identifying User-Agent and >=1.05 seconds between requests.

## Smoke gate

Before resolving the full frame, use exactly the first 200 species in the frozen Step-8C breadth ordering. For each expected observation/photo pair:

- the requested observation must be returned;
- the expected photo ID must still be attached to that observation;
- the photo must expose a URL and licence metadata sufficient for the existing measurement acquisition layer;
- URL strings are recorded but not fetched.

No substitute photo is allowed in the smoke result. Missing/deleted records are recorded as unresolved.

## Full resolver after smoke

If the batch endpoint itself is operational, resolve the complete 42,111 breadth anchors in deterministic batches of <=200. Do not change anchor identities from resolution outcomes. Full resolution therefore requires about 211 requests, rather than one request per species.

## Boundary

This step opens API metadata only. It does not download image bytes, run the flower detector, measure colour, or change the species/anchor selection.
