# RGFCA Step 8D — transport-only recovery after HTTP 422

Date: 2026-09-11 JST

The first Step-8D smoke run (`34495173572`) failed before any image byte or flower-colour outcome was opened. The single URL containing all 200 frozen observation IDs returned HTTP 422 from the iNaturalist API.

A first transport-only recovery (`34495521705`) split the same exact 200 IDs into fixed groups of 50, but this also returned HTTP 422 before any metadata result or image outcome was opened.

## Root cause

Inspection of the current iNaturalist API implementation shows that `GET /v1/observations/:id` explicitly supports multiple comma-separated IDs with `validateMultiIDParam: true`. The observation `show` controller calls `setPerPage(..., {max: 200})`, then rejects the request with HTTP 422 `Too many IDs` when the number of IDs exceeds the effective `per_page` value. Because the failed Step-8D requests omitted `per_page`, the API default was 30. Both 200 IDs and 50 IDs therefore exceeded that default even though the route supports up to 200 when `per_page=200` is supplied.

The scientific IDs and endpoint were not the problem.

## Fixed technical recovery

Use the same exact first 200 species, observation IDs and expected photo IDs from the frozen Step-8C breadth ordering. Request them in one multi-ID call with `per_page=200` explicitly supplied. No anchor identity, species, outcome rule or selection order changes.

The following scientific inputs remain unchanged:

- the same first 200 species in the frozen Step-8C breadth order;
- the same exact observation IDs and expected photo IDs;
- no substitution of failed IDs;
- no image download;
- no flower-colour measurement;
- no change to the 42,111-species universe or tier allocation.

If this corrected request fails, the failure must be diagnosed from the returned HTTP body rather than further arbitrary batch-size tuning.
