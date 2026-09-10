# RGFCA Step 8D — transport-only recovery after HTTP 422

Date: 2026-09-11 JST

The first Step-8D smoke run (`34495173572`) failed before any image byte or flower-colour outcome was opened. The single URL containing all 200 frozen observation IDs returned HTTP 422 from the iNaturalist API.

This is treated as a transport/request-shape failure only. The following scientific inputs remain unchanged:

- the same first 200 species in the frozen Step-8C breadth order;
- the same exact observation IDs and expected photo IDs;
- no substitution of failed IDs;
- no image download;
- no flower-colour measurement;
- no change to the 42,111-species universe or tier allocation.

The technical recovery splits the same ordered 200 observation IDs into four fixed batches of 50 IDs and requests each batch through the same documented multi-ID observation endpoint. Requests are paced at no faster than approximately one request per second. The four responses are combined before the same exact observation/photo identity checks are applied.

This recovery changes only request transport granularity. It cannot select or replace anchors from response content.
