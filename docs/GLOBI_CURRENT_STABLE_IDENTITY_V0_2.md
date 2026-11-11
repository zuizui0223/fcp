# GloBI current stable identity — v0.2

The GloBI route is now split into two distinct source roles.

1. **Namespace discovery** used the versioned review-corpus `datasets.tsv` snapshot already audited in v0.1. It contained 458 namespaces but no geometry/time columns, so that stage correctly terminated as a schema HOLD without opening interaction rows.
2. **Any eventual biological interaction analysis** must use a pinned stable integrated GloBI release or, preferably, the exact original dataset selected through the GloBI source gate.

As checked on 2026-09-30, the official GloBI data page still directs research use to the stable versioned concept DOI `10.5281/zenodo.3950589`. The latest stable integrated association record resolved in the current audit is Zenodo **22761701**, version **v13**, published **2026-09-15**, with primary file `globi_assoc.tar.gz` (reported 889.1 MB; MD5 `de999239cdbc5ba5eb756064c4f3339a`).

The earlier review-corpus record remains useful only to enumerate dataset namespaces. It is not promoted as the primary interaction-response data source.

No interaction row or generalized turnover outcome is opened by this identity update.
