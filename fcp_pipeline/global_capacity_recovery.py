"""Transport-only recovery primitives for the global metadata capacity census.

This module is deliberately biological-outcome agnostic.  It may only be used
when the frozen capacity scan is not evaluable because requests failed with the
explicit iNaturalist HTTP 429 normal-throttling transport error.  Successful
capacity rows are never re-sampled.
"""
from __future__ import annotations

import pandas as pd

RETRYABLE_429 = "HTTPError:HTTP Error 429: normal_throttling"


def retryable_429_mask(audit: pd.DataFrame) -> pd.Series:
    if "request_error" not in audit.columns:
        raise ValueError("capacity audit lacks request_error")
    return audit["request_error"].fillna("").astype(str).eq(RETRYABLE_429)


def nonempty_error_mask(audit: pd.DataFrame) -> pd.Series:
    if "request_error" not in audit.columns:
        raise ValueError("capacity audit lacks request_error")
    return audit["request_error"].fillna("").astype(str).str.len().gt(0)


def frozen_recovery_rows(audit: pd.DataFrame) -> pd.DataFrame:
    """Return every and only explicit 429 row in deterministic global order."""
    required = {"global_row_index", "species", "inat_taxon_id", "request_error"}
    missing = sorted(required - set(audit.columns))
    if missing:
        raise ValueError(f"capacity audit lacks recovery columns: {missing}")
    work = audit.copy()
    work["global_row_index"] = pd.to_numeric(work["global_row_index"], errors="raise").astype(int)
    if work["global_row_index"].duplicated().any() or work["inat_taxon_id"].duplicated().any():
        raise ValueError("capacity audit has duplicate global/taxon identities")
    out = work.loc[retryable_429_mask(work)].copy()
    return out.sort_values("global_row_index", kind="mergesort").reset_index(drop=True)


def merge_transport_recovery(
    original: pd.DataFrame,
    recovered: pd.DataFrame,
) -> pd.DataFrame:
    """Replace only original 429 rows with their exactly-once recovery rows."""
    required = {
        "global_row_index", "species", "inat_taxon_id", "request_error",
        "raw_results", "locally_eligible", "prior_excluded", "wrong_taxon",
        "after_observer_cap", "maximum_span_km", "capped_row_id_sha256",
        "eligible_raw_100", "eligible_raw_80", "eligible_raw_60",
    }
    for name, frame in (("original", original), ("recovered", recovered)):
        missing = sorted(required - set(frame.columns))
        if missing:
            raise ValueError(f"{name} capacity rows lack columns: {missing}")
    base = original.copy()
    base["global_row_index"] = pd.to_numeric(base["global_row_index"], errors="raise").astype(int)
    retry = frozen_recovery_rows(base)
    got = recovered.copy()
    got["global_row_index"] = pd.to_numeric(got["global_row_index"], errors="raise").astype(int)
    if got["global_row_index"].duplicated().any() or got["inat_taxon_id"].duplicated().any():
        raise ValueError("recovery output contains duplicate identities")
    expected_keys = set(zip(retry["global_row_index"].astype(int), retry["inat_taxon_id"].astype(int)))
    got_keys = set(zip(got["global_row_index"].astype(int), got["inat_taxon_id"].astype(int)))
    if got_keys != expected_keys or len(got) != len(retry):
        raise ValueError("recovery rows do not cover exactly the frozen 429 set")
    if len(got):
        original_names = retry.set_index("global_row_index")["species"].astype(str)
        got_names = got.set_index("global_row_index")["species"].astype(str)
        if not original_names.sort_index().equals(got_names.sort_index()):
            raise ValueError("recovery species labels differ from original 429 rows")
    replacement = got.set_index("global_row_index")
    out = base.set_index("global_row_index")
    for idx in replacement.index:
        for column in required - {"global_row_index"}:
            out.at[idx, column] = replacement.at[idx, column]
    out = out.reset_index().sort_values("global_row_index", kind="mergesort").reset_index(drop=True)
    if len(out) != len(base) or out["inat_taxon_id"].nunique() != len(base):
        raise RuntimeError("transport recovery changed the census denominator")
    return out


__all__ = [
    "RETRYABLE_429",
    "frozen_recovery_rows",
    "merge_transport_recovery",
    "nonempty_error_mask",
    "retryable_429_mask",
]
