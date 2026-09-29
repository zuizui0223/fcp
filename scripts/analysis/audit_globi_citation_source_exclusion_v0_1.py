#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from pathlib import Path

RESPONSE_DOIS = {
    "havens_1992": "10.1126/science.257.5073.1107",
    "hadfield_2014": "10.1086/674445",
    "ricciardi_2010": "10.1007/s10641-010-9606-0",
}
MANGAL_NAMESPACE = "globalbioticinteractions/mangal"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def norm(text: str) -> str:
    x = (text or "").strip().lower()
    x = x.replace("https://doi.org/", "").replace("http://doi.org/", "")
    x = x.replace("doi:", "")
    return x


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--citations", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    file_sha = sha256(args.citations)
    with gzip.open(args.citations, "rt", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        headers = list(reader.fieldnames or [])
        rows = list(reader)

    if not rows or not headers:
        raise RuntimeError("empty or headerless GloBI citations file")

    mangal_rows = 0
    matched_dois = {k: 0 for k in RESPONSE_DOIS}
    doi_namespaces = {k: set() for k in RESPONSE_DOIS}
    namespace_fields = []
    citation_like_fields = []

    for h in headers:
        lh = h.lower()
        if "namespace" in lh or "dataset" in lh:
            namespace_fields.append(h)
        if any(tok in lh for tok in ("doi", "citation", "reference", "source")):
            citation_like_fields.append(h)

    for row in rows:
        vals = {k: str(row.get(k) or "") for k in headers}
        joined = "\t".join(vals.values()).lower()
        ns_values = [vals[h].strip() for h in namespace_fields]
        if any(MANGAL_NAMESPACE.lower() in x.lower() for x in ns_values) or MANGAL_NAMESPACE.lower() in joined:
            mangal_rows += 1

        joined_norm = norm(joined)
        for key, doi in RESPONSE_DOIS.items():
            if norm(doi) in joined_norm:
                matched_dois[key] += 1
                for ns in ns_values:
                    if ns:
                        doi_namespaces[key].add(ns)

    all_dois_found = all(v > 0 for v in matched_dois.values())
    status = (
        "GLOBI_CITATION_EXCLUSION_MECHANISM_VERIFIED_FOR_CURRENT_MANGAL_CALIBRATION_DOIS"
        if all_dois_found
        else "HOLD_GLOBI_CITATION_METADATA_DOES_NOT_RESOLVE_ALL_CURRENT_MANGAL_DOIS"
    )

    result = {
        "version": "v0.1",
        "status": status,
        "source": {
            "zenodo_record": 22691479,
            "version": "0.11",
            "file": "citations.tsv.gz",
            "sha256": file_sha,
        },
        "schema": {
            "headers": headers,
            "namespace_fields": namespace_fields,
            "citation_like_fields": citation_like_fields,
        },
        "rows": len(rows),
        "mangal_namespace": MANGAL_NAMESPACE,
        "mangal_rows_detected": mangal_rows,
        "response_doi_match_counts": matched_dois,
        "response_doi_namespaces": {k: sorted(v) for k, v in doi_namespaces.items()},
        "all_current_response_dois_found": all_dois_found,
        "interaction_rows_opened": False,
        "taxon_identities_opened": False,
        "partner_entropy_computed": False,
        "independence_rule_if_pass": (
            "For a Mangal response system, exclude all GloBI interaction records whose dataset namespace is globalbioticinteractions/mangal OR whose citation/source DOI matches the exact response-source DOI family before computing specialization."
            if all_dois_found else
            None
        ),
        "next_gate": (
            "Freeze a row-level FILTER ONLY (namespace/source DOI) on the GloBI interaction archive, verify excluded-row counts and focal-taxon support, then compute no specialization value until that filter receipt is durable."
            if all_dois_found else
            "Remain on HOLD. Do not open GloBI interaction rows to rescue source independence until a new pre-outcome source-mapping contract is frozen."
        ),
        "claim_boundary": "This verifies metadata-level source-exclusion machinery only. It does not establish sufficient focal-taxon coverage and does not authorize partner-specialization computation.",
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
