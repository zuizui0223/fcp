#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

EXPECTED_SHA="6b5191180e9bbfdb7867ff2fdef3e2e305f78cfe6f6405209e455497a06e30c4"

def sha256(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""):
            h.update(chunk)
    return h.hexdigest()

def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--datasets",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()

    observed=sha256(a.datasets)
    if observed!=EXPECTED_SHA:
        raise RuntimeError(f"datasets.tsv SHA mismatch: {observed}")

    with a.datasets.open(newline="",encoding="utf-8") as f:
        rows=list(csv.DictReader(f,delimiter="\t"))
    if not rows or "namespace" not in (rows[0].keys() if rows else []):
        raise RuntimeError("GloBI dataset index missing namespace")

    namespaces=sorted({(r.get("namespace") or "").strip() for r in rows if (r.get("namespace") or "").strip()})
    patterns=[r"mangal",r"mangal\.io",r"mangal-interactions"]
    mangal_like=sorted({ns for ns in namespaces if any(re.search(p,ns,re.I) for p in patterns)})

    result={
      "version":"v0.1",
      "status":"GLOBI_SPECIALIZATION_NAMESPACE_PREFLIGHT_COMPLETE",
      "source":{
        "zenodo_record":22691479,
        "version":"0.11",
        "dataset_index_sha256":EXPECTED_SHA
      },
      "dataset_rows":len(rows),
      "unique_namespaces":len(namespaces),
      "mangal_like_namespaces":mangal_like,
      "mangal_like_namespace_count":len(mangal_like),
      "interaction_rows_opened":False,
      "taxon_identities_opened":False,
      "specialization_values_computed":False,
      "independence_decision":(
        "PARTIAL_NAMESPACE_EXCLUSION_POSSIBLE_SOURCE_DOI_EXCLUSION_STILL_REQUIRED"
        if mangal_like else
        "HOLD_NO_EXPLICIT_MANGAL_NAMESPACE_SOURCE_DOI_EXCLUSION_REQUIRED"
      ),
      "next_gate":"Before any Mangal-response specialization predictor is computed, open only the minimum GloBI citation/source fields needed to exclude the exact Mangal dataset/source DOI family. Do not compute partner entropy in the same step.",
      "claim_boundary":"Dataset namespaces alone cannot prove complete response-source independence because the same underlying study may be indexed under another namespace."
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,indent=2,sort_keys=True))
    return 0

if __name__=="__main__":
    raise SystemExit(main())
