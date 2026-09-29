#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

GEOMETRY_TOKENS=("lat","lon","location","geo","site","coordinate","spatial")
TIME_TOKENS=("date","time","year","temporal","event")


def sha256(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""):
            h.update(chunk)
    return h.hexdigest()


def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--datasets",type=Path,required=True)
    ap.add_argument("--design",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()

    d=json.loads(a.design.read_text())
    expected=d["source"]["dataset_index_sha256"]
    observed=sha256(a.datasets)
    if observed!=expected:
        raise RuntimeError(f"datasets.tsv SHA256 mismatch: {observed}")

    with a.datasets.open(newline="",encoding="utf-8") as f:
        r=csv.DictReader(f,delimiter="\t")
        header=r.fieldnames or []
        rows=list(r)

    lower=[x.lower() for x in header]
    geometry_fields=[h for h,l in zip(header,lower) if any(t in l for t in GEOMETRY_TOKENS)]
    time_fields=[h for h,l in zip(header,lower) if any(t in l for t in TIME_TOKENS)]
    identity_fields=[h for h,l in zip(header,lower) if any(t in l for t in ("namespace","dataset","citation","url","doi","id"))]

    namespace_field=next((h for h in header if "namespace" in h.lower()),None)
    namespaces=[]
    if namespace_field:
        namespaces=sorted({str(row.get(namespace_field,"")).strip() for row in rows if str(row.get(namespace_field,"")).strip()})

    gate=(
        len(rows)>=d["pass_rule"]["minimum_datasets"]
        and (bool(identity_fields) if d["pass_rule"]["require_dataset_identity_field"] else True)
        and (bool(geometry_fields or time_fields) if d["pass_rule"]["require_any_geometry_or_time_field_in_dataset_index"] else True)
    )
    status="PASS_GLOBI_DATASET_INDEX_GEOMETRY_METADATA" if gate else "HOLD_GLOBI_DATASET_INDEX_LACKS_GEOMETRY_TIME_METADATA"

    result={
      "version":"v0.1",
      "status":status,
      "source_record":d["source"]["zenodo_record"],
      "source_version":d["source"]["version"],
      "dataset_index_sha256":observed,
      "dataset_rows":len(rows),
      "header":header,
      "identity_fields":identity_fields,
      "geometry_fields":geometry_fields,
      "time_fields":time_fields,
      "namespace_field":namespace_field,
      "namespace_count":len(namespaces),
      "namespace_set_sha256":hashlib.sha256("\n".join(namespaces).encode()).hexdigest(),
      "interaction_rows_opened":False,
      "taxon_identities_opened":False,
      "biological_outcome_computed":False,
      "full_interaction_outcome_opening_authorized":False,
      "next_gate":(
        "freeze dataset-level geometry/time candidate selection before any interaction row opening"
        if gate else
        "use original-source or dataset-cache metadata only; integrated interaction rows remain sealed"
      )
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,indent=2,sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
