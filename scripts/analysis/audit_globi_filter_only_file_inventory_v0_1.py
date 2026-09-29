#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from pathlib import Path

PREFERRED=("interactions.parquet","interactions.tsv.gz","interactions.tsv")

def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--record-json",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    rec=json.loads(a.record_json.read_text())
    files=rec.get("files") or []
    rows=[]
    for f in files:
        key=str(f.get("key") or f.get("filename") or "")
        rows.append({
          "key":key,
          "size":f.get("size"),
          "checksum":f.get("checksum"),
          "self":(f.get("links") or {}).get("self"),
          "content":(f.get("links") or {}).get("content"),
        })
    by={r["key"]:r for r in rows}
    selected=None
    for name in PREFERRED:
        if name in by:
            selected=by[name]; break
    metadata=rec.get("metadata") or {}
    status="PASS_GLOBI_FILE_INVENTORY_INTERACTION_TABLE_PINNED" if selected else "HOLD_GLOBI_FILE_INVENTORY_NO_INTERACTION_TABLE"
    out={
      "version":"v0.1",
      "status":status,
      "zenodo_record":rec.get("id"),
      "record_version":metadata.get("version"),
      "published":metadata.get("publication_date"),
      "title":metadata.get("title"),
      "file_count":len(rows),
      "files":rows,
      "selection_priority":list(PREFERRED),
      "selected_interaction_table":selected,
      "biological_rows_opened":False,
      "taxon_identities_opened":False,
      "specialization_values_computed":False,
      "next_gate":"Freeze schema-only opening of the pinned interaction table and verify required source/taxon/type fields before any support counts." if selected else "STOP"
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps(out,indent=2,sort_keys=True))
    return 0

if __name__=="__main__":
    raise SystemExit(main())
