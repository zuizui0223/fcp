#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path

def hset(values):
    return hashlib.sha256("\n".join(sorted(set(values))).encode()).hexdigest()

def valid(v):
    return str(v or "").strip() not in {"","NA","N/A","null","None"}

def read_projected(path,allowed):
    with path.open(newline="",encoding="utf-8-sig",errors="replace") as f:
        r=csv.reader(f)
        header=next(r)
        idx={name:header.index(name) for name in allowed if name in header}
        missing=sorted(set(allowed)-set(idx))
        if missing:
            raise RuntimeError(f"{path} missing safe columns {missing}")
        for row in r:
            yield {k:(row[i] if i<len(row) else "") for k,i in idx.items()}

def network_key(row,cols):
    return "|".join(str(row.get(c,"")).strip() for c in cols)

def audit_carni(path,d):
    rows=list(read_projected(path,d["safe_columns"]))
    focal=d["focal_taxon_column"]; keycols=d["local_unit_key"]
    taxa=set(); units=set(); tax_units=defaultdict(set); geo_units=set(); method_units=set()
    for row in rows:
        sp=row[focal].strip()
        key=network_key(row,keycols)
        if not sp or not key: continue
        taxa.add(sp); units.add(key); tax_units[sp].add(key)
        try:
            lat=float(row["decimalLatitude"]); lon=float(row["decimalLongitude"])
            if math.isfinite(lat) and math.isfinite(lon): geo_units.add(key)
        except: pass
        if all(valid(row.get(c)) for c in d["method_columns"]): method_units.add(key)
    repeated=sum(len(v)>=3 for v in tax_units.values())
    return dict(rows_projected=len(rows),focal_taxa=len(taxa),local_units=len(units),
                georeferenced_local_units=len(geo_units),repeated_focal_taxa_ge3=repeated,
                local_units_with_method_metadata=len(method_units),
                focal_taxon_set_sha256=hset(taxa),local_unit_set_sha256=hset(units))

def audit_croppol(field_path,sampling_path,d):
    fields=list(read_projected(field_path,d["safe_field_columns"]))
    sampling=list(read_projected(sampling_path,d["safe_sampling_columns"]))
    focal=d["focal_taxon_column"]; keycols=d["local_unit_key"]
    taxa=set(); units=set(); tax_units=defaultdict(set); geo_units=set()
    for row in fields:
        sp=row[focal].strip(); key=network_key(row,keycols)
        if not sp or not key: continue
        taxa.add(sp); units.add(key); tax_units[sp].add(key)
        try:
            lat=float(row["latitude"]); lon=float(row["longitude"])
            if math.isfinite(lat) and math.isfinite(lon): geo_units.add(key)
        except: pass
    method_by=defaultdict(set)
    for row in sampling:
        key=network_key(row,keycols)
        if not key: continue
        methods=[row.get(c,"").strip() for c in d["method_columns"] if valid(row.get(c))]
        for m in methods: method_by[key].add(m)
    method_units={k for k,v in method_by.items() if v}
    repeated=sum(len(v)>=3 for v in tax_units.values())
    return dict(field_rows_projected=len(fields),sampling_rows_projected=len(sampling),
                focal_taxa=len(taxa),local_units=len(units),georeferenced_local_units=len(geo_units),
                repeated_focal_taxa_ge3=repeated,local_units_with_method_metadata=len(method_units),
                focal_taxon_set_sha256=hset(taxa),local_unit_set_sha256=hset(units))

def decide(x,g):
    tests={
      "temporal_focal_taxa":x["focal_taxa"]>=g["temporal_min_focal_taxa"],
      "spatial_georeferenced_units":x["georeferenced_local_units"]>=g["spatial_min_georeferenced_local_units"],
      "repeated_focal_taxa":x["repeated_focal_taxa_ge3"]>=g["spatial_min_repeated_focal_taxa"],
      "method_metadata_units":x["local_units_with_method_metadata"]>=g["minimum_local_units_with_method_metadata"]
    }
    return tests,all(tests.values())

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--design",type=Path,required=True)
    ap.add_argument("--carni",type=Path,required=True)
    ap.add_argument("--croppol-field",type=Path,required=True)
    ap.add_argument("--croppol-sampling",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args(); d=json.loads(a.design.read_text()); g=d["support_gate"]
    cr=audit_carni(a.carni,d["candidates"]["CARNIDIET"])
    cp=audit_croppol(a.croppol_field,a.croppol_sampling,d["candidates"]["CROPPOL"])
    for x in (cr,cp):
        tests,p=decide(x,g); x["gate_tests"]=tests; x["support_gate_pass"]=p
    n=sum(x["support_gate_pass"] for x in (cr,cp))
    status="PASS_GLOBI_FRESH_METADATA_SUPPORT_CANDIDATE" if n else "HOLD_NO_GLOBI_FRESH_METADATA_SUPPORT_CANDIDATE"
    out={
      "version":"v0.1","status":status,"design":str(a.design),
      "candidates":{"CARNIDIET":cr,"CROPPOL":cp},
      "passing_candidate_count":n,
      "passing_candidates":[k for k,v in [("CARNIDIET",cr),("CROPPOL",cp)] if v["support_gate_pass"]],
      "firewall":{
        "safe_metadata_columns_projected":True,
        "partner_columns_used":False,"partner_values_persisted":False,
        "interaction_weights_used":False,"turnover_outcomes_computed":False,
        "focal_names_persisted":False,"coordinate_values_persisted":False,"date_values_persisted":False
      },
      "next_gate":"For passing candidates, freeze exact local-network reconstruction, phylogeny and external predictor eligibility before any partner column is admitted.",
      "current_fcp_manuscript_changed":False,"chun_el_v0_3_changed":False
    }
    a.out.parent.mkdir(parents=True,exist_ok=True); a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps(out,indent=2,sort_keys=True))
if __name__=="__main__": main()
