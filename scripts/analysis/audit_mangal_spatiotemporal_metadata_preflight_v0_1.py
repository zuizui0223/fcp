#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
DESIGN=ROOT/"data"/"mangal_spatiotemporal_metadata_preflight_design_v0_1.json"

ALLOWED={"dataset","network"}
FORBIDDEN={"node","interaction","trait","environment"}


def fetch_all(base:str, endpoint:str, count:int=1000):
    if endpoint not in ALLOWED or endpoint in FORBIDDEN:
        raise RuntimeError(f"forbidden endpoint: {endpoint}")
    out=[]
    page=0
    while True:
        url=f"{base}/{endpoint}?"+urllib.parse.urlencode({"count":count,"page":page})
        req=urllib.request.Request(url,headers={"User-Agent":"fcp-spatiotemporal-mangal-preflight/0.1"})
        with urllib.request.urlopen(req,timeout=60) as r:
            data=json.loads(r.read().decode("utf-8"))
        if not isinstance(data,list):
            raise RuntimeError(f"unexpected {endpoint} response type: {type(data).__name__}")
        out.extend(data)
        if len(data)<count:
            break
        page+=1
        if page>20:
            raise RuntimeError("pagination safety stop")
    return out


def coord_key(geom):
    if not isinstance(geom,dict):
        return None
    typ=geom.get("type")
    coords=geom.get("coordinates")
    if typ=="Point" and isinstance(coords,list) and len(coords)>=2:
        try:
            return (round(float(coords[0]),6),round(float(coords[1]),6))
        except Exception:
            return None
    return json.dumps(geom,sort_keys=True,separators=(",",":")) if coords is not None else None


def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    d=json.loads(DESIGN.read_text())
    base=d["base_url"]

    datasets=fetch_all(base,"dataset")
    networks=fetch_all(base,"network")

    # Firewall: keep only approved metadata fields in memory/output.
    ds_fields=set(d["allowed_fields"]["dataset"])
    nw_fields=set(d["allowed_fields"]["network"])
    datasets=[{k:v for k,v in x.items() if k in ds_fields} for x in datasets]
    networks=[{k:v for k,v in x.items() if k in nw_fields} for x in networks]

    by=defaultdict(list)
    for n in networks:
        by[n.get("dataset_id")].append(n)

    rows=[]
    gate=d["geometry_gate"]
    for ds in datasets:
        did=ds.get("id")
        ns=[n for n in by.get(did,[]) if n.get("public") is not False]
        geos=[n for n in ns if coord_key(n.get("geom")) is not None]
        coords={coord_key(n.get("geom")) for n in geos}
        dated=[n for n in ns if n.get("date")]
        complete=[n for n in ns if n.get("all_interactions") is True]
        geometry_pass=(
            len(ns)>=gate["min_public_networks"]
            and len(geos)>=gate["min_georeferenced_public_networks"]
            and len(coords)>=gate["min_distinct_locations"]
        )
        rows.append({
            "dataset_id":did,
            "dataset_name":ds.get("name"),
            "reference_id":ds.get("ref_id"),
            "dataset_date":ds.get("date"),
            "public_networks":len(ns),
            "georeferenced_public_networks":len(geos),
            "distinct_locations":len(coords),
            "dated_public_networks":len(dated),
            "all_interactions_true_networks":len(complete),
            "geometry_gate_pass":geometry_pass
        })

    passing=[r for r in rows if r["geometry_gate_pass"]]
    result={
        "version":"v0.1",
        "status":"MANGAL_METADATA_GEOMETRY_PREFLIGHT_COMPLETE",
        "design":"data/mangal_spatiotemporal_metadata_preflight_design_v0_1.json",
        "source":"Mangal API v2",
        "allowed_endpoints_used":["dataset","network"],
        "forbidden_endpoints_used":[],
        "taxon_identities_opened":False,
        "interaction_edges_opened":False,
        "trait_values_opened":False,
        "environment_values_opened":False,
        "datasets_total":len(datasets),
        "networks_total":len(networks),
        "geometry_gate":{
            "criteria":gate,
            "passing_datasets":len(passing),
            "passing_dataset_ids":[r["dataset_id"] for r in passing],
            "passing_dataset_names":[r["dataset_name"] for r in passing]
        },
        "dataset_metadata":rows,
        "full_interaction_spatial_admission_authorized":False,
        "next_gate":d["next_gate_if_geometry_passes"] if passing else "HOLD_NO_GEOMETRY_QUALIFIED_DATASETS",
        "interpretation":"Geometry qualification only. all_interactions is retained as metadata but does not establish standardized effort. No nodes or interaction edges were opened."
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":result["status"],
        "datasets_total":result["datasets_total"],
        "networks_total":result["networks_total"],
        "passing_datasets":result["geometry_gate"]["passing_datasets"],
        "passing_dataset_ids":result["geometry_gate"]["passing_dataset_ids"],
        "passing_dataset_names":result["geometry_gate"]["passing_dataset_names"],
        "next_gate":result["next_gate"]
    },indent=2))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
