#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
DESIGN=ROOT/"data"/"mangal_node_taxonomy_support_design_v0_1.json"
PARENT=ROOT/"results"/"mangal_spatiotemporal_metadata_preflight_v0_1"/"result.json"

ALLOWED={"dataset","network","node"}
FORBIDDEN={"interaction","trait","environment"}


def fetch(base:str, endpoint:str, params:dict, count:int=1000):
    if endpoint not in ALLOWED or endpoint in FORBIDDEN:
        raise RuntimeError(f"forbidden endpoint: {endpoint}")
    out=[]
    page=0
    while True:
        q=dict(params)
        q["count"]=count
        q["page"]=page
        url=f"{base}/{endpoint}?"+urllib.parse.urlencode(q)
        req=urllib.request.Request(url,headers={"User-Agent":"fcp-mangal-node-support/0.1"})
        with urllib.request.urlopen(req,timeout=60) as r:
            data=json.loads(r.read().decode("utf-8"))
        if not isinstance(data,list):
            raise RuntimeError(f"unexpected {endpoint} response type")
        out.extend(data)
        if len(data)<count:
            break
        page+=1
        if page>50:
            raise RuntimeError("pagination safety stop")
    return out


def has_geom(x):
    g=x.get("geom")
    return isinstance(g,dict) and g.get("coordinates") is not None


def hash_ids(ids):
    payload="\n".join(str(x) for x in sorted(ids,key=lambda z:str(z)))
    return hashlib.sha256(payload.encode()).hexdigest()


def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()

    d=json.loads(DESIGN.read_text())
    p=json.loads(PARENT.read_text())
    ids=d["qualified_dataset_ids"]
    parent_ids=[x["dataset_id"] for x in p["geometry_gate"]["passing"]]
    if sorted(ids)!=sorted(parent_ids):
        raise RuntimeError("qualified dataset ids drifted from parent geometry result")

    base="https://mangal.io/api/v2"
    results=[]
    for did in ids:
        networks=fetch(base,"network",{"dataset_id":did})
        networks=[n for n in networks if n.get("public") is not False and has_geom(n)]
        network_ids=[int(n["id"]) for n in networks if n.get("id") is not None]
        occupancy=Counter()
        unique_taxa=set()

        for nid in network_ids:
            nodes=fetch(base,"node",{"network_id":nid})
            seen=set()
            for n in nodes:
                if n.get("node_level")!="taxon":
                    continue
                tax=n.get("taxonomy")
                if not isinstance(tax,dict) or tax.get("rank")!="species":
                    continue
                tid=n.get("taxonomy_id")
                if tid is None:
                    continue
                tid=int(tid)
                seen.add(tid)
                unique_taxa.add(tid)
            for tid in seen:
                occupancy[tid]+=1

        repeated={tid for tid,n in occupancy.items() if n>=2}
        ge3=sum(1 for n in occupancy.values() if n>=3)
        ge5=sum(1 for n in occupancy.values() if n>=5)
        support_pass=len(repeated)>=d["support_definition"]["min_repeated_species"]

        results.append({
            "dataset_id":did,
            "georeferenced_public_networks_audited":len(network_ids),
            "species_level_unique_taxa":len(unique_taxa),
            "species_repeated_ge2_networks":len(repeated),
            "species_repeated_ge3_networks":ge3,
            "species_repeated_ge5_networks":ge5,
            "species_taxonomy_id_set_sha256":hash_ids(unique_taxa),
            "repeated_species_id_set_sha256":hash_ids(repeated),
            "necessary_node_support_pass":support_pass
        })

    passing=[x for x in results if x["necessary_node_support_pass"]]
    out={
        "version":"v0.1",
        "status":"MANGAL_NODE_TAXONOMY_SUPPORT_AUDIT_COMPLETE",
        "design":"data/mangal_node_taxonomy_support_design_v0_1.json",
        "parent_result":"results/mangal_spatiotemporal_metadata_preflight_v0_1/result.json",
        "allowed_endpoints_used":["network","node"],
        "forbidden_endpoints_used":[],
        "taxon_identities_opened":True,
        "taxon_names_persisted":False,
        "interaction_edges_opened":False,
        "trait_values_opened":False,
        "environment_values_opened":False,
        "dataset_results":results,
        "necessary_support_passing_datasets":len(passing),
        "necessary_support_passing_dataset_ids":[x["dataset_id"] for x in passing],
        "full_interaction_spatial_admission_authorized":False,
        "next_gate":(
            d["next_gate_if_pass"] if passing
            else "HOLD_NO_DATASET_WITH_MINIMUM_REPEATED_SPECIES_SUPPORT"
        ),
        "interpretation":"Node/taxonomy support is a necessary ceiling only. No interaction edges, partner profiles, focal-guild roles, or turnover outcomes were opened."
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":out["status"],
        "necessary_support_passing_datasets":out["necessary_support_passing_datasets"],
        "necessary_support_passing_dataset_ids":out["necessary_support_passing_dataset_ids"],
        "dataset_results":results,
        "next_gate":out["next_gate"]
    },indent=2))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
