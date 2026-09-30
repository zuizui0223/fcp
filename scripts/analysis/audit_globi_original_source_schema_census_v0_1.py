#!/usr/bin/env python3
from __future__ import annotations

import argparse
import base64
import json
import os
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

API="https://api.github.com"
UA="fcp-spatiotemporal-globi-schema-census/0.1"

def request_json(url:str,token:str|None):
    headers={"Accept":"application/vnd.github+json","User-Agent":UA,"X-GitHub-Api-Version":"2022-11-28"}
    if token:
        headers["Authorization"]="Bearer "+token
    req=urllib.request.Request(url,headers=headers)
    try:
        with urllib.request.urlopen(req,timeout=60) as r:
            return json.load(r)
    except urllib.error.HTTPError as e:
        if e.code==404:
            return None
        raise

def iter_strings(x,path=""):
    if isinstance(x,dict):
        for k,v in x.items():
            p=f"{path}.{k}" if path else str(k)
            yield str(k),p
            yield from iter_strings(v,p)
    elif isinstance(x,list):
        for i,v in enumerate(x):
            yield from iter_strings(v,f"{path}[{i}]")
    elif isinstance(x,str):
        yield x,path

def norm_tokens(cfg):
    vals=[]
    for text,path in iter_strings(cfg):
        vals.append((str(text).lower(),str(path).lower()))
    return vals

def has_any(vals,tokens):
    for text,path in vals:
        blob=text+" "+path
        if any(tok in blob for tok in tokens):
            return True
    return False

def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--design",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    d=json.loads(a.design.read_text())
    org=d["source_org"]
    token=os.environ.get("GH_TOKEN") or os.environ.get("GITHUB_TOKEN")

    repos=[]
    page=1
    while True:
        url=f"{API}/orgs/{urllib.parse.quote(org)}/repos?type=public&per_page=100&page={page}"
        batch=request_json(url,token)
        if not batch:
            break
        repos.extend(batch)
        if len(batch)<100:
            break
        page+=1
        if page>20:
            raise RuntimeError("repository pagination safety limit exceeded")

    audited=[]
    for repo in repos:
        if repo.get("archived"):
            continue
        full=repo["full_name"]
        branch=repo.get("default_branch") or "main"
        url=f"{API}/repos/{full}/contents/globi.json?ref={urllib.parse.quote(branch,safe='')}"
        item=request_json(url,token)
        if not item or item.get("type")!="file":
            continue
        content=item.get("content","").replace("\n","")
        if not content:
            continue
        raw=base64.b64decode(content)
        try:
            cfg=json.loads(raw.decode("utf-8"))
        except Exception:
            continue
        vals=norm_tokens(cfg)
        source_taxon=has_any(vals,("sourcetaxonname","sourcetaxonid","sourcetaxonspeciesname"))
        target_taxon=has_any(vals,("targettaxonname","targettaxonid","targettaxonspeciesname"))
        lat=has_any(vals,("decimallatitude"," latitude","latitude"))
        lon=has_any(vals,("decimallongitude"," longitude","longitude"))
        temporal=has_any(vals,("eventdate","year_collected","month_collected","day_collected"," date","date_","year","month"))
        context=has_any(vals,("localityname","location_name","locationname","site","plot","study","sample","network"))
        itype=has_any(vals,("interactiontypename","interactiontypeid","interaction type"))
        method=has_any(vals,("collection_method","sampling","method"))
        score=sum([source_taxon,target_taxon,lat,lon,temporal,context,itype,method])
        gate=source_taxon and target_taxon and lat and lon and (temporal or context) and itype
        audited.append({
            "repository":full,
            "default_branch":branch,
            "globi_json_blob_sha":item.get("sha"),
            "source_taxon_schema":source_taxon,
            "target_taxon_schema":target_taxon,
            "latitude_schema":lat,
            "longitude_schema":lon,
            "temporal_schema":temporal,
            "context_schema":context,
            "interaction_type_schema":itype,
            "sampling_method_schema":method,
            "schema_score":score,
            "candidate_gate_pass":gate
        })

    passing=[x for x in audited if x["candidate_gate_pass"]]
    passing=sorted(passing,key=lambda x:(-x["schema_score"],x["repository"]))
    selected=passing[:int(d["ranking"]["max_candidates"])]
    status="PASS_GLOBI_ORIGINAL_SOURCE_SCHEMA_CANDIDATES" if len(selected)>=int(d["pass_rule"]["minimum_schema_candidates"]) else "HOLD_NO_GLOBI_ORIGINAL_SOURCE_SCHEMA_CANDIDATES"

    result={
      "version":"v0.1",
      "status":status,
      "design":str(a.design),
      "source_org":org,
      "github_repositories_seen":len(repos),
      "repositories_with_parseable_globi_json":len(audited),
      "schema_passing_repositories":len(passing),
      "selected_candidate_count":len(selected),
      "selected_candidates":selected,
      "all_passing_repository_names":[x["repository"] for x in passing],
      "firewall":{
        "globi_json_opened":True,
        "readmes_opened":False,
        "source_data_rows_opened":False,
        "interaction_rows_opened":False,
        "taxon_values_opened":False,
        "coordinate_values_opened":False,
        "date_values_opened":False,
        "biological_outcome_computed":False
      },
      "next_gate":d["next_gate_if_pass"] if status.startswith("PASS_") else d["next_gate_if_hold"],
      "current_fcp_manuscript_changed":False,
      "chun_el_v0_3_changed":False
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,indent=2,sort_keys=True))
    return 0

if __name__=="__main__":
    raise SystemExit(main())
