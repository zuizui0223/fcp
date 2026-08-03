#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json,time,urllib.parse,urllib.request
from pathlib import Path

def get_json(url,timeout=45,retries=4):
    last=None
    for i in range(retries):
        try:
            req=urllib.request.Request(url,headers={"Accept":"application/json","User-Agent":"fcp-exploratory-taxon-resolver/1.0"})
            with urllib.request.urlopen(req,timeout=timeout) as r: return json.load(r)
        except Exception as e:
            last=e; time.sleep(2**i)
    raise RuntimeError(url) from last

def main():
    p=argparse.ArgumentParser(); p.add_argument("--input",required=True); p.add_argument("--out",required=True); p.add_argument("--audit",required=True); p.add_argument("--timeout",type=int,default=45); p.add_argument("--retries",type=int,default=4); a=p.parse_args()
    rows=list(csv.DictReader(Path(a.input).open(encoding="utf-8"))); kept=[]; audit=[]
    for row in rows:
        name=row["canonical_name"].strip(); params=urllib.parse.urlencode({"name":name,"kingdom":"Plantae","strict":"true"})
        try: m=get_json("https://api.gbif.org/v1/species/match?"+params,a.timeout,a.retries)
        except Exception as e:
            audit.append({"input_name":name,"accepted":False,"reason":type(e).__name__}); continue
        match=str(m.get("matchType") or ""); rank=str(m.get("rank") or ""); kingdom=str(m.get("kingdom") or ""); family=str(m.get("family") or ""); canonical=str(m.get("canonicalName") or m.get("scientificName") or name).split(" ")[0:2]
        canonical=" ".join(canonical)
        accepted=match in {"EXACT","CONFIDENCE"} and rank=="SPECIES" and kingdom=="Plantae" and bool(family)
        audit.append({"input_name":name,"accepted":accepted,"match_type":match,"rank":rank,"kingdom":kingdom,"family":family,"accepted_name":canonical,"usage_key":m.get("acceptedUsageKey") or m.get("usageKey") or ""})
        if accepted:
            out=dict(row); out["original_candidate_name"]=name; out["canonical_name"]=canonical; out["family"]=family; out["gbif_taxon_validation"]="strict_species_match"; kept.append(out)
        time.sleep(0.05)
    dedup={}
    for row in kept:
        key=row["canonical_name"]
        if key not in dedup or int(row.get("n_supporting_p1_records") or 0)>int(dedup[key].get("n_supporting_p1_records") or 0): dedup[key]=row
    kept=sorted(dedup.values(),key=lambda r:(r["spatial_scale"],r["canonical_name"]))
    if kept:
        with Path(a.out).open("w",newline="",encoding="utf-8") as h:
            w=csv.DictWriter(h,fieldnames=list(kept[0])); w.writeheader(); w.writerows(kept)
    with Path(a.audit).open("w",newline="",encoding="utf-8") as h:
        fields=sorted({k for r in audit for k in r}); w=csv.DictWriter(h,fieldnames=fields); w.writeheader(); w.writerows(audit)
    summary={"requested":len(rows),"strict_species_accepted":len(kept),"within":sum(r["spatial_scale"]=="within_population" for r in kept),"among":sum(r["spatial_scale"]=="among_population" for r in kept),"families":len({r["family"] for r in kept})}
    print(json.dumps(summary,indent=2))
    if len(kept)<20 or len({r["spatial_scale"] for r in kept})<2: raise SystemExit(json.dumps(summary))
if __name__=="__main__": main()
