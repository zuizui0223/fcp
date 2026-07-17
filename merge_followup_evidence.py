#!/usr/bin/env python3
"""Merge targeted OpenAlex evidence into the species ranking before reclassification."""
from __future__ import annotations
import argparse, csv, json
from pathlib import Path


def read(path):
    with Path(path).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def as_int(v):
    try: return int(str(v or 0))
    except ValueError: return 0


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--species", default="data/global_flower_colour_species_ranked.csv")
    ap.add_argument("--followup", action="append", default=[])
    ap.add_argument("--out", default="data/global_flower_colour_species_ranked_resolved.csv")
    ap.add_argument("--qc-out", default="data/followup_merge_qc.json")
    a=ap.parse_args()
    species=read(a.species)
    by_name={r["canonical_name"]:r for r in species}
    before={n:(as_int(r.get("n_works")),as_int(r.get("n_context_matches")),as_int(r.get("max_score"))) for n,r in by_name.items()}
    seen=set(); merged_rows=0; touched=set()
    for path in a.followup:
        p=Path(path)
        if not p.exists(): continue
        for e in read(p):
            name=e.get("canonical_name",""); oid=e.get("openalex_id","")
            if name not in by_name or not oid or (name,oid) in seen: continue
            seen.add((name,oid)); merged_rows+=1; touched.add(name)
            r=by_name[name]; score=as_int(e.get("score")); basis=e.get("evidence_basis","")
            r["n_works"]=str(as_int(r.get("n_works"))+1)
            r["n_context_matches"]=str(as_int(r.get("n_context_matches"))+1)
            if "title_species" in basis or name.lower() in e.get("title","").lower():
                r["n_title_matches"]=str(as_int(r.get("n_title_matches"))+1)
            r["total_score"]=str(as_int(r.get("total_score"))+score)
            if score>as_int(r.get("max_score")):
                r["max_score"]=str(score); r["best_title"]=e.get("title",""); r["best_doi"]=e.get("doi",""); r["best_openalex_id"]=oid; r["best_match_evidence"]=e.get("evidence_snippet","")
    species.sort(key=lambda r:(-as_int(r.get("n_title_matches")),-as_int(r.get("max_score")),-as_int(r.get("total_score")),-as_int(r.get("n_works")),r.get("canonical_name","")))
    for i,r in enumerate(species,1): r["rank"]=str(i)
    with Path(a.out).open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(species[0].keys())); w.writeheader(); w.writerows(species)
    strengthened=sum(1 for n in touched if (as_int(by_name[n].get("n_works")),as_int(by_name[n].get("n_context_matches")),as_int(by_name[n].get("max_score")))!=before[n])
    qc={"followup_files":a.followup,"unique_followup_pairs":len(seen),"merged_rows":merged_rows,"species_touched":len(touched),"species_strengthened":strengthened,"mode":"merge_followup_evidence_v1"}
    Path(a.qc_out).write_text(json.dumps(qc,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(qc,ensure_ascii=False))

if __name__=="__main__": main()
