#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json,re
from collections import Counter,defaultdict
from pathlib import Path
BINOMIAL=re.compile(r"^[A-Z][a-z]{2,} [a-z][a-z-]{2,}$")
TRUTHY={"1","true","yes"}
def truthy(v): return str(v or "").strip().lower() in TRUTHY
def norm(v): return " ".join(str(v or "").split())
def contextual_names(row):
    candidates=[]; seen=set()
    for token in str(row.get("candidate_species") or "").split(";"):
        name=norm(token)
        if BINOMIAL.fullmatch(name) and name not in seen:
            seen.add(name); candidates.append(name)
    title=str(row.get("title") or "").lower()
    hits=[n for n in candidates if n.lower() in title]
    if hits: return hits[:5],"title_exact_mention"
    abstract=str(row.get("abstract") or "")[:1200].lower()
    hits=[n for n in candidates if n.lower() in abstract]
    if hits: return hits[:3],"abstract_first1200_mention"
    return [],"no_contextual_species_mention"
def main():
    p=argparse.ArgumentParser(); p.add_argument("--queue",required=True); p.add_argument("--census"); p.add_argument("--out",required=True); p.add_argument("--qc-out",required=True); a=p.parse_args()
    out=Path(a.out); qcpath=Path(a.qc_out); out.parent.mkdir(parents=True,exist_ok=True)
    diag=Counter(); evidence=defaultdict(list)
    with Path(a.queue).open(newline="",encoding="utf-8") as h:
        reader=csv.DictReader(h); diag["queue_columns"]=len(reader.fieldnames or [])
        for row in reader:
            diag["records_total"]+=1
            if row.get("screen_priority") not in {"P1_high_natural_itv","P1_high_population_itv"}: continue
            diag["records_p1"]+=1
            if any(truthy(row.get(c)) for c in ("cultivated_signal","induced_signal","ontogenetic_signal")):
                diag["records_hard_exclusion"]+=1; continue
            if not truthy(row.get("within_signal")) and not truthy(row.get("among_signal")):
                diag["records_without_direction"]+=1; continue
            names,source=contextual_names(row); diag[source]+=1
            if not names: continue
            diag["records_with_contextual_species"]+=1; diag["contextual_species_mentions"]+=len(names)
            for name in names: evidence[name].append((row,source))
    rows=[]; excluded=Counter()
    for name,pairs in evidence.items():
        records=[r for r,_ in pairs]; within=[r for r in records if truthy(r.get("within_signal"))]; among=[r for r in records if truthy(r.get("among_signal"))]
        if within and among: excluded["mixed_species"]+=1; continue
        if not within and not among: excluded["no_direction_species"]+=1; continue
        directional=within or among
        best=max(directional,key=lambda r:(int(float(r.get("cited_by_count") or 0)),len(r.get("abstract") or "")))
        rows.append({"canonical_name":name,"family":"","spatial_scale":"within_population" if within else "among_population","classification_source":"exploratory_automated","automated_classification_source":"systematic_P1_contextual_binomial_one_direction","candidate_context_source":";".join(sorted({s for _,s in pairs})),"n_supporting_p1_records":len(records),"within_signal_records":len(within),"among_signal_records":len(among),"natural_signal_records":sum(truthy(r.get("natural_signal")) for r in records),"discrete_signal_records":sum(r.get("provisional_variation_form")=="discrete_signal" for r in records),"best_title":best.get("title",""),"best_doi":best.get("doi",""),"best_openalex_id":best.get("record_id",""),"best_match_evidence":(best.get("abstract") or "")[:1200],"review_status":"unreviewed","role":"focal"})
    rows.sort(key=lambda r:(r["spatial_scale"],-r["n_supporting_p1_records"],r["canonical_name"]))
    classes={r["spatial_scale"] for r in rows}; status="complete" if len(rows)>=20 and len(classes)==2 else "insufficient"
    qc={"status":status,"eligible_species_before_gbif_validation":len(rows),"within_species":sum(r["spatial_scale"]=="within_population" for r in rows),"among_species":sum(r["spatial_scale"]=="among_population" for r in rows),"record_diagnostics":dict(diag),"excluded_species":dict(excluded),"semantic_guard":"Exploratory contextual-binomial P1 labels; strict GBIF plant/species validation is required downstream; original 34-species classifications are unchanged."}
    qcpath.write_text(json.dumps(qc,indent=2)+"\n",encoding="utf-8")
    if rows:
        with out.open("w",newline="",encoding="utf-8") as h:
            w=csv.DictWriter(h,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(json.dumps(qc,indent=2))
    if status!="complete": raise SystemExit(json.dumps(qc))
if __name__=="__main__": main()
