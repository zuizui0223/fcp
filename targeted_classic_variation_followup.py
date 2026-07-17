#!/usr/bin/env python3
"""Fourth-pass OpenAlex search using classic flower-colour variation terminology."""
from __future__ import annotations
import argparse, csv, html, json, os, re, time, urllib.parse, urllib.request
from pathlib import Path

OPENALEX="https://api.openalex.org/works"
TERMS='"flower color variation" OR "flower colour variation" OR "floral color polymorphism" OR "floral colour polymorphism" OR "petal color variation" OR "petal colour variation" OR "corolla color variation" OR "corolla colour variation" OR "color morph" OR "colour morph" OR "geographic cline"'
RELEVANT=re.compile(r"\b(?:flower|floral|petal|corolla|labellum|lip)\b.{0,140}\b(?:colou?r|pigment|morph|polymorph|variation|cline)\b|\b(?:colou?r|pigment|morph|polymorph|variation|cline)\b.{0,140}\b(?:flower|floral|petal|corolla|labellum|lip)\b",re.I)
NATURAL=re.compile(r"\b(?:population|wild|geographic|cline|morph frequenc|pollinator|selection|local adaptation|spatial variation)\b",re.I)
NEGATIVE=re.compile(r"\b(?:cultivar|breeding|ornamental|transgenic|mutagenesis|horticultur|floricultur|crop improvement)\b",re.I)
FIELDS=["canonical_name","family","openalex_id","title","year","doi","landing_url","score","evidence_basis","evidence_snippet","query_mode"]

def clean(v): return re.sub(r"\s+"," ",re.sub(r"<[^>]+>"," ",html.unescape(str(v or "")))).strip()
def abstract(inv):
    if not isinstance(inv,dict): return ""
    p=[(i,str(w)) for w,xs in inv.items() if isinstance(xs,list) for i in xs if isinstance(i,int)]
    return clean(" ".join(w for _,w in sorted(p)))
def rows(path):
    with Path(path).open(newline="",encoding="utf-8") as f: return list(csv.DictReader(f))
def request(url,timeout,retries):
    last=None
    for n in range(retries):
        try:
            q=urllib.request.Request(url,headers={"User-Agent":"fcp-classic-variation-followup/1.0","Accept":"application/json"})
            with urllib.request.urlopen(q,timeout=timeout) as r: return json.load(r)
        except Exception as e:
            last=e
            if n+1<retries: time.sleep(2**n)
    raise RuntimeError(str(last))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--deferred",default="data/global_flower_colour_deferred.csv")
    ap.add_argument("--prior",action="append",default=[])
    ap.add_argument("--out",default="data/targeted_classic_variation_followup_works.csv")
    ap.add_argument("--qc-out",default="data/targeted_classic_variation_followup_qc.json")
    ap.add_argument("--max-targets",type=int,default=100)
    ap.add_argument("--batch-size",type=int,default=8)
    ap.add_argument("--per-page",type=int,default=100)
    ap.add_argument("--timeout",type=int,default=30)
    ap.add_argument("--retries",type=int,default=3)
    ap.add_argument("--api-key-env",default="OPENALEX_API_KEY")
    a=ap.parse_args()
    prior=set()
    for path in a.prior:
        if Path(path).exists(): prior.update(r["canonical_name"] for r in rows(path))
    eligible=[]
    for r in rows(a.deferred):
        if r.get("canonical_name") in prior: continue
        cls=r.get("evidence_class","")
        if cls=="possible_polymorphism" or (cls=="insufficient_direct_evidence" and (int(r.get("n_title_matches") or 0)>0 or int(r.get("natural_signal_count") or 0)>0)):
            eligible.append(r)
    targets=eligible[:a.max_targets]
    api=os.environ.get(a.api_key_env,"").strip(); out=[]; errors=[]; seen=set(); requests=0
    for bi in range(0,len(targets),a.batch_size):
        batch=targets[bi:bi+a.batch_size]; genera=sorted({r["canonical_name"].split()[0] for r in batch})
        query=f'({" OR ".join(chr(34)+g+chr(34) for g in genera)}) ({TERMS})'
        params={"search":query,"per-page":min(max(a.per_page,1),200)}
        if api: params["api_key"]=api
        requests+=1
        try: payload=request(OPENALEX+"?"+urllib.parse.urlencode(params),a.timeout,a.retries)
        except RuntimeError as e:
            errors.append({"batch":bi//a.batch_size+1,"genera":genera,"error":str(e)[:500]}); continue
        results=payload.get("results",[])
        if not isinstance(results,list): continue
        for item in results:
            if not isinstance(item,dict): continue
            oid=str(item.get("id") or ""); title=clean(item.get("title")); abs_=abstract(item.get("abstract_inverted_index")); text=f"{title} {abs_}"; low=text.lower()
            for t in batch:
                name=t["canonical_name"]
                if name.lower() not in low or not RELEVANT.search(text) or (name,oid) in seen: continue
                score=14+(10 if name.lower() in title.lower() else 0)+(8 if NATURAL.search(text) else 0)+(5 if re.search(r"polymorph|morph frequenc|geographic cline",text,re.I) else 0)-(10 if NEGATIVE.search(text) else 0)
                if score<14: continue
                seen.add((name,oid)); primary=item.get("primary_location") or {}; landing=str(primary.get("landing_page_url") or oid) if isinstance(primary,dict) else oid
                m=re.search(re.escape(name),abs_,re.I); snippet=abs_[max(0,m.start()-220):min(len(abs_),m.end()+220)] if m else title
                out.append({"canonical_name":name,"family":t.get("family",""),"openalex_id":oid,"title":title,"year":item.get("publication_year") or "","doi":str(item.get("doi") or "").replace("https://doi.org/",""),"landing_url":landing,"score":score,"evidence_basis":"exact_species+classic_variation_term","evidence_snippet":clean(snippet)[:500],"query_mode":"classic_variation"})
    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    with Path(a.out).open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=FIELDS); w.writeheader(); w.writerows(sorted(out,key=lambda r:(-int(r["score"]),r["canonical_name"])))
    qc={"eligible_unresolved_targets":len(targets),"prior_species_excluded":len(prior),"requests_used":requests,"failed_requests":len(errors),"errors":errors,"targets_with_new_evidence":len({r["canonical_name"] for r in out}),"retained_followup_works":len(out),"api_key_present":bool(api),"mode":"targeted_classic_variation_followup_v1"}
    Path(a.qc_out).write_text(json.dumps(qc,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(qc))
if __name__=="__main__": main()
