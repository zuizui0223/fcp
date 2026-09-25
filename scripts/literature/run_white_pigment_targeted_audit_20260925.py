#!/usr/bin/env python3
from __future__ import annotations
import argparse, html, json, os, re, time, urllib.parse, urllib.request
from pathlib import Path
import numpy as np
import pandas as pd

OPENALEX="https://api.openalex.org/works"
MORPHS={"white","yellow_orange","red_pink","blue_purple"}
FLORAL=re.compile(r"\b(flower|floral|petal|corolla|tepal|labellum|bract)\w*\b",re.I)
PIGMENT=re.compile(r"\b(anthocyanin|flavonoid|pigment|pigmentation|MYB|bHLH|WD40|WDR|MBW|DFR|CHS|CHI|F3H|F3.?H|F3.?5.?H|ANS|UFGT|transcriptom|metabolom|expression)\w*\b",re.I)
WHITE=re.compile(r"\b(white[- ]flower|white morph|white[- ]flowered|acyanic|albino|unpigmented|non[- ]pigmented)\w*\b",re.I)
PIGMENTED=re.compile(r"\b(purple|pink|red|blue|magenta|pigmented|colou?r(?:ed)?)\b",re.I)
NATURAL=re.compile(r"\b(wild|natural population|naturally occurring|polymorph|morph|geographic variation|population)\w*\b",re.I)
ARTIFICIAL=re.compile(r"\b(cultivar|breeding|ornamental|transgenic|transformation|mutagenesis|mutant library|horticultur|floricultur)\w*\b",re.I)
ANTHO=re.compile(r"\b(anthocyanin|flavonoid biosynth|flavonoid pathway)\w*\b",re.I)
MYB=re.compile(r"\b(MYB|bHLH|WD40|WDR|MBW)\w*\b",re.I)
STRUCT=re.compile(r"\b(CHS|CHI|F3H|F3.?H|F3.?5.?H|DFR|ANS|UFGT)\b",re.I)
EXPR=re.compile(r"\b(differential expression|expression|transcriptom|transcription|downregulat|upregulat|repress|activat)\w*\b",re.I)
METAB=re.compile(r"\b(HPLC|metabolom|anthocyanin content|flavonoid content|pigment content)\w*\b",re.I)
LOSS=re.compile(r"\b(loss|absent|absence|reduc|suppress|downregulat|inactivat|acyanic|unpigmented|lack)\w*\b",re.I)
TEMP=re.compile(r"\b(heat|temperature|warm|thermal|high-temperature|high temperature)\w*\b",re.I)

Q1='anthocyanin OR flavonoid OR pigment OR pigmentation OR MYB OR bHLH OR WD40 OR DFR OR CHS OR CHI OR F3H OR ANS OR UFGT'
Q2='"white flower" OR "white-flowered" OR "white morph" OR acyanic OR albino OR purple OR pink OR red OR blue transcriptome metabolome expression anthocyanin flavonoid pigment'

def bool_series(s):
    if s.dtype==bool: return s.fillna(False)
    return s.astype(str).str.strip().str.lower().isin({"true","1","yes","y"})

def clean(v):
    return re.sub(r"\s+"," ",re.sub(r"<[^>]+>"," ",html.unescape(str(v or "")))).strip()

def abstract(inv):
    if not isinstance(inv,dict): return ""
    pairs=[(int(i),str(word)) for word,locs in inv.items() if isinstance(locs,list) for i in locs if isinstance(i,int)]
    return clean(" ".join(w for _,w in sorted(pairs)))

def request(url,timeout=45,retries=5):
    last=None
    for n in range(retries):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"fcp-white-pigment-audit/1.0","Accept":"application/json"})
            with urllib.request.urlopen(req,timeout=timeout) as h:
                return json.load(h)
        except Exception as e:
            last=e
            if n+1<retries: time.sleep(min(2**n,16))
    raise RuntimeError(str(last))

def targets(measured,technical,join_key,high_ids):
    d=pd.read_csv(measured,usecols=["photo_id","species","morph","global_classifiable"])
    d=d.loc[bool_series(d.global_classifiable)&d.morph.isin(MORPHS)].copy()
    tech=pd.read_csv(technical,compression="gzip",dtype={"measurement_id":str})
    join=pd.read_csv(join_key,dtype={"measurement_id":str})
    hi=set(pd.read_csv(high_ids,dtype={"measurement_id":str}).measurement_id.astype(str))
    tech=tech.merge(join,on="measurement_id",how="left",validate="one_to_one")
    tech["high_clip"]=tech.measurement_id.astype(str).isin(hi)
    d=d.merge(tech[["photo_id","near_clip_fraction","high_clip"]],on="photo_id",how="left",validate="one_to_one")
    d=d.loc[d.near_clip_fraction.notna() & ~d.high_clip.fillna(False)].copy()
    d["white"]=(d.morph=="white").astype(int)
    q=d.groupby("species").white.agg(["sum","count"]).reset_index()
    q["nonwhite"]=q["count"]-q["sum"]
    q=q.loc[(q["sum"]>=5)&(q.nonwhite>=5)].rename(columns={"sum":"n_white","count":"n_rows"})
    return q.sort_values("species").reset_index(drop=True)

def flags(text):
    return {
      "anthocyanin_pathway":int(bool(ANTHO.search(text))),
      "MYB_MBW_regulation":int(bool(MYB.search(text))),
      "structural_gene":int(bool(STRUCT.search(text))),
      "expression_regulation":int(bool(EXPR.search(text))),
      "metabolite_evidence":int(bool(METAB.search(text))),
      "loss_or_reduction_language":int(bool(LOSS.search(text))),
      "white_contrast":int(bool(WHITE.search(text) and PIGMENTED.search(text))),
      "temperature_link":int(bool(TEMP.search(text))),
      "natural_context":int(bool(NATURAL.search(text))),
      "artificial_context":int(bool(ARTIFICIAL.search(text))),
    }

def auto_tier(f):
    if f["artificial_context"] and not f["natural_context"]: return "excluded_artificial"
    if f["white_contrast"] and f["natural_context"] and (f["MYB_MBW_regulation"] or f["structural_gene"] or f["expression_regulation"]):
        return "candidate_A_pending_manual"
    if f["white_contrast"] and f["natural_context"] and (f["anthocyanin_pathway"] or f["metabolite_evidence"]):
        return "candidate_B_pending_manual"
    return "candidate_C"

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--measured",required=True); p.add_argument("--technical-table",required=True)
    p.add_argument("--join-key",required=True); p.add_argument("--high-clip-ids",required=True)
    p.add_argument("--outdir",required=True); p.add_argument("--batch-size",type=int,default=6)
    p.add_argument("--per-page",type=int,default=100); p.add_argument("--timeout",type=int,default=45)
    p.add_argument("--api-key-env",default="OPENALEX_API_KEY")
    a=p.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    t=targets(a.measured,a.technical_table,a.join_key,a.high_clip_ids)
    t.to_csv(out/"target_species.csv",index=False)
    if len(t)!=281:
        raise SystemExit(f"target panel mismatch: expected 281, got {len(t)}")
    names=t.species.astype(str).tolist(); api=os.environ.get(a.api_key_env,"").strip()
    rows=[]; errors=[]; seen=set(); requests=0
    modes=[("pigment_mechanism",Q1),("white_pigmented_molecular",Q2)]
    for bi in range(0,len(names),a.batch_size):
        batch=names[bi:bi+a.batch_size]
        name_query=" OR ".join(f'"{x}"' for x in batch)
        for mode,terms in modes:
            q=f"({name_query}) ({terms})"
            params={"search":q,"per-page":min(max(a.per_page,1),200)}
            if api: params["api_key"]=api
            requests+=1
            try:
                payload=request(OPENALEX+"?"+urllib.parse.urlencode(params),a.timeout)
            except Exception as e:
                errors.append({"batch":bi//a.batch_size+1,"mode":mode,"error":str(e)[:500]}); continue
            for item in payload.get("results",[]) or []:
                if not isinstance(item,dict): continue
                title=clean(item.get("title")); abs_=abstract(item.get("abstract_inverted_index")); text=clean(title+" "+abs_)
                low=text.lower()
                matched=[sp for sp in batch if sp.lower() in low]
                if not matched or not FLORAL.search(text) or not PIGMENT.search(text): continue
                oid=str(item.get("id") or "")
                loc=item.get("primary_location") or {}
                landing=str(loc.get("landing_page_url") or oid) if isinstance(loc,dict) else oid
                for sp in matched:
                    key=(sp,oid)
                    if key in seen: continue
                    seen.add(key); f=flags(text)
                    m=re.search(re.escape(sp),abs_,re.I)
                    snippet=abs_[max(0,m.start()-220):min(len(abs_),m.end()+500)] if m else title
                    rows.append({
                      "species":sp,"openalex_id":oid,"title":title,"year":item.get("publication_year") or "",
                      "doi":str(item.get("doi") or "").replace("https://doi.org/",""),
                      "landing_url":landing,"query_mode":mode,"auto_tier":auto_tier(f),
                      **f,"evidence_snippet":clean(snippet)[:900]
                    })
            time.sleep(0.12)
    c=pd.DataFrame(rows)
    if len(c):
        # consolidate duplicate species-work retrieved by both fixed query modes
        flagcols=["anthocyanin_pathway","MYB_MBW_regulation","structural_gene","expression_regulation","metabolite_evidence",
                  "loss_or_reduction_language","white_contrast","temperature_link","natural_context","artificial_context"]
        def combine(g):
            r=g.iloc[0].to_dict()
            r["query_mode"]=";".join(sorted(set(g.query_mode.astype(str))))
            for x in flagcols: r[x]=int(g[x].max())
            r["auto_tier"]=auto_tier(r)
            return pd.Series(r)
        c=c.groupby(["species","openalex_id"],as_index=False,group_keys=False).apply(combine,include_groups=False).reset_index(drop=True)
        order={"candidate_A_pending_manual":0,"candidate_B_pending_manual":1,"candidate_C":2,"excluded_artificial":3}
        c["_ord"]=c.auto_tier.map(order).fillna(9)
        c=c.sort_values(["_ord","species","year"],ascending=[True,True,False]).drop(columns="_ord")
    c.to_csv(out/"candidate_works.csv",index=False)
    per=c.groupby("species").agg(
        candidate_works=("openalex_id","nunique"),
        candidate_A=("auto_tier",lambda s:int((s=="candidate_A_pending_manual").sum())),
        candidate_B=("auto_tier",lambda s:int((s=="candidate_B_pending_manual").sum())),
        temperature_candidates=("temperature_link","sum")
    ).reset_index() if len(c) else pd.DataFrame(columns=["species","candidate_works","candidate_A","candidate_B","temperature_candidates"])
    t.merge(per,on="species",how="left").fillna(0).to_csv(out/"species_candidate_summary.csv",index=False)
    summary={
      "schema":"fcp_white_pigment_targeted_audit_v1","status":"complete",
      "target_species":int(len(t)),"requests_used":requests,"failed_requests":len(errors),"errors":errors,
      "candidate_works":int(len(c)),"species_with_any_candidate":int(c.species.nunique()) if len(c) else 0,
      "species_with_candidate_A":int(c.loc[c.auto_tier.eq("candidate_A_pending_manual"),"species"].nunique()) if len(c) else 0,
      "species_with_candidate_B":int(c.loc[c.auto_tier.eq("candidate_B_pending_manual"),"species"].nunique()) if len(c) else 0,
      "species_with_temperature_candidate":int(c.loc[c.temperature_link.eq(1),"species"].nunique()) if len(c) else 0,
      "manual_verification_required":True,
      "hard_nonclaims":["automated tiers are navigation aids, not verified mechanism evidence","not a random sample of all angiosperms","does not infer evolutionary transition direction"]
    }
    (out/"result.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))

if __name__=="__main__":
    main()
