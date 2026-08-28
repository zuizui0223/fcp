#!/usr/bin/env python3
"""Retrieve and deduplicate all OpenAlex OQL v2 title/abstract search blocks.

This search-surface audit records both exact identifier recovery and high-confidence
same-title/version recovery of the historical 34 benchmark sources. The latter prevents
published/preprint DOI differences from being misreported as search misses. No biological
inclusion or spatial classification is performed.
"""
from __future__ import annotations

import argparse, csv, html, json, os, re, time, unicodedata, urllib.parse, urllib.request
from difflib import SequenceMatcher
from pathlib import Path
from typing import Any

USER_AGENT = "fcp-jbi-search-v2/2.0 (https://github.com/zuizui0223/fcp)"
SPACE_RE = re.compile(r"\s+")
TAG_RE = re.compile(r"<[^>]+>")


def clean_text(value: Any) -> str:
    text = html.unescape(str(value or "")); text = TAG_RE.sub(" ", text)
    return SPACE_RE.sub(" ", text).strip()


def normalize_doi(value: Any) -> str:
    doi = clean_text(value).lower()
    doi = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", doi)
    doi = re.sub(r"^doi:\s*", "", doi)
    return doi.strip().rstrip(".")


def normalize_title(value: Any) -> str:
    title = unicodedata.normalize("NFKD", clean_text(value)).lower()
    title = "".join(ch for ch in title if not unicodedata.combining(ch))
    title = re.sub(r"[^a-z0-9\u3040-\u30ff\u3400-\u9fff]+", " ", title)
    return SPACE_RE.sub(" ", title).strip()


def normalize_openalex(value: Any) -> str:
    text = clean_text(value).rstrip("/")
    if "openalex.org/" in text.lower(): return text.rsplit("/", 1)[-1].upper()
    if re.fullmatch(r"W\d+", text, re.I): return text.upper()
    return ""


def stable_key(row: dict[str, Any]) -> str:
    doi = normalize_doi(row.get("doi"))
    if doi: return "doi:" + doi
    return f"title:{normalize_title(row.get('title'))}|year:{row.get('year') or ''}"


def reconstruct_abstract(inverted: Any) -> str:
    if not isinstance(inverted, dict): return ""
    positions=[]
    for word,indexes in inverted.items():
        if isinstance(indexes,list):
            for index in indexes:
                if isinstance(index,int): positions.append((index,str(word)))
    positions.sort(); return clean_text(" ".join(word for _,word in positions))


def request_json(url: str, timeout: int, retries: int) -> dict[str, Any]:
    last=None
    for attempt in range(retries):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":USER_AGENT,"Accept":"application/json"})
            with urllib.request.urlopen(req,timeout=timeout) as response: payload=json.load(response)
            if not isinstance(payload,dict): raise RuntimeError("OpenAlex response was not a JSON object")
            return payload
        except Exception as exc:  # noqa: BLE001
            last=exc
            if attempt+1<retries: time.sleep(2**attempt)
    raise RuntimeError(url) from last


def retrieve_query(query_id: str,oql: str,api_key: str,timeout: int,retries: int):
    cursor="*"; rows=[]; reported=None; pages=0
    while True:
        params={"oql":oql,"per-page":100,"cursor":cursor,
                "select":"id,doi,title,display_name,publication_year,publication_date,type,language,cited_by_count,abstract_inverted_index,primary_location"}
        if api_key: params["api_key"]=api_key
        payload=request_json("https://api.openalex.org/?"+urllib.parse.urlencode(params),timeout,retries); pages+=1
        meta=payload.get("meta") if isinstance(payload.get("meta"),dict) else {}
        if reported is None:
            try: reported=int(meta.get("count"))
            except (TypeError,ValueError): reported=None
        results=payload.get("results")
        if not isinstance(results,list) or not results: break
        for item in results:
            if not isinstance(item,dict): continue
            primary=item.get("primary_location") if isinstance(item.get("primary_location"),dict) else {}
            source=primary.get("source") if isinstance(primary.get("source"),dict) else {}
            rows.append({"query_id":query_id,"record_id":clean_text(item.get("id")),"doi":normalize_doi(item.get("doi")),
                         "title":clean_text(item.get("title") or item.get("display_name")),"abstract":reconstruct_abstract(item.get("abstract_inverted_index")),
                         "year":item.get("publication_year") or "","publication_date":clean_text(item.get("publication_date")),
                         "work_type":clean_text(item.get("type")),"language":clean_text(item.get("language")),
                         "cited_by_count":item.get("cited_by_count") or 0,"journal":clean_text(source.get("display_name")),
                         "url":clean_text(primary.get("landing_page_url") or item.get("id"))})
        nxt=meta.get("next_cursor")
        if not nxt or str(nxt)==cursor: break
        cursor=str(nxt)
        if pages%20==0: print({"query_id":query_id,"pages":pages,"retrieved":len(rows),"reported":reported},flush=True)
        time.sleep(0.02)
    return rows,{"query_id":query_id,"reported_count":reported,"retrieved_count":len(rows),"pages":pages,
                 "truncated":bool(reported is not None and len(rows)<reported)}


def write_csv(path: Path,rows:list[dict[str,Any]],fields:list[str]|None=None):
    path.parent.mkdir(parents=True,exist_ok=True)
    if not rows:
        with path.open("w",newline="",encoding="utf-8") as h: csv.DictWriter(h,fieldnames=fields or ["empty"]).writeheader()
        return
    with path.open("w",newline="",encoding="utf-8") as h:
        w=csv.DictWriter(h,fieldnames=fields or list(rows[0])); w.writeheader(); w.writerows(rows)


def historical_source_key(source_id:str):
    doi=normalize_doi(source_id)
    if doi.startswith("10."): return "doi",doi
    oa=normalize_openalex(source_id)
    if oa: return "openalex",oa
    return "raw",clean_text(source_id).lower()


def resolve_historical_source(source_id:str,api_key:str,timeout:int,retries:int):
    kind,value=historical_source_key(source_id)
    if kind=="doi": params={"filter":"doi:"+value,"per-page":1,"select":"id,doi,title,publication_year"}
    elif kind=="openalex": params={"filter":"openalex:"+value,"per-page":1,"select":"id,doi,title,publication_year"}
    else: return None
    if api_key: params["api_key"]=api_key
    data=request_json("https://api.openalex.org/works?"+urllib.parse.urlencode(params),timeout,retries)
    results=data.get("results")
    return results[0] if isinstance(results,list) and results and isinstance(results[0],dict) else None


def title_version_match(species:str,source_title:str,dedup_rows:list[dict[str,Any]]):
    target=normalize_title(source_title)
    if not target: return None,0.0
    tokens=[x.lower() for x in species.split()[:2]]
    candidates=[]
    for row in dedup_rows:
        nt=normalize_title(row.get("title"))
        if tokens and not all(tok in nt.split() for tok in tokens): continue
        score=SequenceMatcher(None,target,nt).ratio()
        candidates.append((score,row))
    if not candidates: return None,0.0
    score,row=max(candidates,key=lambda x:x[0])
    # >=0.80 is deliberately conservative and requires the focal binomial in the title.
    return (row if score>=0.80 else None),score


def main():
    p=argparse.ArgumentParser(); p.add_argument("--config",default="literature/itv_fcp_search_config_v2.json")
    p.add_argument("--historical-manifest",default="docs/supporting/frozen_classification_manifest.csv")
    p.add_argument("--outdir",required=True); p.add_argument("--api-key-env",default="OPENALEX_API_KEY")
    p.add_argument("--timeout",type=int,default=45); p.add_argument("--retries",type=int,default=3); a=p.parse_args()
    config=json.loads(Path(a.config).read_text(encoding="utf-8")); queries=config.get("queries") or []
    if len(queries)!=15: raise SystemExit(f"Expected 15 v2 query blocks, found {len(queries)}")
    api_key=os.environ.get(a.api_key_env,"").strip(); outdir=Path(a.outdir); outdir.mkdir(parents=True,exist_ok=True)

    raw=[]; logs=[]
    for item in queries:
        rows,log=retrieve_query(str(item["id"]),str(item["oql"]),api_key,a.timeout,a.retries); raw.extend(rows); logs.append(log); print(log,flush=True)
    dedup={}; membership={}
    for row in raw:
        key=stable_key(row); membership.setdefault(key,set()).add(str(row["query_id"])); prev=dedup.get(key)
        if prev is None or len(str(row.get("abstract") or ""))>len(str(prev.get("abstract") or "")): dedup[key]=dict(row)
    dedup_rows=[]
    for key,row in sorted(dedup.items()):
        out=dict(row); out["dedup_key"]=key; out["query_ids"]=";".join(sorted(membership[key])); dedup_rows.append(out)

    with Path(a.historical_manifest).open(newline="",encoding="utf-8") as h: historical=list(csv.DictReader(h))
    if len(historical)!=34: raise SystemExit("Historical manifest must contain 34 rows")
    dois={normalize_doi(r.get("doi")) for r in dedup_rows if normalize_doi(r.get("doi"))}
    oa_ids={normalize_openalex(r.get("record_id")) for r in dedup_rows if normalize_openalex(r.get("record_id"))}
    recovery=[]
    for row in historical:
        species=clean_text(row.get("canonical_name")); source_id=clean_text(row.get("source_id")); kind,value=historical_source_key(source_id)
        exact=(kind=="doi" and value in dois) or (kind=="openalex" and value in oa_ids)
        source_work=resolve_historical_source(source_id,api_key,a.timeout,a.retries)
        source_title=clean_text(source_work.get("title")) if source_work else ""
        version,score=title_version_match(species,source_title,dedup_rows)
        version_hit=bool(version) and not exact
        robust=bool(exact or version)
        recovery.append({"canonical_name":species,"historical_source_id":source_id,"source_kind":kind,
                         "historical_source_openalex_id":normalize_openalex(source_work.get("id")) if source_work else "",
                         "historical_source_title":source_title,"recovered_exact_identifier":int(exact),
                         "title_version_similarity":round(score,4),"recovered_title_version_match":int(version_hit),
                         "matched_v2_title":clean_text(version.get("title")) if version else "",
                         "matched_v2_doi":normalize_doi(version.get("doi")) if version else "",
                         "recovered_in_openalex_v2_direct_queries":int(exact),
                         "recovered_direct_exact_or_version":int(robust)})
        time.sleep(0.02)

    write_csv(outdir/"openalex_oql_v2_raw.csv",raw); write_csv(outdir/"openalex_oql_v2_deduplicated.csv",dedup_rows)
    write_csv(outdir/"openalex_oql_v2_query_log.csv",logs); write_csv(outdir/"openalex_oql_v2_historical34_recovery.csv",recovery)
    truncated=[x for x in logs if x["truncated"]]; exact_count=sum(r["recovered_exact_identifier"] for r in recovery)
    robust_count=sum(r["recovered_direct_exact_or_version"] for r in recovery)
    summary={"status":"complete","protocol_version":config.get("protocol_version"),"query_blocks":len(logs),
             "raw_query_memberships":len(raw),"deduplicated_works":len(dedup_rows),"duplicates_removed":len(raw)-len(dedup_rows),
             "truncated_v2_query_blocks":len(truncated),"truncated_query_ids":[x["query_id"] for x in truncated],
             "historical_34_exact_sources_recovered_by_direct_v2_queries":exact_count,
             "historical_34_direct_exact_or_version_recovered":robust_count,
             "historical_34_exact_identifier_misses":[r["canonical_name"] for r in recovery if not r["recovered_exact_identifier"]],
             "historical_34_robust_direct_misses":[r["canonical_name"] for r in recovery if not r["recovered_direct_exact_or_version"]],
             "crossref_role":config.get("crossref_role"),
             "semantic_guard":"V2 retrieval measures search coverage only; title-version matching is a benchmark metadata diagnostic, not biological evidence."}
    (outdir/"openalex_oql_v2_retrieval_summary.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)
    if truncated: raise SystemExit(f"V2 direct query retrieval unexpectedly truncated: {truncated}")

if __name__=="__main__": main()
