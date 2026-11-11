#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import io
import json
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
DESIGN=ROOT/"data"/"hadfield_calibration_source_recovery_design_v0_1.json"
UA="fcp-hadfield-source-recovery/0.1 (+https://github.com/zuizui0223/fcp)"


def fetch(url:str, accept:str="*/*"):
    req=urllib.request.Request(url,headers={"User-Agent":UA,"Accept":accept,"X-API-Version":"2.1.0"})
    try:
        with urllib.request.urlopen(req,timeout=90) as r:
            return {
                "ok":True,
                "status":getattr(r,"status",200),
                "content_type":r.headers.get("Content-Type"),
                "content_disposition":r.headers.get("Content-Disposition"),
                "body":r.read(),
                "final_url":r.geturl(),
            }
    except urllib.error.HTTPError as e:
        body=e.read()
        return {"ok":False,"status":e.code,"error":str(e),"body":body[:4000]}
    except Exception as e:
        return {"ok":False,"status":None,"error":f"{type(e).__name__}: {e}","body":b""}


def digest(b:bytes):
    return {"bytes":len(b),"sha256":hashlib.sha256(b).hexdigest()}


def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--out-dir",type=Path,required=True)
    a=ap.parse_args()
    a.out_dir.mkdir(parents=True,exist_ok=True)
    d=json.loads(DESIGN.read_text())

    required={int(x["file_id"]):x for x in d["source"]["required_files"]}
    attempts=[]
    recovered={}

    # 1. Public metadata only.
    metadata={}
    for fid,spec in required.items():
        u=f"https://datadryad.org/api/v2/files/{fid}"
        r=fetch(u,"application/json")
        rec={"route":"api_file_metadata","file_id":fid,"url":u,"ok":r["ok"],"status":r["status"]}
        if r["ok"]:
            try:
                j=json.loads(r["body"].decode("utf-8"))
                metadata[fid]=j
                rec["metadata"]={
                    k:j.get(k) for k in ("id","path","size","mimeType","digest","digestType","status") if k in j
                }
            except Exception as e:
                rec["parse_error"]=str(e)
        else:
            rec["error"]=r.get("error")
            rec["body_prefix"]=r.get("body",b"").decode("utf-8","ignore")
        attempts.append(rec)

    # 2. Official API download routes.
    for fid,spec in required.items():
        u=f"https://datadryad.org/api/v2/files/{fid}/download"
        r=fetch(u)
        rec={"route":"api_file_download","file_id":fid,"url":u,"ok":r["ok"],"status":r["status"]}
        if r["ok"]:
            rec.update(digest(r["body"]))
            rec["content_type"]=r.get("content_type")
            rec["content_disposition"]=r.get("content_disposition")
            recovered[fid]=r["body"]
        else:
            rec["error"]=r.get("error")
            rec["body_prefix"]=r.get("body",b"").decode("utf-8","ignore")
        attempts.append(rec)

    # 3. Exact file_stream route only for still-missing authoritative file IDs.
    for fid,spec in required.items():
        if fid in recovered:
            continue
        u=f"https://datadryad.org/downloads/file_stream/{fid}"
        r=fetch(u)
        rec={"route":"landing_file_stream","file_id":fid,"url":u,"ok":r["ok"],"status":r["status"]}
        if r["ok"]:
            rec.update(digest(r["body"]))
            rec["content_type"]=r.get("content_type")
            rec["final_url"]=r.get("final_url")
            # Reject HTML/challenge bodies even if HTTP status is 200.
            ctype=(r.get("content_type") or "").lower()
            if "text/html" in ctype:
                rec["ok"]=False
                rec["rejected"]="html_challenge_or_landing_body"
            else:
                recovered[fid]=r["body"]
        else:
            rec["error"]=r.get("error")
            rec["body_prefix"]=r.get("body",b"").decode("utf-8","ignore")
        attempts.append(rec)

    # 4. Exact DOI full dataset as final authoritative route, no third-party mirrors.
    missing=[fid for fid in required if fid not in recovered]
    if missing:
        doi="doi%3A10.5061%2Fdryad.jf3tj"
        u=f"https://datadryad.org/api/v2/datasets/{doi}/download"
        r=fetch(u)
        rec={"route":"api_dataset_download","url":u,"ok":r["ok"],"status":r["status"]}
        if r["ok"]:
            rec.update(digest(r["body"]))
            rec["content_type"]=r.get("content_type")
            try:
                z=zipfile.ZipFile(io.BytesIO(r["body"]))
                names=z.namelist()
                rec["zip_members"]=names
                byname={Path(n).name:n for n in names}
                for fid in list(missing):
                    path=required[fid]["path"]
                    if path in byname:
                        recovered[fid]=z.read(byname[path])
                rec["required_members_recovered"]=[
                    required[fid]["path"] for fid in missing if fid in recovered
                ]
            except Exception as e:
                rec["zip_parse_error"]=f"{type(e).__name__}: {e}"
        else:
            rec["error"]=r.get("error")
            rec["body_prefix"]=r.get("body",b"").decode("utf-8","ignore")
        attempts.append(rec)

    files=[]
    for fid,spec in required.items():
        b=recovered.get(fid)
        row={
            "file_id":fid,
            "expected_path":spec["path"],
            "role":spec["role"],
            "recovered":b is not None,
            "authoritative_metadata":metadata.get(fid),
        }
        if b is not None:
            row.update(digest(b))
            p=a.out_dir/spec["path"]
            p.write_bytes(b)
        files.append(row)

    all_ok=all(x["recovered"] for x in files)
    out={
        "version":"v0.1",
        "status":"HADFIELD_EXACT_SOURCE_BYTES_RECOVERED" if all_ok else "HOLD_HADFIELD_EXACT_SOURCE_BYTES_INACCESSIBLE",
        "design":"data/hadfield_calibration_source_recovery_design_v0_1.json",
        "required_file_count":len(files),
        "recovered_file_count":sum(int(x["recovered"]) for x in files),
        "files":files,
        "attempts":attempts,
        "third_party_mirror_used":False,
        "substitute_source_used":False,
        "biological_values_analyzed":False,
        "calibration_outcome_opened":False,
        "next_gate":"freeze retrospective calibration estimator before value analysis" if all_ok else "close Hadfield calibration source route; proceed to next finite interaction source",
        "interpretation":"This is an access/provenance result only, not evidence about interaction memory or spatial rewiring."
    }
    (a.out_dir/"result.json").write_text(json.dumps(out,indent=2,sort_keys=True,default=str)+"\n")
    print(json.dumps({
        "status":out["status"],
        "recovered_file_count":out["recovered_file_count"],
        "files":[{k:x.get(k) for k in ("file_id","expected_path","recovered","bytes","sha256")} for x in files],
        "attempts":[{k:x.get(k) for k in ("route","file_id","status","ok","error","rejected")} for x in attempts],
        "next_gate":out["next_gate"]
    },indent=2))
    return 0

if __name__=="__main__":
    raise SystemExit(main())
