#!/usr/bin/env python3
"""Exploratory local colour detection x geography x thermal synthesis (2026-11-27).
Only reads the archived source Table S10. No new photographs or model fits.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
from statistics import mean

SOURCE = Path("manuscripts/climate_colour/Table_S10_common_species.csv.gz")
BLOB = "887991a51009461a640f4cd2ab4b13583f53516b"
N_BOOT = 3999
SEED = 20261127
KEYS = ("geo_detected", "geo_not_detected", "geo_diff",
        "thermal_detected", "thermal_not_detected", "thermal_diff",
        "geo_positive_fraction_detected", "both_positive_fraction_detected")

def load(path):
    raw = path.read_bytes()
    sha = hashlib.sha1(f"blob {len(raw)}\0".encode() + raw).hexdigest()
    if sha != BLOB:
        raise ValueError(f"Source Git blob mismatch: {sha}")
    rows = list(csv.DictReader(io.StringIO(gzip.decompress(raw).decode("utf-8"))))
    assert len(rows) == len({r["inat_taxon_id"] for r in rows}) == 4475
    for r in rows:
        for key in ("geo", "thermal"):
            r[key] = float(r[key]) if r[key] else None
            if r[key] is not None and not math.isfinite(r[key]):
                raise ValueError(f"Non-finite {key}")
        for key in ("repeated", "common", "local5_opportunity", "local5_diff"):
            assert r[key] in ("True", "False")
            r[key] = r[key] == "True"
        assert not r["local5_diff"] or r["local5_opportunity"]
    assert sum(r["repeated"] for r in rows) == 2400
    return rows

def summarize(rows):
    out = {"n":len(rows),"genera":len({r["genus"] for r in rows})}
    for k in ("geo","thermal"):
        v = [r[k] for r in rows if r[k] is not None]
        out[k+"_n"] = len(v)
        out[k+"_mean"] = mean(v) if v else None
        out[k+"_positive"] = sum(x>0 for x in v)
    out["common"] = sum(r["common"] for r in rows)
    out["both_positive"] = sum(r["geo"] is not None and r["thermal"] is not None and r["geo"]>0 and r["thermal"]>0 for r in rows)
    return out

def contrast(rows):
    a = [r for r in rows if r["local5_diff"]]
    b = [r for r in rows if not r["local5_diff"]]
    if not a or not b:
        return None
    g1,g0 = (mean(r["geo"] for r in g) for g in (a,b))
    t1,t0 = (mean(r["thermal"] for r in g) for g in (a,b))
    return (g1,g0,g1-g0,t1,t0,t1-t0,
            mean(r["geo"]>0 for r in a),
            mean(r["geo"]>0 and r["thermal"]>0 for r in a))

def quantile(values, p):
    s=sorted(values)
    x=(len(s)-1)*p
    i=int(x)
    return s[i]+(s[min(i+1,len(s)-1)]-s[i])*(x-i)

def analyze(rows):
    local=[r for r in rows if r["repeated"] and r["local5_opportunity"]]
    common=[r for r in local if r["common"] and r["geo"] is not None and r["thermal"] is not None]
    assert len(local)==415 and len(common)==229
    by_genus={}
    for r in common:
        by_genus.setdefault(r["genus"],[]).append(r)
    groups=list(by_genus.values())
    assert len(groups)==176
    state=SEED
    estimates=contrast(common)
    draws=[[] for _ in KEYS]
    for _ in range(N_BOOT):
        sample=[]
        for _ in groups:
            state=(1664525*state+1013904223)&0xffffffff
            sample.extend(groups[(state*len(groups))//2**32])
        val=contrast(sample)
        if val is not None:
            for dst,v in zip(draws,val):
                dst.append(v)
    stats={k:{"estimate":v,"ci95":[quantile(d,.025),quantile(d,.975)],"valid":len(d)}
           for k,v,d in zip(KEYS,estimates,draws)}
    return {"source_git_blob":BLOB,"source_n":len(rows),"repeated_n":2400,
            "local_n":len(local),"common_n":len(common),"common_genera":len(groups),
            "bootstrap":N_BOOT,"seed":SEED,
            "all_local_detected":summarize([r for r in local if r["local5_diff"]]),
            "all_local_not_detected":summarize([r for r in local if not r["local5_diff"]]),
            "common_detected":summarize([r for r in common if r["local5_diff"]]),
            "common_not_detected":summarize([r for r in common if not r["local5_diff"]]),
            "common_contrasts":stats,
            "scope":"Post hoc descriptive comparison; local detection and geographic estimator reuse source photos; no adjustment for number of local pairs; no causal or inherited-morph inference."}

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source",type=Path,default=SOURCE)
    parser.add_argument("--out",type=Path)
    args=parser.parse_args()
    result=json.dumps(analyze(load(args.source)),indent=2,allow_nan=False)+"\n"
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(result,encoding="utf-8")
    else:
        print(result,end="")
