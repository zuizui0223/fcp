#!/usr/bin/env python3
"""Direct temporal test of symmetric rare-colour increase in FCP photographs.

No imputations, propensity weights, or changes to frozen primary analyses.
The alternative is convergence toward equal WHITE / NONWHITE photo-label
frequencies, NOT a test of all possible forms of frequency-dependent selection.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
S8 = ROOT / "manuscripts/climate_colour/Table_S8_revisits.csv"
S8_BLOB = "450b0df91eb8ba343c48f8d48732223c9191f08a"
SOURCE_SHA = "8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6"
SEED = 20261126
B = 19999
KEY = ["taxon", "cell", "month", "first_year", "later_year",
       "n0", "white0", "n1", "white1"]


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def file_digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(2**20), b""):
            h.update(b)
    return h.hexdigest()


def read_s8():
    raw = S8.read_bytes()
    blob = hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()
    check(blob == S8_BLOB, "frozen S8 Git blob changed")
    with S8.open(newline="") as f:
        return [{k: int(r[k]) for k in KEY} for r in csv.DictReader(f)]


def score(n0, k0, n1, k1):
    p, q = k0/n0, k1/n1
    # Average of (1 - 2*A)*(q - C) over all ORDERED distinct pairs of
    # index-year observer labels A,C. Thus the initial frequency used as
    # rarity predictor is never reused as the baseline comparison outcome.
    return (1 - 2*p)*(q-p) - 2*p*(1-p)/(n0-1)


class XorShift32:
    def __init__(self, seed):
        self.state = seed & 0xffffffff

    def draw(self):
        s = self.state
        s ^= (s << 13) & 0xffffffff
        s ^= s >> 17
        s ^= (s << 5) & 0xffffffff
        self.state = s & 0xffffffff
        return self.state/4294967296


def test(rows):
    """One species-equal conditional randomization test.

    Under the null, observed binary photographic labels are exchangeable
    between the two years within each species-cell-calendar-month window.
    Conditioning on total white labels gives the exact hypergeometric
    allocation across the fixed sample sizes. Independent window shuffles
    are the randomization model; overlapping photographers/individuals
    between distinct windows may limit ecological interpretation.
    """
    check(rows, "no local temporal windows")
    counts = Counter(r["taxon"] for r in rows)
    ns = len(counts)
    obs, naive, expect, variance = 0., 0., 0., 0.
    informative = []
    distributions = []
    for r in rows:
        n0, k0, n1, k1 = (r[k] for k in ["n0","white0","n1","white1"])
        check(n0 >= 2 and n1 >= 2 and 0 <= k0 <= n0 and 0 <= k1 <= n1,
              "ineligible or invalid local year comparison")
        w = 1/(ns*counts[r["taxon"]])
        p, q = k0/n0, k1/n1
        obs += w*score(n0,k0,n1,k1)
        naive += w*(1-2*p)*(q-p)
        K, n = k0+k1, n0+n1
        den = math.comb(n,K)
        cs, first, second = 0., 0., 0.
        cdf, values = [], []
        for k in range(max(0,K-n1), min(n0,K)+1):
            prob = math.comb(n0,k)*math.comb(n1,K-k)/den
            x = w*score(n0,k,n1,K-k)
            cs += prob
            cdf.append(cs)
            values.append(x)
            first += prob*x
            second += prob*x*x
        check(abs(first) < 1e-12, "conditional-null expectation is not zero")
        expect += first
        variance += second-first*first
        distributions.append((cdf,values))
        if len(values)>1:
            informative.append(r)
    rng = XorShift32(SEED)
    sims = []
    exceed = 0
    for _ in range(B):
        s = 0.
        for cdf, values in distributions:
            u = rng.draw()
            for i, prob in enumerate(cdf):
                if u < prob or i == len(cdf)-1:
                    s += values[i]
                    break
        sims.append(s)
        exceed += s >= obs-1e-13
    sims.sort()
    quant = lambda p: sims[math.floor(p*(B-1))]
    return {
        "n_windows":len(rows),
        "n_species":ns,
        "n_informative_windows":len(informative),
        "n_informative_species":len({r["taxon"] for r in informative}),
        "n_initial_photos":sum(r["n0"] for r in rows),
        "n_later_photos":sum(r["n1"] for r in rows),
        "n_initial_mixed_windows":sum(0<r["white0"]<r["n0"] for r in rows),
        "n_one_year_gap":sum(r["later_year"]-r["first_year"]==1 for r in rows),
        "mean_year_gap":sum(r["later_year"]-r["first_year"] for r in rows)/len(rows),
        "naive_score_DIAGNOSTIC_ONLY":naive,
        "primary_score":obs,
        "conditional_null_mean":expect,
        "conditional_null_sd":math.sqrt(max(0.,variance)),
        "null95":[quant(.025),quant(.975)],
        "one_sided_p_for_positive_convergence":(exceed+1)/(B+1),
        "null_replicates":B,
        "seed":SEED,
        "estimand":"Species-equal mean of within-window observer-split rarity x temporal change; dimensionless. Positive = tendency toward a 50:50 white/nonwhite photo-label frequency.",
    }


def construct_all(source):
    import numpy as np
    import pandas as pd

    check(file_digest(source)==SOURCE_SHA,"frozen original 447500-photo source SHA mismatch")
    fields = ["inat_taxon_id","photo_id","observer_id","observed_on",
              "latitude","longitude","positional_accuracy_m",
              "global_classifiable","morph"]
    d = pd.read_csv(source,usecols=fields,low_memory=False)
    check(len(d)==447500 and d.photo_id.nunique()==447500 and
          d.inat_taxon_id.nunique()==4475 and
          d.groupby("inat_taxon_id").size().eq(100).all(),
          "frozen source frame changed")
    classified = d.global_classifiable.astype(str).str.lower().isin(["true","1","yes"])
    check(int(classified.sum())==227107,"classified denominator changed")
    check(d.loc[classified,"morph"].isin(
        ["white","yellow_orange","red_pink","blue_purple"]).all(),
        "source label changed")
    dt = pd.to_datetime(d.observed_on,format="mixed",utc=True,errors="coerce")
    d["month"],d["year"]=dt.dt.month,dt.dt.year
    # Same 10 arcminute grid as archived Table S8 (no climate covariates used).
    valid=classified & d.latitude.between(-90,90) & d.longitude.between(-180,180) & \
        d.month.notna() & d.positional_accuracy_m.between(0,5000) & d.observer_id.notna()
    q=d.loc[valid].copy()
    q["cell"]=(np.clip(np.floor((90-q.latitude)/(1/6)),0,1079).astype(int)*2160
               +np.clip(np.floor((q.longitude+180)/(1/6)),0,2159).astype(int))
    keys=["inat_taxon_id","month","cell","year"]
    q=q.sort_values("photo_id").drop_duplicates(keys+["observer_id"])
    q["white"]=q.morph.eq("white")
    rows=[]
    # Earliest eligible pair of years, chosen without looking at colours.
    for (tax,month,cell), g in q.groupby(["inat_taxon_id","month","cell"],sort=True):
        if len(g)<4:
            continue
        years=[(yr,h) for yr,h in g.groupby("year",sort=True) if len(h)>=2]
        chosen=False
        for i,(yr,h) in enumerate(years):
            for yr2,h2 in years[i+1:]:
                f=h2[~h2.observer_id.isin(h.observer_id)]
                if len(f)>=2:
                    vals=(tax,cell,month,yr,yr2,len(h),h.white.sum(),len(f),f.white.sum())
                    rows.append(dict(zip(KEY,map(int,vals))))
                    chosen=True
                    break
            if chosen:
                break
    # Do not quietly change the previous windows. Reconstruct archived S8
    # exactly inside the old repeated-colour taxon subset.
    freq=pd.read_csv(ROOT/"manuscripts/climate_colour/Table_S6_white_frequencies.csv.gz")
    repeated=set(freq.loc[freq.repeated.astype(str).str.lower().isin(["true","1","yes"]),
                          "inat_taxon_id"].astype(int))
    check(len(repeated)==2400,"old repeated-colour frame changed")
    old={tuple(r[k] for k in KEY) for r in read_s8()}
    reproduced={tuple(r[k] for k in KEY) for r in rows if r["taxon"] in repeated}
    check(old==reproduced,
          f"S8 local-window reconstruction mismatch: old-only={len(old-reproduced)} new-only={len(reproduced-old)}")
    return rows


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source",type=Path)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    archived=read_s8()
    pilot=test(archived)
    check(pilot["n_windows"]==104 and pilot["n_species"]==73 and
          abs(pilot["primary_score"]-0.006992009132420088)<1e-12 and
          abs(pilot["one_sided_p_for_positive_convergence"]-.4218)<1e-12,
          "frozen Table S8 calculation not reproduced")
    out={
        "schema":"fcp_direct_photographic_frequency_convergence_v1",
        "source_S8_git_blob":S8_BLOB,
        "all_species_source_sha256":SOURCE_SHA if args.source else None,
        "archived_repeated_colour_panel_exploratory":pilot,
        "full_4475_registry_primary":None,
        "limitations":[
            "Photo-region WHITE/NONWHITE labels are not genetically validated morphs or individual population censuses.",
            "Year-to-year frequencies can reflect environment, plant turnover, migration, developmental stage, photography, or persistent individuals.",
            "A 10-arcminute geographic cell is not necessarily a biological population.",
            "A positive result tests symmetric convergence toward 0.5 only; cannot identify fitness-frequency dependence or exclude different equilibria.",
            "Null is conditional temporal label exchangeability, not a population-genetic selection model.",
            "No inference that selection is absent from a nonsignificant test."
        ]
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    if args.source:
        rows=construct_all(args.source)
        out["full_4475_registry_primary"]=test(rows)
        out["status"]="HASH_VERIFIED_FULL_SOURCE_TEMPORAL_TEST"
        with (args.output.parent/"fcp_frequency_temporal_windows_all4475.csv").open("w",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=KEY)
            writer.writeheader()
            writer.writerows(rows)
    else:
        out["status"]="ARCHIVED_104_WINDOW_PILOT_ONLY"
    args.output.write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))


if __name__=="__main__":
    main()
