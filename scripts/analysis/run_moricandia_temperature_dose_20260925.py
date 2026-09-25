#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, mannwhitneyu, wilcoxon

SPRING = 14.16666667
MILD = 23.75
HOT = 28.75
TOL = 1e-5

def close(x, y):
    return abs(float(x) - float(y)) <= TOL

def one_sequence(g: pd.DataFrame):
    order={"First":1,"Second":2,"Third":3}
    q=(g.groupby(["Period","Season","weighted.temp"],as_index=False)
         .agg(cyanidin=("Cyanidin","mean"),
              n_flowers=("Cyanidin","size")))
    q["ord"]=q["Period"].map(order)
    q=q.sort_values("ord")
    if len(q)!=3:
        return None
    return tuple(float(x) for x in q["weighted.temp"]), q

def paired_block(pm: pd.DataFrame, ids, p1: str, p2: str, alternative: str):
    w=pm.loc[pm.Individual.isin(ids)].pivot(index="Individual",columns="Period",values="cyanidin")
    w=w.dropna(subset=[p1,p2])
    start=w[p1].astype(float)
    end=w[p2].astype(float)
    delta=end-start
    nonzero=delta[delta.ne(0)]
    wr=wilcoxon(delta.to_numpy(float),zero_method="wilcox",alternative=alternative)
    wr2=wilcoxon(delta.to_numpy(float),zero_method="wilcox",alternative="two-sided")
    if alternative=="less":
        successes=int((nonzero<0).sum())
    elif alternative=="greater":
        successes=int((nonzero>0).sum())
    else:
        successes=max(int((nonzero<0).sum()),int((nonzero>0).sum()))
    sign=binomtest(successes,int(len(nonzero)),0.5,alternative="greater") if len(nonzero) else None
    return {
        "n_individuals":int(len(delta)),
        "mean_start":float(start.mean()),
        "mean_end":float(end.mean()),
        "mean_end_over_start":float(end.mean()/start.mean()),
        "mean_delta":float(delta.mean()),
        "median_delta":float(delta.median()),
        "negative_deltas":int((delta<0).sum()),
        "positive_deltas":int((delta>0).sum()),
        "zero_deltas":int((delta==0).sum()),
        "fraction_negative":float((delta<0).mean()),
        "wilcoxon_one_sided_alternative":alternative,
        "wilcoxon_one_sided_statistic":float(wr.statistic),
        "wilcoxon_one_sided_p":float(wr.pvalue),
        "wilcoxon_two_sided_p":float(wr2.pvalue),
        "sign_test_directional_p":float(sign.pvalue) if sign else np.nan,
        "deltas":{str(int(k) if float(k).is_integer() else k):float(v) for k,v in delta.items()},
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-data",required=True)
    p.add_argument("--outdir",required=True)
    a=p.parse_args()
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)

    raw=pd.read_excel(a.source_data,sheet_name="Figure 3F-J",engine="openpyxl")
    required={"Individual","Period","Season","Cyanidin","weighted.temp"}
    missing=required-set(raw.columns)
    if missing:
        raise SystemExit(f"missing columns {sorted(missing)}")
    d=raw[list(required)].copy()
    for c in ["Individual","Cyanidin","weighted.temp"]:
        d[c]=pd.to_numeric(d[c],errors="coerce")
    d=d.dropna(subset=["Individual","Period","Cyanidin","weighted.temp"]).copy()

    pm=(d.groupby(["Individual","Period","Season","weighted.temp"],as_index=False)
          .agg(cyanidin=("Cyanidin","mean"),n_flowers=("Cyanidin","size")))

    seqs={}
    for ind,g in d.groupby("Individual"):
        z=one_sequence(g)
        if z is not None:
            seqs[float(ind)]=z[0]

    def ids_for(target):
        ids=[]
        for ind,s in seqs.items():
            if len(s)==3 and all(close(x,y) for x,y in zip(s,target)):
                ids.append(ind)
        return ids

    mild_hot=ids_for((SPRING,MILD,HOT))
    mild_spring=ids_for((SPRING,MILD,SPRING))
    spring_hot=ids_for((SPRING,SPRING,HOT))

    if len(mild_hot)<10 or len(mild_spring)<10 or len(spring_hot)<10:
        raise SystemExit(f"unexpected sequence support: mild_hot={len(mild_hot)}, mild_spring={len(mild_spring)}, spring_hot={len(spring_hot)}")

    primary=paired_block(pm,mild_hot,"Second","Third","less")
    reversal=paired_block(pm,mild_spring,"Second","Third","greater")
    spring_to_hot=paired_block(pm,spring_hot,"Second","Third","less")

    def deltas(ids):
        w=pm.loc[pm.Individual.isin(ids)].pivot(index="Individual",columns="Period",values="cyanidin")
        return (w["Third"]-w["Second"]).dropna().astype(float)

    dh=deltas(mild_hot)
    dr=deltas(mild_spring)
    mw=mannwhitneyu(dh.to_numpy(float),dr.to_numpy(float),alternative="two-sided")

    support=bool(primary["mean_delta"]<0 and primary["median_delta"]<0 and primary["wilcoxon_one_sided_p"]<0.05)
    result={
        "schema":"fcp_moricandia_temperature_dose_reanalysis_v1",
        "status":"complete",
        "role":"post_publication_post_preflight_temperature_dose_reanalysis",
        "source_sheet":"Figure 3F-J",
        "temperature_sequences_celsius":{
            "mild_to_hot":[SPRING,MILD,HOT],
            "mild_to_spring":[SPRING,MILD,SPRING],
            "spring_to_hot":[SPRING,SPRING,HOT],
        },
        "source_sequence_counts":{
            "mild_to_hot":len(mild_hot),
            "mild_to_spring":len(mild_spring),
            "spring_to_hot":len(spring_hot),
        },
        "primary_mild_to_hot":primary,
        "reversal_mild_to_spring":reversal,
        "spring_to_hot":spring_to_hot,
        "period2_to_period3_change_mild_hot_vs_mild_spring":{
            "mannwhitney_U":float(mw.statistic),
            "two_sided_p":float(mw.pvalue),
            "mean_change_mild_to_hot":float(dh.mean()),
            "mean_change_mild_to_spring":float(dr.mean()),
        },
        "temperature_dose_support":support,
        "interpretation":"Within the source-data subgroup measured at mild-summer and then hotter-summer conditions under the same stated summer photoperiod, cyanidin declines further at the hotter temperature. Reversal and spring-to-hot sequences provide order/context checks.",
        "hard_nonclaims":[
            "post-publication and post-preflight, not an untouched confirmatory test",
            "does not rescue the failed cross-cohort FCP BIO5 replication",
            "sequential mild-to-hot measurements can still contain period/order carryover",
            "does not prove all white-flower systems are thermally plastic",
        ],
    }
    pm.to_csv(out/"individual_period_cyanidin_means.csv",index=False)
    (out/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))

if __name__=="__main__":
    main()
