#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, wilcoxon

SPRING = 14.16666667
MILD = 23.75
HOT = 28.75
TOL = 1e-5

def close(x, y):
    return abs(float(x) - float(y)) <= TOL

def find_flavonol_column(columns):
    hits=[c for c in columns if re.search(r"flavonol",str(c),re.I)]
    if len(hits)!=1:
        raise RuntimeError(f"expected exactly one flavonol column, found {hits}")
    return hits[0]

def one_sequence(g: pd.DataFrame):
    order={"First":1,"Second":2,"Third":3}
    q=(g.groupby(["Period","Season","weighted.temp"],as_index=False)
         .agg(cyanidin=("Cyanidin","mean"),
              flavonol=("flavonol_value","mean"),
              n_flowers=("Cyanidin","size")))
    q["ord"]=q["Period"].map(order)
    q=q.sort_values("ord")
    if len(q)!=3:
        return None
    return tuple(float(x) for x in q["weighted.temp"]), q

def trait_summary(w, trait):
    mild=w[f"{trait}_mild"].astype(float)
    hot=w[f"{trait}_hot"].astype(float)
    delta=hot-mild
    wr=wilcoxon(delta.to_numpy(float),zero_method="wilcox",alternative="two-sided")
    return {
        "mean_mild":float(mild.mean()),
        "mean_hot":float(hot.mean()),
        "cohort_mean_hot_over_mild":float(hot.mean()/mild.mean()),
        "mean_delta_hot_minus_mild":float(delta.mean()),
        "median_delta_hot_minus_mild":float(delta.median()),
        "individuals_decreasing":int((delta<0).sum()),
        "individuals_increasing":int((delta>0).sum()),
        "two_sided_wilcoxon_p":float(wr.pvalue),
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-data",required=True)
    p.add_argument("--outdir",required=True)
    a=p.parse_args()
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)

    raw=pd.read_excel(a.source_data,sheet_name="Figure 3F-J",engine="openpyxl")
    flav_col=find_flavonol_column(raw.columns)
    required={"Individual","Period","Season","Cyanidin","weighted.temp",flav_col}
    missing=required-set(raw.columns)
    if missing:
        raise SystemExit(f"missing columns {sorted(missing)}")

    d=raw[["Individual","Period","Season","Cyanidin","weighted.temp",flav_col]].copy()
    d=d.rename(columns={flav_col:"flavonol_value"})
    for c in ["Individual","Cyanidin","weighted.temp","flavonol_value"]:
        d[c]=pd.to_numeric(d[c],errors="coerce")
    d=d.dropna(subset=["Individual","Period","Cyanidin","flavonol_value","weighted.temp"]).copy()

    seqs={}
    qs={}
    for ind,g in d.groupby("Individual"):
        z=one_sequence(g)
        if z is not None:
            seqs[float(ind)]=z[0]
            qs[float(ind)]=z[1]

    ids=[ind for ind,s in seqs.items()
         if len(s)==3 and all(close(x,y) for x,y in zip(s,(SPRING,MILD,HOT)))]
    if len(ids)!=15:
        raise SystemExit(f"expected 15 mild->hot individuals, found {len(ids)}")

    rows=[]
    for ind in ids:
        q=qs[ind].set_index("Period")
        rows.append({
            "Individual":ind,
            "cyanidin_mild":float(q.loc["Second","cyanidin"]),
            "cyanidin_hot":float(q.loc["Third","cyanidin"]),
            "flavonol_mild":float(q.loc["Second","flavonol"]),
            "flavonol_hot":float(q.loc["Third","flavonol"]),
        })
    w=pd.DataFrame(rows).sort_values("Individual").reset_index(drop=True)

    if (w[["cyanidin_mild","cyanidin_hot","flavonol_mild","flavonol_hot"]]<=0).any().any():
        bad=w.loc[(w[["cyanidin_mild","cyanidin_hot","flavonol_mild","flavonol_hot"]]<=0).any(axis=1)]
        raise SystemExit(f"nonpositive trait values block log-ratio analysis: {bad.to_dict('records')}")

    w["cyanidin_hot_over_mild"]=w.cyanidin_hot/w.cyanidin_mild
    w["flavonol_hot_over_mild"]=w.flavonol_hot/w.flavonol_mild
    w["log_ratio_cyanidin"]=np.log(w.cyanidin_hot_over_mild)
    w["log_ratio_flavonol"]=np.log(w.flavonol_hot_over_mild)
    w["branch_specificity_B"]=w.log_ratio_cyanidin-w.log_ratio_flavonol

    B=w.branch_specificity_B.to_numpy(float)
    wr1=wilcoxon(B,zero_method="wilcox",alternative="less")
    wr2=wilcoxon(B,zero_method="wilcox",alternative="two-sided")
    nonzero=B[B!=0]
    successes=int((nonzero<0).sum())
    sign=binomtest(successes,len(nonzero),0.5,alternative="greater")

    support=bool(
        len(w)==15 and
        float(np.median(B))<0 and
        int((B<0).sum())>=12 and
        float(wr1.pvalue)<0.05
    )

    result={
        "schema":"fcp_moricandia_anthocyanin_branch_specificity_v1",
        "status":"complete",
        "role":"post_publication_mechanistic_specificity_followup",
        "source_sheet":"Figure 3F-J",
        "flavonol_source_column":str(flav_col),
        "n_individuals":int(len(w)),
        "cyanidin":trait_summary(w,"cyanidin"),
        "flavonol":trait_summary(w,"flavonol"),
        "branch_specificity":{
            "definition":"log(cyanidin_hot/cyanidin_mild) - log(flavonol_hot/flavonol_mild)",
            "mean_B":float(np.mean(B)),
            "median_B":float(np.median(B)),
            "individuals_B_lt_0":int((B<0).sum()),
            "fraction_B_lt_0":float((B<0).mean()),
            "wilcoxon_one_sided_p_B_lt_0":float(wr1.pvalue),
            "wilcoxon_two_sided_p":float(wr2.pvalue),
            "sign_test_directional_p":float(sign.pvalue),
        },
        "verdict":"ANTHOCYANIN_BRANCH_SPECIFIC_HEAT_RESPONSE_SUPPORTED" if support else "ANTHOCYANIN_BRANCH_SPECIFICITY_NOT_SUPPORTED_UNDER_THIS_TEST",
        "interpretation":"A negative branch-specificity contrast means anthocyanin declines proportionally more strongly than UV-absorbing flavonols under hotter summer conditions with the stated summer photoperiod held constant.",
        "hard_nonclaims":[
            "post-publication and not an untouched confirmatory test",
            "does not directly measure metabolic flux with isotope tracing",
            "does not rescue the failed cross-cohort BIO5 replication",
            "does not establish that all natural white-flower systems are thermally plastic"
        ]
    }
    w.to_csv(out/"individual_branch_ratios.csv",index=False)
    (out/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))

if __name__=="__main__":
    main()
