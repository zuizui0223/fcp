#!/usr/bin/env python3
"""Strict finalizer for the frozen RGFCA-v2 synthetic sharedness qualification."""
from __future__ import annotations

import argparse
import json
from math import sqrt
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
MAPPING = ROOT / "docs/supporting/rgfca_sharedness_v2_synthetic_technical_mapping_v1.json"
GATE = ROOT / "docs/supporting/rgfca_sharedness_v2_synthetic_gate_amendment_v1.json"


def wilson(k:int,n:int,z:float=1.959963984540054)->tuple[float,float]:
    if n<=0: return float("nan"),float("nan")
    p=k/n; d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*sqrt((p*(1-p)+z*z/(4*n))/n)/d
    return max(0,c-h),min(1,c+h)


def read_stage(root:Path,stage:str,expected_arms:int)->pd.DataFrame:
    frames=[]
    for arm in range(expected_arms):
        matches=list(root.rglob(f"{stage}_arm_{arm:03d}.csv"))
        if len(matches)!=1: raise RuntimeError(f"{stage} arm {arm} missing/duplicated")
        f=pd.read_csv(matches[0])
        if len(f)!=250 or f["replicate"].tolist()!=list(range(250)) or set(f["stage"].astype(str))!={stage} or set(f["arm_index"].astype(int))!={arm}:
            raise RuntimeError(f"bad census {stage} arm {arm}")
        if not np.isfinite(pd.to_numeric(f["world_score"],errors="raise")).all(): raise RuntimeError("nonfinite score")
        frames.append(f)
    out=pd.concat(frames,ignore_index=True)
    if out.duplicated(["stage","arm_index","replicate"]).any(): raise RuntimeError("duplicate world")
    return out


def main()->int:
    p=argparse.ArgumentParser(); p.add_argument("--input-root",type=Path,required=True); p.add_argument("--output-json",type=Path,required=True); p.add_argument("--output-csv",type=Path,required=True); a=p.parse_args()
    mapping=json.loads(MAPPING.read_text()); gate=json.loads(GATE.read_text())
    cal=read_stage(a.input_root,"calibration",21); ev=read_stage(a.input_root,"evaluation",78)
    rows=[]; q=[]
    for arm,g in cal.groupby("arm_index",sort=True):
        val=float(np.quantile(g["world_score"].to_numpy(float),0.975,method="higher")); q.append(val)
        rows.append({"stage":"calibration","arm_index":int(arm),"arm_id":str(g["arm_id"].iloc[0]),"kind":"nuisance","retention":float(g["retention"].iloc[0]),"n":250,"threshold_or_rate":val,"metric":"q975"})
    threshold=float(max(q))
    summaries=[]
    for arm,g in ev.groupby("arm_index",sort=True):
        score=g["world_score"].to_numpy(float); reject=score>threshold; k=int(reject.sum()); rate=k/len(g); lo,hi=wilson(k,len(g))
        d={"stage":"evaluation","arm_index":int(arm),"arm_id":str(g["arm_id"].iloc[0]),"kind":str(g["kind"].iloc[0]),"scenario_id":str(g["scenario_id"].iloc[0]),"retention":float(g["retention"].iloc[0]),"amplitude":float(g["amplitude"].iloc[0]) if np.isfinite(float(g["amplitude"].iloc[0])) else None,"shared_fraction":float(g["shared_fraction"].iloc[0]),"threshold_sd":float(g["threshold_sd"].iloc[0]),"n":250,"rejections":k,"rate":rate,"wilson95_low":lo,"wilson95_high":hi,"mean_score":float(score.mean())}
        summaries.append(d); rows.append({**d,"metric":"rejection_rate","threshold_or_rate":rate})
    nuisance=[x for x in summaries if x["kind"]=="nuisance"]
    nuisance_pass=all(x["rate"]<=float(mapping["hard_gates"]["every_evaluation_nuisance_rejection_lte"]) for x in nuisance)
    def positive_gate(amp:float,share:float,floor:float):
        m=[x for x in summaries if x["kind"]=="positive" and x["retention"]==0.4 and x["amplitude"]==amp and x["shared_fraction"]==share and x["threshold_sd"] in (0.0,0.25)]
        if len(m)!=2: raise RuntimeError(f"hard positive cell missing amp={amp} share={share}")
        return {"amplitude":amp,"shared_fraction":share,"retention":0.4,"floor":floor,"sd_specific":[{"threshold_sd":x["threshold_sd"],"power":x["rate"],"pass":x["rate"]>=floor} for x in sorted(m,key=lambda z:z["threshold_sd"])],"minimum_power":min(x["rate"] for x in m),"pass":all(x["rate"]>=floor for x in m)}
    gates=[positive_gate(1.0,1.0,0.80),positive_gate(1.0,0.5,0.50),positive_gate(2.0,0.5,0.80)]
    passed=bool(nuisance_pass and all(x["pass"] for x in gates))
    result={"protocol":mapping["protocol"],"status":"complete_rgfca_sharedness_v2_synthetic_qualification","qualification_pass":passed,"global_threshold":threshold,"calibration_worlds":len(cal),"evaluation_worlds":len(ev),"total_worlds":len(cal)+len(ev),"nuisance_gate":{"maximum_rejection_rate":max(x["rate"] for x in nuisance),"ceiling":0.10,"pass":nuisance_pass},"hard_positive_gates":gates,"evaluation_arm_summaries":summaries,"observed_flower_colour_opened":False,"observed_background_colour_opened":False,"image_pixels_opened":False,"ecological_claim_changed":False,"image_acquisition_permitted":False,"next_step_if_pass":"synthetic PASS alone still does not authorize image acquisition; a separate prospective acquisition/measurement firewall is required","next_step_if_fail":"retain FAIL; do not search species count, photo count, axis grid, threshold grid, statistic or power floors post outcome","six_species_used":False,"thirty_four_species_used":False,"claim_ceiling":"Synthetic method qualification only; no observed ecological evidence."}
    a.output_json.parent.mkdir(parents=True,exist_ok=True); a.output_csv.parent.mkdir(parents=True,exist_ok=True)
    a.output_json.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    pd.DataFrame(rows).to_csv(a.output_csv,index=False,lineterminator="\n")
    print(json.dumps({k:result[k] for k in ["status","qualification_pass","global_threshold","calibration_worlds","evaluation_worlds","nuisance_gate","hard_positive_gates"]},indent=2))
    return 0

if __name__=="__main__": raise SystemExit(main())
