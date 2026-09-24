from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]


def load(name: str, rel: str):
    path=ROOT/rel
    spec=importlib.util.spec_from_file_location(name,path)
    assert spec and spec.loader
    mod=importlib.util.module_from_spec(spec)
    sys.modules[name]=mod
    spec.loader.exec_module(mod)
    return mod


MV=load("fcp_v2_mv1_mv3_test", "scripts/analysis/analyze_fcp_v2_mv1_mv3_20260924.py")
MV4=load("fcp_v2_mv4_test", "scripts/analysis/analyze_fcp_v2_mv4_20260924.py")


def synthetic() -> pd.DataFrame:
    rows=[]
    for si in range(4):
        species=f"Species_{si}"
        panel="P" if si<2 else "N"
        white_n=20+si*10
        span=1.0+si
        for i in range(100):
            morph="white" if i<white_n else "red_pink"
            base=np.full(9,0.001,dtype=float)
            if morph=="white":
                base[0]=0.992
            else:
                base[3]=0.60
                base[4]=0.392
            base=base/base.sum()
            rec={
                "measurement_id":f"M{si:02d}{i:03d}",
                "species":species,
                "panel":panel,
                "latitude":float(si)+0.01*i,
                "longitude":span*(i/99.0),
                "heavy_counterfactual":i<20,
                "source_identity_status":"exact_source_sha_match",
                "stratum_roi_unstable_fixed": i in {0,1},
                "base_status":"classified_four_state_morph",
                "base_morph":morph,
            }
            for j,c in enumerate(MV.BIO):
                rec[f"base_fraction_{c}"]=float(base[j])
            for prefix in ("fixed_ev_m0_5","fixed_ev_p0_0","fixed_ev_p0_5","fixed_ev_m1_0","fixed_ev_p1_0"):
                rec[f"{prefix}_status"]="classified_four_state_morph"
                rec[f"{prefix}_morph"]=morph
                for j,c in enumerate(MV.BIO):
                    rec[f"{prefix}_fraction_{c}"]=float(base[j])
            rows.append(rec)
    return pd.DataFrame(rows)


def test_mv1_identity_and_D_identity():
    df=synthetic()
    s=MV.summarize_image_condition(df,"fixed_ev_p0_0")
    assert s["paired_classifiable_rows"]==400
    assert s["coarse_state_flip_probability_equal_image"]==0.0
    assert s["hellinger_mean_equal_image"]==0.0
    d=MV.D_comparison(df,"fixed_ev_p0_0")
    assert d["pooled"]["paired_species"]==4
    assert d["pooled"]["mean_absolute_D_difference"]==0.0


def test_mv3_vectors_construct():
    df=synthetic()
    tab,gates=MV.vector_table(df,"base")
    assert gates["nclass_pass"]==4
    assert len(tab)==4
    w=MV.W_from_vectors(tab)
    assert w is not None and 0.0 <= w <= 1.0


def test_mv4_preserves_100_row_geometry_and_roi_subset():
    df=synthetic()
    MV4.N_PERM=9
    base=MV4.run_condition(df,"base")
    stable=MV4.run_condition(df,"roi_stable_only")
    assert base["eligible_species"]==4
    assert stable["eligible_species"]==4
    assert base["null"].shape==(4,9)
    assert stable["null"].shape==(4,9)
    assert set(base["species_table"]["n_classifiable"])=={100}
    assert set(stable["species_table"]["n_classifiable"])=={98}
