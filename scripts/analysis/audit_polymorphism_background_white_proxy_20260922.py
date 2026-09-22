#!/usr/bin/env python3
"""Recompute the post hoc background-white proxy used in the H2 validity audit."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

MORPHS={"white","yellow_orange","red_pink","blue_purple"}

def main() -> None:
    ap=argparse.ArgumentParser()
    ap.add_argument("--reserve-csv",type=Path,required=True)
    ap.add_argument("--output-json",type=Path,required=True)
    args=ap.parse_args()
    use=["species","morph","background_effective_pixels","background_palette_count_white"]
    df=pd.read_csv(args.reserve_csv,usecols=use)
    x=df[df["morph"].isin(MORPHS)].copy()
    denom=pd.to_numeric(x["background_effective_pixels"],errors="coerce")
    numer=pd.to_numeric(x["background_palette_count_white"],errors="coerce")
    x=x[denom>0].copy()
    x["background_white_fraction"]=(numer[denom>0]/denom[denom>0]).to_numpy()
    x=x[np.isfinite(x["background_white_fraction"])].copy()
    g=(x.assign(is_white=x["morph"].eq("white"))
         .groupby(["species","is_white"],observed=True)["background_white_fraction"]
         .median().unstack())
    g=g.dropna(subset=[False,True])
    diff=g[True]-g[False]
    w=wilcoxon(g[True],g[False],alternative="two-sided",zero_method="wilcox")
    out={
      "definition":"species-paired median background white fraction; white fraction = background_palette_count_white/background_effective_pixels",
      "species_with_both_states":int(len(g)),
      "median_white_classified":float(g[True].median()),
      "median_nonwhite_classified":float(g[False].median()),
      "fraction_species_white_gt_nonwhite":float((diff>0).mean()),
      "median_within_species_difference":float(diff.median()),
      "wilcoxon_statistic":float(w.statistic),
      "wilcoxon_two_sided_p":float(w.pvalue),
    }
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps(out,indent=2,sort_keys=True))
if __name__=="__main__":
    main()
