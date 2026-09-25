#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, re
from pathlib import Path
import numpy as np
import pandas as pd

RX=re.compile(r"(wavelength|reflect|spring|summer|lilac|white|uv|vis|nm)",re.I)

def numeric_profile(s):
    x=pd.to_numeric(s,errors="coerce").dropna()
    if len(x)<20:
        return None
    arr=x.to_numpy(float)
    return {
        "n_numeric":int(len(arr)),
        "min":float(np.nanmin(arr)),
        "max":float(np.nanmax(arr)),
        "monotonic_increasing":bool(np.all(np.diff(arr)>=0)),
        "unique":int(len(np.unique(arr))),
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-data",required=True)
    p.add_argument("--outdir",required=True)
    a=p.parse_args()
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    book=pd.read_excel(a.source_data,sheet_name=None,header=None,engine="openpyxl")
    sheets=[]; hits=[]; profiles=[]
    for name,df in book.items():
        sheets.append({"sheet":str(name),"rows":int(df.shape[0]),"cols":int(df.shape[1])})
        for i,row in df.iterrows():
            vals=["" if pd.isna(x) else str(x) for x in row.tolist()]
            txt=" | ".join(vals)
            if RX.search(txt):
                hits.append({"sheet":str(name),"excel_row":int(i+1),"row_text":txt[:12000]})
        for j in range(df.shape[1]):
            prof=numeric_profile(df.iloc[:,j])
            if prof:
                vals=["" if pd.isna(x) else str(x) for x in df.iloc[:8,j].tolist()]
                profiles.append({"sheet":str(name),"column_index":j,"first_values":" | ".join(vals),**prof})
    pd.DataFrame(sheets).to_csv(out/"sheet_manifest.csv",index=False)
    pd.DataFrame(hits).to_csv(out/"keyword_rows.csv",index=False)
    pd.DataFrame(profiles).to_csv(out/"numeric_column_profiles.csv",index=False)
    result={
      "schema":"fcp_moricandia_spectral_preflight_v1",
      "status":"complete",
      "sheets":sheets,
      "keyword_rows":len(hits),
      "numeric_profiles":len(profiles),
      "rule":"schema-only preflight; no spectral contrast statistic is calculated"
    }
    (out/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
    print(pd.DataFrame(hits).head(120).to_string(index=False))
    print(pd.DataFrame(profiles).head(120).to_string(index=False))
if __name__=="__main__":
    main()
