#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, re
from pathlib import Path
import pandas as pd

KEYWORDS=[
    r"anthocyan", r"flavon", r"treatment", r"spring", r"summer",
    r"mild", r"hot", r"temperature", r"individual", r"plant", r"colour", r"color"
]
RX=re.compile("|".join(KEYWORDS),re.I)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-data",required=True)
    p.add_argument("--outdir",required=True)
    a=p.parse_args()
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    book=pd.read_excel(a.source_data,sheet_name=None,header=None,engine="openpyxl")
    if "Figure 3F-J" not in book:
        raise SystemExit(f"Figure 3F-J sheet not found; sheets={list(book)}")
    d=book["Figure 3F-J"].copy()
    d.to_csv(out/"figure3FJ_raw.csv",index=False,header=False)
    hits=[]
    for i,row in d.iterrows():
        vals=["" if pd.isna(x) else str(x) for x in row.tolist()]
        text=" | ".join(vals)
        if RX.search(text):
            hits.append({"excel_row":int(i+1),"row_text":text[:12000]})
    pd.DataFrame(hits).to_csv(out/"figure3FJ_keyword_rows.csv",index=False)
    nonempty=[]
    for j in range(d.shape[1]):
        vals=[str(x) for x in d.iloc[:,j].dropna().tolist()]
        nonempty.append({
            "column_index":j,
            "n_nonmissing":len(vals),
            "first_values":" | ".join(vals[:12])[:4000]
        })
    pd.DataFrame(nonempty).to_csv(out/"column_preview.csv",index=False)
    result={
        "schema":"fcp_moricandia_temperature_dose_preflight_v1",
        "status":"complete",
        "sheet":"Figure 3F-J",
        "rows":int(d.shape[0]),
        "cols":int(d.shape[1]),
        "keyword_hit_rows":len(hits),
        "role":"schema_only_preflight_before coding temperature-dose model"
    }
    (out/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
    print(pd.DataFrame(hits).head(80).to_string(index=False))

if __name__=="__main__":
    main()
