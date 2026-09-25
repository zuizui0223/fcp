#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, re
from pathlib import Path
import pandas as pd

WHITE = re.compile(r"\b(white|whitish|albino|acyanic|unpigmented|non[- ]pigmented)\b", re.I)
GROUPS = {
    "yellow_orange": re.compile(r"\b(yellow|orange|gold(?:en)?)\b", re.I),
    "red_pink": re.compile(r"\b(red|pink|magenta|scarlet|crimson)\b", re.I),
    "blue_purple": re.compile(r"\b(blue|purple|violet|lavender|lilac)\b", re.I),
    "dark_other": re.compile(r"\b(black|brown|bronze|maroon)\b", re.I),
    "green": re.compile(r"\bgreen\b", re.I),
}
HC_REASONS={
    "known_positive_control",
    "natural_polymorphism_with_direct_and_replicated_support",
    "aggregated_followup_direct_and_natural_support",
}

def classify(text: str) -> tuple[str,list[str]]:
    white=bool(WHITE.search(text))
    groups=[k for k,rx in GROUPS.items() if rx.search(text)]
    if white and groups:
        return "white_plus_nonwhite",groups
    if (not white) and len(groups)>=2:
        return "nonwhite_multicolour",groups
    return "unresolved",groups

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--resolved-queue",default="data/resolved_inputs/global_flower_colour_review_queue_resolved.csv")
    p.add_argument("--outdir",required=True)
    p.add_argument("--min-per-group",type=int,default=20)
    a=p.parse_args()
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(a.resolved_queue,low_memory=False)
    d=d.loc[d.evidence_class.astype(str).eq("natural_polymorphism")].copy()
    d=d.loc[d.canonical_name.astype(str).str.match(r"^[A-Z][A-Za-z-]+ [a-z][A-Za-z-]+$")].copy()
    # The state classification is based only on the literature source text already used
    # to establish natural flower-colour polymorphism; no molecular audit outputs are read.
    rows=[]
    for r in d.itertuples(index=False):
        text=" ".join([
            str(getattr(r,"best_title","") or ""),
            str(getattr(r,"best_match_evidence","") or ""),
        ])
        state,groups=classify(text)
        reason=str(getattr(r,"review_reason","") or "")
        rows.append({
            "canonical_name":str(r.canonical_name),
            "family":str(r.family),
            "review_reason":reason,
            "high_confidence":reason in HC_REASONS,
            "colour_comparison_state":state,
            "white_documented":bool(WHITE.search(text)),
            "nonwhite_groups":";".join(groups),
            "n_nonwhite_groups":len(groups),
            "source_title":str(getattr(r,"best_title","") or ""),
            "source_doi":str(getattr(r,"best_doi","") or ""),
        })
    x=pd.DataFrame(rows).drop_duplicates("canonical_name").sort_values("canonical_name")
    x.to_csv(out/"white_vs_hue_preflight_species.csv",index=False)
    hc=x.loc[x.high_confidence].copy()
    counts=hc.colour_comparison_state.value_counts().to_dict()
    w=int(counts.get("white_plus_nonwhite",0)); h=int(counts.get("nonwhite_multicolour",0))
    status="MOLECULAR_CONTRAST_COVERAGE_PASS" if w>=a.min_per_group and h>=a.min_per_group else "MOLECULAR_CONTRAST_COVERAGE_FAIL"
    result={
      "schema":"fcp_white_vs_hue_molecular_preflight_v1",
      "status":status,
      "source":"resolved natural-FCP literature queue only; molecular audit results not read",
      "natural_fcp_species":int(len(x)),
      "high_confidence_species":int(len(hc)),
      "high_confidence_state_counts":{k:int(v) for k,v in counts.items()},
      "minimum_per_comparison_group":a.min_per_group,
      "white_plus_nonwhite_species":w,
      "nonwhite_multicolour_species":h,
      "unresolved_species":int(counts.get("unresolved",0)),
      "allowed_next_step":"freeze symmetric molecular-evidence search only if coverage passes",
      "hard_nonclaims":[
        "absence of white terminology is not evidence that a species lacks white flowers",
        "unresolved species are never assigned to the nonwhite control group",
        "this preflight contains no molecular-mechanism outcome"
      ]
    }
    (out/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))

if __name__=="__main__":
    main()
