#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
V1 = ROOT / "data/frozen/global_monte_carlo_species_discovery_observation_index_v1.csv.gz"
V2 = ROOT / "data/frozen/global_monte_carlo_species_discovery_v2_observation_index_v1.csv.gz"
OUT = ROOT / "results/rgfca_42111_taxon_cell_conflict_diagnostic_step8e2_20260911"
SEED = 20260911
EXPECTED_PAIRS = 85337


def candidate_key(row: pd.Series) -> str:
    return hashlib.sha256(
        f"{SEED}|{int(row.cell_id)}|{int(row.inat_taxon_id)}|{int(row.observation_id)}|{int(row.photo_id)}".encode("utf-8")
    ).hexdigest()


def pair_order_key(taxon_id: int, cell_id: int) -> str:
    return hashlib.sha256(f"{SEED}|pair|{cell_id}|{taxon_id}".encode("utf-8")).hexdigest()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    parts=[]
    for source,path in (("v1",V1),("v2",V2)):
        x=pd.read_csv(path,usecols=["cell_id","observation_id","photo_id","species","inat_taxon_id"])
        x["discovery_source"]=source
        parts.append(x)
    x=pd.concat(parts,ignore_index=True).dropna().copy()
    for c in ["cell_id","observation_id","photo_id","inat_taxon_id"]:
        x[c]=pd.to_numeric(x[c],errors="raise").astype(int)
    x=x.drop_duplicates(["cell_id","inat_taxon_id","observation_id","photo_id"],keep="first")
    if x[["inat_taxon_id","cell_id"]].drop_duplicates().shape[0] != EXPECTED_PAIRS:
        raise RuntimeError("pair count drift")

    # Independent minimum-hash selection (the failed first attempt).
    x["candidate_hash"]=[candidate_key(r) for r in x.itertuples(index=False)]
    indep=(x.sort_values(["inat_taxon_id","cell_id","candidate_hash"],kind="mergesort")
             .drop_duplicates(["inat_taxon_id","cell_id"],keep="first").copy())
    dup_obs_mask=indep["observation_id"].duplicated(keep=False)
    dup_photo_mask=indep["photo_id"].duplicated(keep=False)
    duplicated_obs_ids=int(indep.loc[dup_obs_mask,"observation_id"].nunique())
    duplicated_photo_ids=int(indep.loc[dup_photo_mask,"photo_id"].nunique())
    affected_pairs=int((dup_obs_mask|dup_photo_mask).sum())

    # Candidate-level cross-pair conflicts.
    photo_pair_n=(x.groupby("photo_id",observed=True)[["inat_taxon_id","cell_id"]]
                    .apply(lambda g: len(g.drop_duplicates())).astype(int))
    obs_pair_n=(x.groupby("observation_id",observed=True)[["inat_taxon_id","cell_id"]]
                  .apply(lambda g: len(g.drop_duplicates())).astype(int))

    # Deterministic constrained greedy matching: pair order and candidate order are
    # both pre-colour hash functions. This is diagnostic only; no frame is frozen here.
    pair_groups={}
    for (taxon_id,cell_id),g in x.groupby(["inat_taxon_id","cell_id"],sort=False,observed=True):
        q=g.sort_values(["candidate_hash","observation_id","photo_id"],kind="mergesort").copy()
        pair_groups[(int(taxon_id),int(cell_id))]=q
    pair_order=sorted(pair_groups,key=lambda k: pair_order_key(k[0],k[1]))
    used_obs=set(); used_photo=set(); matched=[]; unresolved=[]
    for taxon_id,cell_id in pair_order:
        chosen=None
        for r in pair_groups[(taxon_id,cell_id)].itertuples(index=False):
            oid=int(r.observation_id); pid=int(r.photo_id)
            if oid in used_obs or pid in used_photo:
                continue
            chosen=r; break
        if chosen is None:
            unresolved.append((taxon_id,cell_id,len(pair_groups[(taxon_id,cell_id)])))
            continue
        used_obs.add(int(chosen.observation_id)); used_photo.add(int(chosen.photo_id)); matched.append((taxon_id,cell_id))

    unresolved_df=pd.DataFrame(unresolved,columns=["inat_taxon_id","cell_id","candidate_rows"])
    unresolved_df.to_csv(OUT/"greedy_unresolved_pairs.csv",index=False,lineterminator="\n")
    result={
        "analysis":"rgfca_42111_taxon_cell_anchor_conflict_diagnostic_step8e2",
        "status":"complete_metadata_only_conflict_diagnostic",
        "taxon_cell_pairs":EXPECTED_PAIRS,
        "independent_hash_selected_rows":int(len(indep)),
        "duplicated_selected_observation_ids":duplicated_obs_ids,
        "duplicated_selected_photo_ids":duplicated_photo_ids,
        "selected_pairs_affected_by_any_duplicate_id":affected_pairs,
        "candidate_observation_ids_spanning_multiple_pairs":int((obs_pair_n>1).sum()),
        "candidate_photo_ids_spanning_multiple_pairs":int((photo_pair_n>1).sum()),
        "max_pairs_per_observation_id":int(obs_pair_n.max()),
        "max_pairs_per_photo_id":int(photo_pair_n.max()),
        "deterministic_greedy_unique_pairs_matched":int(len(matched)),
        "deterministic_greedy_unresolved_pairs":int(len(unresolved)),
        "image_pixels_opened":False,
        "flower_colour_used":False,
        "boundary":"Diagnostic only. No taxon-cell measurement frame is authorized or frozen by this output."
    }
    (OUT/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    (OUT/"RESULT.md").write_text(
        "# RGFCA Step 8E2 — taxon×cell anchor conflict diagnostic\n\n"
        f"- taxon×cell pairs: **{EXPECTED_PAIRS:,}**\n"
        f"- duplicated selected observation IDs: **{duplicated_obs_ids:,}**\n"
        f"- duplicated selected photo IDs: **{duplicated_photo_ids:,}**\n"
        f"- selected pair rows affected: **{affected_pairs:,}**\n"
        f"- source observation IDs spanning >1 pair: **{int((obs_pair_n>1).sum()):,}**\n"
        f"- source photo IDs spanning >1 pair: **{int((photo_pair_n>1).sum()):,}**\n"
        f"- deterministic unique greedy matched pairs: **{len(matched):,}**\n"
        f"- deterministic greedy unresolved pairs: **{len(unresolved):,}**\n"
        "- image pixels opened: **false**\n",
        encoding="utf-8"
    )
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
