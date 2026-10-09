#!/usr/bin/env python3
"""Retrieve source-frozen 1761 local FCP congener TAXONOMY, not colour outcome.

The 709-tip prior tree only covers 20/1761 exact taxa. Build an explicit
original species/observations taxon-ID ledger from the FIXED source 500km
neighborhood, nested 250km subset, and request only iNaturalist taxon IDs
in <=30-ID groups to recover CURRENT genus/family taxonomy. No source photo
image pixels, colour remeasurement or outcome-selected species replacement.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
import urllib.request
from pathlib import Path

import pandas as pd

from audit_fcp_global42111_local_congener_moisture_null_20261009 import (
    fixed_source_cohort, EXPECTED_ORIGINAL_N
)
from audit_fcp_global42111_local_congeners_20261009 import photo_populations

MAX_TAXA_PER_QUERY=30
MAX_TAXA=1761
WAIT_SECONDS=1.2
TIMEOUT=28
API="https://api.inaturalist.org/v1/taxa/"
USER_AGENT="fcp-source-local-congener-tree-taxonomy-20261010/1.0"
NAME_MATCH="EXACT_CANONICAL_BINOMIAL_CASE_AND_SPACE_ONLY"


def norm(s:str)->str:
    return " ".join(str(s).strip().replace("_"," ").split()).casefold()


def source_manifest(all_data:pd.DataFrame,strict:bool=True)->pd.DataFrame:
    cohort,coverage=photo_populations(all_data,strict=strict)
    largest=fixed_source_cohort(cohort["CLIMATE_ALL"],500)
    inner=fixed_source_cohort(cohort["CLIMATE_ALL"],250)
    if strict and (len(largest)!=1761 or len(inner)!=872):
        raise RuntimeError("Frozen published 500/250km original species identities changed")
    if not set(inner.inat_taxon_id).issubset(set(largest.inat_taxon_id)):
        raise RuntimeError("Original 250km source species not nested within 500km cohort")
    required=["inat_taxon_id","species","genus","genus_cell_id","photo_cell_162"]
    if not set(required).issubset(largest):
        raise RuntimeError("Missing exact original taxon photo manifest identifiers")
    original=largest[required].copy()
    original["in_250km_cohort"]=original.inat_taxon_id.isin(set(inner.inat_taxon_id))
    original["original_scientific_name"]=original.species
    original["original_genus"]=original.genus
    original=original.drop(columns=["species","genus"])
    original["inat_taxon_id"]=pd.to_numeric(original.inat_taxon_id,errors="raise").astype(int)
    original=original.sort_values("inat_taxon_id",kind="stable").reset_index(drop=True)
    if original.inat_taxon_id.duplicated().any() or original.in_250km_cohort.sum()!=872:
        raise RuntimeError("Original source identities duplicated/lost")
    return original


def get_taxa(ids:list[int])->dict:
    if not 1<=len(ids)<=MAX_TAXA_PER_QUERY or len(ids)!=len(set(ids)):
        raise ValueError("Requested more than 30 or duplicate original taxon IDs")
    url=API+",".join(map(str,ids))
    request=urllib.request.Request(url,headers={
        "User-Agent":USER_AGENT,"Accept":"application/json"})
    with urllib.request.urlopen(request,timeout=TIMEOUT) as response:
        data=json.load(response)
    if not isinstance(data,dict) or not isinstance(data.get("results"),list):
        raise RuntimeError("Unexpected iNaturalist public taxonomy envelope")
    if len(data["results"])>MAX_TAXA_PER_QUERY:
        raise RuntimeError("Unbounded taxon search returned too many taxa")
    return data


def interpret(record:dict,taxon:dict|None)->dict:
    item={
        "inat_taxon_id":int(record["inat_taxon_id"]),
        "original_scientific_name":str(record["original_scientific_name"]),
        "original_genus":str(record["original_genus"]),
        "source_genus_cell_id":str(record["genus_cell_id"]),
        "source_photo_cell_162":int(record["photo_cell_162"]),
        "in_fixed_250km_original_source":bool(record["in_250km_cohort"]),
        "current_inat_taxon_name":None,
        "current_inat_rank":None,
        "current_genus":None,
        "current_family":None,
        "taxon_ID_stable":False,
        "taxon_name_exact_match":False,
        "taxon_rank_is_species":False,
        "taxon_genus_matches_original":False,
        "family_ancestor_provenance":"UNKNOWN",
        "status":"NOT_RETURNED",
    }
    if taxon is None:return item
    if not isinstance(taxon,dict) or int(taxon.get("id") or 0)!=item["inat_taxon_id"]:
        raise RuntimeError("Unexpected mismatched public taxonomy ID response")
    item["taxon_ID_stable"]=True
    item["current_inat_taxon_name"]=str(taxon.get("name") or "")
    item["current_inat_rank"]=str(taxon.get("rank") or "")
    item["taxon_name_exact_match"]=norm(item["current_inat_taxon_name"])==norm(item["original_scientific_name"])
    item["taxon_rank_is_species"]=item["current_inat_rank"]=="species"
    ancestors=taxon.get("ancestors")
    if not isinstance(ancestors,list):
        item["status"]="TAXON_RETURNED_WITHOUT_RANKED_ANCESTOR_TAXONOMY"
        return item
    families=[str(z.get("name")) for z in ancestors if isinstance(z,dict) and z.get("rank")=="family" and z.get("name")]
    genera=[str(z.get("name")) for z in ancestors if isinstance(z,dict) and z.get("rank")=="genus" and z.get("name")]
    if len(families)!=1 or len(genera)!=1:
        item["status"]="AMBIGUOUS_OR_UNRESOLVED_FAMILY_GENUS_ANCESTORS"
        return item
    item["current_family"]=families[0]
    item["current_genus"]=genera[0]
    item["taxon_genus_matches_original"]=norm(genera[0])==norm(item["original_genus"])
    item["family_ancestor_provenance"]="I_NATURALIST_RANKED_ANCESTORS"
    item["status"]=(
        "EXACT_SPECIES_FAMILY_GENUS_READY"
        if item["taxon_name_exact_match"] and item["taxon_rank_is_species"] and item["taxon_genus_matches_original"]
        else "CURRENT_TAXONOMY_DRIFT_OR_SUBSPECIES_REVIEW"
    )
    return item


def audit(manifest:pd.DataFrame,client=get_taxa,pause=time.sleep)->tuple[dict,pd.DataFrame]:
    if len(manifest)!=MAX_TAXA or manifest.inat_taxon_id.duplicated().any():
        raise ValueError("Frozen original 1761 taxon manifest changed")
    result=[]
    n_errors=0
    queries=0
    for start in range(0,len(manifest),MAX_TAXA_PER_QUERY):
        batch=manifest.iloc[start:start+MAX_TAXA_PER_QUERY]
        ids=batch.inat_taxon_id.astype(int).tolist()
        queries+=1
        try:
            data=client(ids)
            records={}
            for z in data["results"]:
                k=int(z["id"])
                if k not in ids or k in records:
                    raise RuntimeError("Unrequested or duplicate taxon ID in source response")
                records[k]=z
            piece=[interpret(row,records.get(int(row["inat_taxon_id"]))) for row in batch.to_dict("records")]
        except (RuntimeError,ValueError,TimeoutError,ConnectionError,urllib.error.HTTPError,urllib.error.URLError) as e:
            n_errors+=1
            piece=[{**interpret(row,None),"status":"API_ERROR_UNRESOLVED","error_class":type(e).__name__}
                   for row in batch.to_dict("records")]
        result.extend(piece)
        if start+MAX_TAXA_PER_QUERY<len(manifest):
            pause(WAIT_SECONDS)
    d=pd.DataFrame(result)
    if len(d)!=len(manifest) or set(d.inat_taxon_id)!=set(manifest.inat_taxon_id):
        raise RuntimeError("External taxonomy transport lost immutable source species")
    if int(d.in_fixed_250km_original_source.sum())!=872:
        raise RuntimeError("Original nested 250km taxon cohort changed")
    ready=d.status.eq("EXACT_SPECIES_FAMILY_GENUS_READY")
    summary={
        "schema":"fcp_42111_local_source_taxonomy_family_1761_v1",
        "date_jst":"2026-10-10",
        "status":"CURRENT_PUBLIC_TAXONOMY_TIP_PLACEMENT_PREFLIGHT_ONLY",
        "original_global_photo_species_denominator":42111,
        "original_local_500km_species_denominator":1761,
        "original_local_250km_species_denominator":872,
        "n_500km_ready_exact_original_species_family_genus":int(ready.sum()),
        "n_250km_ready_exact_original_species_family_genus":int((ready&d.in_fixed_250km_original_source).sum()),
        "n_500km_taxonomy_unresolved":int((~ready).sum()),
        "n_250km_taxonomy_unresolved":int((~ready&d.in_fixed_250km_original_source).sum()),
        "n_current_taxonomy_api_batches":queries,
        "n_API_error_batches":n_errors,
        "statuses":{str(k):int(v) for k,v in d.status.value_counts().items()},
        "source_species_selection_fixed_before_taxonomy_query":True,
        "source_photo_colour_used_in_family_assignment":False,
        "original_species_not_replaced_with_current_taxon_synonyms":True,
        "phylogenetic_tree_created":False,
        "phylogenetic_environmental_effect_estimated":False,
        "nonclaims":[
            "Current iNaturalist family/genus metadata is a taxonomic placement, not a measured molecular relationship",
            "Synonyms and name mismatches are unresolved, never silently substituted into original source observations",
            "A V.PhyloMaker2 family/genus graft is not evidence for resolved within-genus diversification",
            "Incomplete taxonomy can induce selective phylogenetic coverage and must be explicitly reported",
        ],
    }
    return summary,d


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-breadth-abiotic",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    a=p.parse_args()
    df=pd.read_csv(a.original_breadth_abiotic,low_memory=False)
    source=source_manifest(df)
    a.outdir.mkdir(parents=True,exist_ok=True)
    source.to_csv(a.outdir/"source_original_fixed_taxa_before_API.csv",index=False)
    summary,records=audit(source)
    records.to_csv(a.outdir/"original_local_species_family_genus.csv",index=False)
    (a.outdir/"result.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    print(json.dumps(summary,indent=2,sort_keys=True),flush=True)
if __name__=="__main__":
    main()
