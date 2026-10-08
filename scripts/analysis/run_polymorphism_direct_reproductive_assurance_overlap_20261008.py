#!/usr/bin/env python3
"""Source-frozen feasibility of DIRECT reproductive assurance evidence across FCP taxa.

Only screens geographic photo-species coverage against published, independently
measured plant mating/autonomous-seed datasets. Does not pool self-compatibility,
outcrossing and autonomous selfing as if one phenotype or causal fitness payoff.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

PHOTO_SHA256 = {
    "discovery": "ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
    "validation": "0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
    "third": "57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186",
}
EXTERNAL_PIN = {
    "source_repo": "zuizui0223/island",
    "commit_sha": "92f007797fbbcb3c6743a4ffb2628fa609329499",
    "rodger_gz_sha256": "fa745c578f3537933fafedc1d36b4ea266348cd83d7f6cbb231c253b0f348d3f",
    "rodger_original_rows": 1528,
    "sources": {
        "rodger": {"doi": "10.1126/sciadv.abd3524", "trait": "autonomous seed/fruit production without pollinators"},
        "goodwillie": {"doi": "10.1146/annurev.ecolsys.36.091704.175539", "trait": "multilocus natural-population genetic outcrossing rate"},
        "razanajatovo": {"doi": "10.1038/ncomms13313", "trait": "experimentally measured autofertility index or self-compatibility index"},
    }
}
PHOTO_STATES = ("white", "yellow_orange", "red_pink", "blue_purple")
MIN_CLASSIFIABLE = 40
MIN_PHOTO_WHITE = 5
MIN_PHOTO_CHROMATIC = 5
MIN_DIRECT_TRAIT_PER_COHORT_FOR_ASSOC = 20
MIN_DIRECT_TRAIT_TOTAL_FOR_ASSOC = 80
MIN_WHITE_COLOURED_PER_COHORT_FOR_ASSOC = 8

SCIENTIFIC_BINOMIAL = re.compile(r"^[A-Z][a-zA-Z-]+\s+[a-z][a-zA-Z-]+$")


def sha(path: Path) -> str:
    digest=hashlib.sha256()
    with path.open("rb") as src:
        for chunk in iter(lambda:src.read(1<<20),b""):
            digest.update(chunk)
    return digest.hexdigest()


def binomial(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    s=re.sub(r"\s+"," ",str(value).replace("_"," ").strip())
    bits=s.split(" ")
    if len(bits)<2: return ""
    first=" ".join(bits[:2])
    return first if SCIENTIFIC_BINOMIAL.fullmatch(first) else ""


def truth(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.lower().str.strip().isin({"true","1","yes","y"})


def read_photos(path: Path,cohort: str) -> pd.DataFrame:
    if sha(path)!=PHOTO_SHA256[cohort]:
        raise ValueError(f"{cohort}: measured photo byte hash has drifted")
    d=pd.read_csv(path,low_memory=False)
    must={"inat_taxon_id","species","photo_id","morph","global_classifiable",
          "latitude","longitude"}
    if not must.issubset(d.columns):
        raise RuntimeError(f"Missing measured photo columns: {sorted(must-set(d.columns))}")
    d=d.loc[truth(d.global_classifiable) & d.morph.isin(PHOTO_STATES), sorted(must)].copy()
    d["latitude"]=pd.to_numeric(d.latitude,errors="coerce")
    d["longitude"]=pd.to_numeric(d.longitude,errors="coerce")
    d=d.loc[d.latitude.between(-90,90) & d.longitude.between(-180,180)].copy()
    counts=d.groupby("inat_taxon_id").size()
    d=d.loc[d.inat_taxon_id.isin(counts[counts>=MIN_CLASSIFIABLE].index)].copy()
    d["canonical_binomial"]=d.species.map(binomial)
    if (d["canonical_binomial"]=="").any():
        raise RuntimeError("Unknown photo scientific binomial")
    if d.photo_id.duplicated().any():
        raise RuntimeError("Repeated photo IDs in eligible source cohort")
    rows=[]
    for taxon,g in d.groupby("inat_taxon_id",sort=True):
        count=g.morph.value_counts()
        white=int(count.get("white",0))
        nc=int(len(g)-white)
        rows.append({
            "cohort":cohort,"inat_taxon_id":int(taxon),
            "species":str(g.species.iloc[0]),
            "canonical_binomial":str(g.canonical_binomial.iloc[0]),
            "genus":str(g.canonical_binomial.iloc[0]).split(" ")[0],
            "n_classifiable":int(len(g)),
            "photo_white":white,"photo_chromatic":nc,
            "white_chromatic_photo_eligible":bool(white>=MIN_PHOTO_WHITE and nc>=MIN_PHOTO_CHROMATIC),
            "two_chromatic_hues_at_least5":bool(sum(count.get(x,0)>=5 for x in PHOTO_STATES[1:])>=2)
        })
    species=pd.DataFrame(rows)
    if species.canonical_binomial.duplicated().any():
        raise RuntimeError("Several taxon IDs have the same source binomial; do not deduplicate silently")
    return species


def file_summary(path:Path) -> tuple[pd.DataFrame,dict]:
    d=pd.read_csv(path,low_memory=False)
    info={"n_rows":int(len(d)), "columns":[str(c) for c in d.columns],
          "header_sha256":hashlib.sha256(",".join(map(str,d.columns)).encode()).hexdigest(),
          "source_file_sha256":sha(path)}
    return d,info


def get_source_names(df:pd.DataFrame,source:str) -> tuple[pd.Series,str]:
    if source=="rodger":
        # Source authorial mapping is fixed by island checkpoint:
        # src/island_v2/rodger_2021_autofertility_checkpoint.py
        #   selection["submitted_species"] = selection["genus.species"]...
        # Do not select raw 'taxon', fuzzy synonyms or max-overlap columns.
        if "genus.species" in df.columns:
            return df["genus.species"].map(binomial),"genus.species"
        return pd.Series([""]*len(df),index=df.index,dtype="string"),"HOLD_MISSING_SOURCE_DEFINED_GENUS_SPECIES"
    if source=="goodwillie":
        for name in ("Genus species","Genus Species","genus species"):
            if name in df.columns:
                return df[name].map(binomial),name
        raise ValueError("Goodwillie source binomial column absent")
    if source=="razanajatovo":
        if "Species" not in df.columns:
            raise ValueError("Razanajatovo original Species column absent")
        return df["Species"].map(binomial),"Species"
    raise ValueError(f"Unrecognized source {source}")


def direct_trait_info(df:pd.DataFrame,source:str) -> dict:
    """Coverage of real measured variables, NEVER infer from a text proxy."""
    if source=="goodwillie":
        needed=["Mean-tm","PopType (0=natural, 1= experimental, 2=seed orchard, 3=agricultural)"]
        if any(v not in df for v in needed):
            raise ValueError("Missing natural-population outcrossing variables")
        natural=pd.to_numeric(df[needed[1]],errors="coerce")==0
        t=pd.to_numeric(df["Mean-tm"],errors="coerce")
        measured=natural & t.between(0,1)
        return {"measured_mask":measured,"measurement_column":"Mean-tm","measurement_type":"genetic_outcrossing_fraction_0_to_1",
                "n_natural_measured_rows":int(measured.sum())}
    if source=="razanajatovo":
        auto=[c for c in df.columns if str(c).startswith("Autofertility_index_")]
        sc=[c for c in df.columns if str(c).startswith("Self-compatibility_index_")]
        if not auto or not sc:
            raise ValueError("Missing direct autofertility vs self-compatibility columns")
        au=df[auto].apply(pd.to_numeric,errors="coerce")
        scv=df[sc].apply(pd.to_numeric,errors="coerce")
        measured=(au.notna() & (au>=0)&(au<=1)).any(axis=1)
        return {"measured_mask":measured,
                "measurement_column":";".join(auto),"measurement_type":"autofertility_index_not_self_compatibility",
                "n_direct_autofertility_rows":int(measured.sum()),
                "n_self_compatibility_only_rows":int((scv.notna().any(axis=1)&~measured).sum())}
    if source=="rodger":
        auto=[c for c in df.columns if str(c).lower().startswith("auto.")]
        if not auto:
            return {"measured_mask":pd.Series(False,index=df.index),
                    "measurement_column":"HOLD_AUTO_COLUMN_MAPPING",
                    "measurement_type":"autonomous_fruit_seed_after_pollinator_exclusion",
                    "status":"HOLD_SOURCE_SCHEMA_NO_AUTO_COLUMNS"}
        au=df[auto].apply(pd.to_numeric,errors="coerce")
        # Island's original checkpoint rejects an entire source row if ANY
        # declared nonmissing exclusion measurement is negative or nonnumeric.
        raw=df[auto].fillna("").astype(str).apply(lambda col:col.str.strip())
        valid_entry=raw.apply(lambda col:~col.str.casefold().isin(("", "na", "nan", "none")))
        valid_numeric=(~valid_entry | (au.notna()&(au>=0))).all(axis=1)
        measured=(au.notna()&(au>=0)).any(axis=1) & valid_numeric
        return {"measured_mask":measured,"measurement_column":";".join(auto),
                "measurement_type":"autonomous_fruit_seed_after_pollinator_exclusion",
                "n_measured_positive_or_zero_rows":int(measured.sum())}
    raise ValueError(source)


def summarize_overlap(photos:pd.DataFrame,trait_dfs:dict[str,pd.DataFrame]) -> tuple[dict,pd.DataFrame]:
    for source in trait_dfs:
        if source not in ("rodger","goodwillie","razanajatovo"):
            raise ValueError("Unexpected external source")
    records=[]
    sources={}
    for source,df in trait_dfs.items():
        taxa,col=get_source_names(df,source)
        info=direct_trait_info(df,source)
        measured=info.pop("measured_mask")
        measured_taxa=set(taxa.loc[measured])-{""}
        named_taxa=set(taxa)-{""}
        joined=photos.loc[photos.canonical_binomial.isin(measured_taxa)].copy()
        total=photos.loc[photos.canonical_binomial.isin(named_taxa)].copy()
        counts=joined.groupby("cohort").size()
        white=joined.loc[joined.white_chromatic_photo_eligible].groupby("cohort").size()
        stats={
            "taxonomic_name_column":col,
            "n_source_rows":len(df),
            "n_exact_species_names":len(named_taxa),
            "n_species_with_direct_measurement":len(measured_taxa),
            "n_fcp_species_any_named_source_overlap":len(total),
            "n_fcp_species_measured_trait_overlap":len(joined),
            "n_fcp_white_chromatic_photo_species_measured_trait_overlap":int(joined.white_chromatic_photo_eligible.sum()),
            "n_measured_overlap_by_cohort":{cohort:int(counts.get(cohort,0)) for cohort in PHOTO_SHA256},
            "n_measured_white_chromatic_overlap_by_cohort":{cohort:int(white.get(cohort,0)) for cohort in PHOTO_SHA256},
            "measurement":info,
            "photo_colour_trait_association_eligible":bool(
                len(joined)>=MIN_DIRECT_TRAIT_TOTAL_FOR_ASSOC and
                all(counts.get(cohort,0)>=MIN_DIRECT_TRAIT_PER_COHORT_FOR_ASSOC
                    and white.get(cohort,0)>=MIN_WHITE_COLOURED_PER_COHORT_FOR_ASSOC
                    for cohort in PHOTO_SHA256)
            ),
        }
        sources[source]=stats
        for _,row in joined.iterrows():
            records.append({
                "source":source,"original_doi":EXTERNAL_PIN["sources"][source]["doi"],
                "cohort":row.cohort,"species":row.canonical_binomial,
                "n_classifiable":row.n_classifiable,
                "n_photo_white":row.photo_white,
                "n_photo_chromatic":row.photo_chromatic,
                "white_chromatic_photo_eligible":bool(row.white_chromatic_photo_eligible),
                "trait_type":info["measurement_type"],
                "exact_species_label_matched":True,
                "source_plant_genotype_linked_to_photo":False,
                "source_anthocyanin_or_fitness_linked_to_photo":False,
            })
    ledger=pd.DataFrame(records)
    result={
        "schema":"fcp_direct_reproductive_assurance_trait_overlap_v1",
        "date_jst":"2026-10-08",
        "status":"complete_source_overlap_feasibility",
        "role":"same_species_independent_direct_mating_trait_coverage_not_causal_selection",
        "island_source_pin":EXTERNAL_PIN,
        "photo_source_pin":PHOTO_SHA256,
        "source_photo_species_total":int(len(photos)),
        "source_photo_species_by_cohort":{c:int((photos.cohort==c).sum()) for c in PHOTO_SHA256},
        "source_trait_panels":sources,
        "minimums":{
            "source_photo_species":MIN_CLASSIFIABLE,
            "white_and_nonwhite_photo_counts_each":5,
            "direct_trait_total":MIN_DIRECT_TRAIT_TOTAL_FOR_ASSOC,
            "direct_trait_per_cohort":MIN_DIRECT_TRAIT_PER_COHORT_FOR_ASSOC,
            "white_chromatic_per_cohort":MIN_WHITE_COLOURED_PER_COHORT_FOR_ASSOC,
        },
        "all_sources_direct_trait_overlap_status":"BIOLOGICAL_ASSOCIATION_FEASIBLE" if
            any(z["photo_colour_trait_association_eligible"] for z in sources.values())
            else "HOLD_NO_SUFFICIENT_INDEPENDENT_SPECIES_COHORT_COVERAGE",
        "hard_nonclaims":[
            "reported self-compatibility is not autonomous selfing; no phenotype proxy replaces measured outcome",
            "population mean outcrossing rate is not autonomous bagging reproduction and is a species-level aggregate",
            "these direct breeding traits are not white-vs-pigmented morph-specific fitness payoffs",
            "case-assembled literature datasets have biased taxonomic/ecological coverage",
            "exact binomial matching without validated synonyms gives conservative lower-bound name overlap",
            "high-depth iNaturalist photographs are not unbiased world angiosperm prevalence of genetic colour polymorphism",
            "photo-derived white states are exposed to clipping and are not verified anthocyanin-loss genotypes",
            "a reproductive trait association would not prove net selective trade-off or balancing selection",
            "Rodger source-defined genus.species field is used; ambiguous raw taxon or unverified synonyms never replace it",
            "no paper manuscript H1/H2 or prior photo geography decision is overwritten",
        ],
        "confirmatory_decisions_changed":False,
    }
    return result,ledger


def main() -> None:
    p=argparse.ArgumentParser()
    for cohort in PHOTO_SHA256:
        p.add_argument("--"+cohort,required=True,type=Path)
    for trait in ("goodwillie","razanajatovo","rodger"):
        p.add_argument("--"+trait,required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    args=p.parse_args()
    if sha(args.rodger)!=EXTERNAL_PIN["rodger_gz_sha256"]:
        raise RuntimeError("Rodger original-source projected CSV gzip SHA256 mismatch")
    photos=pd.concat([read_photos(getattr(args,c),c) for c in PHOTO_SHA256],
                     ignore_index=True)
    if photos.canonical_binomial.duplicated().any():
        raise RuntimeError("One species appears more than once across supposedly disjoint photo cohorts")
    dfs={}
    inventories={}
    for source in ("goodwillie","razanajatovo","rodger"):
        df,desc=file_summary(getattr(args,source))
        dfs[source]=df
        inventories[source]=desc
    res,ledger=summarize_overlap(photos,dfs)
    res["original_external_file_inventory"]=inventories
    args.outdir.mkdir(parents=True,exist_ok=True)
    (args.outdir/"result.json").write_text(
        json.dumps(res,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    ledger.to_csv(args.outdir/"matched_species_direct_trait_opportunity.csv",index=False)
    print(json.dumps(res,indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()
