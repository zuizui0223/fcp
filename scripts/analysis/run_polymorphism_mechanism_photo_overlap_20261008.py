#!/usr/bin/env python3
"""Audit whether source-verified ecological FCP study systems overlap photo cohorts.

The research question is feasibility of a causal mechanism bridge, NOT whether
the eleven purposively chosen literature systems estimate mechanism prevalence.
No trait causation or biological flower-pigment class is inferred from photos.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

COHORT_SHA = {
    "discovery": "ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
    "validation": "0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
    "third": "57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186",
}
MORPHS = ("white", "yellow_orange", "red_pink", "blue_purple")
MIN_CLASSIFIABLE = 40
MIN_PER_MORPH = 5
MIN_LOCAL_PAIRS = 30
EARTH_KM = 6371.0088
MIN_INDEPENDENT_MECHANISM_PHOTO_SYSTEMS = 6
MIN_TRUE_MORPH_CHEMISTRY_DIRECT_LINKS = 3


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for data in iter(lambda:f.read(1<<20),b""):
            h.update(data)
    return h.hexdigest()


def safe_bool(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.casefold().isin(
        ("true", "yes", "1", "y")
    )


def canonical_binomial(s: str) -> str:
    """Do not fuzzy-match scientific names or infer synonyms."""
    terms = str(s).strip().split()
    if len(terms) < 2 or not terms[0] or not terms[1]:
        return ""
    return " ".join(terms[:2])


def registry_check(path: Path) -> dict:
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("schema") != "fcp_morph_maintenance_primary_study_registry_v1":
        raise ValueError("Wrong biological source registry")
    studies = data.get("systems", [])
    if len(studies)<8:
        raise ValueError("Mechanism panel must contain >=8 independent organism systems")
    names=set()
    allowed_prefix="10."
    for s in studies:
        name=s["species"]
        if not name or name!=canonical_binomial(name):
            raise ValueError(f"Not exactly a canonical binomial: {name}")
        if name in names:
            raise ValueError(f"Duplicate biological system: {name}")
        names.add(name)
        if not isinstance(s.get("demonstrated"),list) or not s["demonstrated"]:
            raise ValueError("No primary organism evidence listed")
        if not isinstance(s.get("unresolved"),list) or not s["unresolved"]:
            raise ValueError("Missing scientific inferential limits")
        if not s.get("evidence_role") or not s.get("mechanism"):
            raise ValueError("Unspecified evidence type")
        doi=s.get("source_doi", [])
        if not doi or any(not str(v).startswith(allowed_prefix) for v in doi):
            raise ValueError(f"DOI must identify source: {name}")
    return data


def read_photos(path: Path, cohort: str, *, verify: bool=True) -> pd.DataFrame:
    if verify and sha256(path)!=COHORT_SHA[cohort]:
        raise RuntimeError(f"Frozen {cohort} photos do not match archived bytes")
    d=pd.read_csv(path,low_memory=False)
    required={"species","inat_taxon_id","photo_id","morph","global_classifiable",
              "latitude","longitude","observed_on","observer_id"}
    if not required.issubset(d.columns):
        raise RuntimeError(f"Missing measured source fields {sorted(required-set(d.columns))}")
    d=d.loc[:,sorted(required)].copy()
    d["species_binomial"]=d.species.astype(str).map(canonical_binomial)
    d["global_classifiable"]=safe_bool(d.global_classifiable)
    for v in ("latitude","longitude"):
        d[v]=pd.to_numeric(d[v],errors="coerce")
    d["coordinate_valid"]=d.latitude.between(-90,90)&d.longitude.between(-180,180)
    d["classifiable"]=d.global_classifiable & d.coordinate_valid & d.morph.isin(MORPHS)
    return d


def haversine_matrix(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    la=np.deg2rad(lat.astype(float))
    lo=np.deg2rad(lon.astype(float))
    c=np.cos(la)
    xyz=np.stack([c*np.cos(lo),c*np.sin(lo),np.sin(la)],axis=1)
    return np.arccos(np.clip(xyz@xyz.T,-1.0,1.0))*EARTH_KM


def one_species_photo_summary(g: pd.DataFrame) -> dict:
    """One record for one exact literature species in one photograph cohort."""
    n_source=int(len(g))
    u=g.loc[g.classifiable].copy()
    n=len(u)
    counts=u.morph.value_counts().reindex(MORPHS,fill_value=0)
    white=int(counts["white"])
    nonwhite=int(counts[1:].sum())
    hue_arms=int((counts[1:]>=MIN_PER_MORPH).sum())
    eligible=n>=MIN_CLASSIFIABLE
    photos_regardless_of_fcp=bool(n>0)
    local_pairs=0
    local_white_chromatic=0
    if eligible:
        dist=haversine_matrix(u.latitude.to_numpy(float),
                              u.longitude.to_numpy(float))
        i,j=np.triu_indices(n,k=1)
        near=dist[i,j]<=50.0
        local_pairs=int(near.sum())
        labels=(u.morph.astype(str).to_numpy()=="white")
        local_white_chromatic=int(np.count_nonzero(
            (labels[i]!=labels[j])&near))
    return{
        "observed_in_photo_candidate_cohort":bool(n_source>0),
        "source_photo_rows":n_source,
        "n_classifiable_georeferenced":n,
        "n_photo_white":white,
        "n_photo_nonwhite":nonwhite,
        "n_hue_categories_with_at_least_5_photos":hue_arms,
        "n_distinct_visible_states_at_least_5":int((counts>=MIN_PER_MORPH).sum()),
        "meets_40_photo_eligibility":bool(eligible),
        "meets_5_white_and_5_nonwhite":bool(eligible and white>=5 and nonwhite>=5),
        "meets_5_white_and_one_specific_nonwhite_hue":bool(
            eligible and white>=5 and hue_arms>=1),
        "meets_two_nonwhite_hues_each_5":bool(eligible and hue_arms>=2),
        "photo_neighborhood_50km_pairs":local_pairs,
        "photo_50km_white_nonwhite_pairs":local_white_chromatic,
        "photo_white_state_chemically_confirmed":False,
        "in_pop_fitness_genotype_mapped_to_source_images":False,
    }


def join_panel(registry: dict, dframes: dict[str,pd.DataFrame]) -> tuple[dict,pd.DataFrame]:
    for cohort,d in dframes.items():
        if cohort not in COHORT_SHA:
            raise ValueError("unexpected source cohort")
        if "species_binomial" not in d.columns:
            raise ValueError("source photo input not normalized")
    dct={}
    for cohort,d in dframes.items():
        dct[cohort]={tax:g for tax,g in d.groupby("species_binomial",sort=True)}
    rows=[]
    all_names=set()
    for record in registry["systems"]:
        name=record["species"]
        present=[]
        for cohort in COHORT_SHA:
            group=dct.get(cohort,{}).get(name)
            if group is None:
                continue
            stats=one_species_photo_summary(group)
            if not stats["observed_in_photo_candidate_cohort"]:
                continue
            present.append(cohort)
            rows.append({
                "species":name,"family":record["family"],
                "literature_mechanism":record["mechanism"],
                "literature_evidence_role":record["evidence_role"],
                "primary_study_dois":";".join(record["source_doi"]),
                "cohort":cohort,
                **stats,
            })
        if len(present)>1:
            raise RuntimeError(f"species-disjoint sampling broken: {name} {present}")
        all_names.add(name)
    details=pd.DataFrame(rows)
    matches=sorted(set(details.species)) if len(details) else []
    eligible=details.loc[details.meets_40_photo_eligibility] if len(details) else details
    wmatched=eligible.loc[eligible.meets_5_white_and_5_nonwhite] if len(eligible) else eligible
    photo_true_morph_validated=0 # Not available in measured tables; cannot set using source literature.
    out={
        "schema":"fcp_morph_maintenance_primary_study_photo_overlap_v1",
        "date_jst":"2026-10-08",
        "status":"complete_photo_overlap_diagnostic",
        "role":"prevalence_inference_prohibited_source_mechanism_bridge_feasibility",
        "n_independent_literature_species":len(all_names),
        "n_verified_primary_study_dois":len(set(
            doi for s in registry["systems"] for doi in s["source_doi"])),
        "n_literature_species_any_measured_photo_cohort":len(matches),
        "n_literature_species_ge40_classifiable":len(eligible),
        "n_literature_species_ge5_white_and_ge5_nonwhite":len(wmatched),
        "n_true_morph_genotype_chemistry_photo_links":photo_true_morph_validated,
        "source_species_any_sample":matches,
        "source_species_ge40":sorted(set(eligible.species)) if len(eligible) else [],
        "source_species_white_nonwhite":sorted(set(wmatched.species)) if len(wmatched) else [],
        "by_mechanism_source_systems":{
            s["species"]:{"mechanism":s["mechanism"],"evidence_role":s["evidence_role"],
                          "matched_photo_cohorts":(
                              sorted(details.loc[details.species==s["species"],"cohort"].tolist())
                              if len(details) else [])}
            for s in registry["systems"]
        },
        "direct_biochemical_genotype_fitness_bridge_estimable":bool(
            photo_true_morph_validated>=MIN_TRUE_MORPH_CHEMISTRY_DIRECT_LINKS and
            len(wmatched)>=MIN_INDEPENDENT_MECHANISM_PHOTO_SYSTEMS
        ),
        "minimum_mechanism_matched_species_for_bridge":MIN_INDEPENDENT_MECHANISM_PHOTO_SYSTEMS,
        "minimum_phenotype_biochemical_genotype_validated_links":MIN_TRUE_MORPH_CHEMISTRY_DIRECT_LINKS,
        "hard_nonclaims":[
            "this is a purposive and historically studied sample, never worldwide frequency of selection mechanisms",
            "a photo colour class does not match the published organism colour variant without independently inspecting organ/colour/white clipping",
            "white/coloured photographs do not establish heritable morphs, fitness, plant age or population genotype",
            "absence from 1500 candidate species does not imply absence from natural flowering plants",
            "genetic white pigment loss and environmentally plastic seasonal whitening are distinct mechanisms",
            "photo-local mixture within 50km does not prove a polymorphic interbreeding population",
            "across primary articles, different fitness pathways cannot be summed or pooled as one causal effect",
        ],
        "confirmatory_decisions_changed":False,
    }
    return out,details


def main() -> None:
    parser=argparse.ArgumentParser()
    parser.add_argument("--registry",required=True,type=Path)
    for cohort in COHORT_SHA:
        parser.add_argument("--"+cohort,required=True,type=Path)
    parser.add_argument("--outdir",required=True,type=Path)
    args=parser.parse_args()
    args.outdir.mkdir(parents=True,exist_ok=True)
    reg=registry_check(args.registry)
    cohorts={c:read_photos(getattr(args,c),c) for c in COHORT_SHA}
    out,rows=join_panel(reg,cohorts)
    out["photo_cohort_sha256"]=COHORT_SHA
    out["source_registry_sha256"]=sha256(args.registry)
    (args.outdir/"result.json").write_text(
        json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    rows.to_csv(args.outdir/"study_system_photo_overlap.csv",index=False)
    print(json.dumps(out,indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()
