#!/usr/bin/env python3
"""Metadata-only outcome-blind prospective FCP selection: 2000 new + 730 reserve.

Rebuild 3230 source-frozen species from original P100/P500 allocation, then
remove the ENTIRE previously selected third set (500, including the one with
98 fresh photo records). No colour pixels, flower classifications, locations,
climate, genes, outcomes, historical H2 or prior manuscript statistic is read.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

P100_SHA256 = "1473aad680fe2fa84903c5e11ee104fd0eceef828957f16c2ef1a35f9dd6993c"
P500_SHA256 = "f54d07fb2338a20a3f92808b58f0ebc948ab2d5af3cb01fb26702f861184b7a4"
THIRD_MANIFEST_SHA256 = "16db6a2fc265f0ab20e7dae4e6f9058b907f5c19af3ddab5b6077797dd1def59"
CANDIDATE_3230_SHA256 = "7fc0337074b55a91c0b50773a8d5c3ef82e877074c8b3cacc21f9af58e0ba77e"
P100_ROWS = 3730
P500_ROWS = 500
THIRD_ROWS = 500
UNTOUCHED_ROWS = 2730
PRIMARY_SAMPLE = 2000
RESERVED_SAMPLE = 730
CAPACITY_FLOOR = 100
SALT = "FCP_NEXT2000_SPATIAL_20261008_FIXED_V1"
POOL_COLUMNS = ["inat_taxon_id", "species", "after_observer_cap"]
ALLOWED_PARENT_COLUMNS = set(POOL_COLUMNS) | {"selection_hash", "prospective_rank"}
MANIFEST_COLUMNS = ["prospective_rank", "inat_taxon_id", "species"]
SAMPLED_COLUMNS = ["prospective_rank", "inat_taxon_id", "species", "after_observer_cap", "selection_hash"]
PROHIBITED_TERMS = {"morph","flower_colour","flower_color","white","pigment","colour","color","photos","latitude","longitude","span","bio5","geographic","h2","w"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for part in iter(lambda: f.read(1 << 20), b""):
            h.update(part)
    return h.hexdigest()


def read_parent(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as f:
        rd=csv.DictReader(f)
        names=set(rd.fieldnames or [])
        if not set(POOL_COLUMNS).issubset(names) or names-ALLOWED_PARENT_COLUMNS:
            raise ValueError(f"Parent photo opportunity schema mismatch: {sorted(names)}")
        rows=[]
        for item in rd:
            try:
                tid=int(item["inat_taxon_id"])
                n=int(item["after_observer_cap"])
            except (ValueError,TypeError) as e:
                raise ValueError("Invalid source taxon identity or opportunity depth") from e
            species=str(item["species"]).strip()
            if tid<=0 or len(species.split())<2 or n<CAPACITY_FLOOR:
                raise ValueError(f"Invalid identity or below-U100 capacity: {tid} / {species} / {n}")
            rows.append({"inat_taxon_id":tid,"species":species,"after_observer_cap":n})
    check_unique(rows,"source pool")
    return rows


def read_third_manifest(path: Path) -> list[dict]:
    with path.open(newline="",encoding="utf-8") as f:
        rd=csv.DictReader(f,delimiter="\t")
        if list(rd.fieldnames or [])!=MANIFEST_COLUMNS:
            raise ValueError("Third-cohort manifest schema is not the frozen three-column identity")
        records=[]
        for r in rd:
            records.append({"prospective_rank":int(r["prospective_rank"]),
                            "inat_taxon_id":int(r["inat_taxon_id"]),
                            "species":str(r["species"]).strip()})
    if sorted(r["prospective_rank"] for r in records)!=list(range(1,len(records)+1)):
        raise ValueError("Third-cohort ranks non-unique, missing or out of order")
    check_unique(records,"third cohort")
    return records


def check_unique(rows: list[dict],name: str) -> None:
    for field in ("inat_taxon_id","species"):
        vals=[r[field] for r in rows]
        if len(vals)!=len(set(vals)):
            raise ValueError(f"Duplicate {field} in {name}")


def canonical_pool(path: Path,rows: list[dict]) -> None:
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",encoding="utf-8",newline="") as f:
        wr=csv.DictWriter(f,fieldnames=POOL_COLUMNS,lineterminator="\n")
        wr.writeheader()
        for r in rows: wr.writerow({k:r[k] for k in POOL_COLUMNS})


def allocate(p100: list[dict],p500: list[dict],third: list[dict]) -> tuple[list[dict],list[dict],list[dict]]:
    byid={r["inat_taxon_id"]:r for r in p100}
    p500ids={r["inat_taxon_id"] for r in p500}
    if not p500ids.issubset(byid):
        raise ValueError("Prior P500 contains an unknown P100 taxon")
    for r in p500:
        if byid[r["inat_taxon_id"]]["species"]!=r["species"]:
            raise ValueError("Prior P500 species identity drift")
    pool3230=sorted(
        (r.copy() for r in p100 if r["inat_taxon_id"] not in p500ids),
        key=lambda r:(r["inat_taxon_id"],r["species"]))
    check_unique(pool3230,"reconstructed third candidate pool")
    remaining={r["inat_taxon_id"]:r for r in pool3230}
    thirdids={r["inat_taxon_id"] for r in third}
    if not thirdids.issubset(remaining):
        raise ValueError("Third selected species outside 3230 parent pool")
    for r in third:
        if remaining[r["inat_taxon_id"]]["species"]!=r["species"]:
            raise ValueError("Third selected identity name has drifted")
    unallocated=[r.copy() for r in pool3230 if r["inat_taxon_id"] not in thirdids]
    check_unique(unallocated,"prospective untouched pool")
    # No colour/environment/geography can enter rank. Rank each remaining
    # species by a single precommitted SHA256 identity salt.
    ranked=[]
    for row in unallocated:
        item=row.copy()
        payload=f"{SALT}|{item['inat_taxon_id']}|{item['species']}".encode("utf-8")
        item["selection_hash"]=hashlib.sha256(payload).hexdigest()
        ranked.append(item)
    ranked.sort(key=lambda r:(r["selection_hash"],r["inat_taxon_id"]))
    for i,r in enumerate(ranked,start=1):r["prospective_rank"]=i
    return pool3230,ranked[:PRIMARY_SAMPLE],ranked[PRIMARY_SAMPLE:]


def write_manifest(path: Path,rows: list[dict]) -> None:
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",encoding="utf-8",newline="") as f:
        wr=csv.DictWriter(f,fieldnames=SAMPLED_COLUMNS,delimiter="\t",lineterminator="\n")
        wr.writeheader()
        for r in rows:wr.writerow({k:r[k] for k in SAMPLED_COLUMNS})


def run(p100:Path,p500:Path,third:Path,out:Path,*,check_fingerprints:bool=True) -> dict:
    actual={"p100":sha256(p100),"p500":sha256(p500),"third_manifest":sha256(third)}
    if check_fingerprints:
        expected={"p100":P100_SHA256,"p500":P500_SHA256,"third_manifest":THIRD_MANIFEST_SHA256}
        if actual!=expected: raise RuntimeError(f"Frozen metadata source SHA mismatch: {actual}")
    pp=read_parent(p100);pf=read_parent(p500);tt=read_third_manifest(third)
    if check_fingerprints and [len(pp),len(pf),len(tt)]!=[P100_ROWS,P500_ROWS,THIRD_ROWS]:
        raise RuntimeError(f"Previously allocated species count has drifted: {[len(pp),len(pf),len(tt)]}")
    pool3230,main,reserve=allocate(pp,pf,tt)
    if check_fingerprints and [len(pool3230),len(main),len(reserve)]!=[3230,PRIMARY_SAMPLE,RESERVED_SAMPLE]:
        raise RuntimeError("Prospective denominators are inconsistent with frozen global opportunity")
    out.mkdir(parents=True,exist_ok=True)
    prior=out/"reconstructed_3230_outcome_blind_pool.csv"
    canonical_pool(prior,pool3230)
    if check_fingerprints and sha256(prior)!=CANDIDATE_3230_SHA256:
        raise RuntimeError("Original candidate pool has failed canonical historical byte fingerprint")
    remaining=out/"unused_2730_outcome_blind_pool.csv"
    unused=sorted([r for r in pool3230 if r["inat_taxon_id"] not in
                   {x["inat_taxon_id"] for x in tt}],key=lambda r:(r["inat_taxon_id"],r["species"]))
    canonical_pool(remaining,unused)
    selected=out/"prospective_selected_2000_species.tsv"
    held=out/"prospective_reserve_730_species.tsv"
    write_manifest(selected,main);write_manifest(held,reserve)
    old={r["inat_taxon_id"] for r in pf}|{r["inat_taxon_id"] for r in tt}
    new={r["inat_taxon_id"] for r in main}
    rest={r["inat_taxon_id"] for r in reserve}
    if new&rest or new&old or rest&old or len(new|rest)!=len(unused):
        raise RuntimeError("Broken independence of previously allocated and newly reserved species")
    report={
        "schema":"fcp_next2000_prospective_outcome_blind_allocation_v1",
        "date_jst":"2026-10-08",
        "status":"METADATA_ONLY_COHORT_SELECTION_FROZEN_NO_BIOLOGICAL_OPENING",
        "original_global_discovered_species":42111,
        "original_u100_metadata_capacity_species":4730,
        "excluded_old_discovery_validation_species":1000,
        "excluded_old_p500_selected_even_when_not_measured":500,
        "excluded_old_third_selected_even_when_not_measured":500,
        "previous_observed_high_depth_photo_species":1499,
        "historically_allocated_capacity_species":2000,
        "source_candidate_species_3230":len(pool3230),
        "unused_capacity_species_after_all_allocations":len(unused),
        "new_prospective_sampled_species":len(main),
        "permanently_unopened_reserve_species":len(reserve),
        "potential_total_photographed_species_if_all_new_sample_successful":1499+len(main),
        "primary_photo_depth_target_per_species":100,
        "new_candidate_photo_rows_upper_design_target":len(main)*100,
        "allocation_salt":SALT,
        "ordering":"SHA256(salt|inat_taxon_id|species) ascending; integer taxon_id tie break",
        "source_sha256":actual,
        "historical_3230_canonical_sha256":sha256(prior),
        "unused_2730_canonical_sha256":sha256(remaining),
        "selected_2000_manifest_sha256":sha256(selected),
        "reserve_730_manifest_sha256":sha256(held),
        "selected_2000_first_taxon_ids":[r["inat_taxon_id"] for r in main[:10]],
        "selected_2000_no_prior_taxon_overlap":not bool(new&old),
        "reserve_730_no_prior_or_selected_overlap":not bool(rest&old or rest&new),
        "metadata_only_read_fields":["inat_taxon_id","species","after_observer_cap"],
        "biological_outcome_state":{
            "source_colour_morph_or_palette_labels_read":False,
            "new_photo_pixels_opened":False,
            "new_image_measurement_performed":False,
            "old_prospective_h2_result_consulted_for_species_selection":False,
            "latitude_or_environment_used_for_species_selection":False,
            "species_or_photo_replacement_after_outcome_opening_allowed":False,
            "new_confirmatory_result_available":False
        },
        "next_stage":"independent fresh photo/observation metadata acquisition using prior-ID exclusion and fixed selection; synthetic end-to-end qualification; one authorized colour opening only after durable metadata receipt",
        "nonclaims":[
            "U100 census is a historical 2026-09 opportunity scan; fresh usable >=100 photos are not guaranteed",
            "2000 selected nominal species do not imply 2000 >=40-classifiable measured species",
            "4299 versus 3499 confusion forbidden: 1499 previously photographed + 2000 newly selected = 3499 at most",
            "previous P500 and unused third 500th species are excluded even if their prospective outcome failed",
            "photo-system confirmation cannot establish all-angiosperm colour-polymorphism prevalence, genetic morphs, biochemical costs or net selection"
        ],
        "confirmatory_decisions_changed":False
    }
    (out/"preopening_result.json").write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2,ensure_ascii=False))
    return report


def main()->None:
    p=argparse.ArgumentParser()
    p.add_argument("--p100",required=True,type=Path)
    p.add_argument("--p500",required=True,type=Path)
    p.add_argument("--third-manifest",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    args=p.parse_args()
    run(args.p100,args.p500,args.third_manifest,args.outdir)


if __name__=="__main__":
    main()
