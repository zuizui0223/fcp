#!/usr/bin/env python3
"""FCP all 4,730 U100 species: append every *unmeasured* species, original H9 only.

This code does not alter the frozen 149,900 measured photographs or any ROI-v4,
nine-colour/4-state, H1/D, 50-km or H2/W scientific computation. It selects
ALL unmeasured U100 species by immutable identity, never by observed colour.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import pandas as pd
from scripts.acquisition.expand_existing_fcp_highdepth_20261017 import (
    PER_SPECIES, P100_SHA, P500_SHA, ALLOWED_LICENSES,
    EXCLUSION_FILES, EXCLUSION_TOTAL_N, load_exclusions, sha256, source_rows,
)
from scripts.analysis.append_wave4_to_original_149900_20261017 import OLD_SOURCE

U100 = 4730
ORIGINAL_SPECIES = 1499
WAVE4_SPECIES = 498
REMAINING = U100 - ORIGINAL_SPECIES - WAVE4_SPECIES  # 2733
SHARDS = 4
SOURCE_METADATA_COMMIT = "ca69930986e39e3ea9b2d2f97ab247045ba9db0b"
WAVE4_ARCHIVE_SHA = "acc89e893367c91110bc6bc5d3a983af2cee72ed98b2dbb528f3e9172d7728ab"
SCHEMA = "fcp_u100_all4730_original_h9_only_v1"
LEGACY = ("discovery", "validation", "prospective")


def frozen_old_taxa(paths: dict[str, Path]) -> dict[int, str]:
    if set(paths) != set(LEGACY):
        raise RuntimeError("all 3 historical 149900-photo cohorts required")
    old: dict[int, str] = {}
    ids: set[int] = set()
    for group in LEGACY:
        path = paths[group]
        if sha256(path) != OLD_SOURCE[group]["sha256"]:
            raise RuntimeError("original measured source SHA differs: " + group)
        df = pd.read_csv(path, usecols=["inat_taxon_id", "species", "photo_id"], low_memory=False)
        if len(df) != OLD_SOURCE[group]["rows"]:
            raise RuntimeError("historical photo denominator differs: " + group)
        if df.inat_taxon_id.nunique() != OLD_SOURCE[group]["species"]:
            raise RuntimeError("historical taxon count differs: " + group)
        if not df.groupby("inat_taxon_id").size().eq(PER_SPECIES).all():
            raise RuntimeError("original species do not have 100 frozen photographs")
        for v in df[["inat_taxon_id", "species"]].drop_duplicates().itertuples(index=False):
            key = int(v.inat_taxon_id)
            if key in old:
                raise RuntimeError("historical taxon identity duplicated")
            old[key] = str(v.species).strip()
        photos = set(pd.to_numeric(df.photo_id, errors="raise").astype("int64"))
        if len(photos) != len(df) or photos & ids:
            raise RuntimeError("old high-depth photograph IDs duplicate")
        ids.update(photos)
    if len(old) != ORIGINAL_SPECIES or len(ids) != 149900:
        raise RuntimeError("original 1499 species and 149900 photos not reproduced")
    return old


def validated_wave4(metadata_dir: Path, selection_dir: Path) -> tuple[dict[int, str], dict]:
    report = json.loads((metadata_dir / "metadata_result.json").read_text())
    selection = json.loads((selection_dir / "selection_result.json").read_text())
    allowed = metadata_dir / "wave4_authorized_metadata.csv.gz"
    if (selection.get("wave4_selected_species") != 500
        or selection.get("remaining_candidate_species") != 2730
        or report.get("metadata_gate_pass") is not True
        or report.get("new_exact_100_species") != WAVE4_SPECIES
        or report.get("new_exact_100_photo_rows") != WAVE4_SPECIES * PER_SPECIES
        or report.get("original_photos") != 149900
        or report.get("image_pixels_opened") is not False
        or report.get("selection_sha256") != selection.get("new_selection_manifest_sha256")
        or report.get("outputs", {}).get(allowed.name) != sha256(allowed)):
        raise RuntimeError("new wave4 source photo IDs not authenticated")
    df = pd.read_csv(allowed, low_memory=False)
    if (len(df) != 49800 or df.inat_taxon_id.nunique() != WAVE4_SPECIES
        or df.photo_id.nunique() != 49800 or df.observation_id.nunique() != 49800
        or not df.groupby("inat_taxon_id").size().eq(PER_SPECIES).all()):
        raise RuntimeError("completed wave4 source rows are not 498 × 100")
    return dict(zip(df.inat_taxon_id.astype(int), df.species.astype(str))), report


def construct_all_species_inventory(
    p100: list[dict], old: dict[int, str], wave: dict[int, str],
    *, u100: int = U100, old_n: int = ORIGINAL_SPECIES, wave_n: int = WAVE4_SPECIES,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Deterministic classification by identities, independent of floral phenotype."""
    p100map = {r["inat_taxon_id"]: r["species"] for r in p100}
    if len(p100map) != len(p100) or len(old) != old_n or len(wave) != wave_n:
        raise RuntimeError("U100 identity frame/count changed")
    if set(wave) - set(p100map) or set(old) & set(wave):
        raise RuntimeError("wave4 not inside original 3730 P100 or duplicates old data")
    for tid, name in old.items():
        if tid in p100map and p100map[tid] != name:
            raise RuntimeError("historical measured taxon name differs from parent frame")
    for tid, name in wave.items():
        if p100map[tid] != name:
            raise RuntimeError("wave4 taxon identity/name differs from parent frame")
    union = dict(p100map)
    for tid, name in old.items():
        if tid in union and union[tid] != name:
            raise RuntimeError("duplicate original U100 name mismatch")
        union[tid] = name
    if len(union) != u100 or len(set(union.values())) != u100:
        raise RuntimeError(f"U100 4730 authentic distinct taxa not reproduced, found {len(union)}")
    allrows = []
    for tid, name in sorted(union.items()):
        status = ("original_149900_photo_measured" if tid in old
                  else "wave4_49800_photo_id_authorized" if tid in wave
                  else "unmeasured_original_u100_acquire")
        allrows.append({"inat_taxon_id":tid, "species":name, "u100_status":status,
                        "original_source_photo_rows":PER_SPECIES if tid in old else 0,
                        "wave4_source_photo_rows":PER_SPECIES if tid in wave else 0,
                        "new_candidate_photo_rows":0 if tid in old or tid in wave else PER_SPECIES})
    full = pd.DataFrame(allrows)
    pending = full.loc[full.u100_status == "unmeasured_original_u100_acquire",
                       ["inat_taxon_id", "species"]].reset_index(drop=True)
    if len(pending) != u100-old_n-wave_n:
        raise RuntimeError("ALL unmeasured original U100 taxa not included")
    if full.new_candidate_photo_rows.sum() != len(pending) * PER_SPECIES:
        raise RuntimeError("pending exactly-100 photo accounting drift")
    audit = {
        "schema": SCHEMA,
        "stage": "all_species_identity_inventory_pre_colour",
        "original_U100_distinct_species":u100,
        "original_measured_species":old_n,
        "original_measured_photos":old_n*PER_SPECIES,
        "wave4_authorized_source_species":wave_n,
        "wave4_authorized_source_photo_IDs":wave_n*PER_SPECIES,
        "new_unmeasured_species":int(len(pending)),
        "new_max_metadata_photos":int(len(pending)*PER_SPECIES),
        "full_U100_if_every_species_authorized_photo_rows":u100*PER_SPECIES,
        "pixels_measured_in_this_stage":0,
        "species_selection_conditional_on_flower_colour":False,
        "all_species_included_even_if_prior_H9_attempt_failed":True,
        "old_149900_result_changed":False,
        "status": "IDENTITY_ONLY_NOT_NEW_PIXEL_MEASUREMENT",
    }
    return full, pending, audit


def do_inventory(args: argparse.Namespace) -> dict:
    p100 = source_rows(args.p100, P100_SHA, 3730)
    source_rows(args.p500, P500_SHA, 500)  # cross-check frozen parent frame
    old = frozen_old_taxa({
        "discovery":args.discovery, "validation":args.validation, "prospective":args.prospective
    })
    wave, metadata = validated_wave4(args.wave4_metadata, args.wave4_selection)
    full, pending, audit = construct_all_species_inventory(p100, old, wave)
    args.out.mkdir(parents=True, exist_ok=True)
    a = args.out/"full_u100_4730_species_status.csv"
    b = args.out/"unmeasured_2733_taxa.csv"
    full.to_csv(a, index=False, lineterminator="\n")
    pending.to_csv(b, index=False, lineterminator="\n")
    audit["inputs_sha256"] = {
        "p100":sha256(args.p100),"p500":sha256(args.p500),
        **{g:sha256(getattr(args,g)) for g in LEGACY},
        "wave4_authorized_metadata":sha256(args.wave4_metadata/"wave4_authorized_metadata.csv.gz"),
        "wave4_selected_manifest":sha256(args.wave4_selection/"selected_species_manifest.tsv"),
    }
    audit["outputs_sha256"] = {a.name:sha256(a),b.name:sha256(b)}
    (args.out/"inventory_result.json").write_text(json.dumps(audit,indent=2,sort_keys=True)+"\n")
    return audit


def do_shard(args: argparse.Namespace) -> dict:
    """Unmodified original freeze_h9_metadata; only workload partition differs."""
    from fcp_pipeline.random_photo_h9_pool import freeze_h9_metadata
    from fcp_pipeline.random_photo_pool import InaturalistObservationClient
    from fcp_pipeline.third_cohort_metadata_gate import (
        seal_exact_denominator, validate_fresh_rows,
    )
    if not 0 <= args.shard < SHARDS:
        raise RuntimeError("wrong bulk-acquisition shard index")
    original = json.loads((args.inventory/"inventory_result.json").read_text())
    if original.get("schema")!=SCHEMA or original.get("new_unmeasured_species")!=REMAINING:
        raise RuntimeError("unverified all-4730 original opportunity frame")
    pendingpath=args.inventory/"unmeasured_2733_taxa.csv"
    if sha256(pendingpath)!=original["outputs_sha256"][pendingpath.name]:
        raise RuntimeError("pending identity manifest SHA mismatch")
    pending=pd.read_csv(pendingpath,low_memory=False)
    if len(pending)!=REMAINING or pending.inat_taxon_id.nunique()!=REMAINING:
        raise RuntimeError("all remaining 2733 species missing")
    pending=pending.iloc[args.shard::SHARDS].copy().reset_index(drop=True)
    pending.insert(0,"prospective_rank",range(1,len(pending)+1))
    if len(pending)==0: raise RuntimeError("empty original H9 shard")
    exc_obs,exc_photo,old_ledger=load_exclusions(args.exclusions)
    completed = pd.read_csv(args.wave4_authorized,low_memory=False)
    if (len(completed)!=49800 or completed.inat_taxon_id.nunique()!=498
        or completed.photo_id.nunique()!=49800 or completed.observation_id.nunique()!=49800):
        raise RuntimeError("wave4 authorized source 49800 differs")
    already_obs=set(completed.observation_id.astype("int64"))
    already_photos=set(completed.photo_id.astype("int64"))
    if exc_obs & already_obs or exc_photo & already_photos:
        raise RuntimeError("wave4 and historical photo ID exclusion ledgers overlap")
    exc_obs.update(already_obs)
    exc_photo.update(already_photos)
    if len(exc_obs)!=EXCLUSION_TOTAL_N+49800 or len(exc_photo)!=EXCLUSION_TOTAL_N+49800:
        raise RuntimeError("combined exclusion union changed")
    # No modification of the historical H9 acquisition callable or its arguments.
    client=InaturalistObservationClient(
        user_agent="zuizui0223-fcp-expanded-highdepth/1.0 (github.com/zuizui0223/fcp)",
        request_interval_seconds=1.05,timeout_seconds=45.0,max_retries=0,
    )
    frozen=freeze_h9_metadata(
        client=client, species_frame=pending[["inat_taxon_id","species"]].copy(),
        exclusion_observation_ids=exc_obs, exclusion_photo_ids=exc_photo,
        per_page=200,observer_cap_n=2,fixed_raw_photos=100,
        maximum_positional_accuracy_m=5000,allowed_photo_licenses=ALLOWED_LICENSES,
    )
    audit=frozen.species_audit.merge(
        pending[["inat_taxon_id","prospective_rank"]],on="inat_taxon_id",how="left",validate="one_to_one"
    ).sort_values("prospective_rank",kind="stable")
    obs=frozen.observations
    if obs.empty: raise RuntimeError("no metadata collected from 2733-U100 shard")
    obs=obs.merge(
        pending[["inat_taxon_id","prospective_rank"]],on="inat_taxon_id",how="left",validate="many_to_one"
    ).sort_values(["prospective_rank","h9_selection_order","observation_id","photo_id"],kind="stable")
    validate_fresh_rows(obs,exc_obs,exc_photo)
    successful,authorized=seal_exact_denominator(audit,obs,target_n=100)
    nerr=int(audit.request_error.fillna("").astype(str).str.len().gt(0).sum())
    failed=len(pending)-len(successful)
    verdict="PASS" if nerr / len(pending) <= .05 else "REQUEST_ERROR_CEILING_EXCEEDED"
    args.out.mkdir(parents=True,exist_ok=True)
    files={}
    for key,table,name in [
        ("candidate",obs,"candidate_metadata.csv.gz"),
        ("authorized",authorized,"authorized_metadata.csv.gz"),
        ("species_audit",audit,"species_audit.csv"),
        ("selected",pending,"selected_species.csv"),
    ]:
        p=args.out/name
        if name.endswith(".gz"):
            table.to_csv(p,index=False,lineterminator="\n",
                         compression={"method":"gzip","compresslevel":9,"mtime":0})
        else: table.to_csv(p,index=False,lineterminator="\n")
        files[name]=sha256(p)
    record={
        "schema":SCHEMA,"stage":"historical_H9_shard_metadata_no_pixels",
        "shard":args.shard,"shards":SHARDS,
        "selected_species":len(pending),"complete_100_photo_species":len(successful),
        "below_100_photo_species":failed,"candidate_rows":len(obs),
        "authorized_original_photo_IDs":len(authorized),"request_errors":nerr,
        "metadata_verdict":verdict,"excluded_original_photo_IDs":len(exc_photo),
        "old_frozen_exclusion_audit":old_ledger,"wave4_authorized_sha256":sha256(args.wave4_authorized),
        "all_U100_manifest_sha256":sha256(pendingpath),
        "outputs_sha256":files,"image_pixels_opened":False,
        "same_original_H9_function":True,
        "same_original_H9_params":True,"no_biological_outcome_selection":True,
        "original_149900_replaced":False,
    }
    (args.out/"shard_result.json").write_text(json.dumps(record,indent=2,sort_keys=True,allow_nan=False)+"\n")
    return record



def reconcile_historical_global_seen(
    frames: list[pd.DataFrame], expected_taxa: set[int], *, target_n: int = PER_SPECIES,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Conservatively restore the OLD freeze_h9_metadata global seen IDs.

    Old source 200-result pages were NOT preserved, so never invent substitute
    photos: replay its increasing taxon-ID order using ONLY already selected H9
    candidate metadata from all four source shards. Even partial taxa reserve
    their retained photo and observation IDs, exactly as the original code.
    Taxa with fewer than 100 surviving original IDs remain unresolved.
    """
    required={"inat_taxon_id","observation_id","photo_id"}
    if len(frames)!=SHARDS or any(required-set(g.columns) for g in frames):
        raise RuntimeError("four original H9 candidate tables are mandatory")
    source=pd.concat(frames,ignore_index=True,sort=False)
    if source.empty or not set(source.inat_taxon_id).issubset(expected_taxa):
        raise RuntimeError("original H9 candidate species missing or not in 2733")
    kept_frames=[];details=[];seen_obs:set[int]=set();seen_photo:set[int]=set()
    for tid in sorted(expected_taxa):
        x=source.loc[source.inat_taxon_id.eq(tid)].copy()
        if len(x)>target_n:
            raise RuntimeError("original H9 selected more than 100 photos in one species")
        if x.observation_id.nunique()!=len(x) or x.photo_id.nunique()!=len(x):
            raise RuntimeError("original H9 within-species candidate source duplicates")
        order=(["h9_selection_order"] if "h9_selection_order" in x.columns else [])
        if len(x):x=x.sort_values(order+["observation_id","photo_id"],kind="stable")
        allowed=[];repeated=0
        for row in x.itertuples(index=False):
            oid,pid=int(row.observation_id),int(row.photo_id)
            if oid in seen_obs or pid in seen_photo:
                repeated+=1
                continue
            allowed.append(pid)
            # Original H9 reserves all retained photos, even for taxa
            # whose final count is below the exact 100-photo gate.
            seen_obs.add(oid);seen_photo.add(pid)
        if len(allowed)==target_n:
            accepted=x.loc[x.photo_id.astype("int64").isin(allowed)].copy()
            if len(accepted)!=target_n:
                raise RuntimeError("cannot match unchanged retained original H9 rows")
            kept_frames.append(accepted)
        details.append({"inat_taxon_id":int(tid),
                        "global_selected_original_H9_rows":int(len(x)),
                        "global_cross_shard_rejected_rows":int(repeated),
                        "global_retained_original_H9_rows":len(allowed),
                        "global_complete_100_photo_species":len(allowed)==target_n})
    species=pd.DataFrame(details)
    authorized=(pd.concat(kept_frames,ignore_index=True,sort=False)
                if kept_frames else source.iloc[:0].copy())
    if (authorized.photo_id.nunique()!=len(authorized)
        or authorized.observation_id.nunique()!=len(authorized)
        or int(species.global_complete_100_photo_species.sum())*target_n!=len(authorized)):
        raise RuntimeError("globally duplicated original H9 candidate photographs remain")
    report={
       "original_H9_global_seen_id_reservation_order_restored":True,
       "original_H9_source_query_selection_code_unchanged":True,
       "original_four_shard_inputs_untouched":True,
       "original_200_result_API_pages_not_requeried":True,
       "rejected_photos_not_replaced":True,
       "original_149900_results_unchanged":True,
       "globally_seen_cross_shard_rejected_source_candidate_rows":int(
            species.global_cross_shard_rejected_rows.sum()),
       "original_source_species_affected_by_cross_shard_rejection":int(
            species.global_cross_shard_rejected_rows.gt(0).sum()),
       "post_source_global_100photo_species":int(
            species.global_complete_100_photo_species.sum()),
       "post_source_global_100photo_rows":len(authorized),
       "audit_rule":"all original candidates reserved in ascending taxon id, including incomplete species; if cross-shard identity reuse reduces original selected photos below 100, exclude whole species rather than fabricate another photograph",
    }
    return authorized,species,report



def do_consolidate(args: argparse.Namespace) -> dict:
    inventory=json.loads((args.inventory/"inventory_result.json").read_text())
    if inventory.get("schema")!=SCHEMA or inventory["new_unmeasured_species"]!=REMAINING:
        raise RuntimeError("all species U100 original manifest missing")
    full=pd.read_csv(args.inventory/"full_u100_4730_species_status.csv")
    if len(full)!=U100 or sha256(args.inventory/"full_u100_4730_species_status.csv")!=inventory["outputs_sha256"]["full_u100_4730_species_status.csv"]:
        raise RuntimeError("4730 original identities cannot be audited")
    seen_ids=set(); audits=[];original_complete=[];candidate_frames=[];new_audit=[]
    for shard in range(SHARDS):
        folder=args.parts/f"shard_{shard}"
        d=json.loads((folder/"shard_result.json").read_text())
        if d.get("schema")!=SCHEMA or d.get("shard")!=shard or d.get("shards")!=SHARDS:
            raise RuntimeError("not exactly 4 original H9 shards")
        if d.get("same_original_H9_function") is not True or d.get("image_pixels_opened") is not False:
            raise RuntimeError("shard altered original H9/photo-colour split")
        for name,digest in d["outputs_sha256"].items():
            if sha256(folder/name)!=digest:
                raise RuntimeError("shard content SHA mismatch: "+name)
        selected=pd.read_csv(folder/"selected_species.csv")
        if seen_ids & set(selected.inat_taxon_id):
            raise RuntimeError("selected taxa repeat across shards")
        seen_ids.update(selected.inat_taxon_id)
        au=pd.read_csv(folder/"authorized_metadata.csv.gz",low_memory=False)
        if au.photo_id.nunique()!=len(au) or au.observation_id.nunique()!=len(au):
            raise RuntimeError("original H9 duplicated IDs inside one source shard")
        candidate=pd.read_csv(folder/"candidate_metadata.csv.gz",low_memory=False)
        if (not set(candidate.inat_taxon_id).issubset(set(selected.inat_taxon_id))
            or candidate.photo_id.nunique()!=len(candidate)
            or candidate.observation_id.nunique()!=len(candidate)):
            raise RuntimeError("invalid original H9 candidate species or source IDs")
        authorized_pairs=set(zip(au.inat_taxon_id,au.observation_id,au.photo_id))
        source_pairs=set(zip(candidate.inat_taxon_id,candidate.observation_id,candidate.photo_id))
        if not authorized_pairs.issubset(source_pairs):
            raise RuntimeError("a frozen originally authorized photo is absent from H9 candidates")
        candidate["bulk_shard"]=shard
        candidate_frames.append(candidate)
        if len(au)!=d["authorized_original_photo_IDs"] or au.inat_taxon_id.nunique()!=d["complete_100_photo_species"]:
            raise RuntimeError("all photo rows incomplete for source shard")
        if len(au)>0 and not au.groupby("inat_taxon_id").size().eq(PER_SPECIES).all():
            raise RuntimeError("one or more species have fewer than 100 original photos")
        audit=pd.read_csv(folder/"species_audit.csv")
        if len(audit)!=len(selected):
            raise RuntimeError("lost attempted taxon audit")
        audit["bulk_shard"]=shard
        new_audit.append(audit)
        original_complete.append(au)
        audits.append(d)
    expected=set(pd.read_csv(args.inventory/"unmeasured_2733_taxa.csv").inat_taxon_id)
    if seen_ids != expected or len(seen_ids)!=REMAINING:
        raise RuntimeError("not all remaining 2733 original U100 species attempted")
    originally_full=sum(d["complete_100_photo_species"] for d in audits)
    errs=sum(d["request_errors"] for d in audits)
    original_rows=sum(d["authorized_original_photo_IDs"] for d in audits)
    if original_rows!=originally_full*PER_SPECIES:
        raise RuntimeError("source sharded original 100-photo denominator changed")
    new,global_audit,reconciliation=reconcile_historical_global_seen(candidate_frames,expected)
    old_authorized=pd.concat(original_complete,ignore_index=True,sort=False)
    chosen=set(zip(new.inat_taxon_id,new.observation_id,new.photo_id))
    prior=set(zip(old_authorized.inat_taxon_id,old_authorized.observation_id,old_authorized.photo_id))
    if not chosen.issubset(prior):
        raise RuntimeError("global uniqueness invented a non-authorized original photo")
    good=int(new.inat_taxon_id.nunique())
    rows=len(new)
    audit=pd.concat(new_audit,ignore_index=True)
    if len(audit)!=REMAINING or audit.inat_taxon_id.nunique()!=REMAINING:
        raise RuntimeError("all 2733 original species results not accounted")
    audit=audit.merge(global_audit,on="inat_taxon_id",how="left",validate="one_to_one")
    if len(audit)!=REMAINING or audit.global_complete_100_photo_species.isna().any():
        raise RuntimeError("all 2733 source-global H9 reconciled species not audited")
    if (int(audit.global_complete_100_photo_species.sum())!=good or rows!=good*PER_SPECIES
        or new.photo_id.nunique()!=len(new) or new.observation_id.nunique()!=len(new)):
        raise RuntimeError("global original H9 source duplicate photos remain")
    found=set(new.inat_taxon_id)
    full["final_status"]=full.u100_status
    full.loc[full.inat_taxon_id.isin(found),"final_status"]="new_metadata_100_of_100_photo_IDs"
    failure=set(expected)-found
    full.loc[full.inat_taxon_id.isin(failure),"final_status"]="new_H9_metadata_under_100_or_request_error"
    assert len(full)==4730
    args.out.mkdir(parents=True,exist_ok=True)
    newfile=args.out/"remaining_U100_authorized_metadata.csv.gz"
    new.to_csv(newfile,index=False,lineterminator="\n",
               compression={"method":"gzip","compresslevel":9,"mtime":0})
    afile=args.out/"all_2733_attempted_species_audit.csv"
    audit.to_csv(afile,index=False,lineterminator="\n")
    ffile=args.out/"full_U100_4730_status_after_H9.csv"
    full.to_csv(ffile,index=False,lineterminator="\n")
    result={
       "schema":SCHEMA,"stage":"all_4730_H9_metadata_terminal",
       "U100_identity_species":U100,"old_measured_species":ORIGINAL_SPECIES,
       "wave4_original_photo_metadata_species":WAVE4_SPECIES,
       "new_attempted_species":REMAINING,"new_complete_100_photo_species":good,
       "new_unmet_100_photo_species":len(failure),
       "new_original_photo_IDs":rows,
       "all_original_photo_IDs_if_wave4_pixels_finish":149900+49800+rows,
       "all_species_with_existing_or_full_100_photo_IDs":ORIGINAL_SPECIES+WAVE4_SPECIES+good,
       "species_under_100":len(failure),
       "request_errors":errs,
       "request_error_fraction":errs/REMAINING,
       "historical_source_global_ID_reconciliation":reconciliation,
       "source_shard_full100_species_before_global_ID_reconciliation":originally_full,
       "source_shard_original_photo_rows_before_global_ID_reconciliation":original_rows,
       "metadata_gate_pass":errs/REMAINING <=.05 and good>=int(.6*REMAINING),
       "photo_pixels_opened":False,
       "verified_against_all_4_complete_shards":True,
       "all_2733_attempted_even_if_failed":True,
       "overclaimed_4730_measured":False,
       "old_149900_method_and_results_changed":False,
       "inputs_sha256":{f"shard_{i}":audits[i]["outputs_sha256"] for i in range(SHARDS)},
       "outputs_sha256":{p.name:sha256(p) for p in (newfile,afile,ffile)},
    }
    (args.out/"all4730_metadata_result.json").write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
    return result


def main()->None:
    p=argparse.ArgumentParser()
    sub=p.add_subparsers(dest="stage",required=True)
    a=sub.add_parser("inventory")
    for x in ("p100","p500","discovery","validation","prospective"):
        a.add_argument("--"+x,type=Path,required=True)
    a.add_argument("--wave4-metadata",type=Path,required=True)
    a.add_argument("--wave4-selection",type=Path,required=True)
    a.add_argument("--out",type=Path,required=True)
    b=sub.add_parser("acquire-shard")
    b.add_argument("--inventory",type=Path,required=True)
    b.add_argument("--wave4-authorized",type=Path,required=True)
    b.add_argument("--exclusions",type=Path,nargs=5,required=True)
    b.add_argument("--shard",type=int,required=True)
    b.add_argument("--out",type=Path,required=True)
    c=sub.add_parser("consolidate")
    c.add_argument("--inventory",type=Path,required=True)
    c.add_argument("--parts",type=Path,required=True)
    c.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    res=do_inventory(args) if args.stage=="inventory" else do_shard(args) if args.stage=="acquire-shard" else do_consolidate(args)
    print("FCP_FULL_U100_"+args.stage.upper().replace("-","_")+"_COMPLETE",json.dumps(res,sort_keys=True,allow_nan=False))


if __name__=="__main__":
    main()
