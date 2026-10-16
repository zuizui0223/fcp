#!/usr/bin/env python3
"""Re-evaluate the REAL six-taxon original ROI-v4 photo TECHNICAL pilot.

This examines the genuine 552 new original image measurements preserved in
GitHub Actions #38832829706, but those photos came from a separate LIVE iNat
API query and were NOT verified against the 2026-09-27 dated bulk archive.

Output: descriptive observer-disjoint repeated PHOTO-LABEL detection and
within-ROI colour-mixture warning. Never a verified biological FCP species
count, the registered 1,000-taxon RGFCA P_poly or a genetic polymorphism.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

from fcp_pipeline.photo_polymorphism_screen import (
    MORPHS, PALETTE, POSITIVE, _finalize, sha256,
    observer_split,
)
from fcp_pipeline.rgfca_photo_repeatability import _photo_priority

RUN_ID = 38832829706
SOURCE_ROSTER_SHA = "88aba565e1bdeef48ad3dc5eec0121c52448fda434aca849e542bb615f9d547c"
NEW_ROI_SHA = "2c85859721b2f1d23dbc686aa526163ad69f6d48710aa1281b8123b3938b6047"
CANDIDATE_SHA = "ac2215e8b2a235609f9bc3ce73180438abbb3532a938f8e4554bd131641d9bd0"
DETECTOR_SHA = "f1aaeec4664fe2c178e5cf2bc1f508977bef3e4aa7b40613026cb8ae3de789d5"
ENCODER_SHA = "84ed466ffcc5c1f8d08409bc34a23bb364ab2c15e402cb12d4335a42be0e0951"
DECODER_SHA = "a62f8fa5ea080447c0689418d69e58f1e83e0b7adf9c142e2bd9bcc8045c0b11"
DEPTHS = (20, 32, 40, 44, 48)
N_REPS = 200
SOURCE_TYPE = "iNaturalist_API_created_before_20260927_live_snapshot_NOT_dated_monthly_archive"


def read_csv(path: Path):
    with path.open("r", newline="", encoding="utf-8") as src:
        return list(csv.DictReader(src))


def _draw(taxon: int, species: str, photos: list[dict], repeat: int, depth: int) -> dict:
    if depth not in DEPTHS:
        raise ValueError("Only predeclared 20/32/40/44/48 fold depths for pilot")
    picked = []
    for fold in (0, 1):
        source = sorted((x for x in photos if x["fold"] == fold),
                        key=lambda x: (_photo_priority(repeat,taxon,x["photo_id"]),
                                       int(x["photo_id"])))
        observations = set()
        for row in source:
            if row["observation_id"] in observations:
                continue
            observations.add(row["observation_id"])
            picked.append(row)
            if len(observations) == depth:
                break
    raw_counts = Counter()
    seen_states = defaultdict(set)
    owners = {}
    split_count = [Counter(),Counter()]
    n_valid = 0
    for row in picked:
        if not row["classifiable"]:
            continue
        n_valid += 1
        raw_counts[row["morph"]] += 1
        seen_states[row["observation_id"]].add(row["morph"])
        before = owners.setdefault(row["observation_id"],row["observer_id"])
        if before != row["observer_id"]:
            raise ValueError("Contradictory observation observer: no independent evidence")
        split_count[row["fold"]][row["morph"]] += 1
    result = _finalize(taxon, species, len(photos), len(picked), n_valid,
                       raw_counts, seen_states, owners, split_count)
    n_fold = [sum(v.values()) for v in split_count]
    evaluable = n_valid >= 40 and min(n_fold) >= 20
    return {
        "n_photo":len(picked), "n_classifiable":n_valid,
        "evaluable":bool(evaluable),
        "photo_label_pair_detected":bool(evaluable and
                                        result["photographic_colour_status"] == POSITIVE),
    }


def repeat_descriptive(taxon: int, species: str, photos: list[dict],
                       *, repeat_count: int = N_REPS) -> dict:
    if not 1 <= repeat_count <= 200:
        raise ValueError("Pilot Monte Carlo count must be bounded 1-200")
    good = [x for x in photos if x["classifiable"]]
    counts = Counter(x["morph"] for x in good)
    fold_count = Counter(x["fold"] for x in good)
    rows = {
        "inat_taxon_id":taxon, "species":species,
        "pilot_photos_selected":len(photos),
        "pilot_original_flower_ROI_classifiable":len(good),
        "classifiable_fraction":round(len(good)/len(photos),8) if photos else 0.0,
        "photographic_colour_counts":json.dumps(dict(counts),sort_keys=True),
        "n_classifiable_fold_0":fold_count[0],
        "n_classifiable_fold_1":fold_count[1],
        "unreconciled_date_snapshot":True,
        "genetic_or_population_FCP_proven":False,
    }
    dual = sum(
        (x["fractions"]["white"]>=0.10 and
         x["fractions"]["yellow"]>=0.10) for x in good)
    rows["same_photo_white_and_yellow_both_ge_0_10"] = dual
    rows["same_photo_white_yellow_cooccurrence_fraction"] = (
        round(dual/len(good),8) if good else 0.0)
    rows["within_flower_architecture_QC_required"] = bool(
        (counts.get("white",0)>0 and counts.get("yellow_orange",0)>0)
        or (len(counts)>1))
    rows["photo_colour_structural_confounding_high"] = bool(
        rows["same_photo_white_yellow_cooccurrence_fraction"]>=.50 and
        counts.get("white",0)>0 and counts.get("yellow_orange",0)>0)
    rows["strict_image_FCP_species_confirmed"] = False
    for n in DEPTHS:
        trials = [_draw(taxon,species,photos,k,n)
                  for k in range(1,repeat_count+1)]
        ev = sum(x["evaluable"] for x in trials)
        positive = sum(x["photo_label_pair_detected"] for x in trials)
        rows[f"depth_{n}_repetitions"] = repeat_count
        rows[f"depth_{n}_evaluable_repetitions"] = ev
        rows[f"depth_{n}_positive_photo_label_repetitions"] = positive
        rows[f"depth_{n}_conditional_photo_label_detection"] = (
            "" if ev==0 else round(positive/ev,8))
        rows[f"depth_{n}_median_classifiable_photo_count"] = sorted(
            x["n_classifiable"] for x in trials)[repeat_count//2]
    return rows


def audit(pilot_root: Path, rgfca_roster: Path, outdir: Path) -> dict:
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError("Do not silently overwrite REAL technical photo QC")
    provenance = pilot_root/"source_precolour/rgfca_six_source_precolour_receipt.json"
    selected = pilot_root/"source_precolour/rgfca_six_source_photo_candidates_before_colour.csv"
    measured = pilot_root/"original_roi_results/rgfca_original_new_photo_six_taxa_live_API_technical.csv"
    terminal = pilot_root/"original_roi_results/rgfca_original_new_photo_six_taxa_live_API_technical_receipt.json"
    roster_sha = sha256(rgfca_roster)
    a,b = json.loads(provenance.read_text()),json.loads(terminal.read_text())
    if (roster_sha != SOURCE_ROSTER_SHA or
        a.get("schema")!="fcp_rgfca_1000_enrolled_real_original_API_six_taxon_pilot_v1"
        or a.get("status")!="SELECTED_ORIGINAL_PROVIDER_API_PHOTO_IDS_NO_IMAGE_PIXELS_YET"
        or a.get("source_ROSTER_sha256")!=roster_sha
        or a.get("source_kind") is not None
        or a.get("n_selected_original_photo_IDs")!=552
        or a.get("photo_selection_used_literature_phenotype") is not False
        or a.get("photo_selection_used_climate_soil_or_prior_photo_colours") is not False
        or a.get("NOT_THE_IMMUTABLE_20260927_METADATA_CENSUS") is not True
        or a.get("EXACT_20260927_ARCHIVE_PHOTO_ID_MEMBERSHIP_VERIFIED") is not False
        or sha256(selected)!=CANDIDATE_SHA
        or a.get("selected_metadata_sha256")!=CANDIDATE_SHA
        or b.get("status")!="REAL_NEW_ORIGINAL_PIXELS_TECHNICAL_PILOT_NOT_RGFCA_PPOLY"
        or b.get("source_kind")!=SOURCE_TYPE
        or b.get("exact_dated_20260927_source_membership_admitted") is not False
        or b.get("full_1000_species_ROI_study_complete") is not False
        or b.get("prior_FCP_measured_colours_loaded") is not False
        or b.get("original_pixels_ROI_v4_evaluated")!=552
        or b.get("n_classifiable_original_flower_photos")!=202
        or b.get("selected_real_original_photo_IDs")!=552
        or b.get("new_terminal_photo_csv_sha256")!=NEW_ROI_SHA
        or sha256(measured)!=NEW_ROI_SHA):
        raise ValueError("Real original pilot photo/source provenance mismatch: refuse biology claims")
    method = b.get("model_attestation",{})
    if (method.get("detector_sha256")!=DETECTOR_SHA
        or method.get("efficient_sam_encoder_sha256")!=ENCODER_SHA
        or method.get("efficient_sam_decoder_sha256")!=DECODER_SHA
        or method.get("no_old_measurement_csv_inputs") is not True
        or method.get("location_blind") is not True):
        raise ValueError("Original pixel ROI-v4 model weights not byte identified")
    roster = {r["inat_taxon_id"]:r for r in read_csv(rgfca_roster)}
    choices = read_csv(selected)
    original = {(r["inat_taxon_id"],r["photo_id"]):r for r in choices}
    if len(choices)!=552 or len(original)!=552:
        raise ValueError("Source metadata contains missing/repeated photo IDs")
    per_taxon = defaultdict(list)
    photos = read_csv(measured)
    if len(photos)!=552:
        raise ValueError("Real technical terminal CSV must cover all 552 original photos")
    for row in photos:
        ident,pid = row["inat_taxon_id"], row["photo_id"]
        key = (ident,pid)
        provenance_row = original.get(key)
        if provenance_row is None:
            raise ValueError("Technical ROI terminal not in fixed pre-colour source")
        for field in ("species","observation_id","observer_id","screen_fold"):
            if row[field]!=provenance_row[field]:
                raise ValueError("Original source photo biological attribution changed")
        if (ident not in roster
            or roster[ident]["species"]!=row["species"]
            or roster[ident]["opportunity_stratum"]!=row["source_opportunity_tier"]):
            raise ValueError("Technical pilot taxon not in original random 1000 species")
        if row["global_classifiable"] not in ("true","false"):
            raise ValueError("Missing real ROI flower-classification flag")
        good = row["global_classifiable"]=="true"
        if good != (row["measurement_status"]=="classified_four_state_morph"
                    and row["morph"] in MORPHS):
            raise ValueError("ROI-v4 morphological terminal status inconsistent")
        if not good and row["morph"]!="mixed_uncertain":
            raise ValueError("Unclassifiable ROI-v4 photo must not be used as a flower colour")
        frac = {p:float(row["flower_fraction_"+p]) for p in PALETTE}
        if any(not math.isfinite(v) or v<0 or v>1 for v in frac.values()):
            raise ValueError("Invalid nine-colour original flower ROI")
        if good and abs(sum(frac.values())-1)>1e-7:
            raise ValueError("Original flower ROI does not preserve nine-colour composition")
        fold = int(row["screen_fold"])
        if fold not in (0,1) or fold!=observer_split(int(ident),row["observer_id"]):
            raise ValueError("Photo-to-observer held-out fold changed")
        per_taxon[(int(ident),row["species"])].append({
            "photo_id":pid,"fold":fold,
            "observation_id":row["observation_id"],"observer_id":row["observer_id"],
            "classifiable":good,"morph":row["morph"] if good else None,
            "fractions":frac
        })
    if len(per_taxon)!=6:
        raise ValueError("Original 6-taxon independent imaging pilot lost taxa")
    results = [repeat_descriptive(t,s,pool) for (t,s),pool in sorted(per_taxon.items())]
    for row in results:
        original_counts=b["by_taxon"].get(str(row["inat_taxon_id"]))
        if (original_counts is None
            or original_counts["n_source_photos_in_technical_pilot"]!=row["pilot_photos_selected"]
            or original_counts["n_classifiable"]!=row["pilot_original_flower_ROI_classifiable"]):
            raise ValueError("Real ROI taxon outcome totals disagree with original pixel receipt")
    if sum(x["pilot_photos_selected"] for x in results)!=552 or sum(
            x["pilot_original_flower_ROI_classifiable"] for x in results)!=202:
        raise ValueError("Real source photo/classification totals do not reconcile")
    outdir.mkdir(parents=True,exist_ok=True)
    table = outdir/"real_six_taxa_photo_label_repeatability_TECHNICAL_ONLY.csv"
    with table.open("w",newline="",encoding="utf-8") as f:
        wr=csv.DictWriter(f,fieldnames=list(results[0]),lineterminator="\n")
        wr.writeheader()
        wr.writerows(results)
    technical = {
        "schema":"fcp_real_six_taxon_photo_repeatability_diagnostic_v1",
        "status":"TECHNICAL_PILOT_ONLY_NOT_DATED_SOURCE_PPOLY_OR_CONFIRMED_FCP",
        "source_GitHub_workflow_run":RUN_ID,
        "source_photo_ROI_receipt_sha256":sha256(terminal),
        "live_API_original_photo_source_receipt_sha256":sha256(provenance),
        "RGFCA_balanced_1000_taxon_roster_sha256":roster_sha,
        "selected_photo_metadata_sha256":sha256(selected),
        "new_real_original_photo_ROI_SHA256":sha256(measured),
        "n_species":6,
        "n_new_original_ROI_v4_photo_terminals":552,
        "n_classifiable_flower_ROI":202,
        "ROI_classifiable_fraction":202/552,
        "pilot_repeat_draws_per_species":N_REPS,
        "pilot_depths_per_independent_observer_fold":list(DEPTHS),
        "original_RGFCA_200x250_species_REAL_photo_Ppoly_computed":False,
        "source_20260927_exact_archived_ID_membership_verified":False,
        "independent_new_biological_replicates_from_Monte_Carlo":0,
        "genetic_or_natural_population_FCP_species_confirmed":0,
        "n_technical_label_cooccurrence_high":sum(
            x["photo_colour_structural_confounding_high"] for x in results),
        "manual_FCP_flower_architecture_and_developmental_stage_review_required":True,
        "frozen_original_ROI_v4_used":True,
        "prior_149900_colour_outcomes_used":False,
        "table_SHA256":sha256(table),
        "note":("Repeated observer-disjoint PHOTO LABEL robustness does not "
                "by itself diagnose genetic polymorphism: a single "
                "multi-colour inflorescence, anthesis stage or palette mask can "
                "produce both coarse classes in different photos."),
    }
    (outdir/"real_six_taxa_photo_repeatability_TECHNICAL_ONLY.json").write_text(
        json.dumps(technical,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    return technical


def main(argv=None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--unpacked-original-pilot",required=True,type=Path)
    p.add_argument("--rgfca-1000-roster",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    x=p.parse_args(argv)
    print("FCP_RGFCA_REAL_SIX_PIXELS_TECHNICAL_QC",
          json.dumps(audit(x.unpacked_original_pilot,x.rgfca_1000_roster,x.outdir),
                     sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
