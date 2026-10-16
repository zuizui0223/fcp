#!/usr/bin/env python3
"""RGFCA-style repeated robustness of NEW photographic flower-colour polymorphism.

Primary estimand P_poly(i) = fraction of EVALUABLE balanced repeated photo draws
of species i that reproduce two original ROI-v4 colour states in both fixed,
observer-disjoint halves. It is a resampling STABILITY score conditional on
the same imaged source pool, NOT independent replication, genetic FCP, an
unbiased plant-world prevalence or a significance-test p-value.

Upstream species/photo enrollment is pre-colour and byte frozen. A second
species-specific 80-photo draw (40 per observer-disjoint fold) is repeated
according to the original RGFCA 200 by 250 equal-species schedule. Original
unclassifiable photographs stay in the sampling denominator.
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
    MORPHS, PALETTE, POSITIVE, PROTOCOL, SNAPSHOT, _finalize,
    observer_split, sha256,
)
from scripts.analysis.freeze_fcp_rgfca_photo_repetition import (
    N_SPECIES, REPS, PER_REP, SOURCE_PHOTOS, SOURCE_SHA, FRAME_SHA, SALT,
)


def _read(path: Path):
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def validate_schedule(cohort: Path, schedule: Path, receipt: Path,
                      *, synthetic=False) -> tuple[dict, list[dict]]:
    doc = json.loads(receipt.read_text(encoding="utf-8"))
    if (doc.get("schema") != "fcp_rgfca_photo_polymorphism_balanced_sample_v1"
        or doc.get("status") != "ORIGINAL_SOURCE_BALANCED_REPETITION_PLAN_NO_COLOUR_MEASURED"
        or doc.get("selected_species_csv_sha256") != sha256(cohort)
        or doc.get("realization_schedule_csv_sha256") != sha256(schedule)
        or doc.get("old_149900_flower_colours_used") is not False
        or doc.get("new_flower_images_opened") is not False
        or doc.get("photographic_poly_probabilities_measured") is not False
        or (not synthetic and (
            doc.get("selected_species") != N_SPECIES
            or doc.get("realizations") != REPS
            or doc.get("species_per_realization") != PER_REP
            or doc.get("species_inclusions_each") != 50
            or doc.get("original_source_count_sha256") != SOURCE_SHA
            or doc.get("original_frame_sha256") != FRAME_SHA
            or doc.get("source_photo_opportunities") != SOURCE_PHOTOS))):
        raise ValueError("Cannot load nonoriginal or phenotype-adaptive RGFCA schedule")
    species = _read(cohort)
    assignments = _read(schedule)
    index = {}
    for s in species:
        ident = int(s["inat_taxon_id"])
        if ident in index or not s["species"] or int(s["source_photo_opportunities"]) < 128:
            raise ValueError("Duplicate/ineligible species in RGFCA opportunity frame")
        index[ident] = s
    n_reps = int(doc["realizations"])
    per_rep = int(doc["species_per_realization"])
    n_visits = int(doc["species_inclusions_each"])
    seen, visits, by_rep = set(), Counter(), Counter()
    for r in assignments:
        rep, tid, slot = int(r["realization"]), int(r["inat_taxon_id"]), int(r["slot"])
        if (rep <= 0 or rep > n_reps or slot <= 0 or slot > per_rep
            or tid not in index or index[tid]["species"] != r["species"]
            or (rep, tid) in seen):
            raise ValueError("Changed/repeated original realization photo species")
        seen.add((rep, tid))
        visits[tid] += 1
        by_rep[rep] += 1
    if (len(index) != int(doc["selected_species"])
        or len(assignments) != n_reps * per_rep
        or set(visits.values()) != {n_visits}
        or set(by_rep) != set(range(1, n_reps + 1))
        or any(n != per_rep for n in by_rep.values())):
        raise ValueError("RGFCA species inclusion lost its exact balanced sampling")
    return index, assignments


def load_new_roi(measured: Path, receipt: Path, registry: dict,
                 *, synthetic=False) -> dict[int, list[dict]]:
    evidence = json.loads(receipt.read_text(encoding="utf-8"))
    if (evidence.get("schema") != "fcp_sealed_photo_screen_terminal_v1"
        or evidence.get("status") != "SEALED_PHOTO_SCREEN_UNBIASED_NOT_FINAL"
        or evidence.get("protocol") != PROTOCOL
        or evidence.get("source_date") != SNAPSHOT
        or evidence.get("new_photo_table_sha256") != sha256(measured)
        or evidence.get("historical_measured_flower_colours_reused") is not False
        or evidence.get("photo_ID_selection_pre_colour") is not True
        or evidence.get("source_photo_IDs_sealed") is not True
        or evidence.get("screen_not_entire_original_93million_census") is not True
        or evidence.get("full_source_image_measurement_complete") is not False
        or (not synthetic and
            (evidence.get("synthetic_test_only") is not False
             or evidence.get("original_ROI_v4_method_attested") is not True
             or evidence.get("source_frame_sha256") != FRAME_SHA
             or evidence.get("source_species") != 42111))):
        raise ValueError("No complete original fresh ROI-v4 observer-disjoint PHOTO screen")
    groups: dict[int, list[dict]] = defaultdict(list)
    all_ids = set()
    counted = 0
    with measured.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        required = {"inat_taxon_id","species","photo_id","observation_id",
                    "observer_id","screen_fold","morph","measurement_status",
                    "global_classifiable",*(f"flower_fraction_{x}" for x in PALETTE)}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Missing original nine-colour ROI/observer-fold status")
        for raw in reader:
            counted += 1
            pid, t = str(raw["photo_id"]), raw["inat_taxon_id"]
            if not pid.isdecimal() or pid in all_ids or not t.isdecimal():
                raise ValueError("Duplicate or invalid original measured source photo ID")
            all_ids.add(pid)
            ident = int(t)
            flag = str(raw["global_classifiable"]).lower()
            classified = raw["measurement_status"] == "classified_four_state_morph"
            if (flag not in {"true","false"}
                or (flag=="true") != classified
                or (classified and raw["morph"] not in MORPHS)
                or (not classified and raw["morph"] != "mixed_uncertain")):
                raise ValueError("Original flower-classification ROI terminal not coherent")
            if ident not in registry:
                continue
            if registry[ident]["species"] != raw["species"]:
                raise ValueError("RGFCA species changed after original colour-blind selection")
            observer, obs = raw["observer_id"], raw["observation_id"]
            fold = str(raw["screen_fold"])
            if (not observer or not obs or fold not in {"0","1"}
                or observer_split(ident, observer) != int(fold)):
                raise ValueError("Original observer/photo split was modified")
            if classified:
                v = [float(raw[f"flower_fraction_{p}"]) for p in PALETTE]
                if any(not math.isfinite(x) or x < 0 or x > 1 for x in v) or abs(sum(v)-1)>1e-7:
                    raise ValueError("Original nine-colour flower ROI fractions invalid")
            groups[ident].append({
                "photo_id":pid, "observer_id":observer, "observation_id":obs,
                "fold":int(fold), "morph":raw["morph"] if classified else None,
                "classifiable":classified,
            })
    if counted != evidence.get("n_terminal_measurements"):
        raise ValueError("Whole sealed source photo selection is not reproduced")
    for t, records in groups.items():
        if len(records) > 128:
            raise ValueError("Selected ORIGINAL photo opportunity pool exceeds 128 per species")
    return groups


def _photo_priority(rep: int, taxon: int, photo_id: str) -> str:
    key = "|".join((SALT, "repeated_photo", str(rep), str(taxon), photo_id))
    return hashlib.sha256(key.encode("utf-8")).hexdigest()


def one_draw(taxon: int, species: str, n_source: int, original_rows: list[dict],
             rep: int, *, per_fold: int = 40) -> dict:
    if not 1 <= per_fold <= 64:
        raise ValueError("Primary RGFCA source-photo draw is 40 per observer fold")
    selection = []
    for fold in (0, 1):
        pool = sorted((x for x in original_rows if x["fold"] == fold),
                      key=lambda x: (_photo_priority(rep,taxon,x["photo_id"]),x["photo_id"]))
        seen_obs, chosen = set(), 0
        for row in pool:
            if row["observation_id"] in seen_obs:
                continue
            seen_obs.add(row["observation_id"])
            selection.append(row)
            chosen += 1
            if chosen >= per_fold:
                break
    morph_counts = Counter()
    obs_to_labels: dict[str, set[str]] = defaultdict(set)
    obs_to_observer = {}
    split_photo = [Counter(),Counter()]
    good = 0
    for r in selection:
        if not r["classifiable"]:
            continue
        good += 1
        morph_counts[r["morph"]] += 1
        obs_to_labels[r["observation_id"]].add(r["morph"])
        prior = obs_to_observer.setdefault(r["observation_id"],r["observer_id"])
        if prior != r["observer_id"]:
            raise ValueError("One biological observation belongs to two observers")
        split_photo[r["fold"]][r["morph"]] += 1
    result = _finalize(taxon,species,n_source,len(selection),good,morph_counts,
                       obs_to_labels,obs_to_observer,split_photo)
    per_fold_good = [sum(x.values()) for x in split_photo]
    evaluable = good >= 40 and all(n >= 20 for n in per_fold_good)
    detected = evaluable and result["photographic_colour_status"] == POSITIVE
    return {
        "n_selected":len(selection), "n_classifiable":good,
        "fold_0_classifiable":per_fold_good[0],
        "fold_1_classifiable":per_fold_good[1],
        "evaluable":bool(evaluable), "photo_poly_detected":bool(detected),
        "coarse_morph_status":result["photographic_colour_status"]
    }


def repeat(registry: dict, assignments: list[dict], pool: dict[int,list[dict]],
           outdir: Path, *, synthetic=False) -> dict:
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError("RGFCA photo robustness output must be new, never overwritten")
    grouped = defaultdict(list)
    for row in assignments:
        grouped[int(row["inat_taxon_id"])].append(int(row["realization"]))
    species_rows = []
    wave_counts = {k:Counter() for k in (25,50,100,150,200)}
    for taxon in sorted(registry):
        item = registry[taxon]
        rows = pool.get(taxon, [])
        trials = []
        for rep in grouped[taxon]:
            trial = one_draw(taxon,item["species"],int(item["source_photo_opportunities"]),
                             rows,rep)
            trials.append((rep, trial))
            for cutoff in wave_counts:
                if rep<=cutoff:
                    wave_counts[cutoff]["inclusions"]+=1
                    wave_counts[cutoff]["evaluable"]+=int(trial["evaluable"])
                    wave_counts[cutoff]["positive"]+=int(trial["photo_poly_detected"])
        eligible = sum(v["evaluable"] for _,v in trials)
        positives = sum(v["photo_poly_detected"] for _,v in trials)
        reliability = positives/eligible if eligible else None
        if not rows or eligible < 30:
            tier = "not_evaluable_or_insufficient_repeated_roi_support"
        elif reliability is not None and reliability >= .8:
            tier = "repeatably_photo_polymorphic_phenotype"
        elif positives:
            tier = "photo_variation_sensitive_to_repeated_sampling"
        else:
            tier = "not_detected_in_evaluable_photo_resamples_not_monomorphic"
        species_rows.append({
            "inat_taxon_id":taxon,"species":item["species"],
            "original_source_photo_opportunities":item["source_photo_opportunities"],
            "original_new_roi_pool_photos":len(rows),
            "rgfca_species_inclusions":len(trials),
            "evaluable_40photo_observer_disjoint_draws":eligible,
            "positive_independent_observer_colour_draws":positives,
            "conditional_photo_colour_P_poly":"" if reliability is None else round(reliability,12),
            "repeatability_0_80_status":tier,
            "genetic_polymorphism_established":False,
            "monomorphic_if_not_detected":False,
        })
    if not synthetic and (len(species_rows)!=N_SPECIES or
                          any(int(r["rgfca_species_inclusions"])!=50 for r in species_rows)):
        raise ValueError("Original 1,000 species / 50 balanced inclusions not maintained")
    outdir.mkdir(parents=True, exist_ok=True)
    out_csv = outdir/"species_rgfca_repeated_photo_polymorphism.csv"
    with out_csv.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(species_rows[0]),lineterminator="\n")
        w.writeheader()
        w.writerows(species_rows)
    tiers=Counter(r["repeatability_0_80_status"] for r in species_rows)
    summary={
        "schema":"fcp_photo_rgfca_repeated_observer_disjoint_polymorphism_v1",
        "status":"PHOTO_SCREEN_RESAMPLE_STABILITY_ONLY_NOT_GENETIC_FCP",
        "selected_species":len(registry),
        "balanced_realizations":max(int(r["realization"]) for r in assignments),
        "number_photo_draw_assignments":len(assignments),
        "primary_80_photo_draw_and_40_photo_per_observer_fold":True,
        "n_repeated_phenotypic_polymorphism_species":tiers.get(
            "repeatably_photo_polymorphic_phenotype",0),
        "n_species_unresolved_or_not_evaluable":tiers.get(
            "not_evaluable_or_insufficient_repeated_roi_support",0),
        "detection_tiers":dict(tiers),
        "realization_convergence":{str(n):dict(x) for n,x in wave_counts.items()},
        "same_photos_may_be_reused_between_realizations":True,
        "number_independent_new_biological_worlds":1,
        "conditional_on_selected_1000_species_and_fixed_new_original_photo_pool":True,
        "fcp_global_species_prevalence_estimated":False,
        "non_detected_species_biologically_monomorphic":False,
        "no_old_149900_flower_colour_measurements":True,
        "no_literature_based_FCP_inclusion":True,
        "no_new_white_axis_H2_or_environment_claim":True,
        "species_result_csv_sha256":sha256(out_csv),
    }
    (outdir/"photo_rgfca_repetition_receipt.json").write_text(
        json.dumps(summary,ensure_ascii=False,sort_keys=True,indent=2)+"\n",
        encoding="utf-8")
    return summary


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cohort",type=Path,required=True)
    p.add_argument("--balanced-schedule",type=Path,required=True)
    p.add_argument("--schedule-receipt",type=Path,required=True)
    p.add_argument("--new-original-photo-measurements",type=Path,required=True)
    p.add_argument("--measured-photo-screen-receipt",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    a=p.parse_args(argv)
    registry,assignments=validate_schedule(a.cohort,a.balanced_schedule,a.schedule_receipt)
    photos=load_new_roi(a.new_original_photo_measurements,a.measured_photo_screen_receipt,registry)
    print("FCP_RGFCAPHOTO_REPEATED",json.dumps(
        repeat(registry,assignments,photos,a.outdir),sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
