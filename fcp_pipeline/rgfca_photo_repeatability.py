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
            white_yellow = None
            if classified:
                v = [float(raw[f"flower_fraction_{p}"]) for p in PALETTE]
                if any(not math.isfinite(x) or x < 0 or x > 1 for x in v) or abs(sum(v)-1)>1e-7:
                    raise ValueError("Original nine-colour flower ROI fractions invalid")
                white_yellow = (v[0], v[1])
            groups[ident].append({
                "photo_id":pid, "observer_id":observer, "observation_id":obs,
                "fold":int(fold), "morph":raw["morph"] if classified else None,
                "classifiable":classified, "within_ROI_white_yellow":white_yellow,
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
             rep: int, *, per_fold: int = 60) -> dict:
    if per_fold not in (20, 40, 48, 60, 64):
        raise ValueError("RGFCA photo draw must be one frozen source-only depth 40/80/96/120/128")
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
        # PHOTO IDs are never exported to the biological result CSV.
        # Kept only within the current species to quantify actual overlap
        # among same-source Monte Carlo resamples, without assuming draws
        # are independent biological replicates.
        "_selected_source_photo_IDs":[row["photo_id"] for row in selection],
        "fold_0_classifiable":per_fold_good[0],
        "fold_1_classifiable":per_fold_good[1],
        "evaluable":bool(evaluable), "photo_poly_detected":bool(detected),
        "coarse_morph_status":result["photographic_colour_status"]
    }


def _photo_resample_overlap(series: list[dict]) -> dict:
    """Actual all-pairs intersection / mean draw size, O(R × nphotos).

    For R repeated source-photo draws, photo i appears f_i times.
    The sum of pairwise intersections is Σ_i choose(f_i,2).
    This is a descriptive source/photo sampling statistic, not an
    effective number of independent flowering plants.
    """
    n_draws=len(series)
    photo_use=Counter(
        photo_id for trial in series
        for photo_id in trial["_selected_source_photo_IDs"])
    mean_n=(sum(len(t["_selected_source_photo_IDs"]) for t in series)/n_draws
            if n_draws else 0.0)
    if n_draws<2 or mean_n<=0:
        overlap=None
    else:
        joint=sum(n*(n-1)/2 for n in photo_use.values())
        overlap=(joint/(n_draws*(n_draws-1)/2))/mean_n
    distinct=len({tuple(sorted(x["_selected_source_photo_IDs"])) for x in series})
    if overlap is not None and not -1e-12<=overlap<=1+1e-12:
        raise ValueError("Repeated photographic source ID overlap invalid")
    return {
        "observed_mean_pairwise_source_photo_overlap_share":
            "" if overlap is None else round(overlap,12),
        "distinct_original_source_photo_ID_sets":distinct,
        "mean_selected_original_photo_IDs_per_realization":round(mean_n,12),
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
        depth_trials = {20: [], 40: [], 48: [], 60: [], 64: []}
        for rep in grouped[taxon]:
            trial = one_draw(taxon,item["species"],int(item["source_photo_opportunities"]),
                             rows,rep,per_fold=60)
            trials.append((rep, trial))
            depth_trials[60].append(trial)
            # Frozen source-only draw-depth sensitivity. R=200 species
            # schedules and source-photo opportunities are unchanged.
            for fold_n in (20,40,48,64):
                depth_trials[fold_n].append(
                    one_draw(taxon,item["species"],int(item["source_photo_opportunities"]),
                             rows,rep,per_fold=fold_n))
            for cutoff in wave_counts:
                if rep<=cutoff:
                    wave_counts[cutoff]["inclusions"]+=1
                    wave_counts[cutoff]["evaluable"]+=int(trial["evaluable"])
                    wave_counts[cutoff]["positive"]+=int(trial["photo_poly_detected"])
        eligible = sum(v["evaluable"] for _,v in trials)
        positives = sum(v["photo_poly_detected"] for _,v in trials)
        reliability = positives/eligible if eligible else None
        # The distinct estimands MUST be shown side by side:
        # P_eval=eligible/inclusions; P_poly|eval=positive/eligible;
        # P_joint=positive/inclusions. An impressive conditional rate can
        # otherwise conceal an almost entirely unevaluable source pool.
        evaluation_rate = eligible/len(trials) if trials else None
        unconditional_detection_rate = positives/len(trials) if trials else None
        if (trials and reliability is not None and abs(
                (evaluation_rate*reliability)-unconditional_detection_rate) > 1e-12):
            raise ValueError("FCP resampling outcome denominators inconsistent")
        classifiable_pool = [r for r in rows if r["classifiable"]]
        within_ROI_dual = sum(
            pair is not None and pair[0] >= .10 and pair[1] >= .10
            for pair in (r.get("within_ROI_white_yellow") for r in classifiable_pool))
        dual_fraction = (within_ROI_dual / len(classifiable_pool)
                         if classifiable_pool else 0.0)
        morphs_in_pool = {r["morph"] for r in classifiable_pool}
        # TECHNICAL calibration discovered AFTER the isolated six-tax member
        # pilot, BEFORE the future 1000-taxon complete original-image sample.
        # Photo-hash stability never proves separate biological morphs when
        # one flower's ROI routinely combines e.g. white rays + yellow discs.
        structural_hold = bool(
            len(classifiable_pool) >= 40
            and {"white", "yellow_orange"}.issubset(morphs_in_pool)
            and dual_fraction >= .50)
        depths = {
            str(fold_n*2): {
                "evaluable":sum(v["evaluable"] for v in series),
                "positive":sum(v["photo_poly_detected"] for v in series),
                **_photo_resample_overlap(series),
            }
            for fold_n, series in depth_trials.items()
        }
        primary_overlap = depths["120"]["observed_mean_pairwise_source_photo_overlap_share"]
        # After both six-taxon real-photo pilot waves and BEFORE production
        # original-photo RGFCA outcome measurement, require robust phenotype
        # photo-label detection also at a meaningfully smaller 48+48 draw.
        # In a full 64+64 pool this lower-depth draw shares ~75% of IDs on
        # average across Monte Carlo rounds vs ~94% for 60+60.
        # If every lower-depth draw sees the same photos, the source cannot
        # establish sampling robustness however often the label recurs.
        depth96 = depths["96"]
        eval96 = depth96["evaluable"]
        positive96 = depth96["positive"]
        conditional96 = (positive96/eval96 if eval96 else None)
        overlap96 = depth96["observed_mean_pairwise_source_photo_overlap_share"]
        lower_depth_replicated = bool(
            eval96 >= 30 and conditional96 is not None
            and conditional96 >= 0.80
            and overlap96 != "" and overlap96 <= 0.80
            and depth96["distinct_original_source_photo_ID_sets"] >= 30)
        # Conditional P_poly can approach 1 simply because two source-photo
        # draws reuse >90% of identical images. Report it but never claim
        # independent biological robustness when the effective perturbation
        # was too small.
        high_primary_photo_ID_overlap = (
            primary_overlap != "" and primary_overlap >= 0.90)
        if not rows or eligible < 30:
            tier = "not_evaluable_or_insufficient_repeated_roi_support"
        elif reliability is not None and reliability >= .8 and structural_hold:
            tier = "repeatable_colour_LABELS_WITHIN_ROI_STRUCTURAL_QC_HOLD"
        elif reliability is not None and reliability >= .8 and lower_depth_replicated:
            tier = "repeatably_photo_polymorphic_phenotype"
        elif reliability is not None and reliability >= .8:
            tier = "repeated_photo_labels_but_96photo_perturbation_unresolved_NOT_robust"
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
            "P_evaluable_under_random_photo_draw":"" if evaluation_rate is None else round(evaluation_rate,12),
            "P_detect_both_evaluable_and_multicolour":"" if unconditional_detection_rate is None else round(unconditional_detection_rate,12),
            "P_poly_primary_photo_draw":120,
            "actual_primary_120photo_pairwise_ID_overlap_share":primary_overlap,
            "actual_primary_120photo_distinct_ID_draw_sets":
                depths["120"]["distinct_original_source_photo_ID_sets"],
            "P_poly_96photo_conditional_on_evaluable":"" if conditional96 is None else round(conditional96,12),
            "n_evaluable_independent_observer_photo_redraws_96":eval96,
            "actual_96photo_pairwise_photo_ID_overlap_share":overlap96,
            "source_photo_96depth_provides_meaningful_MC_perturbation":
                lower_depth_replicated,
            "primary_MC_resamples_share_ge90percent_original_photo_IDs":
                high_primary_photo_ID_overlap,
            "draw_depth_sensitivity_counts":json.dumps(depths,sort_keys=True),
            "repeatability_0_80_status":tier,
            "n_classifiable_pool_photos_with_both_white_yellow_ge10":within_ROI_dual,
            "white_yellow_within_single_ROI_fraction":round(dual_fraction,12),
            "structural_organ_colour_review_hold":structural_hold,
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
        "primary_120_photo_draw_and_60_photo_per_observer_fold":True,
        "primary_draw_changed_pre_image_due_96photo_pilot_classifiability":True,
        "sensitivity_photo_draw_sizes":[40,80,96,120,128],
        "post_technical_pilot_96photo_perturbation_qualification":True,
        "robust_photo_label_call_requires_ge30_evaluable_depth96_redraws_and_P_ge080":True,
        "robust_photo_label_call_requires_depth96_photo_ID_overlap_le080":True,
        "source_photo_acquisition_max_per_species":128,
        "probability_estimated_only_conditional_on_classified_original_photos":True,
        "P_eval_evaluable_over_all_balanced_draws_separately_reported":True,
        "P_joint_positive_over_all_balanced_draws_separately_reported":True,
        "P_poly_given_eval_must_not_be_called_overall_detection_probability":True,
        "original_photo_ID_overlap_measured_empirically_each_draw_depth":True,
        "high_photo_ID_reuse_cannot_support_independent_biological_replication":True,
        "n_phenotype_screen_taxa_primary_repeats_share_ge90percent_photo_IDs":sum(
            r["primary_MC_resamples_share_ge90percent_original_photo_IDs"]
            for r in species_rows),
        "n_repeated_phenotypic_polymorphism_species":tiers.get(
            "repeatably_photo_polymorphic_phenotype",0),
        "n_coarse_photo_label_recurring_but_structural_QC_HOLD":tiers.get(
            "repeatable_colour_LABELS_WITHIN_ROI_STRUCTURAL_QC_HOLD",0),
        "post_pilot_intra_ROI_white_yellow_artifact_QC_calibration":True,
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
