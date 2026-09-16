#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import pandas as pd

from fcp_pipeline.random_photo_h9_pool import freeze_h9_metadata
from fcp_pipeline.random_photo_pool import InaturalistObservationClient
from fcp_pipeline.third_cohort_metadata_gate import (
    metadata_gate_verdict,
    seal_exact_denominator,
    validate_fresh_rows,
    validate_selected_manifest,
)

ROOT = Path(__file__).resolve().parents[2]
PROTOCOL = ROOT / "docs" / "POLYMORPHISM_H2_THIRD_COHORT_FRESH_METADATA_PROTOCOL_20260916.md"
SELECTION_RESULT = ROOT / "results" / "polymorphism_h2_third_cohort_selection_20260916" / "result.json"
SELECTION_MANIFEST = ROOT / "results" / "polymorphism_h2_third_cohort_selection_20260916" / "selected_species_manifest.tsv"
EXPECTED_SELECTION_MANIFEST_SHA = "16db6a2fc265f0ab20e7dae4e6f9058b907f5c19af3ddab5b6077797dd1def59"
EXPECTED_CANONICAL_SELECTED_SHA = "4ad1191f39068e0fb2229f84361b1004803e24793190566d21e5c474aef2002a"

EXCLUSION_SOURCES = [
    ROOT / "data" / "frozen" / "random_photo_first_h9_exclusion_ledger_v1.csv",
    ROOT / "data" / "frozen" / "random_photo_first_h9_fresh_metadata_v1.csv",
    ROOT / "results" / "rgfca_42111_breadth_measurement_step8f_20260911" / "rgfca_42111_species_breadth_measured.csv.gz",
    ROOT / "data" / "frozen" / "polymorphism_h2_p500_candidate_metadata_v1.csv.gz",
]
EXPECTED_EXCLUSION_OBS = 178_462
EXPECTED_EXCLUSION_PHOTO = 178_462

OUT = ROOT / "results" / "polymorphism_h2_third_cohort_candidate_metadata_20260916"
METADATA_OUT = ROOT / "data" / "frozen" / "polymorphism_h2_third_cohort_candidate_metadata_v1.csv.gz"
AUTHORIZED_ROWS_OUT = ROOT / "data" / "frozen" / "polymorphism_h2_third_cohort_authorized_metadata_v1.csv.gz"
AUDIT_OUT = OUT / "species_audit.csv"
AUTHORIZED_SPECIES_OUT = OUT / "authorized_species.csv"
RESULT_OUT = OUT / "result.json"

SELECTED_N = 500
TARGET_N = 100
MIN_FULL_SPECIES = 300
OBSERVER_CAP = 2
PER_PAGE = 200
MAX_ACCURACY_M = 5000
REQUEST_INTERVAL_SECONDS = 1.05
REQUEST_TIMEOUT_SECONDS = 45.0
REQUEST_RETRIES = 0
REQUEST_ERROR_CEILING = 0.05
ALLOWED_LICENSES = ("cc0", "cc-by", "cc-by-sa", "cc-by-nc", "cc-by-nc-sa")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def deterministic_gzip_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(
        path,
        index=False,
        lineterminator="\n",
        compression={"method": "gzip", "compresslevel": 9, "mtime": 0},
    )


def load_exclusions() -> tuple[set[int], set[int], list[dict[str, object]]]:
    obs: set[int] = set()
    photo: set[int] = set()
    audit: list[dict[str, object]] = []
    for path in EXCLUSION_SOURCES:
        if not path.exists():
            raise RuntimeError(f"missing exclusion source: {path}")
        header = pd.read_csv(path, nrows=0).columns.tolist()
        if not {"observation_id", "photo_id"}.issubset(header):
            raise RuntimeError(f"exclusion source lacks observation/photo ids: {path}")
        frame = pd.read_csv(path, usecols=["observation_id", "photo_id"], low_memory=False)
        o = pd.to_numeric(frame["observation_id"], errors="coerce").dropna().astype("int64")
        p = pd.to_numeric(frame["photo_id"], errors="coerce").dropna().astype("int64")
        obs.update(o.tolist())
        photo.update(p.tolist())
        audit.append(
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": sha256_file(path),
                "rows": int(len(frame)),
                "unique_observation_ids": int(o.nunique()),
                "unique_photo_ids": int(p.nunique()),
            }
        )
    if len(obs) != EXPECTED_EXCLUSION_OBS or len(photo) != EXPECTED_EXCLUSION_PHOTO:
        raise RuntimeError(
            f"exclusion union drift: obs={len(obs)} photo={len(photo)}; "
            f"expected={EXPECTED_EXCLUSION_OBS}/{EXPECTED_EXCLUSION_PHOTO}"
        )
    return obs, photo, audit


def main() -> int:
    if os.environ.get("GITHUB_ACTIONS") == "true" and os.environ.get("GITHUB_RUN_ATTEMPT", "1") != "1":
        raise RuntimeError("rerun forbidden: third-cohort fresh-metadata gate is one bounded draw")
    if not PROTOCOL.exists():
        raise RuntimeError("frozen metadata protocol is missing")
    if not SELECTION_RESULT.exists() or not SELECTION_MANIFEST.exists():
        raise RuntimeError("frozen third-cohort selection is missing")

    manifest_sha = sha256_file(SELECTION_MANIFEST)
    if manifest_sha != EXPECTED_SELECTION_MANIFEST_SHA:
        raise RuntimeError(f"third-cohort manifest drift: {manifest_sha}")
    selection_receipt = json.loads(SELECTION_RESULT.read_text(encoding="utf-8"))
    if selection_receipt.get("status") != "THIRD_COHORT_SELECTION_FROZEN_PREOUTCOME":
        raise RuntimeError("third-cohort selection receipt status drift")
    if selection_receipt.get("selection", {}).get("canonical_selected_csv_sha256") != EXPECTED_CANONICAL_SELECTED_SHA:
        raise RuntimeError("third-cohort canonical selected CSV fingerprint drift")
    if selection_receipt.get("biological_outcomes_opened") is not False:
        raise RuntimeError("selection receipt no longer certifies preoutcome state")

    selected = pd.read_csv(SELECTION_MANIFEST, sep="\t")
    validate_selected_manifest(selected, expected_n=SELECTED_N)
    selected["inat_taxon_id"] = pd.to_numeric(selected["inat_taxon_id"], errors="raise").astype(int)
    selected["species"] = selected["species"].astype(str)

    exclusion_obs, exclusion_photo, exclusion_audit = load_exclusions()

    client = InaturalistObservationClient(
        user_agent="zuizui0223-fcp-H2-third-cohort-prospective/1.0 (github.com/zuizui0223/fcp)",
        request_interval_seconds=REQUEST_INTERVAL_SECONDS,
        timeout_seconds=REQUEST_TIMEOUT_SECONDS,
        max_retries=REQUEST_RETRIES,
    )
    frozen = freeze_h9_metadata(
        client=client,
        species_frame=selected[["inat_taxon_id", "species"]].copy(),
        exclusion_observation_ids=exclusion_obs,
        exclusion_photo_ids=exclusion_photo,
        per_page=PER_PAGE,
        observer_cap_n=OBSERVER_CAP,
        fixed_raw_photos=TARGET_N,
        maximum_positional_accuracy_m=MAX_ACCURACY_M,
        allowed_photo_licenses=ALLOWED_LICENSES,
    )

    audit = frozen.species_audit.copy()
    if len(audit) != SELECTED_N or audit["inat_taxon_id"].nunique() != SELECTED_N:
        raise RuntimeError(f"candidate audit lost selected species: rows={len(audit)}")
    audit = audit.merge(
        selected[["prospective_rank", "inat_taxon_id"]],
        on="inat_taxon_id",
        how="left",
        validate="one_to_one",
    )
    audit = audit.sort_values("prospective_rank", kind="mergesort").reset_index(drop=True)

    observations = frozen.observations.copy()
    if len(observations):
        observations = observations.merge(
            selected[["prospective_rank", "inat_taxon_id"]],
            on="inat_taxon_id",
            how="left",
            validate="many_to_one",
        )
        observations = observations.sort_values(
            ["prospective_rank", "h9_selection_order", "observation_id", "photo_id"],
            kind="mergesort",
        ).reset_index(drop=True)
    else:
        observations = pd.DataFrame(
            columns=["query_species", "inat_taxon_id", "prospective_rank", "observation_id", "photo_id"]
        )

    validate_fresh_rows(observations, exclusion_obs, exclusion_photo)
    per_species = observations.groupby("inat_taxon_id", sort=False).size() if len(observations) else pd.Series(dtype=int)
    if len(per_species) and int(per_species.max()) > TARGET_N:
        raise RuntimeError("selected metadata exceeds frozen exact-100 target")

    authorized_species, authorized_rows = seal_exact_denominator(audit, observations, target_n=TARGET_N)
    n_full = int(len(authorized_species))
    if len(authorized_rows) != n_full * TARGET_N:
        raise RuntimeError("authorized row count is not exact species x 100 denominator")

    request_errors = int(audit["request_error"].fillna("").astype(str).str.len().gt(0).sum())
    verdict, pixels_may_be_authorized = metadata_gate_verdict(
        request_errors,
        SELECTED_N,
        n_full,
        min_full_species=MIN_FULL_SPECIES,
        error_ceiling=REQUEST_ERROR_CEILING,
    )
    error_fraction = request_errors / SELECTED_N

    OUT.mkdir(parents=True, exist_ok=True)
    deterministic_gzip_csv(observations, METADATA_OUT)
    deterministic_gzip_csv(authorized_rows, AUTHORIZED_ROWS_OUT)
    audit.to_csv(AUDIT_OUT, index=False, lineterminator="\n")
    authorized_species.to_csv(AUTHORIZED_SPECIES_OUT, index=False, lineterminator="\n")

    result = {
        "analysis": "polymorphism_h2_third_cohort_candidate_metadata_acquisition",
        "date_jst": "2026-09-16",
        "status": "complete_fresh_metadata_draw_before_pixels",
        "execution_head": os.environ.get("GITHUB_SHA", "local"),
        "protocol": str(PROTOCOL.relative_to(ROOT)),
        "selection": {
            "selected_species": SELECTED_N,
            "ordered_manifest": str(SELECTION_MANIFEST.relative_to(ROOT)),
            "ordered_manifest_sha256": manifest_sha,
            "canonical_selected_csv_sha256": EXPECTED_CANONICAL_SELECTED_SHA,
            "no_species_replacement": True,
        },
        "decision": {
            "verdict": verdict,
            "pixels_may_be_authorized_in_separate_next_step": pixels_may_be_authorized,
            "biological_pixel_opening_authorized_by_this_gate": False,
            "no_target_relaxation": True,
            "no_species_replacement": True,
            "one_bounded_metadata_draw": True,
        },
        "query": {
            "one_request_per_species": True,
            "per_page": PER_PAGE,
            "order_by": "random",
            "quality_grade": "research",
            "flowering_term_id": 12,
            "flowering_term_value_id": 13,
            "maximum_positional_accuracy_m": MAX_ACCURACY_M,
            "observer_cap": OBSERVER_CAP,
            "target_rows_per_species": TARGET_N,
            "request_interval_seconds": REQUEST_INTERVAL_SECONDS,
            "request_timeout_seconds": REQUEST_TIMEOUT_SECONDS,
            "request_retries": REQUEST_RETRIES,
            "allowed_photo_licenses": list(ALLOWED_LICENSES),
        },
        "transport": {
            "request_attempts": SELECTED_N,
            "request_errors": request_errors,
            "request_error_fraction": error_fraction,
            "request_error_fraction_ceiling": REQUEST_ERROR_CEILING,
        },
        "capacity": {
            "full100_species": n_full,
            "minimum_full100_species_before_separate_pixel_authorization": MIN_FULL_SPECIES,
            "selected_metadata_rows_all_selected_species": int(len(observations)),
            "authorized_species": n_full,
            "authorized_metadata_rows": int(len(authorized_rows)),
            "exact_rows_per_authorized_species": TARGET_N,
        },
        "prior_id_exclusion": {
            "unique_observation_ids": int(len(exclusion_obs)),
            "unique_photo_ids": int(len(exclusion_photo)),
            "sources": exclusion_audit,
        },
        "outcome_firewall": {
            "image_pixels_opened": False,
            "flower_colour_opened": False,
            "morph_opened": False,
            "palette_opened": False,
            "D_opened": False,
            "H2_vectors_opened": False,
            "H2_W_opened": False,
            "structured_null_opened": False,
            "H3_predictors_used": False,
            "P500_recovered_H2_read": False,
            "coordinates_used_only_for_sampling": True,
        },
        "files": {
            "metadata_all_selected": str(METADATA_OUT.relative_to(ROOT)),
            "metadata_all_selected_sha256": sha256_file(METADATA_OUT),
            "authorized_metadata": str(AUTHORIZED_ROWS_OUT.relative_to(ROOT)),
            "authorized_metadata_sha256": sha256_file(AUTHORIZED_ROWS_OUT),
            "species_audit": str(AUDIT_OUT.relative_to(ROOT)),
            "species_audit_sha256": sha256_file(AUDIT_OUT),
            "authorized_species": str(AUTHORIZED_SPECIES_OUT.relative_to(ROOT)),
            "authorized_species_sha256": sha256_file(AUTHORIZED_SPECIES_OUT),
            "protocol_sha256": sha256_file(PROTOCOL),
        },
        "hard_nonclaims": [
            "metadata capacity is not flower-colour polymorphism prevalence",
            "failure to reach 100 fresh rows is not monomorphism",
            "this gate does not test or estimate the white/non-white H2 axis",
            "metadata PASS does not itself authorize biological pixel opening",
        ],
    }
    RESULT_OUT.write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
