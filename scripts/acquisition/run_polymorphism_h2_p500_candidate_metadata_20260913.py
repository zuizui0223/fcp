#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import pandas as pd

from fcp_pipeline.random_photo_h9_pool import freeze_h9_metadata
from fcp_pipeline.random_photo_pool import InaturalistObservationClient

ROOT = Path(__file__).resolve().parents[2]
PROTOCOL = ROOT / "docs" / "POLYMORPHISM_H2_P500_CANDIDATE_METADATA_PROTOCOL_20260913.md"
SELECTION_SCRIPT = ROOT / "scripts" / "analysis" / "select_polymorphism_h2_prospective_u100_20260913.py"
SELECTION_OUT = ROOT / "results" / "polymorphism_h2_prospective_u100_selection_20260913"
P500 = SELECTION_OUT / "p500_frozen_selection.csv"
EXPECTED_P500_SHA = "f54d07fb2338a20a3f92808b58f0ebc948ab0d5af3cb01fb26702f861184b7a4"

EXCLUSION_SOURCES = [
    ROOT / "data" / "frozen" / "random_photo_first_h9_exclusion_ledger_v1.csv",
    ROOT / "data" / "frozen" / "random_photo_first_h9_fresh_metadata_v1.csv",
    ROOT / "results" / "rgfca_42111_breadth_measurement_step8f_20260911" / "rgfca_42111_species_breadth_measured.csv.gz",
]
EXPECTED_EXCLUSION_OBS = 128_464
EXPECTED_EXCLUSION_PHOTO = 128_464

OUT = ROOT / "results" / "polymorphism_h2_p500_candidate_metadata_20260913"
METADATA_OUT = ROOT / "data" / "frozen" / "polymorphism_h2_p500_candidate_metadata_v1.csv.gz"
AUDIT_OUT = OUT / "species_audit.csv"
FULL100_OUT = OUT / "full100_species.csv"
RESULT_OUT = OUT / "result.json"

P500_N = 500
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


def deterministic_gzip_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(
        path,
        index=False,
        lineterminator="\n",
        compression={"method": "gzip", "compresslevel": 9, "mtime": 0},
    )


def main() -> int:
    if not PROTOCOL.exists():
        raise RuntimeError("frozen protocol is missing")

    # Recreate the outcome-blind P500 and require exact byte identity before requests.
    subprocess.run(["python", str(SELECTION_SCRIPT)], cwd=ROOT, check=True)
    p500_sha = sha256_file(P500)
    if p500_sha != EXPECTED_P500_SHA:
        raise RuntimeError(f"P500 selection drift: {p500_sha}")
    p500 = pd.read_csv(P500, usecols=["inat_taxon_id", "species", "prospective_rank", "selection_hash"])
    if len(p500) != P500_N or p500["inat_taxon_id"].nunique() != P500_N or p500["species"].nunique() != P500_N:
        raise RuntimeError("P500 identity fingerprint drift")
    p500["inat_taxon_id"] = pd.to_numeric(p500["inat_taxon_id"], errors="raise").astype(int)
    p500["species"] = p500["species"].astype(str)

    exclusion_obs, exclusion_photo, exclusion_audit = load_exclusions()

    client = InaturalistObservationClient(
        user_agent="zuizui0223-fcp-H2-P500-prospective/1.0 (github.com/zuizui0223/fcp)",
        request_interval_seconds=REQUEST_INTERVAL_SECONDS,
        timeout_seconds=REQUEST_TIMEOUT_SECONDS,
        max_retries=REQUEST_RETRIES,
    )

    frozen = freeze_h9_metadata(
        client=client,
        species_frame=p500[["inat_taxon_id", "species"]].copy(),
        exclusion_observation_ids=exclusion_obs,
        exclusion_photo_ids=exclusion_photo,
        per_page=PER_PAGE,
        observer_cap_n=OBSERVER_CAP,
        fixed_raw_photos=TARGET_N,
        maximum_positional_accuracy_m=MAX_ACCURACY_M,
        allowed_photo_licenses=ALLOWED_LICENSES,
    )

    audit = frozen.species_audit.copy()
    if len(audit) != P500_N or audit["inat_taxon_id"].nunique() != P500_N:
        raise RuntimeError(f"candidate audit lost P500 species: rows={len(audit)}")
    audit = audit.merge(
        p500[["inat_taxon_id", "prospective_rank", "selection_hash"]],
        on="inat_taxon_id",
        how="left",
        validate="one_to_one",
    )
    audit = audit.sort_values("prospective_rank", kind="mergesort").reset_index(drop=True)

    full100 = audit.loc[audit["full_fixed_n"].astype(bool), [
        "prospective_rank", "inat_taxon_id", "species", "after_observer_cap", "retained", "request_error"
    ]].copy()
    full100 = full100.sort_values("prospective_rank", kind="mergesort").reset_index(drop=True)
    full_species = set(full100["species"].astype(str))

    observations = frozen.observations.copy()
    if len(observations):
        observations["full_target_species"] = observations["query_species"].astype(str).isin(full_species)
        observations = observations.merge(
            p500[["inat_taxon_id", "prospective_rank", "selection_hash"]],
            on="inat_taxon_id",
            how="left",
            validate="many_to_one",
        )
        sort_cols = [c for c in ["prospective_rank", "h9_selection_order", "observation_id", "photo_id"] if c in observations.columns]
        observations = observations.sort_values(sort_cols, kind="mergesort").reset_index(drop=True)
    else:
        observations = pd.DataFrame(columns=["query_species", "inat_taxon_id", "full_target_species"])

    # Hard integrity: selected rows are metadata only and never exceed 100 per query species.
    if len(observations):
        per_species = observations.groupby("query_species", sort=False).size()
        if int(per_species.max()) > TARGET_N:
            raise RuntimeError("selected metadata exceeds frozen target")
        if observations["observation_id"].astype(str).duplicated().any():
            raise RuntimeError("duplicate observation IDs in frozen P500 metadata")
        if observations["photo_id"].astype(str).duplicated().any():
            raise RuntimeError("duplicate photo IDs in frozen P500 metadata")
        used_obs = set(pd.to_numeric(observations["observation_id"], errors="raise").astype(int))
        used_photo = set(pd.to_numeric(observations["photo_id"], errors="raise").astype(int))
        if used_obs & exclusion_obs or used_photo & exclusion_photo:
            raise RuntimeError("prior experiment ID leaked into fresh P500 metadata")

    request_errors = int(audit["request_error"].fillna("").astype(str).str.len().gt(0).sum())
    error_fraction = request_errors / P500_N
    n100 = int(len(full100))

    if error_fraction > REQUEST_ERROR_CEILING:
        verdict = "P500_CANDIDATE_METADATA_TRANSPORT_NOT_EVALUABLE"
        pixels_authorized_next = False
    elif n100 < MIN_FULL_SPECIES:
        verdict = "P500_CANDIDATE_METADATA_CAPACITY_NOT_EVALUABLE"
        pixels_authorized_next = False
    else:
        verdict = "P500_CANDIDATE_METADATA_GATE_PASS"
        pixels_authorized_next = True

    OUT.mkdir(parents=True, exist_ok=True)
    deterministic_gzip_csv(observations, METADATA_OUT)
    audit.to_csv(AUDIT_OUT, index=False, lineterminator="\n")
    full100.to_csv(FULL100_OUT, index=False, lineterminator="\n")

    result = {
        "analysis": "polymorphism_h2_p500_candidate_metadata_acquisition",
        "date_jst": "2026-09-13",
        "status": "complete_fresh_metadata_draw_before_pixels",
        "protocol": str(PROTOCOL.relative_to(ROOT)),
        "decision": {
            "verdict": verdict,
            "pixels_may_be_authorized_in_separate_next_step": pixels_authorized_next,
            "no_species_replacement": True,
            "no_target_relaxation": True,
        },
        "P500": {
            "species": P500_N,
            "sha256": p500_sha,
            "legacy_species_reintroduced": False,
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
            "request_attempts": P500_N,
            "request_errors": request_errors,
            "request_error_fraction": error_fraction,
            "request_error_fraction_ceiling": REQUEST_ERROR_CEILING,
        },
        "capacity": {
            "full100_species": n100,
            "minimum_full100_species_before_pixels": MIN_FULL_SPECIES,
            "selected_metadata_rows_all_P500": int(len(observations)),
            "selected_metadata_rows_full100_species": int(observations["full_target_species"].sum()) if len(observations) else 0,
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
            "H2_W_opened": False,
            "H3_predictors_used": False,
            "coordinates_used_only_for_sampling": True,
        },
        "files": {
            "metadata": str(METADATA_OUT.relative_to(ROOT)),
            "metadata_sha256": sha256_file(METADATA_OUT),
            "species_audit": str(AUDIT_OUT.relative_to(ROOT)),
            "species_audit_sha256": sha256_file(AUDIT_OUT),
            "full100_species": str(FULL100_OUT.relative_to(ROOT)),
            "full100_species_sha256": sha256_file(FULL100_OUT),
            "protocol_sha256": sha256_file(PROTOCOL),
        },
        "hard_nonclaims": [
            "metadata capacity is not flower-colour polymorphism prevalence",
            "failure to reach 100 fresh rows is not monomorphism",
            "this step does not test the white/non-white axis",
            "this step does not open image pixels or colour outcomes",
        ],
    }
    RESULT_OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
