#!/usr/bin/env python3
"""One-shot, outcome-blind fresh-metadata gate for FCP v2.

This stage opens no image pixels and no flower-colour outcome. It queries the
already frozen P/N species queues exactly once, excludes every previously used
observation/photo ID supplied on the command line, and freezes the first 200
full-100 species per panel in queue order.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import pandas as pd

from fcp_pipeline.random_photo_h9_pool import freeze_h9_metadata
from fcp_pipeline.random_photo_pool import InaturalistObservationClient

QUEUE_PER_PANEL = 300
TARGET_SPECIES_PER_PANEL = 200
TARGET_ROWS_PER_SPECIES = 100
OBSERVER_CAP = 2
PER_PAGE = 200
MAX_ACCURACY_M = 5000
REQUEST_INTERVAL_SECONDS = 1.05
REQUEST_TIMEOUT_SECONDS = 45.0
REQUEST_RETRIES = 0
REQUEST_ERROR_CEILING = 0.05
EXPECTED_QUEUE_SHA = {
    "P": "6a3c7171988f0400a053ff181cdbfa9a04a3caecc9f6deaf88ee21a9404c3349",
    "N": "f2ee787bda61776ae52983c6b9d03c0232294974870cca7b33868f9bf050bff1",
}
ALLOWED_LICENSES = ("cc0", "cc-by", "cc-by-sa", "cc-by-nc", "cc-by-nc-sa")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--panel-p-queue", type=Path, required=True)
    p.add_argument("--panel-n-queue", type=Path, required=True)
    p.add_argument("--exclusion-source", type=Path, action="append", required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--metadata-out", type=Path, required=True)
    p.add_argument("--authorized-out", type=Path, required=True)
    return p.parse_args()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def canonical_csv_sha(frame: pd.DataFrame) -> str:
    payload = frame.to_csv(index=False, lineterminator="\n").encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def deterministic_gzip_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(
        path,
        index=False,
        lineterminator="\n",
        compression={"method": "gzip", "compresslevel": 9, "mtime": 0},
    )


def read_queue(path: Path, panel: str) -> pd.DataFrame:
    frame = pd.read_csv(path)
    required = {
        "queue_rank",
        "panel",
        "inat_taxon_id",
        "species",
        "after_observer_cap",
        "selection_hash",
    }
    missing = required - set(frame.columns)
    if missing:
        raise RuntimeError(f"panel {panel} queue missing fields: {sorted(missing)}")
    if len(frame) != QUEUE_PER_PANEL:
        raise RuntimeError(f"panel {panel} queue row drift: {len(frame)}")
    if not frame["panel"].astype(str).eq(panel).all():
        raise RuntimeError(f"panel {panel} label drift")
    ranks = pd.to_numeric(frame["queue_rank"], errors="raise").astype(int).tolist()
    if ranks != list(range(1, QUEUE_PER_PANEL + 1)):
        raise RuntimeError(f"panel {panel} queue rank drift")
    frame["inat_taxon_id"] = pd.to_numeric(frame["inat_taxon_id"], errors="raise").astype(int)
    frame["species"] = frame["species"].astype(str)
    if frame["inat_taxon_id"].nunique() != QUEUE_PER_PANEL:
        raise RuntimeError(f"panel {panel} taxon IDs are not unique")
    if frame["species"].nunique() != QUEUE_PER_PANEL:
        raise RuntimeError(f"panel {panel} species are not unique")
    observed_sha = canonical_csv_sha(frame)
    if observed_sha != EXPECTED_QUEUE_SHA[panel]:
        raise RuntimeError(
            f"panel {panel} queue SHA drift: {observed_sha} != {EXPECTED_QUEUE_SHA[panel]}"
        )
    return frame


def load_exclusions(paths: list[Path]) -> tuple[set[int], set[int], list[dict[str, object]]]:
    observation_ids: set[int] = set()
    photo_ids: set[int] = set()
    audit: list[dict[str, object]] = []
    for path in paths:
        if not path.is_file():
            raise RuntimeError(f"missing exclusion source: {path}")
        header = pd.read_csv(path, nrows=0).columns.tolist()
        if not {"observation_id", "photo_id"}.issubset(header):
            raise RuntimeError(f"exclusion source lacks observation/photo IDs: {path}")
        frame = pd.read_csv(path, usecols=["observation_id", "photo_id"], low_memory=False)
        obs = pd.to_numeric(frame["observation_id"], errors="coerce").dropna().astype("int64")
        pho = pd.to_numeric(frame["photo_id"], errors="coerce").dropna().astype("int64")
        observation_ids.update(obs.tolist())
        photo_ids.update(pho.tolist())
        audit.append(
            {
                "path": str(path),
                "sha256": sha256_file(path),
                "rows": int(len(frame)),
                "unique_observation_ids": int(obs.nunique()),
                "unique_photo_ids": int(pho.nunique()),
            }
        )
    return observation_ids, photo_ids, audit


def main() -> int:
    args = parse_args()
    if os.environ.get("GITHUB_ACTIONS") == "true":
        if os.environ.get("GITHUB_RUN_ATTEMPT", "1") != "1":
            raise RuntimeError("rerun forbidden: FCP v2 fresh-metadata gate is one bounded draw")

    panel_p = read_queue(args.panel_p_queue, "P")
    panel_n = read_queue(args.panel_n_queue, "N")
    if set(panel_p["inat_taxon_id"]) & set(panel_n["inat_taxon_id"]):
        raise RuntimeError("P/N queue taxon overlap")
    if set(panel_p["species"]) & set(panel_n["species"]):
        raise RuntimeError("P/N queue species overlap")

    queue = pd.concat([panel_p, panel_n], ignore_index=True)
    exclusion_obs, exclusion_photo, exclusion_audit = load_exclusions(args.exclusion_source)

    client = InaturalistObservationClient(
        user_agent="zuizui0223-fcp-v2-measurement-validity/1.0 (github.com/zuizui0223/fcp)",
        request_interval_seconds=REQUEST_INTERVAL_SECONDS,
        timeout_seconds=REQUEST_TIMEOUT_SECONDS,
        max_retries=REQUEST_RETRIES,
    )
    frozen = freeze_h9_metadata(
        client=client,
        species_frame=queue[["inat_taxon_id", "species"]].copy(),
        exclusion_observation_ids=exclusion_obs,
        exclusion_photo_ids=exclusion_photo,
        per_page=PER_PAGE,
        observer_cap_n=OBSERVER_CAP,
        fixed_raw_photos=TARGET_ROWS_PER_SPECIES,
        maximum_positional_accuracy_m=MAX_ACCURACY_M,
        allowed_photo_licenses=ALLOWED_LICENSES,
    )

    audit = frozen.species_audit.copy()
    audit = audit.merge(
        queue[["panel", "queue_rank", "inat_taxon_id", "species", "selection_hash"]],
        on=["inat_taxon_id", "species"],
        how="left",
        validate="one_to_one",
    )
    if audit["panel"].isna().any() or len(audit) != 2 * QUEUE_PER_PANEL:
        raise RuntimeError("fresh metadata audit lost frozen queue species")
    audit = audit.sort_values(["panel", "queue_rank"], kind="mergesort").reset_index(drop=True)

    observations = frozen.observations.copy()
    if len(observations):
        observations = observations.merge(
            queue[["panel", "queue_rank", "inat_taxon_id", "species", "selection_hash"]],
            on=["inat_taxon_id", "species"],
            how="left",
            validate="many_to_one",
        )
        observations = observations.sort_values(
            ["panel", "queue_rank", "h9_selection_order", "observation_id", "photo_id"],
            kind="mergesort",
        ).reset_index(drop=True)
        obs_ids = pd.to_numeric(observations["observation_id"], errors="raise").astype(int)
        pho_ids = pd.to_numeric(observations["photo_id"], errors="raise").astype(int)
        if obs_ids.duplicated().any() or pho_ids.duplicated().any():
            raise RuntimeError("fresh v2 metadata contains duplicate observation/photo IDs")
        if set(obs_ids) & exclusion_obs or set(pho_ids) & exclusion_photo:
            raise RuntimeError("previously used observation/photo ID leaked into v2 metadata")

    request_errors = int(audit["request_error"].fillna("").astype(str).str.len().gt(0).sum())
    error_fraction = request_errors / (2 * QUEUE_PER_PANEL)

    terminal_frames: list[pd.DataFrame] = []
    panel_summary: dict[str, dict[str, object]] = {}
    panel_pass = True
    for panel in ("P", "N"):
        a = audit.loc[audit["panel"].eq(panel)].sort_values("queue_rank", kind="mergesort")
        full = a.loc[a["full_fixed_n"].astype(bool)].head(TARGET_SPECIES_PER_PANEL).copy()
        full_count_in_queue = int(a["full_fixed_n"].astype(bool).sum())
        panel_ok = len(full) == TARGET_SPECIES_PER_PANEL
        panel_pass = panel_pass and panel_ok
        if panel_ok:
            terminal_frames.append(full)
        panel_summary[panel] = {
            "queue_species": QUEUE_PER_PANEL,
            "full100_species_in_queue": full_count_in_queue,
            "terminal_target": TARGET_SPECIES_PER_PANEL,
            "terminal_species": int(len(full)) if panel_ok else 0,
            "last_terminal_queue_rank": (
                int(full["queue_rank"].max()) if panel_ok else None
            ),
            "pass": bool(panel_ok),
        }

    transport_pass = error_fraction <= REQUEST_ERROR_CEILING
    gate_pass = bool(panel_pass and transport_pass)

    if gate_pass:
        terminal_species = pd.concat(terminal_frames, ignore_index=True)
        terminal_species = terminal_species.sort_values(
            ["panel", "queue_rank"], kind="mergesort"
        ).reset_index(drop=True)
        ids = set(terminal_species["inat_taxon_id"].astype(int))
        authorized = observations.loc[
            observations["inat_taxon_id"].astype(int).isin(ids)
        ].copy()
        authorized = authorized.sort_values(
            ["panel", "queue_rank", "h9_selection_order", "observation_id", "photo_id"],
            kind="mergesort",
        ).reset_index(drop=True)
        counts = authorized.groupby(["panel", "inat_taxon_id"], sort=False).size()
        if len(terminal_species) != 2 * TARGET_SPECIES_PER_PANEL:
            raise RuntimeError("terminal species denominator drift")
        if len(authorized) != 2 * TARGET_SPECIES_PER_PANEL * TARGET_ROWS_PER_SPECIES:
            raise RuntimeError("authorized v2 row denominator drift")
        if len(counts) != 2 * TARGET_SPECIES_PER_PANEL or not counts.eq(
            TARGET_ROWS_PER_SPECIES
        ).all():
            raise RuntimeError("authorized v2 species do not all have exact 100 rows")
        verdict = "FCP_V2_FRESH_METADATA_GATE_PASS"
    else:
        terminal_species = pd.DataFrame(
            columns=["panel", "queue_rank", "inat_taxon_id", "species", "selection_hash"]
        )
        authorized = pd.DataFrame(columns=observations.columns)
        verdict = (
            "FCP_V2_FRESH_METADATA_TRANSPORT_NOT_EVALUABLE"
            if not transport_pass
            else "FCP_V2_FRESH_METADATA_CAPACITY_NOT_EVALUABLE"
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.metadata_out.parent.mkdir(parents=True, exist_ok=True)
    args.authorized_out.parent.mkdir(parents=True, exist_ok=True)

    deterministic_gzip_csv(observations, args.metadata_out)
    deterministic_gzip_csv(authorized, args.authorized_out)
    audit.to_csv(args.output_dir / "species_audit.csv", index=False, lineterminator="\n")
    terminal_species.to_csv(
        args.output_dir / "authorized_species.csv", index=False, lineterminator="\n"
    )

    result = {
        "schema": "fcp_v2_fresh_metadata_gate_v1",
        "date_jst": "2026-09-23",
        "status": "COMPLETE_OUTCOME_BLIND_FRESH_METADATA_DRAW",
        "decision": {
            "verdict": verdict,
            "pixels_may_be_authorized_in_separate_next_step": gate_pass,
            "biological_pixel_opening_authorized_by_this_gate": False,
            "one_bounded_metadata_draw": True,
            "queue_extension_permitted": False,
            "species_replacement_after_pixel_opening": False,
        },
        "queues": {
            "P_sha256": EXPECTED_QUEUE_SHA["P"],
            "N_sha256": EXPECTED_QUEUE_SHA["N"],
            "queue_species_per_panel": QUEUE_PER_PANEL,
            "terminal_target_per_panel": TARGET_SPECIES_PER_PANEL,
        },
        "query": {
            "request_attempts": 2 * QUEUE_PER_PANEL,
            "request_errors": request_errors,
            "request_error_fraction": error_fraction,
            "request_error_fraction_ceiling": REQUEST_ERROR_CEILING,
            "per_page": PER_PAGE,
            "observer_cap": OBSERVER_CAP,
            "target_rows_per_terminal_species": TARGET_ROWS_PER_SPECIES,
            "maximum_positional_accuracy_m": MAX_ACCURACY_M,
            "request_retries": REQUEST_RETRIES,
            "order_by": "random",
            "quality_grade": "research",
            "flowering_term_id": 12,
            "flowering_term_value_id": 13,
            "allowed_photo_licenses": list(ALLOWED_LICENSES),
        },
        "panels": panel_summary,
        "counts": {
            "all_fresh_metadata_rows": int(len(observations)),
            "authorized_species": int(len(terminal_species)),
            "authorized_metadata_rows": int(len(authorized)),
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
            "q_white_projection_opened": False,
            "H2_W_opened": False,
            "spatial_outcome_opened": False,
            "H3_predictors_used": False,
        },
        "files": {
            "metadata_all_queue_species": str(args.metadata_out),
            "metadata_all_queue_species_sha256": sha256_file(args.metadata_out),
            "authorized_metadata": str(args.authorized_out),
            "authorized_metadata_sha256": sha256_file(args.authorized_out),
            "species_audit": str(args.output_dir / "species_audit.csv"),
            "species_audit_sha256": sha256_file(args.output_dir / "species_audit.csv"),
            "authorized_species": str(args.output_dir / "authorized_species.csv"),
            "authorized_species_sha256": sha256_file(
                args.output_dir / "authorized_species.csv"
            ),
        },
        "hard_nonclaims": [
            "metadata capacity is not flower-colour biology",
            "failure to provide 100 fresh rows is not monomorphism",
            "this stage does not test D, H2, spatial organization or H3",
            "metadata PASS does not itself authorize pixel opening",
        ],
    }
    (args.output_dir / "result.json").write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
