#!/usr/bin/env python3
"""One-shot, metadata-only terminalization of FCP v2 Panel P/N.

This stage opens no image pixels and reads no flower-colour outcomes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
from typing import Iterable

import pandas as pd

from fcp_pipeline.random_photo_h9_pool import freeze_h9_metadata
from fcp_pipeline.random_photo_pool import InaturalistObservationClient

EXPECTED_QUEUE_N = 300
TERMINAL_SPECIES_PER_PANEL = 200
ROWS_PER_SPECIES = 100
EXPECTED_PRIOR_IDS = 228_362
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
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def deterministic_gzip_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(
        path,
        index=False,
        lineterminator="\n",
        compression={"method": "gzip", "compresslevel": 9, "mtime": 0},
    )


def validate_queue(frame: pd.DataFrame, panel: str) -> pd.DataFrame:
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
        raise RuntimeError(f"Panel {panel} queue missing columns: {sorted(missing)}")
    x = frame.loc[:, sorted(required)].copy()
    if len(x) != EXPECTED_QUEUE_N:
        raise RuntimeError(f"Panel {panel} queue denominator drift: {len(x)}")
    if set(x["panel"].astype(str)) != {panel}:
        raise RuntimeError(f"Panel {panel} queue panel label drift")
    ranks = pd.to_numeric(x["queue_rank"], errors="raise").astype(int)
    if ranks.tolist() != list(range(1, EXPECTED_QUEUE_N + 1)):
        raise RuntimeError(f"Panel {panel} queue rank drift")
    taxon = pd.to_numeric(x["inat_taxon_id"], errors="raise").astype(int)
    species = x["species"].astype(str).str.strip()
    if taxon.nunique() != EXPECTED_QUEUE_N or species.nunique() != EXPECTED_QUEUE_N:
        raise RuntimeError(f"Panel {panel} queue identity duplication")
    if species.eq("").any():
        raise RuntimeError(f"Panel {panel} queue has blank species")
    x["queue_rank"] = ranks
    x["inat_taxon_id"] = taxon
    x["species"] = species
    return x.sort_values("queue_rank", kind="mergesort").reset_index(drop=True)


def load_exclusion_union(paths: Iterable[Path]) -> tuple[set[int], set[int], list[dict[str, object]]]:
    obs: set[int] = set()
    photo: set[int] = set()
    audit: list[dict[str, object]] = []
    for path in paths:
        if not path.is_file():
            raise RuntimeError(f"missing exclusion source: {path}")
        header = pd.read_csv(path, nrows=0).columns.tolist()
        if not {"observation_id", "photo_id"}.issubset(header):
            raise RuntimeError(f"exclusion source lacks IDs: {path}")
        frame = pd.read_csv(path, usecols=["observation_id", "photo_id"], low_memory=False)
        o = pd.to_numeric(frame["observation_id"], errors="coerce").dropna().astype("int64")
        p = pd.to_numeric(frame["photo_id"], errors="coerce").dropna().astype("int64")
        obs.update(o.tolist())
        photo.update(p.tolist())
        audit.append(
            {
                "path": str(path),
                "sha256": sha256_file(path),
                "rows": int(len(frame)),
                "unique_observation_ids": int(o.nunique()),
                "unique_photo_ids": int(p.nunique()),
            }
        )
    if len(obs) != EXPECTED_PRIOR_IDS or len(photo) != EXPECTED_PRIOR_IDS:
        raise RuntimeError(
            f"prior-ID union drift: obs={len(obs)} photo={len(photo)} "
            f"expected={EXPECTED_PRIOR_IDS}"
        )
    return obs, photo, audit


def terminalize_panel(
    *,
    panel: str,
    queue: pd.DataFrame,
    species_audit: pd.DataFrame,
    observations: pd.DataFrame,
    target_species: int = TERMINAL_SPECIES_PER_PANEL,
    target_rows: int = ROWS_PER_SPECIES,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    audit = species_audit.merge(
        queue[["queue_rank", "inat_taxon_id", "species"]],
        on=["inat_taxon_id", "species"],
        how="left",
        validate="one_to_one",
    )
    if len(audit) != len(queue) or audit["queue_rank"].isna().any():
        raise RuntimeError(f"Panel {panel} audit lost frozen queue identities")
    audit = audit.sort_values("queue_rank", kind="mergesort").reset_index(drop=True)
    full = audit.loc[audit["full_fixed_n"].astype(bool)].copy()
    full = full.sort_values("queue_rank", kind="mergesort").reset_index(drop=True)

    if len(full) < int(target_species):
        terminal_species = full.copy()
        terminal_species["terminal_panel_rank"] = range(1, len(full) + 1)
        terminal_rows = observations.iloc[0:0].copy()
        return audit, terminal_species, terminal_rows

    terminal_species = full.iloc[: int(target_species)].copy()
    terminal_species["terminal_panel_rank"] = range(1, int(target_species) + 1)
    terminal_species.insert(0, "panel", panel)

    terminal_ids = set(
        pd.to_numeric(terminal_species["inat_taxon_id"], errors="raise").astype(int)
    )
    if observations.empty:
        raise RuntimeError(f"Panel {panel} terminal species exist but observations are empty")

    rows = observations.copy()
    rows["inat_taxon_id"] = pd.to_numeric(rows["inat_taxon_id"], errors="raise").astype(int)
    rows = rows.loc[rows["inat_taxon_id"].isin(terminal_ids)].copy()
    rows = rows.merge(
        terminal_species[["queue_rank", "terminal_panel_rank", "inat_taxon_id"]],
        on="inat_taxon_id",
        how="left",
        validate="many_to_one",
    )
    counts = rows.groupby("inat_taxon_id", sort=False).size()
    if len(counts) != int(target_species) or not counts.eq(int(target_rows)).all():
        raise RuntimeError(f"Panel {panel} terminal row denominator is not exact")
    rows.insert(0, "panel", panel)
    rows = rows.sort_values(
        ["terminal_panel_rank", "h9_selection_order", "observation_id", "photo_id"],
        kind="mergesort",
        ignore_index=True,
    )
    return audit, terminal_species, rows


def query_panel(
    *,
    panel: str,
    queue: pd.DataFrame,
    client: object,
    exclusion_obs: set[int],
    exclusion_photo: set[int],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, object]]:
    frozen = freeze_h9_metadata(
        client=client,
        species_frame=queue[["inat_taxon_id", "species"]].copy(),
        exclusion_observation_ids=exclusion_obs,
        exclusion_photo_ids=exclusion_photo,
        per_page=PER_PAGE,
        observer_cap_n=OBSERVER_CAP,
        fixed_raw_photos=ROWS_PER_SPECIES,
        maximum_positional_accuracy_m=MAX_ACCURACY_M,
        allowed_photo_licenses=ALLOWED_LICENSES,
    )
    observations = frozen.observations.copy()
    if len(observations):
        observations = observations.merge(
            queue[["queue_rank", "inat_taxon_id"]],
            on="inat_taxon_id",
            how="left",
            validate="many_to_one",
        )
        observations = observations.sort_values(
            ["queue_rank", "h9_selection_order", "observation_id", "photo_id"],
            kind="mergesort",
            ignore_index=True,
        )
        obs_ids = pd.to_numeric(observations["observation_id"], errors="raise").astype(int)
        photo_ids = pd.to_numeric(observations["photo_id"], errors="raise").astype(int)
        if obs_ids.duplicated().any() or photo_ids.duplicated().any():
            raise RuntimeError(f"Panel {panel} fresh metadata contains duplicate IDs")
        if set(obs_ids) & exclusion_obs or set(photo_ids) & exclusion_photo:
            raise RuntimeError(f"Panel {panel} prior ID leaked into fresh metadata")

    audit, terminal_species, terminal_rows = terminalize_panel(
        panel=panel,
        queue=queue,
        species_audit=frozen.species_audit,
        observations=observations,
    )
    request_errors = int(audit["request_error"].fillna("").astype(str).str.len().gt(0).sum())
    full_species = int(audit["full_fixed_n"].astype(bool).sum())
    error_fraction = request_errors / EXPECTED_QUEUE_N
    passed = (
        error_fraction <= REQUEST_ERROR_CEILING
        and full_species >= TERMINAL_SPECIES_PER_PANEL
        and len(terminal_species) == TERMINAL_SPECIES_PER_PANEL
        and len(terminal_rows) == TERMINAL_SPECIES_PER_PANEL * ROWS_PER_SPECIES
    )
    summary = {
        "panel": panel,
        "queue_species": EXPECTED_QUEUE_N,
        "request_attempts": EXPECTED_QUEUE_N,
        "request_errors": request_errors,
        "request_error_fraction": error_fraction,
        "request_error_fraction_ceiling": REQUEST_ERROR_CEILING,
        "full100_species_in_queue": full_species,
        "terminal_species": int(len(terminal_species)),
        "terminal_rows": int(len(terminal_rows)),
        "target_species": TERMINAL_SPECIES_PER_PANEL,
        "target_rows_per_species": ROWS_PER_SPECIES,
        "pass": bool(passed),
    }
    return audit, observations, terminal_species, terminal_rows, summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel-p-queue", type=Path, required=True)
    parser.add_argument("--panel-n-queue", type=Path, required=True)
    parser.add_argument("--exclusion-source", type=Path, action="append", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    if os.environ.get("GITHUB_ACTIONS") == "true" and os.environ.get("GITHUB_RUN_ATTEMPT", "1") != "1":
        raise RuntimeError("rerun forbidden: FCP v2 fresh-metadata draw is one-shot")

    p_queue = validate_queue(pd.read_csv(args.panel_p_queue), "P")
    n_queue = validate_queue(pd.read_csv(args.panel_n_queue), "N")
    if set(p_queue["inat_taxon_id"]) & set(n_queue["inat_taxon_id"]):
        raise RuntimeError("Panel P/N queue taxon overlap")
    if set(p_queue["species"]) & set(n_queue["species"]):
        raise RuntimeError("Panel P/N queue species overlap")

    exclusion_obs, exclusion_photo, exclusion_audit = load_exclusion_union(args.exclusion_source)

    client = InaturalistObservationClient(
        user_agent="zuizui0223-fcp-v2-measurement-validity/1.0 (github.com/zuizui0223/fcp)",
        request_interval_seconds=REQUEST_INTERVAL_SECONDS,
        timeout_seconds=REQUEST_TIMEOUT_SECONDS,
        max_retries=REQUEST_RETRIES,
    )

    panel_outputs = {}
    for panel, queue in (("P", p_queue), ("N", n_queue)):
        audit, all_rows, terminal_species, terminal_rows, summary = query_panel(
            panel=panel,
            queue=queue,
            client=client,
            exclusion_obs=exclusion_obs,
            exclusion_photo=exclusion_photo,
        )
        panel_outputs[panel] = {
            "audit": audit,
            "all_rows": all_rows,
            "terminal_species": terminal_species,
            "terminal_rows": terminal_rows,
            "summary": summary,
        }

    p_terminal = panel_outputs["P"]["terminal_species"]
    n_terminal = panel_outputs["N"]["terminal_species"]
    p_rows = panel_outputs["P"]["terminal_rows"]
    n_rows = panel_outputs["N"]["terminal_rows"]

    passed = bool(panel_outputs["P"]["summary"]["pass"] and panel_outputs["N"]["summary"]["pass"])
    if passed:
        combined_species = pd.concat([p_terminal, n_terminal], ignore_index=True)
        combined_rows = pd.concat([p_rows, n_rows], ignore_index=True)
        if len(combined_species) != 400 or len(combined_rows) != 40_000:
            raise RuntimeError("combined terminal denominator drift")
        if combined_species["inat_taxon_id"].nunique() != 400:
            raise RuntimeError("combined terminal taxon overlap")
        if combined_species["species"].nunique() != 400:
            raise RuntimeError("combined terminal species overlap")
        if combined_rows["observation_id"].nunique() != 40_000:
            raise RuntimeError("combined terminal observation overlap")
        if combined_rows["photo_id"].nunique() != 40_000:
            raise RuntimeError("combined terminal photo overlap")
    else:
        combined_species = pd.concat([p_terminal, n_terminal], ignore_index=True)
        combined_rows = pd.concat([p_rows, n_rows], ignore_index=True)

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    for panel in ("P", "N"):
        x = panel_outputs[panel]
        x["audit"].to_csv(out / f"panel_{panel}_species_audit.csv", index=False, lineterminator="\n")
        deterministic_gzip_csv(x["all_rows"], out / f"panel_{panel}_all_fresh_metadata.csv.gz")
        x["terminal_species"].to_csv(
            out / f"panel_{panel}_terminal_species.csv", index=False, lineterminator="\n"
        )
        deterministic_gzip_csv(
            x["terminal_rows"], out / f"panel_{panel}_terminal_metadata.csv.gz"
        )

    combined_species.to_csv(
        out / "terminal_species_manifest.tsv",
        sep="\t",
        index=False,
        lineterminator="\n",
    )
    deterministic_gzip_csv(combined_rows, out / "authorized_metadata_400x100.csv.gz")

    result = {
        "schema": "fcp_v2_fresh_metadata_terminalization_v1",
        "date_jst": "2026-09-23",
        "status": (
            "FCP_V2_FRESH_METADATA_PASS"
            if passed
            else "FCP_V2_UNDERIDENTIFIED_PREPIXEL"
        ),
        "panels": {
            panel: panel_outputs[panel]["summary"]
            for panel in ("P", "N")
        },
        "terminal": {
            "species": int(len(combined_species)),
            "rows": int(len(combined_rows)),
            "target_species": 400,
            "target_rows": 40_000,
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
            "q_white_or_W_opened": False,
            "spatial_colour_outcome_opened": False,
            "H3_predictors_opened": False,
        },
        "no_rescue": {
            "queue_extension_allowed": False,
            "queue_reordering_allowed": False,
            "species_replacement_after_pixel_opening_allowed": False,
            "metadata_rerun_allowed": False,
        },
        "files": {},
    }
    for path in sorted(out.iterdir()):
        if path.is_file() and path.name != "result.json":
            result["files"][path.name] = {
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
            }
    (out / "result.json").write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
