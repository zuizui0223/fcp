"""Metadata-only admission and strict blinded measurement for RGFCA replication.

The reserve is the entire complement of the original hash-selected 500 taxa.
No colour outcome is consulted in selection. Historical inputs are read from
their exact Git objects, avoiding checkout newline transformations on Windows.
"""
from __future__ import annotations

import hashlib
import io
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from fcp_pipeline.global_measurement_budget import select_measurement_rows

ROOT = Path(__file__).resolve().parents[1]
SOURCE_COMMIT = "29584f3ad7ae0cd99a1d8f43459252af38f615da"
PROTOCOL = "rgfca-reserve-species-replication-v1"
CONTRACT = ROOT / "docs/supporting/rgfca_reserve_replication_contract_v1.json"
AUDIT = ROOT / "docs/supporting/rgfca_reserve_metadata_audit_v1.json"
PREFIX = "rgfca_reserve_replication"
CANDIDATE = "data/frozen/global_monte_carlo_candidate_photos_v1.csv"
DISCOVERY = "data/derived/global_monte_carlo_measured_photos_v1.csv"
PRIOR = (
    "data/derived/random_photo_first_measured_photos_v1.csv",
    "data/frozen/random_photo_first_h9_fresh_metadata_v1.csv",
    "data/frozen/random_photo_first_h9_exclusion_ledger_v1.csv",
)
ID_FIELDS = ["inat_taxon_id", "photo_id", "observation_id"]
COLOURS = ["colour_white", "colour_yellow_orange", "colour_red_pink", "colour_blue_purple"]


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def git_bytes(path: str) -> bytes:
    return subprocess.check_output(["git", "show", f"{SOURCE_COMMIT}:{path}"], cwd=ROOT)


def id_digest(values) -> str:
    return sha(("\n".join(map(str, sorted(int(v) for v in values))) + "\n").encode())


def reserve_rows(candidate: pd.DataFrame, discovery_ids: pd.DataFrame, *, target=100, budget=500) -> pd.DataFrame:
    """Require the original budget exactly; return its unselected complement."""
    for field in ("photo_id", "observation_id"):
        if candidate[field].isna().any() or candidate[field].duplicated().any():
            raise ValueError(f"candidate {field} is missing or duplicated")
        if discovery_ids[field].isna().any() or discovery_ids[field].duplicated().any():
            raise ValueError(f"discovery {field} is missing or duplicated")
    selected = select_measurement_rows(candidate, target_photos_per_species=target, maximum_species=budget, seed=20260918)
    a = selected[ID_FIELDS].sort_values(ID_FIELDS).reset_index(drop=True).astype("int64")
    b = discovery_ids[ID_FIELDS].sort_values(ID_FIELDS).reset_index(drop=True).astype("int64")
    if not a.equals(b):
        raise ValueError("discovery IDs differ from the original outcome-blind species budget")
    reserve = candidate.loc[~candidate.inat_taxon_id.isin(selected.inat_taxon_id)].copy()
    if len(reserve) != budget * target or reserve.inat_taxon_id.nunique() != budget:
        raise ValueError("entire reserve must contain exactly the other 500 taxa / 50,000 rows")
    return reserve.sort_values(["inat_taxon_id", "photo_id"], kind="stable").reset_index(drop=True)


def overlap_counts(reserve: pd.DataFrame, prior: pd.DataFrame) -> dict:
    return {field: len(set(reserve[field].astype(int)) & set(prior[field].dropna().astype(int)))
            for field in ("photo_id", "observation_id")}


def geometry_audit(reserve: pd.DataFrame) -> pd.DataFrame:
    from fcp_pipeline.global_g3 import _upper_triangle_distances_km

    numeric = reserve[["latitude", "longitude", "positional_accuracy_m", "observer_id"]].to_numpy(float)
    if not np.isfinite(numeric).all():
        raise ValueError("missing or nonfinite required metadata")
    if not reserve.latitude.between(-90, 90).all() or not reserve.longitude.between(-180, 180).all():
        raise ValueError("coordinate out of range")
    if not reserve.positional_accuracy_m.between(0, 5000).all():
        raise ValueError("positional accuracy exceeds the fixed 5-km ceiling")
    if int(reserve.groupby(["inat_taxon_id", "observer_id"]).size().max()) > 2:
        raise ValueError("observer cap exceeds two photos per species")
    dates = pd.to_datetime(reserve.observed_on, errors="raise")
    if dates.isna().any():
        raise ValueError("missing observation date")
    frame = reserve.assign(quarter=dates.dt.quarter, year=dates.dt.year)
    rows = []
    for taxon, g in frame.groupby("inat_taxon_id", sort=True):
        dist = _upper_triangle_distances_km(g.latitude.to_numpy(), g.longitude.to_numpy())
        row = {"inat_taxon_id": int(taxon), "species": str(g.species.iloc[0]), "raw_photos": len(g),
               "observers": g.observer_id.nunique(), "quarters": g.quarter.nunique(),
               "years": g.year.nunique(), "maximum_span_km": float(dist.max()),
               "median_pair_distance_km": float(np.median(dist))}
        for step in (0.25, 0.5):
            cells = np.floor(g[["latitude", "longitude"]].to_numpy() / step).astype(int)
            _, counts = np.unique(cells, axis=0, return_counts=True)
            label = str(step).replace(".", "p")
            row[f"cells_{label}_degree"] = len(counts)
            row[f"maximum_cell_photos_{label}_degree"] = int(counts.max())
        u, v = np.triu_indices(len(g), k=1)
        for cap in (100, 250, 500):
            keep = dist <= cap
            row[f"pairs_within_{cap}_km"] = int(keep.sum())
            row[f"photos_with_neighbor_within_{cap}_km"] = int(len(np.unique(np.r_[u[keep], v[keep]])))
        rows.append(row)
    # Metadata summaries only: millimetre-scale kilometre rounding avoids
    # platform trig/CSV differences without altering coordinates or inference.
    return pd.DataFrame(rows).round(6)


def audit_reserve() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    inputs = {path: git_bytes(path) for path in (CANDIDATE, DISCOVERY, *PRIOR)}
    if sha(inputs[CANDIDATE]) != "f1319461d8883f3094cec8ac0e5fc247ff902b464f34146af275575b94edc9d2":
        raise ValueError("original candidate hash drift")
    if sha(inputs[DISCOVERY]) != "ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4":
        raise ValueError("original discovery hash drift")
    candidate = pd.read_csv(io.BytesIO(inputs[CANDIDATE]))
    discovery = pd.read_csv(io.BytesIO(inputs[DISCOVERY]), usecols=ID_FIELDS)
    reserve = reserve_rows(candidate, discovery)
    overlaps = {DISCOVERY: overlap_counts(reserve, discovery)}
    for path in PRIOR:
        previous = pd.read_csv(io.BytesIO(inputs[path]), usecols=["photo_id", "observation_id"])
        overlaps[path] = overlap_counts(reserve, previous)
    if any(value for counts in overlaps.values() for value in counts.values()):
        raise ValueError("reserve overlaps an already measured or prior-experiment observation/photo")
    geometry = geometry_audit(reserve)
    report = {
        "protocol": PROTOCOL, "status": "metadata_only_reserve_admission_passed",
        "source_commit": SOURCE_COMMIT, "rows": len(reserve), "species": reserve.inat_taxon_id.nunique(),
        "photos_per_species": 100, "maximum_positional_accuracy_m": float(reserve.positional_accuracy_m.max()),
        "maximum_observer_photos_per_species": int(reserve.groupby(["inat_taxon_id", "observer_id"]).size().max()),
        "selection": "entire complement of the original fixed hash-ranked 500-taxon measurement budget; no replacement",
        "ids_sha256": {f: id_digest(reserve[f].unique()) for f in ID_FIELDS},
        "source_sha256": {p: sha(b) for p, b in inputs.items()}, "overlap_counts": overlaps,
        "discovery_taxon_overlap": int(len(set(reserve.inat_taxon_id) & set(discovery.inat_taxon_id))),
        "geometry_summary": {col: {"min": float(geometry[col].min()), "median": float(geometry[col].median()),
                                   "max": float(geometry[col].max())}
                             for col in geometry.select_dtypes(include="number") if col != "inat_taxon_id"},
        "candidate_pixels_opened_by_audit": False, "colour_fields_parsed": False,
        "measurement_authorized": False,
        "independence_ceiling": "New observation and photo IDs, species-disjoint from RGFCA discovery, same iNaturalist sampling frame and measurement model. Shared observers/regions and systematic measurement error remain possible. This audit cannot prove nobody outside the recorded workflows ever viewed these public photos.",
    }
    return reserve, geometry, report


def verify_audit(reserve: pd.DataFrame, report: dict) -> None:
    if report.get("status") != "metadata_only_reserve_admission_passed" or report.get("protocol") != PROTOCOL:
        raise ValueError("reserve audit not passed")
    if len(reserve) != 50000 or reserve.inat_taxon_id.nunique() != 500:
        raise ValueError("reserve denominator changed")
    if report["ids_sha256"] != {f: id_digest(reserve[f].unique()) for f in ID_FIELDS}:
        raise ValueError("reserve identities changed")


def build_firewall(reserve: pd.DataFrame, report: dict, out: Path, contract: dict) -> dict:
    # Imported only after materializing the exact validated infrastructure.
    from fcp_pipeline.photo_first_measurement_execution import WORKER_FIELDS, ACQUISITION_FIELDS, semantic_shard, compute_partition
    from scripts.analysis.build_global_monte_carlo_measurement_firewall import measurement_id, measurement_batch

    verify_audit(reserve, report)
    if contract.get("protocol") != PROTOCOL or contract.get("status") != "frozen_before_reserve_pixels":
        raise ValueError("replication contract not frozen")
    m = reserve.copy()
    m.insert(0, "measurement_id", [measurement_id(v, salt=PROTOCOL) for v in m.photo_id])
    m.insert(1, "measurement_batch", [measurement_batch(v) for v in m.measurement_id])
    if m.measurement_id.nunique() != 50000:
        raise ValueError("replication blinded ID collision")
    worker = pd.DataFrame({"measurement_id": m.measurement_id, "image_filename": m.measurement_id + ".jpg", "photo_license": m.photo_license})
    acquire = worker.assign(photo_url_large=m.photo_url_large)[list(ACQUISITION_FIELDS)]
    if tuple(worker.columns) != WORKER_FIELDS:
        raise ValueError("blind worker interface drift")
    (out / "sealed_keys").mkdir(parents=True, exist_ok=True)
    join = out / "sealed_keys/metadata_join_key.csv"
    m.drop(columns="photo_url_large").to_csv(join, index=False, lineterminator="\n")
    files = {"sealed_keys/metadata_join_key.csv": sha(join.read_bytes())}
    assignments = []
    for batch in range(2):
        keep = m.measurement_batch.eq(batch)
        for suffix, frame in (("worker_packet/measurement_manifest.csv", worker), ("sealed_keys/acquisition_key.csv", acquire)):
            path = out / f"batch_{batch}" / suffix
            path.parent.mkdir(parents=True, exist_ok=True)
            frame.loc[keep].to_csv(path, index=False, lineterminator="\n")
            files[path.relative_to(out).as_posix()] = sha(path.read_bytes())
        for mid in m.loc[keep, "measurement_id"]:
            assignments.append((mid, batch, semantic_shard(mid), compute_partition(mid)))
    assignments = pd.DataFrame(assignments, columns=["measurement_id", "measurement_batch", "semantic_shard", "compute_partition"])
    assignments.to_csv(out / "partition_assignments.csv", index=False, lineterminator="\n")
    files["partition_assignments.csv"] = sha((out / "partition_assignments.csv").read_bytes())
    result = {"protocol": PROTOCOL, "status": "replication_measurement_firewall_frozen_before_pixels",
              "frozen_rows": 50000, "frozen_species": 500, "terminal_partitions": 256,
              "worker_fields": list(WORKER_FIELDS), "files_sha256": files,
              "reserve_ids_sha256": report["ids_sha256"], "candidate_pixels_opened": False,
              "coordinate_colour_join_opened": False, "contract_sha256": sha(CONTRACT.read_bytes())}
    (out / "measurement_firewall_manifest.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8", newline="\n")
    return result


def palette_vectors(joined: pd.DataFrame) -> pd.DataFrame:
    """Identical nine-to-four probability mapping to the discovery measurement."""
    joined = joined.copy()
    groups = (("white",), ("yellow", "orange", "bronze"), ("red", "pink", "magenta"), ("blue", "purple"))
    for name, palette in zip(COLOURS, groups):
        joined[name] = sum(pd.to_numeric(joined[f"flower_fraction_{p}"], errors="raise") for p in palette)
    joined["global_classifiable"] = joined.morph.isin({"white", "yellow_orange", "red_pink", "blue_purple"}) & joined.measurement_status.eq("classified_four_state_morph")
    values = joined.loc[joined.global_classifiable, COLOURS].to_numpy(float)
    if not np.isfinite(values).all() or (values < 0).any() or not np.allclose(values.sum(axis=1), 1, rtol=0, atol=1e-8):
        raise ValueError("classifiable replication vectors are invalid")
    from fcp_pipeline.photo_first_measurement import REFERENCE_RGB
    counts = joined.loc[joined.global_classifiable, [f"background_palette_count_{p}" for p in REFERENCE_RGB]].apply(pd.to_numeric, errors="raise").to_numpy(float)
    pixels = pd.to_numeric(joined.loc[joined.global_classifiable, "background_effective_pixels"], errors="raise").to_numpy(float)
    if (not np.isfinite(counts).all() or (counts < 0).any() or not np.equal(counts, np.floor(counts)).all()
            or not np.isfinite(pixels).all() or (pixels < 100).any() or not np.equal(counts.sum(axis=1), pixels).all()):
        raise ValueError("matched background counts are incomplete or invalid")
    return joined


def reassemble(firewall_dir: Path, results_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    from fcp_pipeline.photo_first_measurement_execution import reassemble_complete_measurement, validate_terminal_partition_results

    firewall = json.loads((firewall_dir / "measurement_firewall_manifest.json").read_text())
    if firewall.get("protocol") != PROTOCOL or firewall.get("status") != "replication_measurement_firewall_frozen_before_pixels":
        raise ValueError("replication firewall identity mismatch")
    if firewall["contract_sha256"] != sha(CONTRACT.read_bytes()):
        raise ValueError("replication contract changed after firewall creation")
    for path, expected in firewall["files_sha256"].items():
        if sha((firewall_dir / path).read_bytes()) != expected:
            raise ValueError(f"firewall file changed: {path}")
    assignment = pd.read_csv(firewall_dir / "partition_assignments.csv")
    frames = []
    receipt_hashes = {}
    expected_names = {f"b{b}_partition_s{s:02d}_p{p:02d}.{ext}" for b in range(2) for s in range(32) for p in range(4) for ext in ("csv", "json")}
    if {p.name for p in results_dir.iterdir()} != expected_names:
        raise ValueError("not_evaluable_incomplete_or_unexpected_replication_partitions")
    for batch in range(2):
        for semantic in range(32):
            for compute in range(4):
                stem = f"b{batch}_partition_s{semantic:02d}_p{compute:02d}"
                path = results_dir / f"{stem}.csv"
                receipt = json.loads((results_dir / f"{stem}.json").read_text())
                if receipt.get("status") != "complete_random_photo_first_terminal_partition":
                    raise ValueError(f"nonterminal partition: {stem}")
                if receipt.get("semantic_shard") != semantic or receipt.get("compute_partition") != compute:
                    raise ValueError(f"partition receipt identity changed: {stem}")
                if any(receipt.get(k) is not False for k in ("source_urls_present", "species_present", "coordinates_present")):
                    raise ValueError("terminal receipt leaked context")
                frame = pd.read_csv(path, dtype={"measurement_id": str}).fillna("")
                expected = assignment.loc[(assignment.measurement_batch == batch) & (assignment.semantic_shard == semantic) & (assignment.compute_partition == compute), "measurement_id"]
                validate_terminal_partition_results(frame, expected.tolist())
                if receipt.get("terminal_rows") != len(frame):
                    raise ValueError(f"partition receipt row mismatch: {stem}")
                frames.append(frame)
                receipt_hashes[stem] = {"csv": sha(path.read_bytes()), "json": sha((results_dir / f"{stem}.json").read_bytes())}
    workers = pd.concat([pd.read_csv(firewall_dir / f"batch_{b}/worker_packet/measurement_manifest.csv") for b in range(2)], ignore_index=True)
    metadata = pd.read_csv(firewall_dir / "sealed_keys/metadata_join_key.csv")
    if len(workers) != 50000 or metadata.inat_taxon_id.nunique() != 500:
        raise ValueError("replication denominator drift")
    joined = palette_vectors(reassemble_complete_measurement(frames, workers, metadata, expected_partition_receipts=256).joined_photos)
    if not joined.groupby("inat_taxon_id").size().eq(100).all():
        raise ValueError("replication per-species denominator drift")
    support = joined.groupby(["inat_taxon_id", "species"]).agg(raw_photos=("photo_id", "size"), classifiable_photos=("global_classifiable", "sum")).reset_index()
    support["evaluable"] = support.classifiable_photos >= 40
    result = {"protocol": PROTOCOL, "status": "complete_reserve_location_blind_measurement_and_join",
              "rows": len(joined), "species": len(support), "classifiable_photos": int(joined.global_classifiable.sum()),
              "evaluable_species": int(support.evaluable.sum()), "measurement_gate_pass": bool(support.evaluable.sum() >= 250),
              "terminal_partitions": 256, "inference_run": False,
              "coordinate_colour_join_opened_after_complete_measurement": True,
              "firewall_sha256": sha((firewall_dir / "measurement_firewall_manifest.json").read_bytes()),
              "partition_hashes": receipt_hashes}
    return joined, support, result
