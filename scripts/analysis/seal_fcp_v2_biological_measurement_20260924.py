#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path

import pandas as pd

EXPECTED_ROWS = 40_000
EXPECTED_PARTITIONS = 256
EXPECTED_HEAVY = 8_000
BIO_MORPHS = {"white", "yellow_orange", "red_pink", "blue_purple"}


def canonical_sha(frame: pd.DataFrame) -> str:
    if frame["measurement_id"].astype(str).nunique() != len(frame):
        raise RuntimeError("biological table IDs are not unique")
    ordered = frame.sort_values("measurement_id", kind="mergesort").reset_index(drop=True)
    buf = io.StringIO()
    ordered.to_csv(buf, index=False, lineterminator="\n")
    return hashlib.sha256(buf.getvalue().encode("utf-8")).hexdigest()


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--partition-results-dir", type=Path, required=True)
    p.add_argument("--join-dir", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()

    files = sorted(args.partition_results_dir.rglob("biological_b*_s*_p*.csv"))
    if len(files) != EXPECTED_PARTITIONS:
        raise RuntimeError(f"biological partition census incomplete: {len(files)}")
    parts = [pd.read_csv(path, dtype={"measurement_id": str}).fillna("") for path in files]
    biological = pd.concat(parts, ignore_index=True, sort=False)
    if len(biological) != EXPECTED_ROWS:
        raise RuntimeError(f"biological row census drift: {len(biological)}")
    if biological["measurement_id"].nunique() != EXPECTED_ROWS:
        raise RuntimeError("biological measurement IDs are not unique")

    technical_summary = json.loads((args.join_dir / "technical_summary.json").read_text())
    if technical_summary.get("status") != "RESPONSE_BLIND_TECHNICAL_SEAL_FROZEN":
        raise RuntimeError("technical seal status drift")
    if technical_summary.get("biological_outcomes_opened") is not False:
        raise RuntimeError("technical seal was not response blind at handoff")

    join = pd.read_csv(args.join_dir / "sealed_join_key.csv", dtype={"measurement_id": str})
    strata = pd.read_csv(args.join_dir / "technical_strata.csv", dtype={"measurement_id": str})
    if len(join) != EXPECTED_ROWS or join["measurement_id"].nunique() != EXPECTED_ROWS:
        raise RuntimeError("sealed biological join key census drift")
    if len(strata) != EXPECTED_ROWS or strata["measurement_id"].nunique() != EXPECTED_ROWS:
        raise RuntimeError("technical strata census drift")

    joined = join.merge(
        strata,
        on="measurement_id",
        how="inner",
        validate="one_to_one",
    ).merge(
        biological,
        on="measurement_id",
        how="inner",
        validate="one_to_one",
        suffixes=("", "_bio"),
    )
    if len(joined) != EXPECTED_ROWS:
        raise RuntimeError("post-measurement biological join lost rows")

    heavy = joined["heavy_counterfactual"].fillna("").astype(str).str.casefold().isin({"true", "1", "yes"})
    if int(heavy.sum()) != EXPECTED_HEAVY:
        raise RuntimeError("heavy-counterfactual denominator drift after biological opening")

    source = joined["source_identity_status"].astype(str)
    exact = source.eq("exact_source_sha_match")
    pass_b_failed = source.eq("pass_b_acquisition_failed")
    drift = source.eq("source_byte_drift")
    measurement_failed = source.eq("exact_source_sha_match_measurement_failed")
    if int((exact | pass_b_failed | drift | measurement_failed).sum()) != EXPECTED_ROWS:
        raise RuntimeError("unknown biological source identity terminal state")

    base_morph = joined["base_morph"].astype(str)
    base_status = joined["base_status"].astype(str)
    classified = base_status.eq("classified_four_state_morph") & base_morph.isin(BIO_MORPHS)
    mixed = base_morph.eq("mixed_uncertain")

    digest = canonical_sha(joined)
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    joined.to_csv(
        out / "biological_table.csv.gz",
        index=False,
        compression={"method": "gzip", "compresslevel": 9, "mtime": 0},
        lineterminator="\n",
    )

    summary = {
        "schema": "fcp_v2_biological_pass_seal_v1",
        "status": "BIOLOGICAL_PASS_B_SEALED",
        "rows": EXPECTED_ROWS,
        "partition_receipts": EXPECTED_PARTITIONS,
        "technical_table_sha256": technical_summary["technical_table_sha256"],
        "biological_table_sha256": digest,
        "source_identity": {
            "exact_sha_match_rows": int(exact.sum()),
            "pass_b_acquisition_failed_rows": int(pass_b_failed.sum()),
            "source_byte_drift_rows": int(drift.sum()),
            "exact_sha_match_measurement_failed_rows": int(measurement_failed.sum()),
            "replacement_rows": 0,
        },
        "base_measurement": {
            "classified_four_state_rows": int(classified.sum()),
            "mixed_uncertain_rows": int(mixed.sum()),
            "other_or_missing_rows": int(EXPECTED_ROWS - classified.sum() - mixed.sum()),
        },
        "heavy_counterfactual_rows": int(heavy.sum()),
        "biological_outcomes_opened": True,
        "species_opened_after_complete_partition_measurement": True,
        "coordinates_opened_after_complete_partition_measurement": True,
        "q_white_or_W_computed": False,
        "MV1_computed": False,
        "MV2_computed": False,
        "MV3_computed": False,
        "MV4_computed": False,
        "no_rescue": {
            "source_drift_replacement_allowed": False,
            "species_replacement_allowed": False,
            "counterfactual_retuning_allowed": False,
        },
    }
    (out / "biological_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
