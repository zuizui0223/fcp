"""Build the preregistered angiosperm background-sampling manifest.

This script expands data/background_sampling_strata.csv into one row per planned
species draw.  Species names are intentionally left blank until a frozen taxonomic
pool is supplied.  The manifest is therefore independent of whether a species is
already known to be colour-polymorphic.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

FIELDS = [
    "sample_id",
    "stratum_id",
    "broad_clade",
    "family",
    "expected_pollination_axis",
    "species",
    "taxon_status",
    "pollination_mode",
    "pollination_evidence",
    "display_structure",
    "display_opportunity",
    "outcome_class",
    "within_population_coexistence",
    "among_population_differentiation",
    "outcome_evidence",
    "sampling_effort_proxy",
    "range_size_proxy",
    "environmental_heterogeneity_proxy",
    "connectivity_proxy",
    "evidence_status",
    "notes",
]


def build_manifest(strata_path: Path, output_path: Path) -> int:
    rows: list[dict[str, str]] = []
    with strata_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"stratum_id", "broad_clade", "family", "expected_pollination_axis", "target_n"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"missing required strata columns: {sorted(missing)}")
        for stratum in reader:
            target_n = int(stratum["target_n"])
            if target_n <= 0:
                raise ValueError(f"target_n must be positive for {stratum['stratum_id']}")
            for draw in range(1, target_n + 1):
                row = {field: "" for field in FIELDS}
                row.update(
                    sample_id=f"{stratum['stratum_id']}-{draw:03d}",
                    stratum_id=stratum["stratum_id"],
                    broad_clade=stratum["broad_clade"],
                    family=stratum["family"],
                    expected_pollination_axis=stratum["expected_pollination_axis"],
                    taxon_status="unassigned",
                    evidence_status="pending_random_draw",
                )
                rows.append(row)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--strata", default="data/background_sampling_strata.csv")
    parser.add_argument("--out", default="data/background_sample_manifest.csv")
    args = parser.parse_args()
    n = build_manifest(Path(args.strata), Path(args.out))
    print(f"wrote {n} preregistered species-draw slots to {args.out}")


if __name__ == "__main__":
    main()
