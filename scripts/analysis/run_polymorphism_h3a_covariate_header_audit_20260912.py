#!/usr/bin/env python3
"""Header-only audit for H3a sampling-opportunity covariates.

This script reads zero data rows. It does not read morph/fine_state values, compute D,
or evaluate any outcome-covariate association.
"""
from __future__ import annotations

import json
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "polymorphism_h3a_covariate_header_audit_20260912"
OUT.mkdir(parents=True, exist_ok=True)
FILES = {
    "discovery": ROOT / "data" / "derived" / "global_monte_carlo_measured_photos_v1.csv",
    "reserve": ROOT / "data" / "derived" / "rgfca_reserve_replication_measured_photos_v1.csv",
}


def main() -> None:
    cohorts = {}
    for cohort, path in FILES.items():
        header = pd.read_csv(path, nrows=0)
        cols = list(map(str, header.columns))
        cohorts[cohort] = {
            "path": str(path.relative_to(ROOT)),
            "columns": cols,
            "n_columns": len(cols),
        }
    result = {
        "analysis": "polymorphism_h3a_covariate_header_audit",
        "date_jst": "2026-09-12",
        "rows_read": 0,
        "D_computed": False,
        "outcome_values_read": False,
        "association_computed": False,
        "cohorts": cohorts,
    }
    (OUT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
