#!/usr/bin/env python3
"""Verify the SF empirical transport against its frozen v0.6 receipt."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


NUMERIC_KEYS = (
    "species_conditioned_mean_rho",
    "species_conditioned_p_upper",
    "spread_spatial_rho",
    "spread_spatial_p_upper",
    "naive_pooled_rho",
    "naive_pooled_p",
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--observed", type=Path, required=True)
    parser.add_argument("--expected", type=Path, required=True)
    args = parser.parse_args()

    observed = json.loads(args.observed.read_text(encoding="utf-8"))
    expected = json.loads(args.expected.read_text(encoding="utf-8"))

    if observed["analysis"]["taxa"] != expected["fixture"]["taxa"]:
        raise RuntimeError("taxon count drift")
    if observed["analysis"]["observations"] != expected["fixture"]["observations"]:
        raise RuntimeError("observation count drift")
    if (
        observed["analysis"]["permutations_per_taxon"]
        != expected["analysis"]["permutations_per_taxon"]
    ):
        raise RuntimeError("permutation count drift")

    for key in NUMERIC_KEYS:
        a = float(observed["analysis"][key])
        b = float(expected["analysis"][key])
        if not math.isclose(a, b, rel_tol=0.0, abs_tol=1e-12):
            raise RuntimeError(f"{key} drift: observed={a!r}, expected={b!r}")

    print(
        json.dumps(
            {
                "status": "PASS",
                "schema": expected["schema"],
                "taxa": observed["analysis"]["taxa"],
                "observations": observed["analysis"]["observations"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
