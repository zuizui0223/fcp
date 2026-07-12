"""Generate a phase-grid CSV for the minimal floral colour polymorphism model."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from model import classify_phase


def frange(start: float, stop: float, n: int) -> list[float]:
    if n < 2:
        raise ValueError("n must be >= 2")
    step = (stop - start) / (n - 1)
    return [start + i * step for i in range(n)]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--b-min", type=float, default=-2.0)
    parser.add_argument("--b-max", type=float, default=2.0)
    parser.add_argument("--d-min", type=float, default=-2.0)
    parser.add_argument("--d-max", type=float, default=2.0)
    parser.add_argument("--n", type=int, default=101)
    parser.add_argument("--out", type=Path, default=Path("results/phase_grid.csv"))
    args = parser.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "balancing_strength", "directional_bias", "invasion_A_into_B",
            "invasion_B_into_A", "phase", "equilibrium_p_A"
        ])
        for b in frange(args.b_min, args.b_max, args.n):
            for d in frange(args.d_min, args.d_max, args.n):
                result = classify_phase(b, d)
                writer.writerow([
                    result.balancing_strength,
                    result.directional_bias,
                    result.invasion_A_into_B,
                    result.invasion_B_into_A,
                    result.phase,
                    "" if result.equilibrium_p_A is None else result.equilibrium_p_A,
                ])


if __name__ == "__main__":
    main()
