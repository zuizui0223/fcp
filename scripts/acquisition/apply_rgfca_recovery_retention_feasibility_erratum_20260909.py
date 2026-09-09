#!/usr/bin/env python3
"""Materialize the frozen fail-closed retention-feasibility erratum at runtime.

The parent recovery runner is preserved byte-for-byte as the record of the
first failed execution. This patcher makes two exact source substitutions and
fails if the expected parent source has drifted. It changes no query, page,
identity, pixel, colour, M1-M3, or biological threshold rule; it only makes the
already-frozen <=2 photos/observer retention rule executable.
"""
from __future__ import annotations

import argparse
from pathlib import Path

ERRATUM_COMMIT = "d98dacc1d4f42c47d482b51d93c75465ced3c47f"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    return p.parse_args()


def replace_once(text: str, old: str, new: str, label: str) -> str:
    n = text.count(old)
    if n != 1:
        raise RuntimeError(f"expected exactly one {label} patch target, found {n}")
    return text.replace(old, new, 1)


def main() -> int:
    args = parse_args()
    source = args.source.read_text(encoding="utf-8")

    old_group = '''                n_observer = int(g["observer_id"].nunique())
                qualified = (
                    n_photo >= qmod.MIN_PHOTOS_PER_GROUP
                    and n_obs >= qmod.MIN_PHOTOS_PER_GROUP
                    and n_observer >= qmod.MIN_OBSERVERS_PER_GROUP
                )
                local_groups.append({
'''
    new_group = '''                n_observer = int(g["observer_id"].nunique())
                capped_capacity = int(
                    g.groupby("observer_id")["photo_id"].nunique().clip(upper=2).sum()
                )
                qualified = (
                    n_photo >= qmod.MIN_PHOTOS_PER_GROUP
                    and n_obs >= qmod.MIN_PHOTOS_PER_GROUP
                    and n_observer >= qmod.MIN_OBSERVERS_PER_GROUP
                    and capped_capacity >= qmod.MIN_PHOTOS_PER_GROUP
                )
                local_groups.append({
'''
    source = replace_once(source, old_group, new_group, "retention-feasibility")

    old_field = '''                    "n_observers": n_observer,
                    "centroid_x_m": float(g["x_m"].mean()),
'''
    new_field = '''                    "n_observers": n_observer,
                    "max_two_per_observer_photo_capacity": capped_capacity,
                    "centroid_x_m": float(g["x_m"].mean()),
'''
    source = replace_once(source, old_field, new_field, "group-audit")

    old_manifest = '''        "protocol_commit": PROTOCOL_COMMIT,
        "status": "complete_metadata_only_recovery_shard",
'''
    new_manifest = f'''        "protocol_commit": PROTOCOL_COMMIT,
        "execution_erratum_commit": "{ERRATUM_COMMIT}",
        "status": "complete_metadata_only_recovery_shard",
'''
    source = replace_once(source, old_manifest, new_manifest, "manifest-lineage")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(source, encoding="utf-8")
    print({
        "source": str(args.source),
        "output": str(args.output),
        "execution_erratum_commit": ERRATUM_COMMIT,
        "candidate_pixels_opened": False,
        "flower_colour_used": False,
        "M1_M3_scores_used": False,
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
