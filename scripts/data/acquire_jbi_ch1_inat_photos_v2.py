#!/usr/bin/env python3
"""Run Chapter 1 iNaturalist acquisition with frozen repository→iNat taxon overrides."""

from __future__ import annotations

import json
from pathlib import Path
import sys

import scripts.data.acquire_jbi_ch1_inat_photos as core


def config_path_from_argv() -> Path:
    for index, arg in enumerate(sys.argv[1:], start=1):
        if arg == "--config" and index + 1 < len(sys.argv):
            return Path(sys.argv[index + 1])
        if arg.startswith("--config="):
            return Path(arg.split("=", 1)[1])
    return Path("docs/supporting/jbi_ch1_inat_acquisition_v1.json")


config = json.loads(config_path_from_argv().read_text(encoding="utf-8"))
overrides = config.get("taxon_overrides", {})
_original_exact_taxon = core.exact_taxon


def exact_taxon_with_override(base_url: str, name: str, *, pause: float):
    override = overrides.get(name)
    if override:
        required = {"id", "name", "rank"}
        missing = required - set(override)
        if missing:
            raise RuntimeError(f"taxon override for {name!r} missing {sorted(missing)}")
        print(
            f"Taxon override: {name} -> {override['name']} (iNat {override['id']}, {override['rank']})",
            flush=True,
        )
        return {
            "id": int(override["id"]),
            "name": str(override["name"]),
            "rank": str(override["rank"]),
        }
    return _original_exact_taxon(base_url, name, pause=pause)


core.exact_taxon = exact_taxon_with_override


if __name__ == "__main__":
    raise SystemExit(core.main())
