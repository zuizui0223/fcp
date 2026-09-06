#!/usr/bin/env python3
"""Technical-only recovery wrapper for the frozen scale-moderation replication.

The authorized first execution stopped before reading the biological pool because
Python 3.12 dataclasses require a dynamically imported module to be registered in
sys.modules.  This wrapper changes only that import plumbing, then calls the exact
frozen replication main function unchanged.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TARGET = ROOT / "scripts/analysis/run_global_rgfca_worldclim_scale_moderation_replication.py"
WORLDCLIM = ROOT / "scripts/analysis/run_global_rgfca_worldclim_edge_secondary.py"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def fixed_load_worldclim_module():
    return load_module("rgfca_worldclim_secondary", WORLDCLIM)


def main() -> int:
    replication = load_module("rgfca_scale_replication", TARGET)
    replication.load_worldclim_module = fixed_load_worldclim_module
    return int(replication.main())


if __name__ == "__main__":
    raise SystemExit(main())
