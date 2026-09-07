#!/usr/bin/env python3
"""Stage-separated RGFCA reserve audit, blind measurement and final reassembly."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

from fcp_pipeline.rgfca_reserve_replication import (
    ROOT, AUDIT, CONTRACT, PREFIX, PROTOCOL, audit_reserve, build_firewall, reassemble, sha,
)


def write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8", newline="\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["audit", "firewall", "reassemble"])
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--firewall-dir", type=Path)
    parser.add_argument("--results-dir", type=Path)
    args = parser.parse_args()
    if args.stage == "audit":
        _, geometry, report = audit_reserve()
        out = ROOT / "data/frozen/rgfca_reserve_geometry_audit_v1.csv"
        geometry.to_csv(out, index=False, lineterminator="\n")
        report["geometry_table_sha256"] = sha(out.read_bytes())
        if AUDIT.exists() and json.loads(AUDIT.read_text(encoding="utf-8")) != report:
            raise RuntimeError("frozen reserve audit changed; do not overwrite")
        write_json(AUDIT, report)
        print(json.dumps(report, indent=2))
    elif args.stage == "firewall":
        if args.output_dir is None:
            parser.error("firewall requires --output-dir")
        auth_path = ROOT / "docs/supporting/rgfca_reserve_measurement_authorization_v1.json"
        auth = json.loads(auth_path.read_text(encoding="utf-8"))
        if auth.get("status") != "authorize_one_reserve_measurement_after_preflight":
            raise RuntimeError("reserve pixel measurement is not authorized")
        if os.environ.get("GITHUB_RUN_ATTEMPT") != "1":
            raise RuntimeError("a rerun cannot silently reopen reserve pixels")
        if os.environ.get("GITHUB_REF_NAME") != auth["branch"]:
            raise RuntimeError("measurement authorization branch mismatch")
        if sha(CONTRACT.read_bytes()) != auth["contract_sha256"] or sha(AUDIT.read_bytes()) != auth["metadata_audit_sha256"]:
            raise RuntimeError("authorized contract or metadata audit hash drift")
        for path, expected in auth["code_sha256"].items():
            if sha((ROOT / path).read_bytes()) != expected:
                raise RuntimeError(f"authorized code hash drift: {path}")
        if (ROOT / f"data/derived/{PREFIX}_measured_photos_v1.csv").exists():
            raise RuntimeError("reserve measurement already exists; no repeat measurement")
        reserve, _, _ = audit_reserve()
        frozen = json.loads(AUDIT.read_text(encoding="utf-8"))
        contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
        result = build_firewall(reserve, frozen, args.output_dir, contract)
        print(json.dumps(result, indent=2))
    else:
        if args.firewall_dir is None or args.results_dir is None:
            parser.error("reassemble requires --firewall-dir and --results-dir")
        joined, support, result = reassemble(args.firewall_dir, args.results_dir)
        paths = [ROOT / f"data/derived/{PREFIX}_{suffix}_v1.csv" for suffix in ("measured_photos", "species_support")]
        manifest = ROOT / f"docs/supporting/{PREFIX}_measurement_result_v1.json"
        if any(p.exists() for p in [*paths, manifest]):
            raise RuntimeError("replication result already exists; do not overwrite")
        for path, frame in zip(paths, (joined, support)):
            frame.to_csv(path, index=False, lineterminator="\n")
        result["files_sha256"] = {p.relative_to(ROOT).as_posix(): sha(p.read_bytes()) for p in paths}
        result["github_run_id"] = os.environ.get("GITHUB_RUN_ID")
        write_json(manifest, result)
        print(json.dumps({k: v for k, v in result.items() if k != "partition_hashes"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
