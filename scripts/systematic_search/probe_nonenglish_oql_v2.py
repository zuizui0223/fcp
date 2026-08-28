#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json, os, time, urllib.parse, urllib.request
from pathlib import Path
from typing import Any

USER_AGENT = "fcp-jbi-nonenglish-probe/1.0 (https://github.com/zuizui0223/fcp)"


def request(url: str, timeout: int, retries: int) -> dict[str, Any]:
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                data = json.load(r)
            if not isinstance(data, dict): raise RuntimeError("non-object response")
            return data
        except Exception as e:  # noqa: BLE001
            last = e
            if attempt + 1 < retries: time.sleep(2 ** attempt)
    raise RuntimeError(url) from last


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--config", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--summary-out", required=True)
    p.add_argument("--api-key-env", default="OPENALEX_API_KEY")
    p.add_argument("--timeout", type=int, default=30)
    p.add_argument("--retries", type=int, default=3)
    a = p.parse_args()
    cfg = json.loads(Path(a.config).read_text(encoding="utf-8"))
    api_key = os.environ.get(a.api_key_env, "").strip()
    rows = []
    for item in cfg.get("queries", []):
        params = {"oql": item["oql"], "per-page": 1}
        if api_key: params["api_key"] = api_key
        data = request("https://api.openalex.org/?" + urllib.parse.urlencode(params), a.timeout, a.retries)
        meta = data.get("meta") if isinstance(data.get("meta"), dict) else {}
        count = int(meta.get("count") or 0)
        row = {"query_id": item["id"], "language": item["language"], "count": count, "oql": item["oql"]}
        rows.append(row); print(row, flush=True); time.sleep(0.05)
    if len(rows) != 5: raise SystemExit(f"Expected 5 probes, got {len(rows)}")
    out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    summary = {
        "status": "complete",
        "counts": {r["query_id"]: r["count"] for r in rows},
        "total_membership_upper_bound": sum(r["count"] for r in rows),
        "guard": "Counts are scope diagnostics only; broader blocks require relevance review before adoption."
    }
    Path(a.summary_out).write_text(json.dumps(summary, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)

if __name__ == "__main__": main()
