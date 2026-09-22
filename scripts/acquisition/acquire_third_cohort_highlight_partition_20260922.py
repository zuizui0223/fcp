#!/usr/bin/env python3
from __future__ import annotations

import argparse
from io import BytesIO
import hashlib
from pathlib import Path
import time
from urllib.request import Request, urlopen

import pandas as pd
from PIL import Image

USER_AGENT = "fcp-third-cohort-highlight-validity/1.0 (github.com/zuizui0223/fcp)"


def _download(url: str, retries: int = 4, timeout: float = 45.0) -> tuple[bytes | None, str]:
    last = ""
    for attempt in range(retries + 1):
        try:
            req = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "image/*"})
            with urlopen(req, timeout=timeout) as response:
                payload = response.read()
            if not payload:
                raise RuntimeError("empty_image_payload")
            with Image.open(BytesIO(payload)) as im:
                im.verify()
            return payload, ""
        except Exception as exc:
            last = f"{type(exc).__name__}:{str(exc)[:240]}"
            if attempt < retries:
                time.sleep(min(16.0, 2.0 ** attempt))
    return None, last or "unknown_acquisition_failure"


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--worker-manifest", type=Path, required=True)
    p.add_argument("--acquisition-key", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()

    worker = pd.read_csv(args.worker_manifest, dtype=str).fillna("")
    key = pd.read_csv(args.acquisition_key, dtype=str).fillna("")
    if list(worker.columns) != ["measurement_id", "image_filename", "photo_license"]:
        raise RuntimeError("worker manifest schema drift")
    if list(key.columns) != ["measurement_id", "image_filename", "photo_url_large", "photo_license"]:
        raise RuntimeError("acquisition key schema drift")
    if worker["measurement_id"].nunique() != len(worker):
        raise RuntimeError("worker IDs are not unique")
    merged = worker.merge(
        key, on=["measurement_id", "image_filename", "photo_license"],
        how="inner", validate="one_to_one"
    )
    if len(merged) != len(worker):
        raise RuntimeError("worker/acquisition key mismatch")

    out = args.output_dir
    images = out / "images"
    images.mkdir(parents=True, exist_ok=True)
    receipts = []
    success = []
    for i, row in enumerate(merged.itertuples(index=False), start=1):
        payload, error = _download(str(row.photo_url_large))
        if payload is None:
            receipts.append({
                "measurement_id": row.measurement_id,
                "acquisition_status": "image_acquisition_failed",
                "reacquired_image_sha256": "",
                "image_bytes": 0,
                "failure_reason": error,
            })
        else:
            digest = hashlib.sha256(payload).hexdigest()
            (images / str(row.image_filename)).write_bytes(payload)
            success.append({
                "measurement_id": row.measurement_id,
                "image_filename": row.image_filename,
                "photo_license": row.photo_license,
            })
            receipts.append({
                "measurement_id": row.measurement_id,
                "acquisition_status": "acquired_and_decode_verified",
                "reacquired_image_sha256": digest,
                "image_bytes": len(payload),
                "failure_reason": "",
            })
        if i % 25 == 0 or i == len(merged):
            print(f"acquired_or_failed={i}/{len(merged)}", flush=True)

    pd.DataFrame(success, columns=["measurement_id","image_filename","photo_license"]).to_csv(
        out / "successful_worker_manifest.csv", index=False, lineterminator="\n"
    )
    pd.DataFrame(receipts).to_csv(
        out / "acquisition_receipt.csv", index=False, lineterminator="\n"
    )


if __name__ == "__main__":
    main()
