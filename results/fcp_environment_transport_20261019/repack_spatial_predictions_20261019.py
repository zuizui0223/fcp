#!/usr/bin/env python3
"""Repack an existing NPZ with ZIP_LZMA, preserving exact ordered NPY bytes.

After a normal scientific rerun, invoke explicitly (same-path output is safe):
  python repack_spatial_predictions_20261019.py spatial_oof_predictions.npz \
      spatial_oof_predictions.npz

The original archive must be preserved separately for a release audit. This
helper prints packaging evidence; the release receipt is updated separately.
No fitting, scoring, array serialization, or receipt history is changed here.
"""

import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import tempfile
import zipfile

import numpy as np


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def repack(source, destination):
    source, destination = Path(source).resolve(), Path(destination).resolve()
    old_sha, old_bytes = sha256(source), source.stat().st_size
    handle, temporary_name = tempfile.mkstemp(
        prefix=destination.name + ".", suffix=".npz", dir=destination.parent
    )
    os.close(handle)
    temporary = Path(temporary_name)
    try:
        with zipfile.ZipFile(source) as src, zipfile.ZipFile(
            temporary, "w", compression=zipfile.ZIP_LZMA, allowZip64=True
        ) as dst:
            names = src.namelist()
            require(len(names) == len(set(names)), "Duplicate archive members")
            require(all(name.endswith(".npy") for name in names), "Non-NPY member")
            dst.comment = src.comment
            for info in src.infolist():
                copied = copy.copy(info)
                copied.compress_type = zipfile.ZIP_LZMA
                dst.writestr(copied, src.read(info), compress_type=zipfile.ZIP_LZMA)

        members, arrays = {}, {}
        with zipfile.ZipFile(source) as src, zipfile.ZipFile(temporary) as dst:
            require(src.namelist() == dst.namelist(), "Member order changed")
            require(dst.testzip() is None, "Archive CRC verification failed")
            for name in names:
                original = src.read(name)
                require(original == dst.read(name), "Member bytes changed: " + name)
                members[name] = hashlib.sha256(original).hexdigest()
        with np.load(source, allow_pickle=False) as src, np.load(
            temporary, allow_pickle=False
        ) as dst:
            require(src.files == dst.files, "Array names/order changed")
            for name in src.files:
                original, packed = src[name], dst[name]
                require(original.dtype == packed.dtype, "Array dtype changed")
                require(original.shape == packed.shape, "Array shape changed")
                require(original.tobytes() == packed.tobytes(), "Array bytes changed")
                equal = (
                    np.array_equal(original, packed, equal_nan=True)
                    if original.dtype.kind in "fc"
                    else np.array_equal(original, packed)
                )
                require(equal, "Array values changed")
                arrays[name] = {
                    "dtype": str(original.dtype),
                    "shape": list(original.shape),
                    "values": int(original.size),
                }
        new_bytes, new_sha = temporary.stat().st_size, sha256(temporary)
        require(new_bytes < 12_000_000, "Archive exceeds release packaging target")
        require(sha256(source) == old_sha, "Source changed during repacking")
        report = {
            "compression": "ZIP_LZMA (ZIP method 14)",
            "old_sha256": old_sha,
            "old_bytes": old_bytes,
            "new_sha256": new_sha,
            "new_bytes": new_bytes,
            "base64_characters": 4 * ((new_bytes + 2) // 3),
            "ordered_members_identical": True,
            "member_bytes_identical": True,
            "member_sha256": members,
            "np_load_allow_pickle_false_passed": True,
            "array_dtypes_shapes_bytes_values_identical": True,
            "arrays": arrays,
        }
        os.replace(temporary, destination)
        return report
    finally:
        temporary.unlink(missing_ok=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    print(json.dumps(repack(args.source, args.destination), indent=2, sort_keys=True))
