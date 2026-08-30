#!/usr/bin/env python3
"""Compatibility wrapper for the frozen direct Florence colour extractor.

The measurement, prompt, box-selection, and CIELAB rules remain those in
``extract_jbi_ch1_direct_florence_colour.py``.  This wrapper only makes model loading
robust to the AutoModel class rename across Transformers releases.  Candidate classes
are tried in a fixed order; the first class that accepts the frozen model is recorded in
every output row.
"""

from __future__ import annotations

import platform
from typing import Any

import torch
import transformers
from transformers import AutoProcessor

from scripts.data import extract_jbi_ch1_direct_florence_colour as base


AUTO_MODEL_CLASS_ORDER = (
    "AutoModelForCausalLM",
    "AutoModelForImageTextToText",
    "AutoModelForVision2Seq",
    "AutoModelForSeq2SeqLM",
)


def robust_load_model(model_id: str) -> tuple[Any, Any, dict[str, Any]]:
    processor = AutoProcessor.from_pretrained(model_id, trust_remote_code=True)
    failures: list[str] = []
    model = None
    selected_class = None
    selected_dtype_keyword = None
    for class_name in AUTO_MODEL_CLASS_ORDER:
        model_class = getattr(transformers, class_name, None)
        if model_class is None:
            failures.append(f"{class_name}: unavailable")
            continue
        for dtype_keyword in ("dtype", "torch_dtype"):
            kwargs = {"trust_remote_code": True, dtype_keyword: torch.float32}
            try:
                model = model_class.from_pretrained(model_id, **kwargs)
                selected_class = class_name
                selected_dtype_keyword = dtype_keyword
                break
            except Exception as exc:
                failures.append(f"{class_name}/{dtype_keyword}: {type(exc).__name__}: {exc}")
        if model is not None:
            break
    if model is None or selected_class is None:
        raise RuntimeError(
            "no fixed AutoModel candidate could load the Florence model; " + " | ".join(failures)
        )
    model.eval()
    config = getattr(model, "config", None)
    metadata = {
        "model_id": model_id,
        "model_commit_hash": getattr(config, "_commit_hash", None),
        "model_architectures": getattr(config, "architectures", None),
        "processor_class": type(processor).__name__,
        "model_class": type(model).__name__,
        "auto_model_class": selected_class,
        "dtype_keyword": selected_dtype_keyword,
        "auto_model_candidate_order": list(AUTO_MODEL_CLASS_ORDER),
        "torch_version": torch.__version__,
        "transformers_version": transformers.__version__,
        "python_version": platform.python_version(),
    }
    return processor, model, metadata


base.load_model = robust_load_model


if __name__ == "__main__":
    base.main()
