#!/usr/bin/env python3
"""Path-independent entrypoint for the symmetric direct Florence colour extractor.

Running a file by path places ``scripts/data`` rather than the repository root at
``sys.path[0]``.  This entrypoint inserts the resolved repository root before importing
the frozen v2 compatibility wrapper.  It changes no model, prompt, ROI, colour, split,
or inferential rule.
"""

from __future__ import annotations

import sys
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

# Importing v2 installs the fixed-order robust model loader into the base implementation.
from scripts.data import extract_jbi_ch1_direct_florence_colour_v2 as _compat  # noqa: E402,F401
from scripts.data import extract_jbi_ch1_direct_florence_colour as _base  # noqa: E402


if __name__ == "__main__":
    _base.main()
