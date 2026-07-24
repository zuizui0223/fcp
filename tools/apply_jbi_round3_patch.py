#!/usr/bin/env python3
from pathlib import Path
parts = Path(__file__).resolve().parent / "jbi_patch_parts"
source = "".join((parts / f"part{i:02d}.pyfrag").read_text(encoding="utf-8") for i in range(5))
exec(compile(source, str(parts / "combined.py"), "exec"))
