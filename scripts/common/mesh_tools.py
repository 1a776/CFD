#!/usr/bin/env python3

"""Mesh dictionary helpers used by the study driver."""

from __future__ import annotations

import re
from pathlib import Path

from advection_tools import mesh_resolution


BLOCK_PATTERN = re.compile(
    r"(hex\s+\(0\s+1\s+2\s+3\s+4\s+5\s+6\s+7\)\s+)"
    r"\(\s*\d+\s+\d+\s+1\s*\)"
)


def patch_block_mesh_resolution(case: Path, resolution: int) -> None:
    """Set the structured quad mesh to N by N by 1 cells."""
    path = case / "system" / "blockMeshDict"
    text = path.read_text(encoding="utf-8")
    updated, count = BLOCK_PATTERN.subn(
        rf"\g<1>({resolution} {resolution} 1)", text, count=1
    )
    if count != 1:
        raise RuntimeError(f"Cannot patch blockMeshDict for {case}")
    path.write_text(updated, encoding="utf-8")


__all__ = ["mesh_resolution", "patch_block_mesh_resolution"]
