#!/usr/bin/env python3

"""Planned helpers for the solid-rotation advection example."""

from __future__ import annotations

from pathlib import Path


def write_case_initial_field(case: Path, nx: int, ny: int) -> Path:
    """Placeholder for LeVeque-style slotted-disk/cone/hump data."""
    raise NotImplementedError(
        "solid_rotation initial field generation has not been implemented yet"
    )
