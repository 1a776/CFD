#!/usr/bin/env python3

"""Sine-wave initial and exact fields for the advection equation."""

from __future__ import annotations

from pathlib import Path

from advection_tools import exact_values, write_initial_field


def write_case_initial_field(case: Path, nx: int, ny: int) -> Path:
    """Write phi(x,y,0)=sin(2*pi*(x+y)) into 0.orig/T."""
    return write_initial_field(case, nx, ny)


__all__ = ["exact_values", "write_case_initial_field"]
