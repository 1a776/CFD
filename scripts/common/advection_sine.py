#!/usr/bin/env python3

"""Sine-wave initial and exact fields for the advection equation."""

from __future__ import annotations

import math
from pathlib import Path

from advection_tools import exact_values, write_initial_field
from foam_fields import write_scalar_field


def exact_values_at_centres(
    centres: list[tuple[float, float, float]],
    time_value: float,
    velocity: tuple[float, float, float] = (1.0, 1.0, 0.0),
) -> list[float]:
    """Evaluate the translated sine wave at arbitrary cell centres."""
    phase_speed = float(velocity[0]) + float(velocity[1])
    return [
        math.sin(2.0 * math.pi * (x + y - phase_speed * time_value))
        for x, y, _ in centres
    ]


def write_case_initial_field(
    case: Path,
    nx: int,
    ny: int,
    velocity: tuple[float, float, float] = (1.0, 1.0, 0.0),
) -> Path:
    """Write phi(x,y,0)=sin(2*pi*(x+y)) into 0.orig/T."""
    return write_initial_field(case, nx, ny, velocity)


def write_case_initial_field_from_centres(
    case: Path,
    centres: list[tuple[float, float, float]],
    velocity: tuple[float, float, float] = (1.0, 1.0, 0.0),
) -> Path:
    """Write the sine-wave initial field using actual OpenFOAM cell centres."""
    values = exact_values_at_centres(centres, 0.0, velocity)
    return write_scalar_field(
        case / "0.orig" / "T",
        values,
        {
            "xMin": "cyclic",
            "xMax": "cyclic",
            "yMin": "cyclic",
            "yMax": "cyclic",
            "zMin": "empty",
            "zMax": "empty",
        },
        object_name="T",
        location="0",
    )


__all__ = [
    "exact_values",
    "exact_values_at_centres",
    "write_case_initial_field",
    "write_case_initial_field_from_centres",
]
