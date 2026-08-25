#!/usr/bin/env python3

"""Solid-rotation velocity and profile fields for advection tests."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

from foam_fields import write_scalar_field, write_vector_field


SCALAR_BOUNDARY = {
    "xMin": {"type": "fixedValue", "value": "uniform 0"},
    "xMax": {"type": "fixedValue", "value": "uniform 0"},
    "yMin": {"type": "fixedValue", "value": "uniform 0"},
    "yMax": {"type": "fixedValue", "value": "uniform 0"},
    "zMin": "empty",
    "zMax": "empty",
}

VELOCITY_BOUNDARY = {
    "xMin": "zeroGradient",
    "xMax": "zeroGradient",
    "yMin": "zeroGradient",
    "yMax": "zeroGradient",
    "zMin": "empty",
    "zMax": "empty",
}


def _profile_parameters(config: dict[str, Any]) -> dict[str, float]:
    def point(name: str, default: tuple[float, float]) -> tuple[float, float]:
        value = config.get(name, default)
        if not isinstance(value, list | tuple) or len(value) != 2:
            raise RuntimeError(f"initialProfile.{name} must contain two numbers")
        return (float(value[0]), float(value[1]))

    return {
        "radius": float(config.get("radius", 0.15)),
        "disk_center": point("diskCenter", (0.5, 0.75)),
        "cone_center": point("coneCenter", (0.5, 0.25)),
        "hump_center": point("humpCenter", (0.25, 0.5)),
        "slot_half_width": float(config.get("slotHalfWidth", 0.025)),
        "slot_top_y": float(config.get("slotTopY", 0.85)),
    }


def solid_rotation_velocity(
    x: float,
    y: float,
    center: tuple[float, float] = (0.5, 0.5),
    angular_velocity: float = 1.0,
) -> tuple[float, float, float]:
    """Evaluate u=(omega*(yc-y), omega*(x-xc), 0)."""
    xc, yc = center
    return (
        angular_velocity * (yc - y),
        angular_velocity * (x - xc),
        0.0,
    )


def slotted_disk(
    x: float,
    y: float,
    radius: float = 0.15,
    center: tuple[float, float] = (0.5, 0.75),
    slot_half_width: float = 0.025,
    slot_top_y: float = 0.85,
) -> float:
    """Evaluate the slotted disk."""
    xc, yc = center
    distance = math.hypot(x - xc, y - yc)
    outside_slot = abs(x - xc) >= slot_half_width or y >= slot_top_y
    return 1.0 if distance <= radius and outside_slot else 0.0


def cone(
    x: float,
    y: float,
    radius: float = 0.15,
    center: tuple[float, float] = (0.5, 0.25),
) -> float:
    """Evaluate the cone."""
    xc, yc = center
    distance = math.hypot(x - xc, y - yc)
    if distance > radius:
        return 0.0
    return 1.0 - distance / radius


def cosine_hump(
    x: float,
    y: float,
    radius: float = 0.15,
    center: tuple[float, float] = (0.25, 0.5),
) -> float:
    """Evaluate the smooth cosine hump."""
    xc, yc = center
    distance = math.hypot(x - xc, y - yc)
    if distance > radius:
        return 0.0
    return 0.25 * (1.0 + math.cos(math.pi * distance / radius))


def solid_rotation_profile(x: float, y: float, config: dict[str, Any] | None = None) -> float:
    """Evaluate the slotted-disk + cone + cosine-hump initial profile."""
    params = _profile_parameters(config or {})
    return (
        slotted_disk(
            x,
            y,
            params["radius"],
            params["disk_center"],
            params["slot_half_width"],
            params["slot_top_y"],
        )
        + cone(x, y, params["radius"], params["cone_center"])
        + cosine_hump(x, y, params["radius"], params["hump_center"])
    )


def _cell_centres(
    nx: int,
    ny: int,
    domain: tuple[float, float, float, float],
) -> list[tuple[float, float]]:
    xmin, xmax, ymin, ymax = domain
    dx = (xmax - xmin) / nx
    dy = (ymax - ymin) / ny
    return [
        (xmin + (i + 0.5) * dx, ymin + (j + 0.5) * dy)
        for j in range(ny)
        for i in range(nx)
    ]


def write_case_initial_field(
    case: Path,
    nx: int,
    ny: int,
    domain: tuple[float, float, float, float] = (0.0, 1.0, 0.0, 1.0),
    profile_config: dict[str, Any] | None = None,
) -> Path:
    """Write the solid-rotation composite profile into 0.orig/T."""
    values = [
        solid_rotation_profile(x, y, profile_config)
        for x, y in _cell_centres(nx, ny, domain)
    ]
    return write_scalar_field(case / "0.orig" / "T", values, SCALAR_BOUNDARY)


def write_case_initial_field_from_centres(
    case: Path,
    centres: list[tuple[float, float, float]],
    profile_config: dict[str, Any] | None = None,
) -> Path:
    """Write the profile at the actual cell centres of an unstructured mesh."""
    values = [
        solid_rotation_profile(x, y, profile_config)
        for x, y, _ in centres
    ]
    return write_scalar_field(case / "0.orig" / "T", values, SCALAR_BOUNDARY)


def write_case_velocity_field(
    case: Path,
    nx: int,
    ny: int,
    domain: tuple[float, float, float, float] = (0.0, 1.0, 0.0, 1.0),
    velocity_config: dict[str, Any] | None = None,
) -> Path:
    """Write the cell-centred solid-rotation velocity into 0.orig/U."""
    velocity_config = velocity_config or {}
    center_value = velocity_config.get("center", [0.5, 0.5])
    if not isinstance(center_value, list | tuple) or len(center_value) != 2:
        raise RuntimeError("velocityModel.center must contain two numbers")
    center = (float(center_value[0]), float(center_value[1]))
    angular_velocity = float(velocity_config.get("angularVelocity", 1.0))
    values = [
        solid_rotation_velocity(x, y, center, angular_velocity)
        for x, y in _cell_centres(nx, ny, domain)
    ]
    return write_vector_field(case / "0.orig" / "U", values, VELOCITY_BOUNDARY)


def write_case_velocity_field_from_centres(
    case: Path,
    centres: list[tuple[float, float, float]],
    velocity_config: dict[str, Any] | None = None,
) -> Path:
    """Write the rotation velocity at the actual cell centres."""
    velocity_config = velocity_config or {}
    center_value = velocity_config.get("center", [0.5, 0.5])
    if not isinstance(center_value, list | tuple) or len(center_value) != 2:
        raise RuntimeError("velocityModel.center must contain two numbers")
    center = (float(center_value[0]), float(center_value[1]))
    angular_velocity = float(velocity_config.get("angularVelocity", 1.0))
    values = [
        solid_rotation_velocity(x, y, center, angular_velocity)
        for x, y, _ in centres
    ]
    return write_vector_field(case / "0.orig" / "U", values, VELOCITY_BOUNDARY)


__all__ = [
    "solid_rotation_velocity",
    "solid_rotation_profile",
    "write_case_initial_field",
    "write_case_initial_field_from_centres",
    "write_case_velocity_field",
    "write_case_velocity_field_from_centres",
]
