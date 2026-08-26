#!/usr/bin/env python3

"""Helpers for the diffusion-equation benchmark cases."""

from __future__ import annotations

import math
import re
from pathlib import Path

from foam_fields import write_scalar_field


def discontinuous_exact_value(x: float, y: float, time_value: float) -> float:
    """Exact solution of the discontinuous diffusion test for mu=1."""
    if time_value <= 0.0:
        return 1.0 if abs(x) <= 1.0 and abs(y) <= 1.0 else 0.0
    denominator = 2.0 * math.sqrt(time_value)
    x_part = -math.erf((-1.0 - x) / denominator) + math.erf((1.0 - x) / denominator)
    y_part = -math.erf((-1.0 - y) / denominator) + math.erf((1.0 - y) / denominator)
    return 0.25 * x_part * y_part


def gaussian_exact_value(
    x: float,
    y: float,
    time_value: float,
    mu: float = 1.0,
) -> float:
    """Exact solution of the smooth Gaussian diffusion test."""
    denominator = 1.0 + 4.0 * time_value
    return math.exp(-mu * (x * x + y * y) / denominator) / denominator


def structured_centres(
    nx: int,
    ny: int,
    domain: tuple[float, float, float, float],
) -> list[tuple[float, float]]:
    """Cell centres for a structured quad mesh."""
    xmin, xmax, ymin, ymax = domain
    dx = (xmax - xmin) / nx
    dy = (ymax - ymin) / ny
    return [
        (xmin + (i + 0.5) * dx, ymin + (j + 0.5) * dy)
        for j in range(ny)
        for i in range(nx)
    ]


def discontinuous_values(
    nx: int,
    ny: int,
    domain: tuple[float, float, float, float],
    time_value: float,
) -> list[float]:
    """Cell-centred exact values for the discontinuous diffusion case."""
    return [
        discontinuous_exact_value(x, y, time_value)
        for x, y in structured_centres(nx, ny, domain)
    ]


def gaussian_values(
    nx: int,
    ny: int,
    domain: tuple[float, float, float, float],
    time_value: float,
    mu: float = 1.0,
) -> list[float]:
    """Cell-centred exact values for the smooth Gaussian diffusion case."""
    return [
        gaussian_exact_value(x, y, time_value, mu)
        for x, y in structured_centres(nx, ny, domain)
    ]


def write_discontinuous_initial_field(
    case: Path,
    nx: int,
    ny: int,
    domain: tuple[float, float, float, float],
    field_name: str = "phi",
) -> Path:
    """Write the discontinuous initial condition to 0.orig/<field_name>."""
    boundary_types = {
        "xMin": "zeroGradient",
        "xMax": "zeroGradient",
        "yMin": "zeroGradient",
        "yMax": "zeroGradient",
        "zMin": "empty",
        "zMax": "empty",
    }
    return write_scalar_field(
        case / "0.orig" / field_name,
        discontinuous_values(nx, ny, domain, 0.0),
        boundary_types,
        object_name=field_name,
        location="0",
    )


def write_discontinuous_initial_field_from_centres(
    case: Path,
    centres: list[tuple[float, float, float]] | list[tuple[float, float]],
    field_name: str = "phi",
) -> Path:
    """Write the discontinuous initial condition using arbitrary cell centres."""
    boundary_types = {
        "xMin": "zeroGradient",
        "xMax": "zeroGradient",
        "yMin": "zeroGradient",
        "yMax": "zeroGradient",
        "zMin": "empty",
        "zMax": "empty",
    }
    values = [
        1.0 if abs(float(x)) <= 1.0 and abs(float(y)) <= 1.0 else 0.0
        for x, y, *_ in centres
    ]
    return write_scalar_field(
        case / "0.orig" / field_name,
        values,
        boundary_types,
        object_name=field_name,
        location="0",
    )


def _gaussian_boundary_block(patch_name: str, mu: float = 1.0) -> str:
    return f"""        type codedFixedValue;
        value uniform 0;
        name {patch_name}GaussianDirichlet;
        code
        #{{ 
            const scalar t = this->db().time().value();
            const scalar denom = 1.0 + 4.0*t;
            const scalar invDenom = 1.0/denom;
            const scalar mu = {mu:.16g};
            const vectorField& Cf = patch().Cf();
            scalarField values(patch().size(), 0.0);
            forAll(Cf, i)
            {{
                const scalar x = Cf[i].x();
                const scalar y = Cf[i].y();
                values[i] = invDenom*exp(-mu*(x*x + y*y)*invDenom);
            }}
            operator==(values);
        #}};"""


def write_gaussian_initial_field(
    case: Path,
    nx: int,
    ny: int,
    domain: tuple[float, float, float, float],
    field_name: str = "phi",
    mu: float = 1.0,
) -> Path:
    """Write the Gaussian diffusion initial condition to 0.orig/<field_name>."""
    boundary_types = {
        "xMin": {"__raw__": _gaussian_boundary_block("xMin", mu)},
        "xMax": {"__raw__": _gaussian_boundary_block("xMax", mu)},
        "yMin": {"__raw__": _gaussian_boundary_block("yMin", mu)},
        "yMax": {"__raw__": _gaussian_boundary_block("yMax", mu)},
        "zMin": "empty",
        "zMax": "empty",
    }
    return write_scalar_field(
        case / "0.orig" / field_name,
        gaussian_values(nx, ny, domain, 0.0, mu),
        boundary_types,
        object_name=field_name,
        location="0",
    )


def write_gaussian_initial_field_from_centres(
    case: Path,
    centres: list[tuple[float, float, float]] | list[tuple[float, float]],
    field_name: str = "phi",
    mu: float = 1.0,
) -> Path:
    """Write the Gaussian diffusion initial condition using arbitrary cell centres."""
    boundary_types = {
        "xMin": {"__raw__": _gaussian_boundary_block("xMin", mu)},
        "xMax": {"__raw__": _gaussian_boundary_block("xMax", mu)},
        "yMin": {"__raw__": _gaussian_boundary_block("yMin", mu)},
        "yMax": {"__raw__": _gaussian_boundary_block("yMax", mu)},
        "zMin": "empty",
        "zMax": "empty",
    }
    values = [
        gaussian_exact_value(float(x), float(y), 0.0, mu)
        for x, y, *_ in centres
    ]
    return write_scalar_field(
        case / "0.orig" / field_name,
        values,
        boundary_types,
        object_name=field_name,
        location="0",
    )


def parse_diffusion_solver_log(path: Path) -> list[dict[str, float]]:
    """Extract one record per completed explicit diffusion time step."""
    if not path.exists():
        raise RuntimeError(f"Solver log not found: {path}")

    time_pattern = re.compile(
        r"Time = ([^\s]+)\s+step = (\d+)\s+deltaT = ([^\s]+)"
    )
    rmin_pattern = re.compile(r"Rphi min\s*=\s*([^\s]+)")
    rmax_pattern = re.compile(r"Rphi max\s*=\s*([^\s]+)")
    phimin_pattern = re.compile(r"phi min\s*=\s*([^\s]+)")
    phimax_pattern = re.compile(r"phi max\s*=\s*([^\s]+)")
    mass_pattern = re.compile(r"mass\s*=\s*([^\s]+)")

    records: list[dict[str, float]] = []
    current: dict[str, float] | None = None
    for line in path.read_text(encoding="utf-8").splitlines():
        match = time_pattern.search(line)
        if match:
            current = {
                "time": float(match.group(1)),
                "step": float(match.group(2)),
                "deltaT": float(match.group(3)),
            }
            continue
        if current is None:
            continue
        for pattern, name in (
            (rmin_pattern, "RphiMin"),
            (rmax_pattern, "RphiMax"),
            (phimin_pattern, "minT"),
            (phimax_pattern, "maxT"),
        ):
            match = pattern.search(line)
            if match:
                current[name] = float(match.group(1))
                if name == "maxT":
                    current["amplitude"] = current["maxT"] - current.get("minT", current["maxT"])
                break
        match = mass_pattern.search(line)
        if match:
            current["mass"] = float(match.group(1))
            records.append(current)
            current = None

    if not records:
        raise RuntimeError(f"No completed diffusion time-step records found in {path}")
    return records
