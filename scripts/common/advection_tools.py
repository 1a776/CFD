#!/usr/bin/env python3

"""Shared utilities for the student linear-advection cases."""

from __future__ import annotations

import math
import re
from pathlib import Path


FIELD_PATTERN = re.compile(
    r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;",
    re.S,
)
BLOCK_PATTERN = re.compile(
    r"hex\s+\(0\s+1\s+2\s+3\s+4\s+5\s+6\s+7\)\s+"
    r"\(\s*(\d+)\s+(\d+)\s+1\s*\)"
)


def read_scalar_field(path: Path) -> list[float]:
    """Read a nonuniform OpenFOAM scalar field."""
    text = path.read_text(encoding="utf-8")
    match = FIELD_PATTERN.search(text)
    if not match:
        raise RuntimeError(f"Cannot parse scalar field: {path}")

    count = int(match.group(1))
    values = [float(item) for item in re.findall(r"[-+0-9.eE]+", match.group(2))]
    if len(values) != count:
        raise RuntimeError(
            f"{path}: expected {count} scalar values, got {len(values)}"
        )
    return values


def numeric_times(case: Path) -> list[tuple[float, Path]]:
    """Return numeric OpenFOAM time directories in ascending order."""
    times: list[tuple[float, Path]] = []
    for path in case.iterdir():
        if not path.is_dir():
            continue
        try:
            times.append((float(path.name), path))
        except ValueError:
            continue
    return sorted(times)


def latest_time(case: Path) -> tuple[float, Path]:
    times = numeric_times(case)
    if not times:
        raise RuntimeError(f"No numeric time directories found in {case}")
    return times[-1]


def mesh_resolution(case: Path) -> tuple[int, int]:
    """Read the x/y cell counts from system/blockMeshDict."""
    path = case / "system" / "blockMeshDict"
    match = BLOCK_PATTERN.search(path.read_text(encoding="utf-8"))
    if not match:
        raise RuntimeError(f"Cannot read mesh resolution from {path}")
    return int(match.group(1)), int(match.group(2))


def control_value(case: Path, name: str, default: float) -> float:
    path = case / "system" / "controlDict"
    match = re.search(rf"^\s*{re.escape(name)}\s+([^;]+);", path.read_text(encoding="utf-8"), re.M)
    if not match:
        return default
    return float(match.group(1).strip())


def exact_values(nx: int, ny: int, time_value: float) -> list[float]:
    """Return sin(2*pi*(x+y-2*t)) at cell centres."""
    values: list[float] = []
    for j in range(ny):
        y = (j + 0.5) / ny
        for i in range(nx):
            x = (i + 0.5) / nx
            values.append(math.sin(2.0 * math.pi * (x + y - 2.0 * time_value)))
    return values


def write_initial_field(case: Path, nx: int, ny: int) -> Path:
    """Write T(x,y,0)=sin(2*pi*(x+y)) to case/0.orig/T."""
    output = case / "0.orig" / "T"
    output.parent.mkdir(parents=True, exist_ok=True)

    values = exact_values(nx, ny, 0.0)
    body = "\n".join(f"    {value:.16e}" for value in values)
    output.write_text(
        f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /    O peration     |
    \\\\  /    A nd           |
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       volScalarField;
    location    \"0\";
    object      T;
}}

dimensions      [0 0 0 0 0 0 0];

internalField   nonuniform List<scalar>
{len(values)}
(
{body}
);

boundaryField
{{
    xMin {{ type cyclic; }}
    xMax {{ type cyclic; }}
    yMin {{ type cyclic; }}
    yMax {{ type cyclic; }}
    zMin {{ type empty; }}
    zMax {{ type empty; }}
}}
""",
        encoding="utf-8",
    )
    return output


def parse_solver_log(path: Path) -> list[dict[str, float]]:
    """Extract one record per completed time step from a solver log."""
    if not path.exists():
        raise RuntimeError(f"Solver log not found: {path}")

    time_pattern = re.compile(
        r"Time = ([^\s]+)\s+step = (\d+)\s+"
        r"deltaT = ([^\s]+)\s+maxCo = ([^\s]+)"
    )
    residual_pattern = re.compile(r"residual integral\s*=\s*([^\s]+)")
    min_pattern = re.compile(r"T min\s*=\s*([^\s]+)")
    max_pattern = re.compile(r"T max\s*=\s*([^\s]+)")

    records: list[dict[str, float]] = []
    current: dict[str, float] | None = None
    for line in path.read_text(encoding="utf-8").splitlines():
        match = time_pattern.search(line)
        if match:
            current = {
                "time": float(match.group(1)),
                "step": float(match.group(2)),
                "deltaT": float(match.group(3)),
                "maxCo": float(match.group(4)),
            }
            continue
        if current is None:
            continue

        match = residual_pattern.search(line)
        if match:
            current["residualIntegral"] = float(match.group(1))
            continue
        match = min_pattern.search(line)
        if match:
            current["minT"] = float(match.group(1))
            continue
        match = max_pattern.search(line)
        if match:
            current["maxT"] = float(match.group(1))
            current["amplitude"] = 0.5 * (
                current["maxT"] - current.get("minT", current["maxT"])
            )
            records.append(current)
            current = None

    if not records:
        raise RuntimeError(f"No completed time-step records found in {path}")
    return records


def normalized_errors(
    numerical: list[float], exact: list[float], cell_volume: float
) -> tuple[float, float, float]:
    """Compute volume-weighted normalized L1, L2 and Linf errors."""
    if len(numerical) != len(exact):
        raise RuntimeError("Numerical and exact fields have different sizes")
    error = [value - reference for value, reference in zip(numerical, exact)]
    denominator_l1 = sum(abs(reference) for reference in exact) * cell_volume
    denominator_l2 = sum(reference * reference for reference in exact) * cell_volume
    denominator_linf = max(abs(reference) for reference in exact)
    if denominator_l1 <= 0.0 or denominator_l2 <= 0.0 or denominator_linf <= 0.0:
        raise RuntimeError("Exact field norm is zero")

    l1 = sum(abs(value) for value in error) * cell_volume / denominator_l1
    l2 = math.sqrt(sum(value * value for value in error) * cell_volume / denominator_l2)
    linf = max(abs(value) for value in error) / denominator_linf
    return l1, l2, linf
