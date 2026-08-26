#!/usr/bin/env python3

"""Small OpenFOAM field readers and writers used by case preparation."""

from __future__ import annotations

import re
import json
from pathlib import Path
from typing import Any


_SCALAR_FIELD_PATTERN = re.compile(
    r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;",
    re.S,
)
_VECTOR_FIELD_PATTERN = re.compile(
    r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;",
    re.S,
)
_VECTOR_PATTERN = re.compile(
    r"\(\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\)"
)
_UNIFORM_VECTOR_PATTERN = re.compile(
    r"(internalField\s+uniform\s+)"
    r"\(\s*[-+0-9.eE]+\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s*\)"
)


def read_scalar_field(path: Path) -> list[float]:
    """Read an OpenFOAM nonuniform scalar field."""
    text = path.read_text(encoding="utf-8")
    match = _SCALAR_FIELD_PATTERN.search(text)
    if not match:
        raise RuntimeError(f"Cannot parse nonuniform scalar field: {path}")

    count = int(match.group(1))
    values = [float(item) for item in re.findall(r"[-+0-9.eE]+", match.group(2))]
    if len(values) != count:
        raise RuntimeError(f"{path}: expected {count} scalar values, got {len(values)}")
    return values


def read_vector_field(path: Path) -> list[tuple[float, float, float]]:
    """Read an OpenFOAM nonuniform vector field."""
    text = path.read_text(encoding="utf-8")
    match = _VECTOR_FIELD_PATTERN.search(text)
    if not match:
        raise RuntimeError(f"Cannot parse nonuniform vector field: {path}")

    count = int(match.group(1))
    values = [
        (float(x), float(y), float(z))
        for x, y, z in _VECTOR_PATTERN.findall(match.group(2))
    ]
    if len(values) != count:
        raise RuntimeError(f"{path}: expected {count} vector values, got {len(values)}")
    return values


def _boundary_entry(patch: str, spec: str | dict[str, Any]) -> str:
    """Format one OpenFOAM boundaryField entry."""
    if isinstance(spec, dict):
        raw_body = spec.get("__raw__")
        if raw_body is not None:
            body = str(raw_body).rstrip()
        else:
            lines = [f"        {key} {value};" for key, value in spec.items()]
            body = "\n".join(lines)
    else:
        body = f"        type {spec};"
    return f"    {patch}\n    {{\n{body}\n    }}"


def write_scalar_field(
    path: Path,
    values: list[float],
    boundary_types: dict[str, str | dict[str, Any]],
    object_name: str = "T",
    location: str = "0",
) -> Path:
    """Write a nonuniform OpenFOAM volScalarField."""
    body = "\n".join(f"    {value:.16e}" for value in values)
    boundary_body = "\n".join(
        _boundary_entry(patch, patch_type)
        for patch, patch_type in boundary_types.items()
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /    O peration     | M anipulation
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       volScalarField;
    location    "{location}";
    object      {object_name};
}}

dimensions      [0 0 0 0 0 0 0];

internalField   nonuniform List<scalar>
{len(values)}
(
{body}
);

boundaryField
{{
{boundary_body}
}}
""",
        encoding="utf-8",
    )
    return path


def write_vector_field(
    path: Path,
    values: list[tuple[float, float, float]],
    boundary_types: dict[str, str | dict[str, Any]],
    object_name: str = "U",
    location: str = "0",
) -> Path:
    """Write a nonuniform OpenFOAM volVectorField."""
    body = "\n".join(
        f"    ({x:.16e} {y:.16e} {z:.16e})" for x, y, z in values
    )
    boundary_body = "\n".join(
        _boundary_entry(patch, patch_type)
        for patch, patch_type in boundary_types.items()
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /    O peration     | M anipulation
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       volVectorField;
    location    "{location}";
    object      {object_name};
}}

dimensions      [0 1 -1 0 0 0 0];

internalField   nonuniform List<vector>
{len(values)}
(
{body}
);

boundaryField
{{
{boundary_body}
}}
""",
        encoding="utf-8",
    )
    return path


def patch_uniform_vector_field(
    path: Path,
    value: tuple[float, float, float] | list[float],
) -> Path:
    """Replace a volVectorField's uniform internal value in-place."""
    if len(value) != 3:
        raise ValueError("A vector field value must contain exactly three components")

    vector_text = "({:g} {:g} {:g})".format(
        float(value[0]),
        float(value[1]),
        float(value[2]),
    )
    text = path.read_text(encoding="utf-8")
    updated, count = _UNIFORM_VECTOR_PATTERN.subn(
        rf"\g<1>{vector_text}",
        text,
        count=1,
    )
    if count != 1:
        raise RuntimeError(f"Cannot find uniform internalField in vector field: {path}")
    path.write_text(updated, encoding="utf-8")
    return path


def read_cell_geometry(case: Path) -> tuple[list[tuple[float, float, float]], list[float]]:
    """Read cell centres C and cell volumes Vc written by foamPostProcess."""
    centres = read_vector_field(case / "constant" / "C")
    volumes = read_scalar_field(case / "constant" / "Vc")
    if len(centres) != len(volumes):
        raise RuntimeError(
            f"Cell geometry size mismatch in {case}: "
            f"C={len(centres)}, Vc={len(volumes)}"
        )
    return centres, volumes


def read_tri_geometry_metadata(path: Path) -> dict[str, object]:
    """Read the auxiliary geometry JSON emitted by the Gmsh helper."""
    data = json.loads(path.read_text(encoding="utf-8"))
    if "nodes" not in data or "triangles" not in data:
        raise RuntimeError(f"Missing nodes/triangles in {path}")
    return data


__all__ = [
    "read_cell_geometry",
    "read_scalar_field",
    "read_vector_field",
    "patch_uniform_vector_field",
    "read_tri_geometry_metadata",
    "write_scalar_field",
    "write_vector_field",
]
