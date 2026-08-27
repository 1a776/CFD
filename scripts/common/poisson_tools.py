#!/usr/bin/env python3

"""Field writers for the manufactured Poisson benchmark."""

from __future__ import annotations

import math
from pathlib import Path

from foam_fields import write_scalar_field


def poisson_exact_value(x: float, y: float) -> float:
    """Return phi_exact = cos(pi*x)*cos(pi*y)."""
    return math.cos(math.pi * x) * math.cos(math.pi * y)


def poisson_source_value(x: float, y: float) -> float:
    """Return omega = laplacian(phi_exact)."""
    return -2.0 * math.pi**2 * poisson_exact_value(x, y)


def _dirichlet_boundary_block(patch_name: str) -> str:
    """Return a codedFixedValue block evaluated at boundary face centres."""
    return f"""        // Dirichlet condition from the manufactured exact solution:
        //     phi_boundary(x,y) = cos(pi*x)*cos(pi*y)
        type codedFixedValue;
        value uniform 0;
        name {patch_name}PoissonDirichlet;
        code
        #{{ 
            const vectorField& Cf = patch().Cf();
            scalarField values(patch().size(), 0.0);
            forAll(Cf, i)
            {{
                const scalar x = Cf[i].x();
                const scalar y = Cf[i].y();
                values[i] = cos(constant::mathematical::pi*x)
                          * cos(constant::mathematical::pi*y);
            }}
            operator==(values);
        #}};"""


def _boundary_types() -> dict[str, str | dict[str, str]]:
    return {
        "xMin": {"__raw__": _dirichlet_boundary_block("xMin")},
        "xMax": {"__raw__": _dirichlet_boundary_block("xMax")},
        "yMin": {"__raw__": _dirichlet_boundary_block("yMin")},
        "yMax": {"__raw__": _dirichlet_boundary_block("yMax")},
        "zMin": "empty",
        "zMax": "empty",
    }


def write_poisson_solution_field(
    case: Path,
    centres: list[tuple[float, float, float]] | list[tuple[float, float]],
    field_name: str = "phi",
) -> Path:
    """Write a zero initial guess with exact Dirichlet boundary values."""
    values = [0.0 for _ in centres]
    return write_scalar_field(
        case / "0.orig" / field_name,
        values,
        _boundary_types(),
        object_name=field_name,
        location="0",
        dimensions="[0 0 0 0 0 0 0]",
        header_comment=(
            "// Unknown field for the Poisson equation.\n"
            "// The internalField is a zero initial guess for the linear solver.\n"
            "// The boundaryField is the exact Dirichlet condition:\n"
            "//     phi = cos(pi*x)*cos(pi*y) on the boundary.\n"
        ),
    )


def write_poisson_source_field(
    case: Path,
    centres: list[tuple[float, float, float]] | list[tuple[float, float]],
    field_name: str = "omega",
) -> Path:
    """Write omega at cell centres for the manufactured solution."""
    values = [
        poisson_source_value(float(x), float(y))
        for x, y, *_ in centres
    ]
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
        values,
        boundary_types,
        object_name=field_name,
        location="0",
        dimensions="[0 -2 0 0 0 0 0]",
        header_comment=(
            "// Known source field for the manufactured Poisson solution.\n"
            "//     omega(x,y) = -2*pi^2*cos(pi*x)*cos(pi*y)\n"
            "// so that laplacian(phi_exact) = omega.\n"
        ),
    )


def write_poisson_fields(
    case: Path,
    centres: list[tuple[float, float, float]] | list[tuple[float, float]],
    solution_name: str = "phi",
    source_name: str = "omega",
) -> tuple[Path, Path]:
    """Write both fields needed by poissonFoamStudent."""
    return (
        write_poisson_solution_field(case, centres, solution_name),
        write_poisson_source_field(case, centres, source_name),
    )


__all__ = [
    "poisson_exact_value",
    "poisson_source_value",
    "write_poisson_fields",
    "write_poisson_solution_field",
    "write_poisson_source_field",
]
