#!/usr/bin/env python3

"""Configuration loader for one linear-advection study case."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from paths import CONFIG_DIR, PROJECT_ROOT, solver_case_dir, solver_cases_dir


DEFAULT_RESOLUTIONS = (10, 20, 40, 80)


@dataclass(frozen=True)
class CaseConfig:
    """Small typed view of one JSON case configuration."""

    path: Path
    solver_family: str
    case_name: str
    description: str
    equation: str
    problem: str
    mesh_type: str
    mesh_backend: str
    gmsh_python: str | None
    scheme_name: str
    div_scheme: str
    laplacian_scheme: str
    sn_grad_scheme: str
    scalar_field: str
    source_field: str
    n_non_orthogonal_correctors: int
    linear_solver: str
    linear_tolerance: float
    diffusivity: float
    diffusion_co: float
    advection_diffusion_co: float
    max_delta_t: float | None
    solver: str
    resolutions: tuple[int, ...]
    end_time: float
    max_co: float
    velocity: tuple[float, float, float]
    domain: tuple[float, float, float, float]
    velocity_model: dict[str, Any]
    initial_profile: dict[str, Any]
    boundary_condition: str
    postprocess: dict[str, Any]
    thickness: float
    implemented: bool
    grad_t_scheme: str | None

    @property
    def case_root(self) -> Path:
        return solver_cases_dir(self.solver_family) / self.case_name

    def case_dir(self, resolution: int) -> Path:
        return solver_case_dir(self.solver_family, self.case_name, resolution)

    def require_implemented(self) -> None:
        if not self.implemented:
            raise RuntimeError(
                f"Config is marked as planned but not implemented yet: {self.path}"
            )


def parse_resolutions(value: str) -> list[int]:
    """Parse command-line resolution lists such as 10,20,40,80."""
    try:
        values = sorted({int(item.strip()) for item in value.split(",") if item.strip()})
    except ValueError as exc:
        raise argparse.ArgumentTypeError("resolutions must be integers") from exc
    if not values or any(item <= 0 for item in values):
        raise argparse.ArgumentTypeError("resolutions must be positive integers")
    return values


def resolve_config_path(value: str | Path) -> Path:
    """Accept either a JSON path or a config stem under scripts/configs."""
    path = Path(value).expanduser()
    candidates: list[Path]
    if path.suffix == ".json":
        candidates = [path, PROJECT_ROOT / path]
    else:
        candidates = [
            path,
            CONFIG_DIR / f"{path.name}.json",
            CONFIG_DIR / path.name,
        ]
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    raise FileNotFoundError(f"Cannot find config: {value}")


def _tuple3(value: Any, name: str) -> tuple[float, float, float]:
    if not isinstance(value, list | tuple) or len(value) != 3:
        raise RuntimeError(f"{name} must be a list with three numbers")
    return (float(value[0]), float(value[1]), float(value[2]))


def _tuple4(value: Any, name: str) -> tuple[float, float, float, float]:
    if not isinstance(value, list | tuple) or len(value) != 4:
        raise RuntimeError(f"{name} must be a list with four numbers")
    return (float(value[0]), float(value[1]), float(value[2]), float(value[3]))


def _dict_value(value: Any, name: str) -> dict[str, Any]:
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise RuntimeError(f"{name} must be an object")
    return dict(value)


def load_config(value: str | Path) -> CaseConfig:
    path = resolve_config_path(value)
    data = json.loads(path.read_text(encoding="utf-8"))
    resolutions = tuple(int(item) for item in data.get("resolutions", DEFAULT_RESOLUTIONS))
    if not resolutions:
        raise RuntimeError(f"Config has no resolutions: {path}")

    return CaseConfig(
        path=path,
        solver_family=str(data.get("solverFamily", "01_advection_equation")),
        case_name=str(data["caseName"]),
        description=str(data.get("description", "")),
        equation=str(data.get("equation", str(data.get("problem", "sine_wave_advection")).split("_")[-1])),
        problem=str(data.get("problem", "sine_wave_advection")),
        mesh_type=str(data.get("meshType", "quad")),
        mesh_backend=str(
            data.get("meshBackend", "blockMesh" if data.get("meshType", "quad") == "quad" else "")
        ),
        gmsh_python=(
            str(data["gmshPython"])
            if data.get("gmshPython") is not None
            else None
        ),
        scheme_name=str(data.get("schemeName", "")),
        div_scheme=str(data.get("divScheme", "")),
        laplacian_scheme=str(data.get("laplacianScheme", "Gauss linear corrected")),
        sn_grad_scheme=str(data.get("snGradScheme", "corrected")),
        scalar_field=str(data.get("scalarField", "T")),
        source_field=str(data.get("sourceField", "omega")),
        n_non_orthogonal_correctors=int(
            data.get("nNonOrthogonalCorrectors", 0)
        ),
        linear_solver=str(data.get("linearSolver", "GAMG")),
        linear_tolerance=float(data.get("linearTolerance", 1e-12)),
        diffusivity=float(data.get("mu", 1.0)),
        diffusion_co=float(data.get("diffusionCo", data.get("maxCo", 0.5))),
        advection_diffusion_co=float(
            data.get("advectionDiffusionCo", data.get("diffusionCo", data.get("maxCo", 0.5)))
        ),
        max_delta_t=(
            float(data["maxDeltaT"])
            if data.get("maxDeltaT") is not None
            else None
        ),
        solver=str(data.get("solver", "explicitAdvectionFoamStudent")),
        resolutions=resolutions,
        end_time=float(data.get("endTime", 1.0)),
        max_co=float(data.get("maxCo", 0.2)),
        velocity=_tuple3(data.get("velocity", [1.0, 1.0, 0.0]), "velocity"),
        domain=_tuple4(data.get("domain", [0.0, 1.0, 0.0, 1.0]), "domain"),
        velocity_model=_dict_value(data.get("velocityModel"), "velocityModel"),
        initial_profile=_dict_value(data.get("initialProfile"), "initialProfile"),
        boundary_condition=str(data.get("boundaryCondition", "periodicXY")),
        postprocess=_dict_value(data.get("postprocess"), "postprocess"),
        thickness=float(data.get("thickness", 0.1)),
        implemented=bool(data.get("implemented", True)),
        grad_t_scheme=data.get("gradTScheme"),
    )
