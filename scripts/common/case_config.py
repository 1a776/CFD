#!/usr/bin/env python3

"""Configuration loader for one linear-advection study case."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from paths import CASES_DIR, CONFIG_DIR, PROJECT_ROOT


DEFAULT_RESOLUTIONS = (10, 20, 40, 80)


@dataclass(frozen=True)
class CaseConfig:
    """Small typed view of one JSON case configuration."""

    path: Path
    case_name: str
    description: str
    problem: str
    mesh_type: str
    scheme_name: str
    div_scheme: str
    solver: str
    template_resolution: int
    resolutions: tuple[int, ...]
    end_time: float
    max_co: float
    velocity: tuple[float, float, float]
    thickness: float
    implemented: bool
    grad_t_scheme: str | None

    @property
    def case_root(self) -> Path:
        return CASES_DIR / self.case_name

    @property
    def template_case(self) -> Path:
        return self.case_root / f"N{self.template_resolution}"

    def case_dir(self, resolution: int) -> Path:
        return self.case_root / f"N{resolution}"

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


def load_config(value: str | Path) -> CaseConfig:
    path = resolve_config_path(value)
    data = json.loads(path.read_text(encoding="utf-8"))
    resolutions = tuple(int(item) for item in data.get("resolutions", DEFAULT_RESOLUTIONS))
    if not resolutions:
        raise RuntimeError(f"Config has no resolutions: {path}")

    return CaseConfig(
        path=path,
        case_name=str(data["caseName"]),
        description=str(data.get("description", "")),
        problem=str(data.get("problem", "sine_wave_advection")),
        mesh_type=str(data.get("meshType", "quad")),
        scheme_name=str(data.get("schemeName", "")),
        div_scheme=str(data["divScheme"]),
        solver=str(data.get("solver", "explicitAdvectionFoamStudent")),
        template_resolution=int(data.get("templateResolution", 20)),
        resolutions=resolutions,
        end_time=float(data.get("endTime", 1.0)),
        max_co=float(data.get("maxCo", 0.2)),
        velocity=_tuple3(data.get("velocity", [1.0, 1.0, 0.0]), "velocity"),
        thickness=float(data.get("thickness", 0.1)),
        implemented=bool(data.get("implemented", True)),
        grad_t_scheme=data.get("gradTScheme"),
    )
