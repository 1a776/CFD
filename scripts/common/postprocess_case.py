#!/usr/bin/env python3

"""Post-process one structured sine-wave advection case."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
MPLCONFIGDIR = PROJECT_ROOT / "build" / "matplotlib"
MPLCONFIGDIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np

from advection_sine import exact_values_at_centres
from advection_tools import (
    control_value,
    exact_values,
    latest_time,
    mesh_resolution,
    normalized_errors,
    parse_solver_log,
    read_scalar_field,
)
from foam_fields import read_cell_geometry, read_tri_geometry_metadata

THICKNESS = 0.1
DEFAULT_VELOCITY = (1.0, 1.0, 0.0)


def default_case(default_case_name: str | None) -> Path:
    environment_case = os.environ.get("VIBEFLOW_CASE_DIR")
    if environment_case:
        return Path(environment_case).expanduser().resolve()
    if default_case_name:
        return PROJECT_ROOT / "cases" / default_case_name / "N20"
    raise SystemExit("Please provide --case-dir or set VIBEFLOW_CASE_DIR")


def _case_velocity(case: Path) -> tuple[float, float, float]:
    metadata_path = case / "metadata.json"
    if not metadata_path.exists():
        return DEFAULT_VELOCITY
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    velocity = metadata.get("velocity", list(DEFAULT_VELOCITY))
    if not isinstance(velocity, list | tuple) or len(velocity) != 3:
        return DEFAULT_VELOCITY
    return (float(velocity[0]), float(velocity[1]), float(velocity[2]))


def write_time_history(records: list[dict[str, float]], output: Path) -> None:
    fields = [
        "time",
        "step",
        "deltaT",
        "maxCo",
        "residualIntegral",
        "minT",
        "maxT",
        "amplitude",
    ]
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(
            {field: record.get(field, "") for field in fields}
            for record in records
        )


def write_field_data(
    initial: np.ndarray,
    numerical: np.ndarray,
    exact: np.ndarray,
    data_output: Path,
    error_output: Path,
) -> None:
    ny, nx = numerical.shape
    with data_output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            ["i", "j", "x", "y", "initial", "numerical", "exact", "error", "absError"]
        )
        for j in range(ny):
            y = (j + 0.5) / ny
            for i in range(nx):
                x = (i + 0.5) / nx
                error = float(numerical[j, i] - exact[j, i])
                writer.writerow(
                    [
                        i,
                        j,
                        f"{x:.16e}",
                        f"{y:.16e}",
                        f"{initial[j, i]:.16e}",
                        f"{numerical[j, i]:.16e}",
                        f"{exact[j, i]:.16e}",
                        f"{error:.16e}",
                        f"{abs(error):.16e}",
                    ]
                )

    with error_output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["i", "j", "x", "y", "numerical", "exact", "error", "absError"])
        for j in range(ny):
            y = (j + 0.5) / ny
            for i in range(nx):
                x = (i + 0.5) / nx
                error = float(numerical[j, i] - exact[j, i])
                writer.writerow(
                    [
                        i,
                        j,
                        f"{x:.16e}",
                        f"{y:.16e}",
                        f"{numerical[j, i]:.16e}",
                        f"{exact[j, i]:.16e}",
                        f"{error:.16e}",
                        f"{abs(error):.16e}",
                    ]
                )


def write_tri_field_data(
    centres: list[tuple[float, float, float]],
    volumes: list[float],
    initial: list[float],
    numerical: list[float],
    exact: list[float],
    data_output: Path,
    error_output: Path,
) -> None:
    """Write cell-centred data for an unstructured triangular mesh."""
    if not (
        len(centres)
        == len(volumes)
        == len(initial)
        == len(numerical)
        == len(exact)
    ):
        raise RuntimeError("Triangular field and geometry sizes do not match")

    with data_output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "cell",
                "x",
                "y",
                "z",
                "volume",
                "initial",
                "numerical",
                "exact",
                "error",
                "absError",
            ]
        )
        for cell, ((x, y, z), volume, initial_value, numerical_value, exact_value) in enumerate(
            zip(centres, volumes, initial, numerical, exact)
        ):
            error = numerical_value - exact_value
            writer.writerow(
                [
                    cell,
                    f"{x:.16e}",
                    f"{y:.16e}",
                    f"{z:.16e}",
                    f"{volume:.16e}",
                    f"{initial_value:.16e}",
                    f"{numerical_value:.16e}",
                    f"{exact_value:.16e}",
                    f"{error:.16e}",
                    f"{abs(error):.16e}",
                ]
            )

    with error_output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "cell",
                "x",
                "y",
                "z",
                "volume",
                "numerical",
                "exact",
                "error",
                "absError",
            ]
        )
        for cell, ((x, y, z), volume, numerical_value, exact_value) in enumerate(
            zip(centres, volumes, numerical, exact)
        ):
            error = numerical_value - exact_value
            writer.writerow(
                [
                    cell,
                    f"{x:.16e}",
                    f"{y:.16e}",
                    f"{z:.16e}",
                    f"{volume:.16e}",
                    f"{numerical_value:.16e}",
                    f"{exact_value:.16e}",
                    f"{error:.16e}",
                    f"{abs(error):.16e}",
                ]
            )


def plot_field_comparison(
    initial: np.ndarray,
    numerical: np.ndarray,
    exact: np.ndarray,
    output: Path,
    final_time: float,
) -> None:
    error = numerical - exact
    maximum = max(
        float(np.max(np.abs(initial))),
        float(np.max(np.abs(numerical))),
        float(np.max(np.abs(exact))),
    )
    error_maximum = max(float(np.max(np.abs(error))), 1.0e-16)
    ny, nx = numerical.shape
    x_edges = np.linspace(0.0, 1.0, nx + 1)
    y_edges = np.linspace(0.0, 1.0, ny + 1)

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
    panels = [
        (axes[0, 0], initial, "Initial field T(x,y,0)", -maximum, maximum),
        (axes[0, 1], numerical, f"Numerical field, t={final_time:g}", -maximum, maximum),
        (axes[1, 0], exact, f"Exact field, t={final_time:g}", -maximum, maximum),
        (axes[1, 1], error, "Numerical minus exact", -error_maximum, error_maximum),
    ]

    for axis, values, title, lower, upper in panels:
        image = axis.pcolormesh(
            x_edges,
            y_edges,
            values,
            shading="auto",
            cmap="coolwarm",
            vmin=lower,
            vmax=upper,
        )
        axis.set_aspect("equal")
        axis.set_xlabel("x")
        axis.set_ylabel("y")
        axis.set_title(title)
        fig.colorbar(image, ax=axis, shrink=0.86)

    fig.suptitle("Periodic linear-advection field comparison")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def _triangulation_from_metadata(
    metadata: dict[str, object],
) -> mtri.Triangulation:
    """Build a Matplotlib triangulation from Gmsh node/connectivity metadata."""
    nodes = np.asarray(metadata["nodes"], dtype=float)
    triangles = np.asarray(metadata["triangles"], dtype=int)
    if nodes.ndim != 2 or nodes.shape[1] != 2:
        raise RuntimeError("Triangular mesh metadata has invalid node coordinates")
    if triangles.ndim != 2 or triangles.shape[1] != 3:
        raise RuntimeError("Triangular mesh metadata has invalid connectivity")
    return mtri.Triangulation(nodes[:, 0], nodes[:, 1], triangles)


def _match_triangles_to_cells(
    cell_centres: list[tuple[float, float, float]],
    triangulation: mtri.Triangulation,
    tolerance: float = 1.0e-8,
) -> list[int]:
    """Return the OpenFOAM cell index associated with each plotted triangle."""
    triangle_nodes = triangulation.triangles
    triangle_centres = np.column_stack(
        (
            np.mean(triangulation.x[triangle_nodes], axis=1),
            np.mean(triangulation.y[triangle_nodes], axis=1),
        )
    )
    cells = np.asarray([(x, y) for x, y, _ in cell_centres], dtype=float)
    if len(cells) != len(triangle_centres):
        raise RuntimeError(
            "The number of 2-D Gmsh triangles does not match OpenFOAM cells: "
            f"triangles={len(triangle_centres)}, cells={len(cells)}"
        )

    # The transfinite Gmsh mesh has the same geometric centroids as the
    # imported OpenFOAM prism cells. Quantization avoids relying on numbering.
    buckets: dict[tuple[int, int], list[int]] = {}
    for cell_index, (x, y) in enumerate(cells):
        key = (round(float(x) / tolerance), round(float(y) / tolerance))
        buckets.setdefault(key, []).append(cell_index)

    triangle_to_cell: list[int] = []
    used: set[int] = set()
    for triangle_index, (x, y) in enumerate(triangle_centres):
        key = (round(float(x) / tolerance), round(float(y) / tolerance))
        candidates = buckets.get(key, [])
        if not candidates:
            distance = np.sum((cells - np.asarray([x, y])) ** 2, axis=1)
            candidates = [int(np.argmin(distance))]
        unused = [candidate for candidate in candidates if candidate not in used]
        if not unused:
            raise RuntimeError(
                "Could not construct a one-to-one triangle-to-cell mapping "
                f"for triangle {triangle_index}"
            )
        selected = min(
            unused,
            key=lambda candidate: float(np.sum((cells[candidate] - [x, y]) ** 2)),
        )
        used.add(selected)
        triangle_to_cell.append(selected)

    if len(used) != len(cells):
        raise RuntimeError("Some OpenFOAM cells were not matched to Gmsh triangles")
    return triangle_to_cell


def plot_tri_field_comparison(
    triangulation: mtri.Triangulation,
    initial: np.ndarray,
    numerical: np.ndarray,
    exact: np.ndarray,
    triangle_to_cell: list[int],
    output: Path,
    final_time: float,
) -> None:
    """Plot cell-centred fields on the actual triangular mesh."""
    error = numerical - exact
    initial_faces = initial[triangle_to_cell]
    numerical_faces = numerical[triangle_to_cell]
    exact_faces = exact[triangle_to_cell]
    error_faces = error[triangle_to_cell]
    maximum = max(
        float(np.max(np.abs(initial_faces))),
        float(np.max(np.abs(numerical_faces))),
        float(np.max(np.abs(exact_faces))),
    )
    error_maximum = max(float(np.max(np.abs(error_faces))), 1.0e-16)
    panels = [
        (initial_faces, "Initial field T(x,y,0)", -maximum, maximum),
        (numerical_faces, f"Numerical field, t={final_time:g}", -maximum, maximum),
        (exact_faces, f"Exact field, t={final_time:g}", -maximum, maximum),
        (error_faces, "Numerical minus exact", -error_maximum, error_maximum),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
    for axis, (values, title, lower, upper) in zip(axes.flat, panels):
        image = axis.tripcolor(
            triangulation,
            facecolors=values,
            shading="flat",
            cmap="coolwarm",
            vmin=lower,
            vmax=upper,
        )
        axis.triplot(triangulation, color="black", linewidth=0.25, alpha=0.25)
        axis.set_aspect("equal")
        axis.set_xlabel("x")
        axis.set_ylabel("y")
        axis.set_title(title)
        fig.colorbar(image, ax=axis, shrink=0.86)

    fig.suptitle("Periodic linear-advection field comparison on triangular cells")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_diagonal_profile(
    numerical: np.ndarray, exact: np.ndarray, output: Path, final_time: float
) -> None:
    count = min(numerical.shape)
    coordinate = (np.arange(count) + 0.5) / max(numerical.shape)
    index = np.arange(count)

    fig, axis = plt.subplots(figsize=(8, 4.8), constrained_layout=True)
    axis.plot(coordinate, numerical[index, index], "o-", markersize=3.5, label="numerical")
    axis.plot(coordinate, exact[index, index], "--", linewidth=1.5, label="exact")
    axis.set_xlabel("diagonal coordinate, x = y")
    axis.set_ylabel("T")
    axis.set_title(f"Diagonal profile at t={final_time:g}")
    axis.grid(True, alpha=0.25)
    axis.legend()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_tri_diagonal_profile(
    centres: list[tuple[float, float, float]],
    numerical: np.ndarray,
    exact: np.ndarray,
    output: Path,
    final_time: float,
    resolution: int,
) -> None:
    """Plot cells nearest to x=y as a practical diagonal slice."""
    coordinates = np.asarray([(x, y) for x, y, _ in centres], dtype=float)
    distance = np.abs(coordinates[:, 0] - coordinates[:, 1])
    band = max(0.35 / max(resolution, 1), 1.0e-10)
    selected = np.flatnonzero(distance <= band)
    if len(selected) < 2:
        selected = np.argsort(distance)[: min(2 * max(resolution, 1), len(distance))]
    selected = selected[np.argsort(np.sum(coordinates[selected], axis=1))]
    diagonal_coordinate = 0.5 * np.sum(coordinates[selected], axis=1)

    fig, axis = plt.subplots(figsize=(8, 4.8), constrained_layout=True)
    axis.plot(
        diagonal_coordinate,
        numerical[selected],
        "o",
        markersize=3.5,
        label="numerical cell centres",
    )
    axis.plot(
        diagonal_coordinate,
        exact[selected],
        "--",
        linewidth=1.5,
        label="exact at same centres",
    )
    axis.set_xlabel("diagonal coordinate, cells nearest x = y")
    axis.set_ylabel("T")
    axis.set_title(f"Approximate diagonal profile at t={final_time:g}")
    axis.grid(True, alpha=0.25)
    axis.legend()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_amplitude_history(records: list[dict[str, float]], output: Path) -> None:
    fig, axis = plt.subplots(figsize=(8, 4.8), constrained_layout=True)
    axis.plot(
        [record["time"] for record in records],
        [record["amplitude"] for record in records],
        "o-",
        markersize=2.5,
        label="numerical amplitude",
    )
    axis.axhline(1.0, linestyle="--", color="black", label="exact amplitude")
    axis.set_xlabel("time")
    axis.set_ylabel("amplitude, (max(T) - min(T))/2")
    axis.set_title("Numerical dissipation seen through wave amplitude")
    axis.grid(True, alpha=0.25)
    axis.legend()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_cfl_history(
    records: list[dict[str, float]], output: Path, target_co: float
) -> None:
    fig, axis = plt.subplots(figsize=(8, 4.8), constrained_layout=True)
    axis.plot(
        [record["time"] for record in records],
        [record["maxCo"] for record in records],
        "o-",
        markersize=2.5,
        label="actual maxCo",
    )
    axis.axhline(target_co, linestyle="--", color="black", label="target maxCo")
    axis.set_xlabel("time")
    axis.set_ylabel("maximum Courant number")
    axis.set_title("CFL stability monitor")
    axis.grid(True, alpha=0.25)
    axis.legend()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def _read_mesh_type(case: Path) -> str:
    metadata_path = case / "metadata.json"
    if not metadata_path.exists():
        return "quad"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    return str(metadata.get("meshType", "quad"))


def _postprocess_quad_case(
    case: Path, output_root: Path, target_co: float
) -> dict[str, object]:
    case = case.resolve()
    velocity = _case_velocity(case)
    nx, ny = mesh_resolution(case)
    final_time, final_dir = latest_time(case)
    initial_values = read_scalar_field(case / "0" / "T")
    numerical_values = read_scalar_field(final_dir / "T")
    if len(initial_values) != nx * ny or len(numerical_values) != nx * ny:
        raise RuntimeError(f"Field size does not match mesh in {case}")

    exact_values_at_final = exact_values(nx, ny, final_time, velocity)
    initial = np.asarray(initial_values, dtype=float).reshape(ny, nx)
    numerical = np.asarray(numerical_values, dtype=float).reshape(ny, nx)
    exact = np.asarray(exact_values_at_final, dtype=float).reshape(ny, nx)
    records = parse_solver_log(case / "log.explicitAdvectionFoamStudent")

    cell_volume = THICKNESS / (nx * ny)
    l1, l2, linf = normalized_errors(
        numerical_values, exact_values_at_final, cell_volume
    )
    initial_mass = float(np.sum(initial) * cell_volume)
    final_mass = float(np.sum(numerical) * cell_volume)
    mass_scale = float(np.sum(np.abs(initial)) * cell_volume)
    mass_change = final_mass - initial_mass
    solver_log = (case / "log.explicitAdvectionFoamStudent").read_text(encoding="utf-8")
    mesh_log_path = case / "log.checkMesh"
    mesh_log = mesh_log_path.read_text(encoding="utf-8") if mesh_log_path.exists() else ""
    configured_end_time = control_value(case, "endTime", 1.0)

    data_dir = output_root / "data" / "cases" / case.parent.name / case.name
    figure_dir = output_root / "figures" / "cases" / case.parent.name / case.name
    data_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    write_time_history(records, data_dir / "time_history.csv")
    write_field_data(
        initial,
        numerical,
        exact,
        data_dir / "field_data.csv",
        data_dir / "error_field.csv",
    )

    summary: dict[str, object] = {
        "case": str(case),
        "caseName": case.parent.name,
        "mesh": "quad",
        "problem": "sine",
        "resolution": nx,
        "nominalH": 1.0 / nx,
        "nCells": nx * ny,
        "finalTime": final_time,
        "finalDirectory": str(final_dir),
        "finalTimeError": abs(final_time - configured_end_time),
        "normalizedL1": l1,
        "normalizedL2": l2,
        "normalizedLinf": linf,
        "initialMass": initial_mass,
        "finalMass": final_mass,
        "massChange": mass_change,
        "normalizedMassError": abs(mass_change) / mass_scale if mass_scale else 0.0,
        "maxAbsResidualIntegral": max(
            abs(record.get("residualIntegral", 0.0)) for record in records
        ),
        "minCo": min(record["maxCo"] for record in records),
        "maxCo": max(record["maxCo"] for record in records),
        "targetCo": target_co,
        "initialAmplitude": 0.5 * (float(initial.max()) - float(initial.min())),
        "finalAmplitude": 0.5 * (float(numerical.max()) - float(numerical.min())),
        "minFinal": float(numerical.min()),
        "maxFinal": float(numerical.max()),
        "maxAbsFinal": float(np.max(np.abs(numerical))),
        "timeSteps": len(records),
        "meshOK": "Mesh OK." in mesh_log,
        "boundedBelowInitial": float(numerical.min()) >= float(initial.min()) - 1.0e-12,
        "boundedAboveInitial": float(numerical.max()) <= float(initial.max()) + 1.0e-12,
        "solverEnded": "End" in solver_log and "Stage 5 time loop completed" in solver_log,
        "solverFatal": "FOAM FATAL ERROR" in solver_log,
    }
    (data_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    plot_field_comparison(
        initial, numerical, exact, figure_dir / "field_comparison.png", final_time
    )
    plot_diagonal_profile(
        numerical, exact, figure_dir / "diagonal_profile.png", final_time
    )
    plot_amplitude_history(records, figure_dir / "amplitude_history.png")
    plot_cfl_history(records, figure_dir / "cfl_history.png", target_co)

    print(f"case={case}")
    print(f"resolution={nx}")
    print(f"finalTime={final_time:.16g}")
    print(f"normalizedL1={l1:.16e}")
    print(f"normalizedL2={l2:.16e}")
    print(f"normalizedLinf={linf:.16e}")
    print(f"data={data_dir}")
    print(f"figures={figure_dir}")
    return summary


def _postprocess_tri_case(
    case: Path, output_root: Path, target_co: float
) -> dict[str, object]:
    """Post-process a triangular-prism case using real cell geometry."""
    case = case.resolve()
    velocity = _case_velocity(case)
    centres, volumes = read_cell_geometry(case)
    geometry_metadata = read_tri_geometry_metadata(case / "mesh" / "mesh_geometry.json")
    triangulation = _triangulation_from_metadata(geometry_metadata)
    triangle_to_cell = _match_triangles_to_cells(centres, triangulation)

    final_time, final_dir = latest_time(case)
    initial_values = read_scalar_field(case / "0" / "T")
    numerical_values = read_scalar_field(final_dir / "T")
    if len(initial_values) != len(centres) or len(numerical_values) != len(centres):
        raise RuntimeError(
            f"Field size does not match triangular mesh in {case}: "
            f"cells={len(centres)}, initial={len(initial_values)}, "
            f"final={len(numerical_values)}"
        )

    exact_values_at_final = exact_values_at_centres(centres, final_time, velocity)
    initial = np.asarray(initial_values, dtype=float)
    numerical = np.asarray(numerical_values, dtype=float)
    exact = np.asarray(exact_values_at_final, dtype=float)
    volumes_array = np.asarray(volumes, dtype=float)
    records = parse_solver_log(case / "log.explicitAdvectionFoamStudent")

    l1, l2, linf = normalized_errors(
        numerical_values, exact_values_at_final, volumes
    )
    initial_mass = float(np.dot(initial, volumes_array))
    final_mass = float(np.dot(numerical, volumes_array))
    mass_scale = float(np.dot(np.abs(initial), volumes_array))
    mass_change = final_mass - initial_mass
    solver_log = (case / "log.explicitAdvectionFoamStudent").read_text(encoding="utf-8")
    mesh_log_path = case / "log.checkMesh"
    mesh_log = mesh_log_path.read_text(encoding="utf-8") if mesh_log_path.exists() else ""
    configured_end_time = control_value(case, "endTime", 1.0)
    case_metadata = json.loads((case / "metadata.json").read_text(encoding="utf-8"))
    resolution = int(case_metadata.get("resolution", geometry_metadata["resolution"]))

    data_dir = output_root / "data" / "cases" / case.parent.name / case.name
    figure_dir = output_root / "figures" / "cases" / case.parent.name / case.name
    data_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    write_time_history(records, data_dir / "time_history.csv")
    write_tri_field_data(
        centres,
        volumes,
        initial_values,
        numerical_values,
        exact_values_at_final,
        data_dir / "field_data.csv",
        data_dir / "error_field.csv",
    )

    summary: dict[str, object] = {
        "case": str(case),
        "caseName": case.parent.name,
        "mesh": "tri",
        "problem": "sine",
        "resolution": resolution,
        "nominalH": 1.0 / resolution,
        "nCells": len(centres),
        "finalTime": final_time,
        "finalDirectory": str(final_dir),
        "finalTimeError": abs(final_time - configured_end_time),
        "normalizedL1": l1,
        "normalizedL2": l2,
        "normalizedLinf": linf,
        "initialMass": initial_mass,
        "finalMass": final_mass,
        "massChange": mass_change,
        "normalizedMassError": abs(mass_change) / mass_scale if mass_scale else 0.0,
        "maxAbsResidualIntegral": max(
            abs(record.get("residualIntegral", 0.0)) for record in records
        ),
        "minCo": min(record["maxCo"] for record in records),
        "maxCo": max(record["maxCo"] for record in records),
        "targetCo": target_co,
        "initialAmplitude": 0.5 * (float(initial.max()) - float(initial.min())),
        "finalAmplitude": 0.5 * (float(numerical.max()) - float(numerical.min())),
        "minFinal": float(numerical.min()),
        "maxFinal": float(numerical.max()),
        "maxAbsFinal": float(np.max(np.abs(numerical))),
        "timeSteps": len(records),
        "meshOK": "Mesh OK." in mesh_log,
        "boundedBelowInitial": float(numerical.min()) >= float(initial.min()) - 1.0e-12,
        "boundedAboveInitial": float(numerical.max()) <= float(initial.max()) + 1.0e-12,
        "solverEnded": "End" in solver_log and "Stage 5 time loop completed" in solver_log,
        "solverFatal": "FOAM FATAL ERROR" in solver_log,
    }
    (data_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    plot_tri_field_comparison(
        triangulation,
        initial,
        numerical,
        exact,
        triangle_to_cell,
        figure_dir / "field_comparison.png",
        final_time,
    )
    plot_tri_diagonal_profile(
        centres,
        numerical,
        exact,
        figure_dir / "diagonal_profile.png",
        final_time,
        resolution,
    )
    plot_amplitude_history(records, figure_dir / "amplitude_history.png")
    plot_cfl_history(records, figure_dir / "cfl_history.png", target_co)

    print(f"case={case}")
    print(f"resolution={resolution}")
    print(f"mesh=tri")
    print(f"nCells={len(centres)}")
    print(f"finalTime={final_time:.16g}")
    print(f"normalizedL1={l1:.16e}")
    print(f"normalizedL2={l2:.16e}")
    print(f"normalizedLinf={linf:.16e}")
    print(f"data={data_dir}")
    print(f"figures={figure_dir}")
    return summary


def postprocess_case(case: Path, output_root: Path, target_co: float) -> dict[str, object]:
    mesh_type = _read_mesh_type(case)
    if mesh_type == "quad":
        return _postprocess_quad_case(case, output_root, target_co)
    if mesh_type == "tri":
        return _postprocess_tri_case(case, output_root, target_co)
    raise NotImplementedError(f"Unsupported mesh type in metadata.json: {mesh_type}")


def main(default_case_name: str | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-dir", type=Path, default=None)
    parser.add_argument("--output-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--target-co", type=float, default=0.2)
    args = parser.parse_args()
    postprocess_case(
        args.case_dir.resolve() if args.case_dir else default_case(default_case_name),
        args.output_root.resolve(),
        args.target_co,
    )


if __name__ == "__main__":
    main()
