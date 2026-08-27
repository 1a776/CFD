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
from advection_diffusion_tools import (
    exact_value as advection_diffusion_exact_value,
    parse_advection_diffusion_solver_log,
    rotating_peak_exact_value,
    rotating_peak_structured_values,
    rotating_peak_values_at_centres,
    structured_values as advection_diffusion_structured_values,
)
from advection_tools import (
    control_value,
    exact_values,
    latest_time,
    mesh_resolution,
    normalized_errors,
    parse_solver_log,
    read_scalar_field,
)
from diffusion_tools import (
    discontinuous_exact_value,
    discontinuous_values,
    gaussian_exact_value,
    gaussian_values,
    parse_diffusion_solver_log,
    structured_centres,
)
from poisson_tools import poisson_exact_value
from foam_fields import read_cell_geometry, read_tri_geometry_metadata
from paths import solver_data_dir, solver_figure_dir

THICKNESS = 0.1
DEFAULT_VELOCITY = (1.0, 1.0, 0.0)


def default_case(default_case_name: str | None) -> Path:
    environment_case = os.environ.get("VIBEFLOW_CASE_DIR")
    if environment_case:
        return Path(environment_case).expanduser().resolve()
    if default_case_name:
        return PROJECT_ROOT / "cases" / "01_advection_equation" / default_case_name / "N20"
    raise SystemExit("Please provide --case-dir or set VIBEFLOW_CASE_DIR")


def _case_namespace(case: Path) -> tuple[str, str]:
    """Return solver family and case name from the namespaced case path."""
    metadata_path = case / "metadata.json"
    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        family = metadata.get("solverFamily")
        name = metadata.get("caseName")
        if family and name:
            return str(family), str(name)
    resolved = case.resolve()
    return resolved.parent.parent.name, resolved.parent.name


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


def write_structured_field_data(
    centres: list[tuple[float, float]],
    initial: list[float],
    numerical: list[float],
    exact: list[float],
    data_output: Path,
    error_output: Path,
) -> None:
    """Write structured cell-centred data with physical coordinates."""
    if not (len(centres) == len(initial) == len(numerical) == len(exact)):
        raise RuntimeError("Structured field data sizes do not match")
    with data_output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["cell", "x", "y", "initial", "numerical", "exact", "error", "absError"])
        for cell, ((x, y), initial_value, numerical_value, exact_value) in enumerate(
            zip(centres, initial, numerical, exact)
        ):
            error = numerical_value - exact_value
            writer.writerow(
                [
                    cell,
                    f"{x:.16e}",
                    f"{y:.16e}",
                    f"{initial_value:.16e}",
                    f"{numerical_value:.16e}",
                    f"{exact_value:.16e}",
                    f"{error:.16e}",
                    f"{abs(error):.16e}",
                ]
            )
    with error_output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["cell", "x", "y", "numerical", "exact", "error", "absError"])
        for cell, ((x, y), numerical_value, exact_value) in enumerate(
            zip(centres, numerical, exact)
        ):
            error = numerical_value - exact_value
            writer.writerow(
                [
                    cell,
                    f"{x:.16e}",
                    f"{y:.16e}",
                    f"{numerical_value:.16e}",
                    f"{exact_value:.16e}",
                    f"{error:.16e}",
                    f"{abs(error):.16e}",
                ]
            )


def write_rotation_field_data(
    initial: np.ndarray,
    numerical: np.ndarray,
    data_output: Path,
) -> None:
    """Write structured solid-rotation final field data."""
    ny, nx = numerical.shape
    with data_output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["i", "j", "x", "y", "initial", "numerical", "change", "absChange"])
        for j in range(ny):
            y = (j + 0.5) / ny
            for i in range(nx):
                x = (i + 0.5) / nx
                change = float(numerical[j, i] - initial[j, i])
                writer.writerow(
                    [
                        i,
                        j,
                        f"{x:.16e}",
                        f"{y:.16e}",
                        f"{initial[j, i]:.16e}",
                        f"{numerical[j, i]:.16e}",
                        f"{change:.16e}",
                        f"{abs(change):.16e}",
                    ]
                )


def write_tri_rotation_field_data(
    centres: list[tuple[float, float, float]],
    volumes: list[float],
    initial: list[float],
    numerical: list[float],
    data_output: Path,
) -> None:
    """Write initial/final data for solid rotation on triangular cells."""
    if not (
        len(centres)
        == len(volumes)
        == len(initial)
        == len(numerical)
    ):
        raise RuntimeError("Triangular rotation field and geometry sizes do not match")

    with data_output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            ["cell", "x", "y", "z", "volume", "initial", "numerical", "change", "absChange"]
        )
        for cell, ((x, y, z), volume, initial_value, numerical_value) in enumerate(
            zip(centres, volumes, initial, numerical)
        ):
            change = numerical_value - initial_value
            writer.writerow(
                [
                    cell,
                    f"{x:.16e}",
                    f"{y:.16e}",
                    f"{z:.16e}",
                    f"{volume:.16e}",
                    f"{initial_value:.16e}",
                    f"{numerical_value:.16e}",
                    f"{change:.16e}",
                    f"{abs(change):.16e}",
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


def plot_advection_diffusion_field_comparison(
    initial: np.ndarray,
    numerical: np.ndarray,
    exact: np.ndarray,
    output: Path,
    final_time: float,
    domain: tuple[float, float, float, float],
) -> None:
    """Plot the third-problem scalar field with phi-specific labels."""
    xmin, xmax, ymin, ymax = domain
    error = numerical - exact
    initial_maximum = max(float(np.max(np.abs(initial))), 1.0e-16)
    numerical_maximum = max(float(np.max(np.abs(numerical))), 1.0e-40)
    exact_maximum = max(float(np.max(np.abs(exact))), 1.0e-40)
    error_maximum = max(float(np.max(np.abs(error))), 1.0e-16)
    ny, nx = numerical.shape
    x_edges = np.linspace(xmin, xmax, nx + 1)
    y_edges = np.linspace(ymin, ymax, ny + 1)
    panels = [
        (initial, "Initial field phi(x,y,0)", -initial_maximum, initial_maximum),
        (numerical, f"Numerical field, t={final_time:g}", -numerical_maximum, numerical_maximum),
        (exact, f"Exact field, t={final_time:g}", -exact_maximum, exact_maximum),
        (error, "Numerical minus exact", -error_maximum, error_maximum),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
    for axis, (values, title, lower, upper) in zip(axes.flat, panels):
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
    fig.suptitle("Periodic sine-wave advection-diffusion field comparison")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_advection_diffusion_diagonal_profile(
    numerical: np.ndarray,
    exact: np.ndarray,
    output: Path,
    final_time: float,
) -> None:
    """Plot phi along x=y for the structured advection-diffusion case."""
    count = min(numerical.shape)
    coordinate = (np.arange(count) + 0.5) / max(numerical.shape)
    index = np.arange(count)
    fig, axis = plt.subplots(figsize=(8, 4.8), constrained_layout=True)
    axis.plot(coordinate, numerical[index, index], "o-", markersize=3.5, label="numerical")
    axis.plot(coordinate, exact[index, index], "--", linewidth=1.5, label="exact")
    axis.set_xlabel("diagonal coordinate, x = y")
    axis.set_ylabel("phi")
    axis.set_title(f"Diagonal phi profile at t={final_time:g}")
    axis.grid(True, alpha=0.25)
    axis.legend()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_advection_diffusion_amplitude_history(
    records: list[dict[str, float]],
    output: Path,
    mu: float,
) -> None:
    """Compare numerical decay with exp(-8*pi^2*mu*t)."""
    times = [record["time"] for record in records]
    numerical_amplitude = [record["amplitude"] for record in records]
    exact_amplitude = [math.exp(-8.0 * math.pi * math.pi * mu * time) for time in times]
    fig, axis = plt.subplots(figsize=(8, 4.8), constrained_layout=True)
    axis.plot(times, numerical_amplitude, "o-", markersize=2.5, label="numerical amplitude")
    axis.plot(times, exact_amplitude, "--", linewidth=1.5, label="exact exp(-8*pi^2*mu*t)")
    axis.set_xlabel("time")
    axis.set_ylabel("amplitude, (max(phi)-min(phi))/2")
    axis.set_title("Diffusive amplitude decay")
    axis.grid(True, alpha=0.25)
    axis.legend()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_advection_diffusion_stability_history(
    records: list[dict[str, float]],
    output: Path,
    target: float,
) -> None:
    """Plot the combined explicit stability number logged by the solver."""
    fig, axis = plt.subplots(figsize=(8, 4.8), constrained_layout=True)
    axis.plot(
        [record["time"] for record in records],
        [record["maxCo"] for record in records],
        "o-",
        markersize=2.5,
        label="actual advectionDiffusionCo",
    )
    axis.axhline(target, linestyle="--", color="black", label="target safety coefficient")
    axis.set_xlabel("time")
    axis.set_ylabel("combined explicit stability number")
    axis.set_title("Advection-diffusion explicit stability monitor")
    axis.grid(True, alpha=0.25)
    axis.legend()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def _rotating_peak_parameters(
    metadata: dict[str, object],
) -> tuple[float, tuple[float, float]]:
    profile = metadata.get("initialProfile", {})
    if not isinstance(profile, dict):
        profile = {}
    centre_raw = profile.get("center", [0.0, 0.5])
    if not isinstance(centre_raw, list | tuple) or len(centre_raw) != 2:
        centre_raw = [0.0, 0.5]
    return (
        float(profile.get("diffusionStartTime", math.pi / 2.0)),
        (float(centre_raw[0]), float(centre_raw[1])),
    )


def plot_rotating_peak_midline_profile(
    centres: list[tuple[float, float]],
    numerical_values: list[float],
    exact_values_at_final: list[float],
    output: Path,
    final_time: float,
    centre_y: float,
) -> None:
    """Plot the horizontal line closest to the final peak centre."""
    rows: dict[float, list[tuple[float, float, float]]] = {}
    for (x, y), numerical, exact in zip(
        centres, numerical_values, exact_values_at_final
    ):
        rows.setdefault(round(float(y), 12), []).append((float(x), numerical, exact))
    y_key = min(rows, key=lambda value: abs(value - centre_y))
    points = sorted(rows[y_key], key=lambda item: item[0])

    fig, axis = plt.subplots(figsize=(8, 4.8), constrained_layout=True)
    axis.plot(
        [item[0] for item in points],
        [item[1] for item in points],
        "o-",
        markersize=3.0,
        label="numerical",
    )
    axis.plot(
        [item[0] for item in points],
        [item[2] for item in points],
        "--",
        linewidth=1.5,
        label="exact",
    )
    axis.set_xlabel(f"x at cell centres nearest y={centre_y:g}")
    axis.set_ylabel("phi")
    axis.set_title(f"Rotating peak horizontal profile at tau={final_time:g}")
    axis.grid(True, alpha=0.25)
    axis.legend()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_rotating_peak_structured_final_contour(
    numerical: np.ndarray,
    exact: np.ndarray,
    output: Path,
    final_time: float,
    domain: tuple[float, float, float, float],
) -> None:
    """Plot final numerical, exact and error contours for a structured mesh."""
    xmin, xmax, ymin, ymax = domain
    error = numerical - exact
    maximum = max(float(np.max(numerical)), float(np.max(exact)), 1.0e-40)
    error_maximum = max(float(np.max(np.abs(error))), 1.0e-40)
    ny, nx = numerical.shape
    x = np.linspace(xmin + 0.5 * (xmax - xmin) / nx, xmax - 0.5 * (xmax - xmin) / nx, nx)
    y = np.linspace(ymin + 0.5 * (ymax - ymin) / ny, ymax - 0.5 * (ymax - ymin) / ny, ny)
    x_grid, y_grid = np.meshgrid(x, y)
    positive_levels = np.linspace(0.0, maximum, 18)
    error_levels = np.linspace(-error_maximum, error_maximum, 18)

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.2), constrained_layout=True)
    panels = [
        (axes[0], numerical, "Numerical final contours", positive_levels, "viridis"),
        (axes[1], exact, "Exact final contours", positive_levels, "viridis"),
        (axes[2], error, "Numerical minus exact", error_levels, "coolwarm"),
    ]
    for axis, values, title, levels, cmap in panels:
        filled = axis.contourf(x_grid, y_grid, values, levels=levels, cmap=cmap)
        axis.contour(x_grid, y_grid, values, levels=levels, colors="black", linewidths=0.35, alpha=0.55)
        axis.set_aspect("equal")
        axis.set_xlabel("x")
        axis.set_ylabel("y")
        axis.set_title(title)
        fig.colorbar(filled, ax=axis, shrink=0.86)
    fig.suptitle(f"Rotating sharp-peak contours, t={final_time:g}")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_rotating_peak_tri_final_contour(
    triangulation: mtri.Triangulation,
    numerical: np.ndarray,
    exact: np.ndarray,
    triangle_to_cell: list[int],
    output: Path,
    final_time: float,
) -> None:
    """Plot final numerical, exact and error contours on triangular cells."""
    numerical_faces = numerical[triangle_to_cell]
    exact_faces = exact[triangle_to_cell]
    error_faces = numerical_faces - exact_faces
    maximum = max(float(np.max(numerical_faces)), float(np.max(exact_faces)), 1.0e-40)
    error_maximum = max(float(np.max(np.abs(error_faces))), 1.0e-40)
    positive_levels = np.linspace(0.0, maximum, 18)
    error_levels = np.linspace(-error_maximum, error_maximum, 18)

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.2), constrained_layout=True)
    panels = [
        (axes[0], numerical_faces, "Numerical final contours", positive_levels, "viridis"),
        (axes[1], exact_faces, "Exact final contours", positive_levels, "viridis"),
        (axes[2], error_faces, "Numerical minus exact", error_levels, "coolwarm"),
    ]
    for axis, values, title, levels, cmap in panels:
        filled = axis.tripcolor(
            triangulation,
            facecolors=values,
            shading="flat",
            cmap=cmap,
            vmin=float(levels[0]),
            vmax=float(levels[-1]),
        )
        node_values = np.zeros(len(triangulation.x), dtype=float)
        node_counts = np.zeros(len(triangulation.x), dtype=float)
        for triangle_index, nodes in enumerate(triangulation.triangles):
            node_values[nodes] += values[triangle_index]
            node_counts[nodes] += 1.0
        node_values /= np.maximum(node_counts, 1.0)
        axis.tricontour(
            triangulation,
            node_values,
            levels=levels,
            colors="black",
            linewidths=0.35,
            alpha=0.55,
        )
        axis.triplot(triangulation, color="black", linewidth=0.16, alpha=0.18)
        axis.set_aspect("equal")
        axis.set_xlabel("x")
        axis.set_ylabel("y")
        axis.set_title(title)
        fig.colorbar(filled, ax=axis, shrink=0.86)
    fig.suptitle(f"Rotating sharp-peak contours on triangular cells, t={final_time:g}")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_structured_field_comparison(
    initial: np.ndarray,
    numerical: np.ndarray,
    exact: np.ndarray,
    output: Path,
    final_time: float,
    domain: tuple[float, float, float, float],
    title: str,
) -> None:
    """Plot initial, numerical, exact and error fields on a rectangular domain."""
    xmin, xmax, ymin, ymax = domain
    error = numerical - exact
    maximum = max(
        float(np.max(np.abs(initial))),
        float(np.max(np.abs(numerical))),
        float(np.max(np.abs(exact))),
        1.0e-16,
    )
    error_maximum = max(float(np.max(np.abs(error))), 1.0e-16)
    ny, nx = numerical.shape
    x_edges = np.linspace(xmin, xmax, nx + 1)
    y_edges = np.linspace(ymin, ymax, ny + 1)

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
    panels = [
        (axes[0, 0], initial, "Initial field", 0.0, maximum, "viridis"),
        (axes[0, 1], numerical, f"Numerical field, t={final_time:g}", 0.0, maximum, "viridis"),
        (axes[1, 0], exact, f"Exact field, t={final_time:g}", 0.0, maximum, "viridis"),
        (axes[1, 1], error, "Numerical minus exact", -error_maximum, error_maximum, "coolwarm"),
    ]
    for axis, values, panel_title, lower, upper, cmap in panels:
        image = axis.pcolormesh(
            x_edges,
            y_edges,
            values,
            shading="auto",
            cmap=cmap,
            vmin=lower,
            vmax=upper,
        )
        axis.set_aspect("equal")
        axis.set_xlabel("x")
        axis.set_ylabel("y")
        axis.set_title(panel_title)
        fig.colorbar(image, ax=axis, shrink=0.86)
    fig.suptitle(title)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_structured_midline_profile(
    centres: list[tuple[float, float]],
    numerical: np.ndarray,
    exact: np.ndarray,
    output: Path,
    final_time: float,
    domain: tuple[float, float, float, float],
) -> None:
    """Plot the y=0 midline profile for the discontinuous diffusion case."""
    xmin, xmax, ymin, ymax = domain
    ny, nx = numerical.shape
    dy = (ymax - ymin) / ny
    mid_j = min(range(ny), key=lambda j: abs(ymin + (j + 0.5) * dy))
    x_values = [centres[mid_j * nx + i][0] for i in range(nx)]

    fig, axis = plt.subplots(figsize=(8, 4.8), constrained_layout=True)
    axis.plot(x_values, numerical[mid_j, :], "o-", markersize=3.0, label="numerical")
    axis.plot(x_values, exact[mid_j, :], "--", linewidth=1.5, label="exact")
    axis.set_xlabel("x at cell centres nearest y=0")
    axis.set_ylabel("phi")
    axis.set_title(f"Midline profile at t={final_time:g}")
    axis.grid(True, alpha=0.25)
    axis.legend()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_rotation_field_comparison(
    initial: np.ndarray,
    numerical: np.ndarray,
    output: Path,
    final_time: float,
) -> None:
    """Plot initial, final and cycle-change fields for solid rotation."""
    change = numerical - initial
    maximum = max(float(np.max(initial)), float(np.max(numerical)), 1.0e-16)
    change_maximum = max(float(np.max(np.abs(change))), 1.0e-16)
    ny, nx = numerical.shape
    x_edges = np.linspace(0.0, 1.0, nx + 1)
    y_edges = np.linspace(0.0, 1.0, ny + 1)

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.2), constrained_layout=True)
    panels = [
        (axes[0], initial, "Initial profile", 0.0, maximum, "viridis"),
        (axes[1], numerical, f"Final profile, t={final_time:g}", 0.0, maximum, "viridis"),
        (axes[2], change, "Final minus initial", -change_maximum, change_maximum, "coolwarm"),
    ]
    for axis, values, title, lower, upper, cmap in panels:
        image = axis.pcolormesh(
            x_edges,
            y_edges,
            values,
            shading="auto",
            cmap=cmap,
            vmin=lower,
            vmax=upper,
        )
        axis.set_aspect("equal")
        axis.set_xlabel("x")
        axis.set_ylabel("y")
        axis.set_title(title)
        fig.colorbar(image, ax=axis, shrink=0.82)
    fig.suptitle("Solid-rotation advection profile")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_rotation_final_contour(
    numerical: np.ndarray,
    output: Path,
    final_time: float,
) -> None:
    """Plot the final contour requested by the solid-rotation benchmark."""
    ny, nx = numerical.shape
    x = (np.arange(nx) + 0.5) / nx
    y = (np.arange(ny) + 0.5) / ny
    x_grid, y_grid = np.meshgrid(x, y)

    fig, axis = plt.subplots(figsize=(6.0, 5.4), constrained_layout=True)
    levels = np.linspace(0.05, max(float(np.max(numerical)), 0.05), 16)
    filled = axis.contourf(x_grid, y_grid, numerical, levels=levels, cmap="viridis")
    axis.contour(x_grid, y_grid, numerical, levels=levels, colors="black", linewidths=0.35, alpha=0.55)
    axis.set_aspect("equal")
    axis.set_xlabel("x")
    axis.set_ylabel("y")
    axis.set_title(f"Solid rotation final contours, t={final_time:g}")
    fig.colorbar(filled, ax=axis, shrink=0.86)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_tri_rotation_field_comparison(
    triangulation: mtri.Triangulation,
    initial: np.ndarray,
    numerical: np.ndarray,
    triangle_to_cell: list[int],
    output: Path,
    final_time: float,
) -> None:
    """Plot initial, final and cycle difference on a triangular mesh."""
    initial_faces = initial[triangle_to_cell]
    numerical_faces = numerical[triangle_to_cell]
    change_faces = (numerical - initial)[triangle_to_cell]
    maximum = max(
        float(np.max(np.abs(initial_faces))),
        float(np.max(np.abs(numerical_faces))),
        1.0e-16,
    )
    change_maximum = max(float(np.max(np.abs(change_faces))), 1.0e-16)
    panels = [
        (initial_faces, "Initial profile", -maximum, maximum, "viridis"),
        (numerical_faces, f"Final profile, t={final_time:g}", -maximum, maximum, "viridis"),
        (change_faces, "Final minus initial", -change_maximum, change_maximum, "coolwarm"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.2), constrained_layout=True)
    for axis, (values, title, lower, upper, cmap) in zip(axes, panels):
        image = axis.tripcolor(
            triangulation,
            facecolors=values,
            shading="flat",
            cmap=cmap,
            vmin=lower,
            vmax=upper,
        )
        axis.triplot(triangulation, color="black", linewidth=0.25, alpha=0.25)
        axis.set_aspect("equal")
        axis.set_xlabel("x")
        axis.set_ylabel("y")
        axis.set_title(title)
        fig.colorbar(image, ax=axis, shrink=0.82)
    fig.suptitle("Solid-rotation advection on triangular cells")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_tri_rotation_final_contour(
    triangulation: mtri.Triangulation,
    numerical: np.ndarray,
    triangle_to_cell: list[int],
    output: Path,
    final_time: float,
) -> None:
    """Plot final solid-rotation contours on the actual triangular mesh."""
    values = numerical[triangle_to_cell]
    maximum = max(float(np.max(values)), 0.05)
    levels = np.linspace(0.05, maximum, 16)
    node_values = np.zeros(len(triangulation.x), dtype=float)
    node_counts = np.zeros(len(triangulation.x), dtype=float)
    for triangle_index, nodes in enumerate(triangulation.triangles):
        node_values[nodes] += values[triangle_index]
        node_counts[nodes] += 1.0
    node_values /= np.maximum(node_counts, 1.0)

    fig, axis = plt.subplots(figsize=(6.0, 5.4), constrained_layout=True)
    filled = axis.tripcolor(
        triangulation,
        facecolors=values,
        shading="flat",
        cmap="viridis",
        vmin=0.0,
        vmax=maximum,
    )
    axis.tricontour(
        triangulation,
        node_values,
        levels=levels,
        colors="black",
        linewidths=0.35,
        alpha=0.55,
    )
    axis.triplot(triangulation, color="black", linewidth=0.18, alpha=0.2)
    axis.set_aspect("equal")
    axis.set_xlabel("x")
    axis.set_ylabel("y")
    axis.set_title(f"Solid rotation final contours, t={final_time:g}")
    fig.colorbar(filled, ax=axis, shrink=0.86)
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
    field_name: str = "T",
    figure_title: str = "Periodic linear-advection field comparison on triangular cells",
    independent_scales: bool = False,
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
        1.0e-16,
    )
    error_maximum = max(float(np.max(np.abs(error_faces))), 1.0e-16)
    initial_maximum = max(float(np.max(np.abs(initial_faces))), 1.0e-16)
    numerical_maximum = max(float(np.max(np.abs(numerical_faces))), 1.0e-40)
    exact_maximum = max(float(np.max(np.abs(exact_faces))), 1.0e-40)
    panels = [
        (
            initial_faces,
            f"Initial field {field_name}(x,y,0)",
            -initial_maximum if independent_scales else -maximum,
            initial_maximum if independent_scales else maximum,
        ),
        (
            numerical_faces,
            f"Numerical field, t={final_time:g}",
            -numerical_maximum if independent_scales else -maximum,
            numerical_maximum if independent_scales else maximum,
        ),
        (
            exact_faces,
            f"Exact field, t={final_time:g}",
            -exact_maximum if independent_scales else -maximum,
            exact_maximum if independent_scales else maximum,
        ),
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

    fig.suptitle(figure_title)
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
    field_name: str = "T",
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
    axis.set_ylabel(field_name)
    axis.set_title(f"Approximate diagonal {field_name} profile at t={final_time:g}")
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


def _read_problem(case: Path) -> str:
    metadata_path = case / "metadata.json"
    if not metadata_path.exists():
        return "sine_wave_advection"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    return str(metadata.get("problem", "sine_wave_advection"))


def _postprocess_poisson_case(
    case: Path, output_root: Path, target_co: float
) -> dict[str, object]:
    """Post-process the steady manufactured Poisson benchmark on either mesh."""
    case = case.resolve()
    metadata = json.loads((case / "metadata.json").read_text(encoding="utf-8"))
    field_name = str(metadata.get("scalarField", "phi"))
    domain_raw = metadata.get("domain", [0.0, 1.0, 0.0, 1.0])
    if not isinstance(domain_raw, list | tuple) or len(domain_raw) != 4:
        raise RuntimeError(f"Invalid Poisson domain in metadata.json: {case}")
    domain = tuple(float(value) for value in domain_raw)
    mesh_type = str(metadata.get("meshType", "quad"))

    # A steady OpenFOAM solve writes the converged field to time 0.
    # Keep the source field from 0.orig so the initial guess remains traceable.
    _, final_dir = latest_time(case)
    initial_values = read_scalar_field(case / "0.orig" / field_name)
    numerical_values = read_scalar_field(final_dir / field_name)

    if mesh_type == "quad":
        nx, ny = mesh_resolution(case)
        centres2d = structured_centres(nx, ny, domain)
        centres = [(x, y, 0.0) for x, y in centres2d]
        volumes = [
            (domain[1] - domain[0]) / nx
            * (domain[3] - domain[2]) / ny
            * THICKNESS
        ] * (nx * ny)
    elif mesh_type == "tri":
        centres, volumes = read_cell_geometry(case)
        geometry_metadata = read_tri_geometry_metadata(case / "mesh" / "mesh_geometry.json")
        triangulation = _triangulation_from_metadata(geometry_metadata)
        triangle_to_cell = _match_triangles_to_cells(centres, triangulation)
    else:
        raise NotImplementedError(f"Unsupported Poisson mesh type: {mesh_type}")

    if len(initial_values) != len(centres) or len(numerical_values) != len(centres):
        raise RuntimeError(
            f"Poisson field size does not match mesh in {case}: "
            f"cells={len(centres)}, initial={len(initial_values)}, "
            f"final={len(numerical_values)}"
        )

    exact_values = [
        poisson_exact_value(float(x), float(y)) for x, y, *_ in centres
    ]
    l1, l2, linf = normalized_errors(numerical_values, exact_values, volumes)
    initial = np.asarray(initial_values, dtype=float)
    numerical = np.asarray(numerical_values, dtype=float)
    exact = np.asarray(exact_values, dtype=float)
    volumes_array = np.asarray(volumes, dtype=float)
    initial_mass = float(np.dot(initial, volumes_array))
    final_mass = float(np.dot(numerical, volumes_array))
    mass_scale = float(np.dot(np.abs(exact), volumes_array))
    solver_log_path = case / "log.poissonFoamStudent"
    solver_log = solver_log_path.read_text(encoding="utf-8") if solver_log_path.exists() else ""
    mesh_log_path = case / "log.checkMesh"
    mesh_log = mesh_log_path.read_text(encoding="utf-8") if mesh_log_path.exists() else ""
    resolution = int(metadata.get("resolution", 0))
    if resolution <= 0:
        resolution = int(round(1.0 / max(float(metadata.get("nominalH", 1.0)), 1.0e-30)))
    final_time = 0.0

    solver_family, case_name = _case_namespace(case)
    data_dir = solver_data_dir(solver_family, case_name, resolution)
    figure_dir = solver_figure_dir(solver_family, case_name, resolution)
    data_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    # Steady cases have no physical time history; retain a traceable one-row
    # record so downstream tools can use the same file layout as other studies.
    write_time_history(
        [
            {
                "time": final_time,
                "step": 0.0,
                "deltaT": 0.0,
                "maxCo": 0.0,
                "residualIntegral": 0.0,
                "minT": float(numerical.min()),
                "maxT": float(numerical.max()),
                "amplitude": 0.5 * float(numerical.max() - numerical.min()),
            }
        ],
        data_dir / "time_history.csv",
    )
    if mesh_type == "quad":
        initial_grid = initial.reshape(ny, nx)
        numerical_grid = numerical.reshape(ny, nx)
        exact_grid = exact.reshape(ny, nx)
        write_structured_field_data(
            [(x, y) for x, y in centres2d],
            initial_values,
            numerical_values,
            exact_values,
            data_dir / "field_data.csv",
            data_dir / "error_field.csv",
        )
        plot_structured_field_comparison(
            initial_grid,
            numerical_grid,
            exact_grid,
            figure_dir / "field_comparison.png",
            final_time,
            domain,
            "Poisson manufactured solution on quadrilateral cells",
        )
    else:
        write_tri_field_data(
            centres,
            volumes,
            initial_values,
            numerical_values,
            exact_values,
            data_dir / "field_data.csv",
            data_dir / "error_field.csv",
        )
        plot_tri_field_comparison(
            triangulation,
            initial,
            numerical,
            exact,
            triangle_to_cell,
            figure_dir / "field_comparison.png",
            final_time,
            field_name="phi",
            figure_title="Poisson manufactured solution on triangular cells",
        )

    summary: dict[str, object] = {
        "case": str(case),
        "solverFamily": solver_family,
        "caseName": case_name,
        "mesh": mesh_type,
        "problem": "poisson_manufactured",
        "resolution": resolution,
        "nominalH": 1.0 / resolution,
        "nCells": len(centres),
        "finalTime": final_time,
        "finalDirectory": str(final_dir),
        "finalTimeError": 0.0,
        "normalizedL1": l1,
        "normalizedL2": l2,
        "normalizedLinf": linf,
        "initialMass": initial_mass,
        "finalMass": final_mass,
        "massChange": final_mass - initial_mass,
        "normalizedMassError": abs(final_mass - initial_mass) / mass_scale if mass_scale else 0.0,
        "maxAbsResidualIntegral": 0.0,
        "minCo": 0.0,
        "maxCo": 0.0,
        "targetCo": 0.0,
        "initialAmplitude": 0.5 * float(initial.max() - initial.min()),
        "finalAmplitude": 0.5 * float(numerical.max() - numerical.min()),
        "minFinal": float(numerical.min()),
        "maxFinal": float(numerical.max()),
        "maxAbsFinal": float(np.max(np.abs(numerical))),
        "timeSteps": 1,
        "meshOK": "Mesh OK." in mesh_log,
        "boundedBelowInitial": "",
        "boundedAboveInitial": "",
        "solverEnded": "End" in solver_log,
        "solverFatal": "FOAM FATAL ERROR" in solver_log,
        "linearSolver": metadata.get("linearSolver", ""),
        "linearTolerance": metadata.get("linearTolerance", ""),
        "nNonOrthogonalCorrectors": metadata.get("nNonOrthogonalCorrectors", ""),
    }
    (data_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    print(f"case={case}")
    print(f"resolution={resolution}")
    print(f"mesh={mesh_type}")
    print(f"nCells={len(centres)}")
    print(f"normalizedL1={l1:.16e}")
    print(f"normalizedL2={l2:.16e}")
    print(f"normalizedLinf={linf:.16e}")
    print(f"data={data_dir}")
    print(f"figures={figure_dir}")
    return summary


def _postprocess_quad_advection_case(
    case: Path, output_root: Path, target_co: float
) -> dict[str, object]:
    problem = _read_problem(case)
    if _read_problem(case) == "diffusion_discontinuity":
        return _postprocess_quad_diffusion_discontinuity_case(case, output_root, target_co)
    if _read_problem(case) == "poisson_manufactured":
        return _postprocess_poisson_case(case, output_root, target_co)
    if _read_problem(case) == "solid_rotation_advection":
        return _postprocess_quad_solid_rotation_case(case, output_root, target_co)

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

    solver_family, case_name = _case_namespace(case)
    data_dir = solver_data_dir(solver_family, case_name, int(nx))
    figure_dir = solver_figure_dir(solver_family, case_name, int(nx))
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
        "solverFamily": solver_family,
        "caseName": case_name,
        "mesh": "quad",
        "problem": problem,
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
        "solverEnded": (
            "End" in solver_log
            if problem == "diffusion_discontinuity"
            else "End" in solver_log and "Stage 5 time loop completed" in solver_log
        ),
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


def _postprocess_quad_diffusion_discontinuity_case(
    case: Path, output_root: Path, target_co: float
) -> dict[str, object]:
    """Post-process the structured discontinuous diffusion benchmark."""
    case = case.resolve()
    metadata = json.loads((case / "metadata.json").read_text(encoding="utf-8"))
    field_name = str(metadata.get("scalarField", "phi"))
    domain = tuple(float(value) for value in metadata.get("domain", [-5.0, 5.0, -5.0, 5.0]))
    if len(domain) != 4:
        raise RuntimeError(f"Invalid domain in metadata.json: {case}")
    domain4 = (domain[0], domain[1], domain[2], domain[3])

    nx, ny = mesh_resolution(case)
    final_time, final_dir = latest_time(case)
    initial_values = read_scalar_field(case / "0" / field_name)
    numerical_values = read_scalar_field(final_dir / field_name)
    if len(initial_values) != nx * ny or len(numerical_values) != nx * ny:
        raise RuntimeError(f"Field size does not match mesh in {case}")

    centres = structured_centres(nx, ny, domain4)
    exact_values_at_final = discontinuous_values(nx, ny, domain4, final_time)
    initial = np.asarray(initial_values, dtype=float).reshape(ny, nx)
    numerical = np.asarray(numerical_values, dtype=float).reshape(ny, nx)
    exact = np.asarray(exact_values_at_final, dtype=float).reshape(ny, nx)
    records = parse_diffusion_solver_log(case / "log.explicitDiffusionFoamStudent")

    dx = (domain4[1] - domain4[0]) / nx
    dy = (domain4[3] - domain4[2]) / ny
    cell_volume = dx * dy * THICKNESS
    l1, l2, linf = normalized_errors(
        numerical_values, exact_values_at_final, cell_volume
    )
    initial_mass = float(np.sum(initial) * cell_volume)
    final_mass = float(np.sum(numerical) * cell_volume)
    mass_scale = float(np.sum(np.abs(initial)) * cell_volume)
    mass_change = final_mass - initial_mass
    solver_log = (case / "log.explicitDiffusionFoamStudent").read_text(encoding="utf-8")
    mesh_log_path = case / "log.checkMesh"
    mesh_log = mesh_log_path.read_text(encoding="utf-8") if mesh_log_path.exists() else ""
    configured_end_time = control_value(case, "endTime", 0.2)

    solver_family, case_name = _case_namespace(case)
    data_dir = solver_data_dir(solver_family, case_name, int(nx))
    figure_dir = solver_figure_dir(solver_family, case_name, int(nx))
    data_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    time_history_records = [
        {
            "time": record.get("time", 0.0),
            "step": record.get("step", 0.0),
            "deltaT": record.get("deltaT", 0.0),
            "maxCo": target_co,
            "residualIntegral": "",
            "minT": record.get("minT", ""),
            "maxT": record.get("maxT", ""),
            "amplitude": record.get("amplitude", ""),
        }
        for record in records
    ]
    write_time_history(time_history_records, data_dir / "time_history.csv")
    write_structured_field_data(
        centres,
        initial_values,
        numerical_values,
        exact_values_at_final,
        data_dir / "field_data.csv",
        data_dir / "error_field.csv",
    )

    summary: dict[str, object] = {
        "case": str(case),
        "solverFamily": solver_family,
        "caseName": case_name,
        "mesh": "quad",
        "problem": "diffusion_discontinuity",
        "resolution": nx,
        "nominalH": max(dx, dy),
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
        "maxAbsResidualIntegral": "",
        "minCo": target_co,
        "maxCo": target_co,
        "targetCo": target_co,
        "initialAmplitude": float(initial.max()) - float(initial.min()),
        "finalAmplitude": float(numerical.max()) - float(numerical.min()),
        "minFinal": float(numerical.min()),
        "maxFinal": float(numerical.max()),
        "maxAbsFinal": float(np.max(np.abs(numerical))),
        "timeSteps": len(records),
        "meshOK": "Mesh OK." in mesh_log,
        "boundedBelowInitial": float(numerical.min()) >= float(initial.min()) - 1.0e-12,
        "boundedAboveInitial": float(numerical.max()) <= float(initial.max()) + 1.0e-12,
        "solverEnded": "End" in solver_log,
        "solverFatal": "FOAM FATAL ERROR" in solver_log,
    }
    (data_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    plot_structured_field_comparison(
        initial,
        numerical,
        exact,
        figure_dir / "field_comparison.png",
        final_time,
        domain4,
        "Discontinuous diffusion field comparison",
    )
    plot_structured_midline_profile(
        centres,
        numerical,
        exact,
        figure_dir / "midline_profile.png",
        final_time,
        domain4,
    )
    plot_cfl_history(time_history_records, figure_dir / "diffusion_step_history.png", target_co)

    print(f"case={case}")
    print(f"resolution={nx}")
    print(f"finalTime={final_time:.16g}")
    print(f"normalizedL1={l1:.16e}")
    print(f"normalizedL2={l2:.16e}")
    print(f"normalizedLinf={linf:.16e}")
    print(f"data={data_dir}")
    print(f"figures={figure_dir}")
    return summary


def _postprocess_quad_solid_rotation_case(
    case: Path, output_root: Path, target_co: float
) -> dict[str, object]:
    """Post-process the structured solid-rotation benchmark."""
    case = case.resolve()
    nx, ny = mesh_resolution(case)
    final_time, final_dir = latest_time(case)
    initial_values = read_scalar_field(case / "0" / "T")
    numerical_values = read_scalar_field(final_dir / "T")
    if len(initial_values) != nx * ny or len(numerical_values) != nx * ny:
        raise RuntimeError(f"Field size does not match mesh in {case}")

    initial = np.asarray(initial_values, dtype=float).reshape(ny, nx)
    numerical = np.asarray(numerical_values, dtype=float).reshape(ny, nx)
    records = parse_solver_log(case / "log.explicitAdvectionFoamStudent")
    cell_volume = THICKNESS / (nx * ny)
    initial_mass = float(np.sum(initial) * cell_volume)
    final_mass = float(np.sum(numerical) * cell_volume)
    mass_scale = float(np.sum(np.abs(initial)) * cell_volume)
    cycle_l1 = (
        float(np.sum(np.abs(numerical - initial)) * cell_volume) / mass_scale
        if mass_scale
        else 0.0
    )
    mass_change = final_mass - initial_mass
    solver_log = (case / "log.explicitDiffusionFoamStudent").read_text(encoding="utf-8")
    mesh_log_path = case / "log.checkMesh"
    mesh_log = mesh_log_path.read_text(encoding="utf-8") if mesh_log_path.exists() else ""
    configured_end_time = control_value(case, "endTime", 2.0 * math.pi)
    metadata = json.loads((case / "metadata.json").read_text(encoding="utf-8"))
    resolution = int(metadata["resolution"])

    solver_family, case_name = _case_namespace(case)
    data_dir = solver_data_dir(solver_family, case_name, int(nx))
    figure_dir = solver_figure_dir(solver_family, case_name, int(nx))
    data_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    write_time_history(records, data_dir / "time_history.csv")
    write_rotation_field_data(initial, numerical, data_dir / "field_data.csv")

    summary: dict[str, object] = {
        "case": str(case),
        "solverFamily": solver_family,
        "caseName": case_name,
        "mesh": "quad",
        "problem": "solid_rotation_advection",
        "resolution": nx,
        "nominalH": 1.0 / nx,
        "nCells": nx * ny,
        "finalTime": final_time,
        "finalDirectory": str(final_dir),
        "finalTimeError": abs(final_time - configured_end_time),
        "cycleL1AgainstInitial": cycle_l1,
        "initialMass": initial_mass,
        "finalMass": final_mass,
        "massChange": mass_change,
        "normalizedMassError": abs(mass_change) / mass_scale if mass_scale else 0.0,
        "minCo": min(record["maxCo"] for record in records),
        "maxCo": max(record["maxCo"] for record in records),
        "targetCo": target_co,
        "minInitial": float(initial.min()),
        "maxInitial": float(initial.max()),
        "minFinal": float(numerical.min()),
        "maxFinal": float(numerical.max()),
        "timeSteps": len(records),
        "meshOK": "Mesh OK." in mesh_log,
        "solverEnded": (
            "End" in solver_log
            if problem == "diffusion_discontinuity"
            else "End" in solver_log and "Stage 5 time loop completed" in solver_log
        ),
        "solverFatal": "FOAM FATAL ERROR" in solver_log,
    }
    (data_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    plot_rotation_field_comparison(
        initial,
        numerical,
        figure_dir / "field_comparison.png",
        final_time,
    )
    plot_rotation_final_contour(
        numerical,
        figure_dir / "contour_final.png",
        final_time,
    )
    plot_cfl_history(records, figure_dir / "cfl_history.png", target_co)

    print(f"case={case}")
    print(f"resolution={nx}")
    print(f"finalTime={final_time:.16g}")
    print(f"cycleL1AgainstInitial={cycle_l1:.16e}")
    print(f"data={data_dir}")
    print(f"figures={figure_dir}")
    return summary


def _postprocess_tri_solid_rotation_case(
    case: Path, output_root: Path, target_co: float
) -> dict[str, object]:
    """Post-process solid rotation on the actual triangular mesh."""
    case = case.resolve()
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

    initial = np.asarray(initial_values, dtype=float)
    numerical = np.asarray(numerical_values, dtype=float)
    volumes_array = np.asarray(volumes, dtype=float)
    records = parse_solver_log(case / "log.explicitAdvectionFoamStudent")
    initial_mass = float(np.dot(initial, volumes_array))
    final_mass = float(np.dot(numerical, volumes_array))
    mass_scale = float(np.dot(np.abs(initial), volumes_array))
    mass_change = final_mass - initial_mass
    cycle_l1 = (
        float(np.dot(np.abs(numerical - initial), volumes_array)) / mass_scale
        if mass_scale
        else 0.0
    )
    solver_log = (case / "log.explicitDiffusionFoamStudent").read_text(encoding="utf-8")
    mesh_log_path = case / "log.checkMesh"
    mesh_log = mesh_log_path.read_text(encoding="utf-8") if mesh_log_path.exists() else ""
    configured_end_time = control_value(case, "endTime", 2.0 * math.pi)
    metadata = json.loads((case / "metadata.json").read_text(encoding="utf-8"))
    resolution = int(metadata["resolution"])

    solver_family, case_name = _case_namespace(case)
    data_dir = solver_data_dir(solver_family, case_name, resolution)
    figure_dir = solver_figure_dir(solver_family, case_name, resolution)
    data_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    write_time_history(records, data_dir / "time_history.csv")
    write_tri_rotation_field_data(
        centres,
        volumes,
        initial_values,
        numerical_values,
        data_dir / "field_data.csv",
    )

    summary: dict[str, object] = {
        "case": str(case),
        "solverFamily": solver_family,
        "caseName": case_name,
        "mesh": "tri",
        "problem": "solid_rotation_advection",
        "resolution": resolution,
        "nominalH": 1.0 / resolution,
        "nCells": len(centres),
        "finalTime": final_time,
        "finalDirectory": str(final_dir),
        "finalTimeError": abs(final_time - configured_end_time),
        "cycleL1AgainstInitial": cycle_l1,
        "initialMass": initial_mass,
        "finalMass": final_mass,
        "massChange": mass_change,
        "normalizedMassError": abs(mass_change) / mass_scale if mass_scale else 0.0,
        "minCo": min(record["maxCo"] for record in records),
        "maxCo": max(record["maxCo"] for record in records),
        "targetCo": target_co,
        "minInitial": float(initial.min()),
        "maxInitial": float(initial.max()),
        "minFinal": float(numerical.min()),
        "maxFinal": float(numerical.max()),
        "timeSteps": len(records),
        "meshOK": "Mesh OK." in mesh_log,
        "solverEnded": "End" in solver_log,
        "solverFatal": "FOAM FATAL ERROR" in solver_log,
    }
    (data_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    plot_tri_rotation_field_comparison(
        triangulation,
        initial,
        numerical,
        triangle_to_cell,
        figure_dir / "field_comparison.png",
        final_time,
    )
    plot_tri_rotation_final_contour(
        triangulation,
        numerical,
        triangle_to_cell,
        figure_dir / "contour_final.png",
        final_time,
    )
    plot_cfl_history(records, figure_dir / "cfl_history.png", target_co)

    print(f"case={case}")
    print(f"resolution={summary['resolution']}")
    print("mesh=tri")
    print(f"nCells={len(centres)}")
    print(f"finalTime={final_time:.16g}")
    print(f"cycleL1AgainstInitial={cycle_l1:.16e}")
    print(f"data={data_dir}")
    print(f"figures={figure_dir}")
    return summary


def _postprocess_tri_advection_diffusion_case(
    case: Path, output_root: Path, target_co: float
) -> dict[str, object]:
    """Post-process the third problem's first example on triangular cells."""
    case = case.resolve()
    metadata = json.loads((case / "metadata.json").read_text(encoding="utf-8"))
    field_name = str(metadata.get("scalarField", "phi"))
    velocity = _case_velocity(case)
    mu = float(metadata.get("mu", 1.0))
    resolution = int(metadata["resolution"])

    centres, volumes = read_cell_geometry(case)
    geometry_metadata = read_tri_geometry_metadata(case / "mesh" / "mesh_geometry.json")
    triangulation = _triangulation_from_metadata(geometry_metadata)
    triangle_to_cell = _match_triangles_to_cells(centres, triangulation)
    final_time, final_dir = latest_time(case)
    initial_values = read_scalar_field(case / "0" / field_name)
    numerical_values = read_scalar_field(final_dir / field_name)
    if len(initial_values) != len(centres) or len(numerical_values) != len(centres):
        raise RuntimeError(
            f"Field size does not match triangular mesh in {case}: "
            f"cells={len(centres)}, initial={len(initial_values)}, "
            f"final={len(numerical_values)}"
        )

    exact_values_at_final = [
        advection_diffusion_exact_value(float(x), float(y), final_time, velocity, mu)
        for x, y, _ in centres
    ]
    records = parse_advection_diffusion_solver_log(
        case / "log.explicitAdvectionDiffusionFoamStudent"
    )
    initial = np.asarray(initial_values, dtype=float)
    numerical = np.asarray(numerical_values, dtype=float)
    exact = np.asarray(exact_values_at_final, dtype=float)
    volumes_array = np.asarray(volumes, dtype=float)

    l1, l2, linf = normalized_errors(numerical_values, exact_values_at_final, volumes)
    total_volume = float(np.sum(volumes_array))
    absolute_l1 = float(np.dot(np.abs(numerical - exact), volumes_array) / total_volume)
    absolute_l2 = float(
        math.sqrt(np.dot((numerical - exact) ** 2, volumes_array) / total_volume)
    )
    absolute_linf = float(np.max(np.abs(numerical - exact)))
    exact_amplitude = math.exp(-8.0 * math.pi * math.pi * mu * final_time)
    initial_mass = float(np.dot(initial, volumes_array))
    final_mass = float(np.dot(numerical, volumes_array))
    mass_scale = float(np.dot(np.abs(initial), volumes_array))
    mass_change = final_mass - initial_mass
    solver_log = (case / "log.explicitAdvectionDiffusionFoamStudent").read_text(
        encoding="utf-8"
    )
    mesh_log_path = case / "log.checkMesh"
    mesh_log = mesh_log_path.read_text(encoding="utf-8") if mesh_log_path.exists() else ""
    configured_end_time = control_value(case, "endTime", 1.0)

    solver_family, case_name = _case_namespace(case)
    data_dir = solver_data_dir(solver_family, case_name, resolution)
    figure_dir = solver_figure_dir(solver_family, case_name, resolution)
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
        "solverFamily": solver_family,
        "caseName": case_name,
        "mesh": "tri",
        "problem": "sine_wave_advection_diffusion",
        "resolution": resolution,
        "nominalH": 1.0 / resolution,
        "nCells": len(centres),
        "finalTime": final_time,
        "finalDirectory": str(final_dir),
        "finalTimeError": abs(final_time - configured_end_time),
        "normalizedL1": l1,
        "normalizedL2": l2,
        "normalizedLinf": linf,
        "absoluteL1": absolute_l1,
        "absoluteL2": absolute_l2,
        "absoluteLinf": absolute_linf,
        "exactAmplitudeAtFinal": exact_amplitude,
        "initialMass": initial_mass,
        "finalMass": final_mass,
        "massChange": mass_change,
        "normalizedMassError": abs(mass_change) / mass_scale if mass_scale else 0.0,
        "maxAbsResidualIntegral": "",
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
        "solverEnded": "End" in solver_log,
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
        field_name="phi",
        figure_title="Periodic sine-wave advection-diffusion on triangular cells",
        independent_scales=True,
    )
    plot_tri_diagonal_profile(
        centres,
        numerical,
        exact,
        figure_dir / "diagonal_profile.png",
        final_time,
        resolution,
        field_name="phi",
    )
    plot_advection_diffusion_amplitude_history(
        records, figure_dir / "amplitude_history.png", mu
    )
    plot_advection_diffusion_stability_history(
        records, figure_dir / "advection_diffusion_stability_history.png", target_co
    )

    print(f"case={case}")
    print(f"resolution={resolution}")
    print("mesh=tri")
    print(f"nCells={len(centres)}")
    print(f"finalTime={final_time:.16g}")
    print(f"normalizedL1={l1:.16e}")
    print(f"normalizedL2={l2:.16e}")
    print(f"normalizedLinf={linf:.16e}")
    print(f"data={data_dir}")
    print(f"figures={figure_dir}")
    return summary


def _postprocess_tri_rotating_peak_advection_diffusion_case(
    case: Path, output_root: Path, target_co: float
) -> dict[str, object]:
    """Post-process the third problem's rotating sharp-peak case on triangles."""
    case = case.resolve()
    metadata = json.loads((case / "metadata.json").read_text(encoding="utf-8"))
    field_name = str(metadata.get("scalarField", "phi"))
    mu = float(metadata.get("mu", 1.0e-3))
    resolution = int(metadata["resolution"])
    diffusion_start_time, initial_center = _rotating_peak_parameters(metadata)

    centres, volumes = read_cell_geometry(case)
    geometry_metadata = read_tri_geometry_metadata(case / "mesh" / "mesh_geometry.json")
    triangulation = _triangulation_from_metadata(geometry_metadata)
    triangle_to_cell = _match_triangles_to_cells(centres, triangulation)
    final_time, final_dir = latest_time(case)
    initial_values = read_scalar_field(case / "0" / field_name)
    numerical_values = read_scalar_field(final_dir / field_name)
    exact_values_at_final = rotating_peak_values_at_centres(
        centres,
        final_time,
        mu,
        diffusion_start_time,
        initial_center,
    )
    records = parse_advection_diffusion_solver_log(
        case / "log.explicitAdvectionDiffusionFoamStudent"
    )

    initial = np.asarray(initial_values, dtype=float)
    numerical = np.asarray(numerical_values, dtype=float)
    exact = np.asarray(exact_values_at_final, dtype=float)
    volumes_array = np.asarray(volumes, dtype=float)
    l1, l2, linf = normalized_errors(numerical_values, exact_values_at_final, volumes)
    initial_mass = float(np.dot(initial, volumes_array))
    final_mass = float(np.dot(numerical, volumes_array))
    mass_scale = float(np.dot(np.abs(initial), volumes_array))
    mass_change = final_mass - initial_mass
    solver_log = (case / "log.explicitAdvectionDiffusionFoamStudent").read_text(
        encoding="utf-8"
    )
    mesh_log_path = case / "log.checkMesh"
    mesh_log = mesh_log_path.read_text(encoding="utf-8") if mesh_log_path.exists() else ""
    configured_end_time = control_value(case, "endTime", 2.0 * math.pi)

    solver_family, case_name = _case_namespace(case)
    data_dir = solver_data_dir(solver_family, case_name, resolution)
    figure_dir = solver_figure_dir(solver_family, case_name, resolution)
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

    final_center = (
        initial_center[0] * math.cos(final_time) - initial_center[1] * math.sin(final_time),
        initial_center[0] * math.sin(final_time) + initial_center[1] * math.cos(final_time),
    )
    summary: dict[str, object] = {
        "case": str(case),
        "solverFamily": solver_family,
        "caseName": case_name,
        "mesh": "tri",
        "problem": "rotating_peak_advection_diffusion",
        "resolution": resolution,
        "nominalH": 2.0 / resolution,
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
        "maxAbsResidualIntegral": "",
        "minCo": min(record["maxCo"] for record in records),
        "maxCo": max(record["maxCo"] for record in records),
        "targetCo": target_co,
        "initialAmplitude": float(initial.max()) - float(initial.min()),
        "finalAmplitude": float(numerical.max()) - float(numerical.min()),
        "minFinal": float(numerical.min()),
        "maxFinal": float(numerical.max()),
        "maxAbsFinal": float(np.max(np.abs(numerical))),
        "timeSteps": len(records),
        "meshOK": "Mesh OK." in mesh_log,
        "boundedBelowInitial": float(numerical.min()) >= -1.0e-12,
        "boundedAboveInitial": float(numerical.max()) <= float(initial.max()) + 1.0e-12,
        "solverEnded": "End" in solver_log,
        "solverFatal": "FOAM FATAL ERROR" in solver_log,
        "diffusionStartTime": diffusion_start_time,
        "initialCenter": list(initial_center),
        "finalCenterExact": list(final_center),
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
        field_name="phi",
        figure_title="Rotating sharp-peak advection-diffusion on triangular cells",
        independent_scales=True,
    )
    plot_rotating_peak_tri_final_contour(
        triangulation,
        numerical,
        exact,
        triangle_to_cell,
        figure_dir / "contour_final.png",
        final_time,
    )
    plot_tri_diagonal_profile(
        centres,
        numerical,
        exact,
        figure_dir / "diagonal_profile.png",
        final_time,
        resolution,
        field_name="phi",
    )
    plot_advection_diffusion_stability_history(
        records, figure_dir / "advection_diffusion_stability_history.png", target_co
    )

    print(f"case={case}")
    print(f"resolution={resolution}")
    print("mesh=tri")
    print(f"nCells={len(centres)}")
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
    problem = _read_problem(case)
    if problem == "poisson_manufactured":
        return _postprocess_poisson_case(case, output_root, target_co)
    if problem == "solid_rotation_advection":
        return _postprocess_tri_solid_rotation_case(case, output_root, target_co)
    if problem == "sine_wave_advection_diffusion":
        return _postprocess_tri_advection_diffusion_case(case, output_root, target_co)
    if problem == "rotating_peak_advection_diffusion":
        return _postprocess_tri_rotating_peak_advection_diffusion_case(
            case, output_root, target_co
        )
    case = case.resolve()
    field_name = str(json.loads((case / "metadata.json").read_text(encoding="utf-8")).get("scalarField", "T"))
    centres, volumes = read_cell_geometry(case)
    geometry_metadata = read_tri_geometry_metadata(case / "mesh" / "mesh_geometry.json")
    triangulation = _triangulation_from_metadata(geometry_metadata)
    triangle_to_cell = _match_triangles_to_cells(centres, triangulation)

    final_time, final_dir = latest_time(case)
    initial_values = read_scalar_field(case / "0" / field_name)
    numerical_values = read_scalar_field(final_dir / field_name)
    if len(initial_values) != len(centres) or len(numerical_values) != len(centres):
        raise RuntimeError(
            f"Field size does not match triangular mesh in {case}: "
            f"cells={len(centres)}, initial={len(initial_values)}, "
            f"final={len(numerical_values)}"
        )

    if problem == "diffusion_discontinuity":
        exact_values_at_final = [
            discontinuous_exact_value(float(x), float(y), final_time)
            for x, y, *_ in centres
        ]
        records = parse_diffusion_solver_log(case / "log.explicitDiffusionFoamStudent")
    elif problem == "diffusion_gaussian":
        mu = float(json.loads((case / "metadata.json").read_text(encoding="utf-8")).get("mu", 1.0))
        exact_values_at_final = [
            gaussian_exact_value(float(x), float(y), final_time, mu)
            for x, y, *_ in centres
        ]
        records = parse_diffusion_solver_log(case / "log.explicitDiffusionFoamStudent")
    else:
        velocity = _case_velocity(case)
        exact_values_at_final = exact_values_at_centres(centres, final_time, velocity)
        records = parse_solver_log(case / "log.explicitAdvectionFoamStudent")

    initial = np.asarray(initial_values, dtype=float)
    numerical = np.asarray(numerical_values, dtype=float)
    exact = np.asarray(exact_values_at_final, dtype=float)
    volumes_array = np.asarray(volumes, dtype=float)

    l1, l2, linf = normalized_errors(
        numerical_values, exact_values_at_final, volumes
    )
    initial_mass = float(np.dot(initial, volumes_array))
    final_mass = float(np.dot(numerical, volumes_array))
    mass_scale = float(np.dot(np.abs(initial), volumes_array))
    mass_change = final_mass - initial_mass
    solver_log = (case / "log.explicitDiffusionFoamStudent").read_text(encoding="utf-8")
    mesh_log_path = case / "log.checkMesh"
    mesh_log = mesh_log_path.read_text(encoding="utf-8") if mesh_log_path.exists() else ""
    configured_end_time = control_value(case, "endTime", 1.0)
    case_metadata = json.loads((case / "metadata.json").read_text(encoding="utf-8"))
    resolution = int(case_metadata.get("resolution", geometry_metadata["resolution"]))

    solver_family, case_name = _case_namespace(case)
    data_dir = solver_data_dir(solver_family, case_name, resolution)
    figure_dir = solver_figure_dir(solver_family, case_name, resolution)
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
        "solverFamily": solver_family,
        "caseName": case_name,
        "mesh": "tri",
        "problem": problem,
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
        "maxAbsResidualIntegral": (
            max(abs(record.get("residualIntegral", 0.0)) for record in records)
            if records and "residualIntegral" in records[0]
            else ""
        ),
        "minCo": (
            min(record["maxCo"] for record in records)
            if records and "maxCo" in records[0]
            else target_co
        ),
        "maxCo": (
            max(record["maxCo"] for record in records)
            if records and "maxCo" in records[0]
            else target_co
        ),
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
        "solverEnded": "End" in solver_log,
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
    if records and "maxCo" in records[0]:
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


def _postprocess_quad_case(
    case: Path, output_root: Path, target_co: float
) -> dict[str, object]:
    problem = _read_problem(case)
    if problem == "poisson_manufactured":
        return _postprocess_poisson_case(case, output_root, target_co)
    if problem == "diffusion_discontinuity":
        return _postprocess_quad_diffusion_discontinuity_case(case, output_root, target_co)
    if problem == "diffusion_gaussian":
        return _postprocess_quad_diffusion_gaussian_case(case, output_root, target_co)
    if problem == "solid_rotation_advection":
        return _postprocess_quad_solid_rotation_case(case, output_root, target_co)
    if problem == "sine_wave_advection_diffusion":
        return _postprocess_quad_advection_diffusion_case(case, output_root, target_co)
    if problem == "rotating_peak_advection_diffusion":
        return _postprocess_quad_rotating_peak_advection_diffusion_case(
            case, output_root, target_co
        )
    return _postprocess_quad_advection_case(case, output_root, target_co)


def _postprocess_quad_advection_diffusion_case(
    case: Path, output_root: Path, target_co: float
) -> dict[str, object]:
    """Post-process the third problem's first periodic sine-wave example."""
    case = case.resolve()
    metadata = json.loads((case / "metadata.json").read_text(encoding="utf-8"))
    field_name = str(metadata.get("scalarField", "phi"))
    domain_raw = metadata.get("domain", [0.0, 1.0, 0.0, 1.0])
    domain = tuple(float(value) for value in domain_raw)
    if len(domain) != 4:
        raise RuntimeError(f"Invalid domain in metadata.json: {case}")
    domain4 = (domain[0], domain[1], domain[2], domain[3])
    velocity = _case_velocity(case)
    mu = float(metadata.get("mu", 1.0))

    nx, ny = mesh_resolution(case)
    final_time, final_dir = latest_time(case)
    initial_values = read_scalar_field(case / "0" / field_name)
    numerical_values = read_scalar_field(final_dir / field_name)
    if len(initial_values) != nx * ny or len(numerical_values) != nx * ny:
        raise RuntimeError(f"Field size does not match mesh in {case}")

    exact_values_at_final = advection_diffusion_structured_values(
        nx, ny, domain4, final_time, velocity, mu
    )
    initial = np.asarray(initial_values, dtype=float).reshape(ny, nx)
    numerical = np.asarray(numerical_values, dtype=float).reshape(ny, nx)
    exact = np.asarray(exact_values_at_final, dtype=float).reshape(ny, nx)
    records = parse_advection_diffusion_solver_log(
        case / "log.explicitAdvectionDiffusionFoamStudent"
    )

    dx = (domain4[1] - domain4[0]) / nx
    dy = (domain4[3] - domain4[2]) / ny
    cell_volume = dx * dy * THICKNESS
    l1, l2, linf = normalized_errors(numerical_values, exact_values_at_final, cell_volume)
    absolute_l1 = float(np.sum(np.abs(numerical - exact)) * cell_volume / (dx * dy * nx * ny * THICKNESS))
    absolute_l2 = float(np.sqrt(np.sum((numerical - exact) ** 2) * cell_volume / (dx * dy * nx * ny * THICKNESS)))
    absolute_linf = float(np.max(np.abs(numerical - exact)))
    exact_amplitude = math.exp(-8.0 * math.pi * math.pi * mu * final_time)
    initial_mass = float(np.sum(initial) * cell_volume)
    final_mass = float(np.sum(numerical) * cell_volume)
    mass_scale = float(np.sum(np.abs(initial)) * cell_volume)
    mass_change = final_mass - initial_mass
    solver_log = (case / "log.explicitAdvectionDiffusionFoamStudent").read_text(
        encoding="utf-8"
    )
    mesh_log_path = case / "log.checkMesh"
    mesh_log = mesh_log_path.read_text(encoding="utf-8") if mesh_log_path.exists() else ""
    configured_end_time = control_value(case, "endTime", 1.0)

    solver_family, case_name = _case_namespace(case)
    data_dir = solver_data_dir(solver_family, case_name, nx)
    figure_dir = solver_figure_dir(solver_family, case_name, nx)
    data_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    write_time_history(records, data_dir / "time_history.csv")
    write_structured_field_data(
        [
            (
                domain4[0] + (i + 0.5) * dx,
                domain4[2] + (j + 0.5) * dy,
            )
            for j in range(ny)
            for i in range(nx)
        ],
        initial_values,
        numerical_values,
        exact_values_at_final,
        data_dir / "field_data.csv",
        data_dir / "error_field.csv",
    )

    summary: dict[str, object] = {
        "case": str(case),
        "solverFamily": solver_family,
        "caseName": case_name,
        "mesh": "quad",
        "problem": "sine_wave_advection_diffusion",
        "resolution": nx,
        "nominalH": max(dx, dy),
        "nCells": nx * ny,
        "finalTime": final_time,
        "finalDirectory": str(final_dir),
        "finalTimeError": abs(final_time - configured_end_time),
        "normalizedL1": l1,
        "normalizedL2": l2,
        "normalizedLinf": linf,
        "absoluteL1": absolute_l1,
        "absoluteL2": absolute_l2,
        "absoluteLinf": absolute_linf,
        "exactAmplitudeAtFinal": exact_amplitude,
        "initialMass": initial_mass,
        "finalMass": final_mass,
        "massChange": mass_change,
        "normalizedMassError": abs(mass_change) / mass_scale if mass_scale else 0.0,
        "maxAbsResidualIntegral": "",
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
        "solverEnded": "End" in solver_log,
        "solverFatal": "FOAM FATAL ERROR" in solver_log,
    }
    (data_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    plot_advection_diffusion_field_comparison(
        initial, numerical, exact, figure_dir / "field_comparison.png", final_time, domain4
    )
    plot_advection_diffusion_diagonal_profile(
        numerical, exact, figure_dir / "diagonal_profile.png", final_time
    )
    plot_advection_diffusion_amplitude_history(
        records, figure_dir / "amplitude_history.png", mu
    )
    plot_advection_diffusion_stability_history(
        records, figure_dir / "advection_diffusion_stability_history.png", target_co
    )

    print(f"case={case}")
    print(f"resolution={nx}")
    print(f"finalTime={final_time:.16g}")
    print(f"normalizedL1={l1:.16e}")
    print(f"normalizedL2={l2:.16e}")
    print(f"normalizedLinf={linf:.16e}")
    print(f"data={data_dir}")
    print(f"figures={figure_dir}")
    return summary


def _postprocess_quad_rotating_peak_advection_diffusion_case(
    case: Path, output_root: Path, target_co: float
) -> dict[str, object]:
    """Post-process the structured rotating sharp-peak advection-diffusion case."""
    case = case.resolve()
    metadata = json.loads((case / "metadata.json").read_text(encoding="utf-8"))
    field_name = str(metadata.get("scalarField", "phi"))
    domain_raw = metadata.get("domain", [-1.0, 1.0, -1.0, 1.0])
    domain = tuple(float(value) for value in domain_raw)
    if len(domain) != 4:
        raise RuntimeError(f"Invalid domain in metadata.json: {case}")
    domain4 = (domain[0], domain[1], domain[2], domain[3])
    mu = float(metadata.get("mu", 1.0e-3))
    diffusion_start_time, initial_center = _rotating_peak_parameters(metadata)

    nx, ny = mesh_resolution(case)
    final_time, final_dir = latest_time(case)
    initial_values = read_scalar_field(case / "0" / field_name)
    numerical_values = read_scalar_field(final_dir / field_name)
    if len(initial_values) != nx * ny or len(numerical_values) != nx * ny:
        raise RuntimeError(f"Field size does not match mesh in {case}")

    centres = structured_centres(nx, ny, domain4)
    exact_values_at_final = rotating_peak_structured_values(
        nx,
        ny,
        domain4,
        final_time,
        mu,
        diffusion_start_time,
        initial_center,
    )
    initial = np.asarray(initial_values, dtype=float).reshape(ny, nx)
    numerical = np.asarray(numerical_values, dtype=float).reshape(ny, nx)
    exact = np.asarray(exact_values_at_final, dtype=float).reshape(ny, nx)
    records = parse_advection_diffusion_solver_log(
        case / "log.explicitAdvectionDiffusionFoamStudent"
    )

    dx = (domain4[1] - domain4[0]) / nx
    dy = (domain4[3] - domain4[2]) / ny
    cell_volume = dx * dy * THICKNESS
    l1, l2, linf = normalized_errors(
        numerical_values, exact_values_at_final, cell_volume
    )
    initial_mass = float(np.sum(initial) * cell_volume)
    final_mass = float(np.sum(numerical) * cell_volume)
    mass_scale = float(np.sum(np.abs(initial)) * cell_volume)
    mass_change = final_mass - initial_mass
    solver_log = (case / "log.explicitAdvectionDiffusionFoamStudent").read_text(
        encoding="utf-8"
    )
    mesh_log_path = case / "log.checkMesh"
    mesh_log = mesh_log_path.read_text(encoding="utf-8") if mesh_log_path.exists() else ""
    configured_end_time = control_value(case, "endTime", 2.0 * math.pi)
    final_center = (
        initial_center[0] * math.cos(final_time) - initial_center[1] * math.sin(final_time),
        initial_center[0] * math.sin(final_time) + initial_center[1] * math.cos(final_time),
    )

    solver_family, case_name = _case_namespace(case)
    data_dir = solver_data_dir(solver_family, case_name, nx)
    figure_dir = solver_figure_dir(solver_family, case_name, nx)
    data_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    write_time_history(records, data_dir / "time_history.csv")
    write_structured_field_data(
        centres,
        initial_values,
        numerical_values,
        exact_values_at_final,
        data_dir / "field_data.csv",
        data_dir / "error_field.csv",
    )

    summary: dict[str, object] = {
        "case": str(case),
        "solverFamily": solver_family,
        "caseName": case_name,
        "mesh": "quad",
        "problem": "rotating_peak_advection_diffusion",
        "resolution": nx,
        "nominalH": max(dx, dy),
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
        "maxAbsResidualIntegral": "",
        "minCo": min(record["maxCo"] for record in records),
        "maxCo": max(record["maxCo"] for record in records),
        "targetCo": target_co,
        "initialAmplitude": float(initial.max()) - float(initial.min()),
        "finalAmplitude": float(numerical.max()) - float(numerical.min()),
        "minFinal": float(numerical.min()),
        "maxFinal": float(numerical.max()),
        "maxAbsFinal": float(np.max(np.abs(numerical))),
        "timeSteps": len(records),
        "meshOK": "Mesh OK." in mesh_log,
        "boundedBelowInitial": float(numerical.min()) >= -1.0e-12,
        "boundedAboveInitial": float(numerical.max()) <= float(initial.max()) + 1.0e-12,
        "solverEnded": "End" in solver_log,
        "solverFatal": "FOAM FATAL ERROR" in solver_log,
        "diffusionStartTime": diffusion_start_time,
        "initialCenter": list(initial_center),
        "finalCenterExact": list(final_center),
    }
    (data_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    plot_advection_diffusion_field_comparison(
        initial,
        numerical,
        exact,
        figure_dir / "field_comparison.png",
        final_time,
        domain4,
    )
    plot_rotating_peak_structured_final_contour(
        numerical,
        exact,
        figure_dir / "contour_final.png",
        final_time,
        domain4,
    )
    plot_rotating_peak_midline_profile(
        centres,
        numerical_values,
        exact_values_at_final,
        figure_dir / "peak_profile.png",
        final_time,
        final_center[1],
    )
    plot_advection_diffusion_stability_history(
        records, figure_dir / "advection_diffusion_stability_history.png", target_co
    )

    print(f"case={case}")
    print(f"resolution={nx}")
    print(f"finalTime={final_time:.16g}")
    print(f"normalizedL1={l1:.16e}")
    print(f"normalizedL2={l2:.16e}")
    print(f"normalizedLinf={linf:.16e}")
    print(f"data={data_dir}")
    print(f"figures={figure_dir}")
    return summary


def _postprocess_quad_diffusion_gaussian_case(
    case: Path, output_root: Path, target_co: float
) -> dict[str, object]:
    """Post-process the structured Gaussian diffusion benchmark."""
    case = case.resolve()
    metadata = json.loads((case / "metadata.json").read_text(encoding="utf-8"))
    field_name = str(metadata.get("scalarField", "phi"))
    domain = tuple(float(value) for value in metadata.get("domain", [-5.0, 5.0, -5.0, 5.0]))
    if len(domain) != 4:
        raise RuntimeError(f"Invalid domain in metadata.json: {case}")
    domain4 = (domain[0], domain[1], domain[2], domain[3])
    mu = float(metadata.get("mu", 1.0))

    nx, ny = mesh_resolution(case)
    final_time, final_dir = latest_time(case)
    initial_values = read_scalar_field(case / "0" / field_name)
    numerical_values = read_scalar_field(final_dir / field_name)
    if len(initial_values) != nx * ny or len(numerical_values) != nx * ny:
        raise RuntimeError(f"Field size does not match mesh in {case}")

    centres = structured_centres(nx, ny, domain4)
    exact_values_at_final = gaussian_values(nx, ny, domain4, final_time, mu)
    initial = np.asarray(initial_values, dtype=float).reshape(ny, nx)
    numerical = np.asarray(numerical_values, dtype=float).reshape(ny, nx)
    exact = np.asarray(exact_values_at_final, dtype=float).reshape(ny, nx)
    records = parse_diffusion_solver_log(case / "log.explicitDiffusionFoamStudent")

    dx = (domain4[1] - domain4[0]) / nx
    dy = (domain4[3] - domain4[2]) / ny
    cell_volume = dx * dy * THICKNESS
    l1, l2, linf = normalized_errors(
        numerical_values, exact_values_at_final, cell_volume
    )
    initial_mass = float(np.sum(initial) * cell_volume)
    final_mass = float(np.sum(numerical) * cell_volume)
    mass_scale = float(np.sum(np.abs(initial)) * cell_volume)
    mass_change = final_mass - initial_mass
    solver_log = (case / "log.explicitDiffusionFoamStudent").read_text(encoding="utf-8")
    mesh_log_path = case / "log.checkMesh"
    mesh_log = mesh_log_path.read_text(encoding="utf-8") if mesh_log_path.exists() else ""
    configured_end_time = control_value(case, "endTime", 0.2)

    solver_family, case_name = _case_namespace(case)
    data_dir = solver_data_dir(solver_family, case_name, int(nx))
    figure_dir = solver_figure_dir(solver_family, case_name, int(nx))
    data_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    time_history_records = [
        {
            "time": record.get("time", 0.0),
            "step": record.get("step", 0.0),
            "deltaT": record.get("deltaT", 0.0),
            "maxCo": target_co,
            "residualIntegral": "",
            "minT": record.get("minT", ""),
            "maxT": record.get("maxT", ""),
            "amplitude": record.get("amplitude", ""),
        }
        for record in records
    ]
    write_time_history(time_history_records, data_dir / "time_history.csv")
    write_structured_field_data(
        centres,
        initial_values,
        numerical_values,
        exact_values_at_final,
        data_dir / "field_data.csv",
        data_dir / "error_field.csv",
    )

    summary: dict[str, object] = {
        "case": str(case),
        "solverFamily": solver_family,
        "caseName": case_name,
        "mesh": "quad",
        "problem": "diffusion_gaussian",
        "resolution": nx,
        "nominalH": max(dx, dy),
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
        "maxAbsResidualIntegral": "",
        "minCo": target_co,
        "maxCo": target_co,
        "targetCo": target_co,
        "initialAmplitude": float(initial.max()) - float(initial.min()),
        "finalAmplitude": float(numerical.max()) - float(numerical.min()),
        "minFinal": float(numerical.min()),
        "maxFinal": float(numerical.max()),
        "maxAbsFinal": float(np.max(np.abs(numerical))),
        "timeSteps": len(records),
        "meshOK": "Mesh OK." in mesh_log,
        "boundedBelowInitial": float(numerical.min()) >= float(initial.min()) - 1.0e-12,
        "boundedAboveInitial": float(numerical.max()) <= float(initial.max()) + 1.0e-12,
        "solverEnded": "End" in solver_log,
        "solverFatal": "FOAM FATAL ERROR" in solver_log,
    }
    (data_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    plot_structured_field_comparison(
        initial,
        numerical,
        exact,
        figure_dir / "field_comparison.png",
        final_time,
        domain4,
        "Gaussian diffusion field comparison",
    )
    plot_structured_midline_profile(
        centres,
        numerical,
        exact,
        figure_dir / "midline_profile.png",
        final_time,
        domain4,
    )
    plot_cfl_history(time_history_records, figure_dir / "diffusion_step_history.png", target_co)

    print(f"case={case}")
    print(f"resolution={nx}")
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
