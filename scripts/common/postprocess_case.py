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

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from advection_tools import (
    control_value,
    exact_values,
    latest_time,
    mesh_resolution,
    normalized_errors,
    parse_solver_log,
    read_scalar_field,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]
THICKNESS = 0.1


def default_case(default_case_name: str | None) -> Path:
    environment_case = os.environ.get("VIBEFLOW_CASE_DIR")
    if environment_case:
        return Path(environment_case).expanduser().resolve()
    if default_case_name:
        return PROJECT_ROOT / "cases" / default_case_name / "N20"
    raise SystemExit("Please provide --case-dir or set VIBEFLOW_CASE_DIR")


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


def postprocess_case(case: Path, output_root: Path, target_co: float) -> dict[str, object]:
    case = case.resolve()
    nx, ny = mesh_resolution(case)
    final_time, final_dir = latest_time(case)
    initial_values = read_scalar_field(case / "0" / "T")
    numerical_values = read_scalar_field(final_dir / "T")
    if len(initial_values) != nx * ny or len(numerical_values) != nx * ny:
        raise RuntimeError(f"Field size does not match mesh in {case}")

    exact_values_at_final = exact_values(nx, ny, final_time)
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
