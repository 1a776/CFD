#!/usr/bin/env python3

"""Post-process the fifth-problem equilateral triangular cavity."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/student_project_matplotlib")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
from scipy.ndimage import maximum_filter, minimum_filter
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = PROJECT_ROOT / "data" / "05_navier_stokes_equation" / "cases"
FIGURE_ROOT = PROJECT_ROOT / "figures" / "05_navier_stokes_equation" / "cases"

REFERENCE_U = {
    100: [(0.0, 1.0), (-0.0571, 0.81766), (-0.0714, 0.77171),
          (-0.0857, 0.72676), (-0.1, 0.68320), (-0.4429, 0.16395),
          (-0.8714, -0.16653), (-1.1857, -0.29476), (-1.3143, -0.27349),
          (-1.5, -0.18897), (-1.5857, -0.14445), (-1.6429, -0.11669),
          (-1.6714, -0.10376), (-1.6857, -0.09758), (-1.9429, -0.02046),
          (-1.9571, -0.01798), (-2.1143, -0.00032), (-2.1857, 0.00312),
          (-2.2, 0.00355), (-2.3, 0.00481), (-3.0, 0.0)],
    200: [(0.0, 1.0), (-0.0571, 0.78999), (-0.0714, 0.73771),
          (-0.0857, 0.68813), (-0.1, 0.64201), (-0.4429, 0.23983),
          (-0.8714, -0.10758), (-1.1857, -0.33718), (-1.3143, -0.37775),
          (-1.5, -0.29791), (-1.5857, -0.22425), (-1.6429, -0.17469),
          (-1.6714, -0.15138), (-1.6857, -0.14028), (-1.9429, -0.01231),
          (-1.9571, -0.00880), (-2.1143, 0.01304), (-2.1857, 0.01543),
          (-2.2, 0.01553), (-2.3, 0.01371), (-3.0, 0.0)],
    500: [(0.0, 1.0), (-0.0571, 0.73838), (-0.0714, 0.67965),
          (-0.0857, 0.62977), (-0.1, 0.58927), (-0.4429, 0.28014),
          (-0.8714, -0.05682), (-1.1857, -0.26938), (-1.3143, -0.35935),
          (-1.5, -0.45972), (-1.5857, -0.41095), (-1.6429, -0.33494),
          (-1.6714, -0.28993), (-1.6857, -0.26698), (-1.9429, -0.01293),
          (-1.9571, -0.00764), (-2.1143, 0.02439), (-2.1857, 0.02737),
          (-2.2, 0.02732), (-2.3, 0.02277), (-3.0, 0.0)],
}
REFERENCE_V = {
    100: [(1.1547, 0.0), (1.0854, -0.14663), (1.0623, -0.18510),
          (1.0277, -0.23363), (0.9815, -0.28188), (0.9122, -0.32187),
          (0.8545, -0.33011), (0.0, 0.10458), (-0.4734, 0.18995),
          (-0.4965, 0.18981), (-0.6813, 0.17807), (-0.7852, 0.16183),
          (-0.8545, 0.14537), (-0.8891, 0.13501), (-0.9122, 0.12723),
          (-1.1547, 0.0)],
    200: [(1.1547, 0.0), (1.0854, -0.20538), (1.0623, -0.25872),
          (1.0277, -0.32360), (0.9815, -0.38138), (0.9122, -0.41260),
          (0.8545, -0.40093), (0.0, 0.09033), (-0.4734, 0.26513),
          (-0.4965, 0.26566), (-0.6813, 0.23822), (-0.7852, 0.20468),
          (-0.8545, 0.17678), (-0.8891, 0.16088), (-0.9122, 0.14946),
          (-1.1547, 0.0)],
    500: [(1.1547, 0.0), (1.0854, -0.33354), (1.0623, -0.40554),
          (1.0277, -0.47331), (0.9815, -0.49887), (0.9122, -0.45508),
          (0.8545, -0.39707), (0.0, 0.05136), (-0.4734, 0.28742),
          (-0.4965, 0.29968), (-0.6813, 0.35731), (-0.7852, 0.32956),
          (-0.8545, 0.28510), (-0.8891, 0.25700), (-0.9122, 0.23645),
          (-1.1547, 0.0)],
}
REFERENCE_PRIMARY = {
    100: (-0.2482, -1.3669, 0.3315, -0.6445),
    200: (-0.2624, -1.2518, 0.2030, -0.7266),
    500: (-0.2774, -1.1791, 0.1319, -0.7793),
}


def _final_time(case: Path) -> str:
    candidates = []
    for path in case.iterdir():
        if not path.is_dir():
            continue
        try:
            value = float(path.name)
        except ValueError:
            continue
        if value > 0.0 and (path / "U").exists():
            candidates.append((value, path.name))
    if not candidates:
        raise RuntimeError(f"No positive time directory with U found in {case}; run the solver before post-processing.")
    return max(candidates)[1]


def _cell_centres(case: Path, time_name: str) -> np.ndarray:
    """Read centres from the final time, falling back to constant/C."""
    final_path = case / time_name / "C"
    if not final_path.exists():
        final_path = case / "constant" / "C"
    return _read_field(final_path, True)[:, :2]


def _read_field(path: Path, vector: bool) -> np.ndarray:
    text = path.read_text(encoding="utf-8")
    kind = "vector" if vector else "scalar"
    match = re.search(
        rf"internalField\s+nonuniform\s+List<{kind}>\s+(\d+)\s*\((.*?)\)\s*;",
        text,
        re.S,
    )
    if not match:
        raise RuntimeError(f"Cannot parse {path}")
    count = int(match.group(1))
    if vector:
        values = re.findall(
            r"\(\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\)",
            match.group(2),
        )
        result = np.asarray(values, dtype=float)
    else:
        result = np.asarray(re.findall(r"[-+0-9.eE]+", match.group(2)), dtype=float)[:, None]
    if len(result) != count:
        raise RuntimeError(f"{path}: expected {count} values, found {len(result)}")
    return result


def _grid(case: Path, field: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    centres = _cell_centres(case, _final_time(case))
    values = field[:, :2]
    x = np.linspace(-math.sqrt(3.0), math.sqrt(3.0), 240)
    y = np.linspace(-3.0, 0.0, 260)
    xx, yy = np.meshgrid(x, y)
    mask = (yy >= math.sqrt(3.0) * xx - 3.0 - 1e-9) & (
        yy >= -math.sqrt(3.0) * xx - 3.0 - 1e-9
    )
    ux = griddata(centres, values[:, 0], (xx, yy), method="linear")
    uy = griddata(centres, values[:, 1], (xx, yy), method="linear")
    ux[~mask] = np.nan
    uy[~mask] = np.nan
    return xx, yy, np.stack((ux, uy), axis=0)


def _vorticity(xx: np.ndarray, yy: np.ndarray, ux: np.ndarray, uy: np.ndarray) -> np.ndarray:
    dy = float(yy[1, 0] - yy[0, 0])
    dx = float(xx[0, 1] - xx[0, 0])
    omega = np.full_like(ux, np.nan)
    omega[1:-1, 1:-1] = (
        (uy[1:-1, 2:] - uy[1:-1, :-2]) / (2.0 * dx)
        - (ux[2:, 1:-1] - ux[:-2, 1:-1]) / (2.0 * dy)
    )
    return omega


def _streamfunction(
    xx: np.ndarray, yy: np.ndarray, ux: np.ndarray, uy: np.ndarray, omega: np.ndarray
) -> np.ndarray:
    dy = float(yy[1, 0] - yy[0, 0])
    dx = float(xx[0, 1] - xx[0, 0])
    interior = np.isfinite(ux) & np.isfinite(uy)
    boundary = np.zeros_like(interior, dtype=bool)
    for j in range(ux.shape[0]):
        for i in range(ux.shape[1]):
            if not interior[j, i]:
                continue
            for dj, di in ((0, -1), (0, 1), (-1, 0), (1, 0)):
                nj, ni = j + dj, i + di
                if not (0 <= nj < interior.shape[0] and 0 <= ni < interior.shape[1]) or not interior[nj, ni]:
                    boundary[j, i] = True
                    break
    unknown = interior & ~boundary
    if not np.any(boundary) or not np.any(unknown):
        raise RuntimeError("Unable to build a valid streamfunction Poisson system for the triangular cavity.")
    ids = -np.ones_like(unknown, dtype=int)
    ids[unknown] = np.arange(np.count_nonzero(unknown))
    matrix = lil_matrix((np.count_nonzero(unknown), np.count_nonzero(unknown)))
    rhs = np.zeros(np.count_nonzero(unknown))
    for j, i in zip(*np.nonzero(unknown)):
        row = ids[j, i]
        matrix[row, row] = -2.0 / dx**2 - 2.0 / dy**2
        rhs[row] = -float(omega[j, i])
        for dj, di, coeff in ((0, -1, 1.0 / dx**2), (0, 1, 1.0 / dx**2),
                               (-1, 0, 1.0 / dy**2), (1, 0, 1.0 / dy**2)):
            nj, ni = j + dj, i + di
            if not (0 <= nj < unknown.shape[0] and 0 <= ni < unknown.shape[1]):
                continue
            if unknown[nj, ni]:
                matrix[row, ids[nj, ni]] = coeff
            elif boundary[nj, ni]:
                rhs[row] -= 0.0 * coeff
    psi = np.full_like(ux, np.nan)
    psi[boundary] = 0.0
    psi[unknown] = spsolve(matrix.tocsr(), rhs)
    if not np.any(np.isfinite(psi)):
        raise RuntimeError("Streamfunction solve returned only NaN values.")
    return psi


def postprocess(case: Path) -> Path:
    metadata = json.loads((case / "metadata.json").read_text(encoding="utf-8"))
    time_name = _final_time(case)
    field = _read_field(case / time_name / "U", True)
    xx, yy, grid = _grid(case, field)
    omega = _vorticity(xx, yy, grid[0], grid[1])
    psi = _streamfunction(xx, yy, grid[0], grid[1], omega)
    re_number = int(metadata["reynolds"])
    case_name = str(metadata["caseName"])
    data_dir = DATA_ROOT / case_name
    figure_dir = FIGURE_ROOT / case_name
    data_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    centreline_x = np.full(160, 0.0)
    centreline_y = np.linspace(-3.0, 0.0, 160)
    horizontal = griddata(
        _cell_centres(case, time_name),
        field[:, 0],
        (centreline_x, centreline_y),
        method="linear",
    )
    if np.any(~np.isfinite(horizontal)):
        fallback = griddata(
            _cell_centres(case, time_name),
            field[:, 0],
            (centreline_x, centreline_y),
            method="nearest",
        )
        horizontal = np.where(np.isfinite(horizontal), horizontal, fallback)
    horizontal[0] = 0.0
    horizontal[-1] = 1.0
    # At y=-1, the triangle inequalities give
    #     -2/sqrt(3) <= x <= 2/sqrt(3).
    # Sampling the enclosing square would silently trigger nearest-neighbour
    # fallback outside the physical domain and contaminate v(x,-1).
    line_half_width = 2.0 / math.sqrt(3.0)
    sample_x = np.linspace(-line_half_width + 1e-5, line_half_width - 1e-5, 220)
    horizontal_y = np.full_like(sample_x, -1.0)
    vertical = griddata(
        _cell_centres(case, time_name),
        field[:, 1],
        (sample_x, horizontal_y),
        method="linear",
    )
    if np.any(~np.isfinite(vertical)):
        fallback = griddata(
            _cell_centres(case, time_name),
            field[:, 1],
            (sample_x, horizontal_y),
            method="nearest",
        )
        vertical = np.where(np.isfinite(vertical), vertical, fallback)

    with (data_dir / "u_centerline.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["y", "u"])
        writer.writerows((float(y), float(u)) for y, u in zip(centreline_y, horizontal))
    with (data_dir / "v_horizontal.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["x", "v"])
        writer.writerows((float(x), float(v)) for x, v in zip(sample_x, vertical))

    finite_psi = np.isfinite(psi)
    min_index = np.nanargmin(psi)
    main_y, main_x = np.unravel_index(min_index, psi.shape)
    local_min = finite_psi & (psi == minimum_filter(np.nan_to_num(psi, nan=np.inf), size=9))
    local_max = finite_psi & (psi == maximum_filter(np.nan_to_num(psi, nan=-np.inf), size=9))
    vortices = []
    for label, selection in (("minimum", local_min), ("maximum", local_max)):
        for j, i in zip(*np.nonzero(selection)):
            if abs(psi[j, i]) < 1e-5:
                continue
            vortices.append({"type": label, "x": float(xx[j, i]), "y": float(yy[j, i]), "psi": float(psi[j, i])})
    primary = REFERENCE_PRIMARY.get(re_number)
    summary = {
        "caseName": case_name,
        "solver": metadata.get("solver"),
        "reynolds": re_number,
        "viscosity": metadata.get("viscosity"),
        "resolution": metadata.get("resolution"),
        "cells": int(len(field)),
        "finalTime": time_name,
        "mainVortex": {
            "x": float(xx[main_y, main_x]),
            "y": float(yy[main_y, main_x]),
            "psi": float(psi[main_y, main_x]),
            "omega": float(omega[main_y, main_x]),
        },
        "vortices": vortices,
        "referencePrimary": (
            {"psi": primary[0], "omega": primary[1], "x": primary[2], "y": primary[3]}
            if primary else None
        ),
        "uMin": float(np.nanmin(field[:, 0])),
        "uMax": float(np.nanmax(field[:, 0])),
        "vMin": float(np.nanmin(field[:, 1])),
        "vMax": float(np.nanmax(field[:, 1])),
        "horizontalSampleRange": [-line_half_width + 1e-5, line_half_width - 1e-5],
        "dataFiles": ["u_centerline.csv", "v_horizontal.csv", "summary.json"],
        "figureFiles": ["field_streamlines.png", "centerline_comparison.png", "streamfunction_vortices.png"],
    }
    (data_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    speed = np.hypot(grid[0], grid[1])
    fig, axes = plt.subplots(1, 2, figsize=(13, 6), constrained_layout=True)
    axes[0].pcolormesh(xx, yy, speed, shading="auto", cmap="viridis")
    axes[0].set_title(f"|U|, Re={re_number}")
    axes[1].streamplot(xx[0], yy[:, 0], grid[0], grid[1], density=1.4, color=speed, cmap="plasma")
    axes[1].set_title("Streamlines")
    for axis in axes:
        axis.set_aspect("equal")
        axis.set_xlabel("x")
        axis.set_ylabel("y")
    fig.savefig(figure_dir / "field_streamlines.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    axes[0].plot(horizontal, centreline_y, label="OpenFOAM")
    axes[0].plot([u for _, u in REFERENCE_U[re_number]], [y for y, _ in REFERENCE_U[re_number]], "k--", label="Kohno-Bathe")
    axes[0].set_xlabel("u(0,y)")
    axes[0].set_ylabel("y")
    axes[1].plot(sample_x, vertical, label="OpenFOAM")
    axes[1].plot([x for x, _ in REFERENCE_V[re_number]], [v for _, v in REFERENCE_V[re_number]], "k--", label="Kohno-Bathe")
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("v(x,-1)")
    for axis in axes:
        axis.grid(True, alpha=0.3)
        axis.legend()
    fig.savefig(figure_dir / "centerline_comparison.png", dpi=180)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(7, 6), constrained_layout=True)
    levels = np.linspace(np.nanmin(psi), np.nanmax(psi), 31)
    axis.contourf(xx, yy, psi, levels=levels, cmap="RdBu_r")
    axis.contour(xx, yy, psi, levels=levels[::3], colors="k", linewidths=0.35)
    axis.plot([-math.sqrt(3), math.sqrt(3), 0, -math.sqrt(3)], [0, 0, -3, 0], "k-", linewidth=1.2)
    axis.plot(xx[main_y, main_x], yy[main_y, main_x], "ko", label="primary vortex")
    axis.set_aspect("equal")
    axis.set_title(f"Streamfunction and vortices, Re={re_number}")
    axis.set_xlabel("x")
    axis.set_ylabel("y")
    axis.legend()
    fig.savefig(figure_dir / "streamfunction_vortices.png", dpi=180)
    plt.close(fig)
    return data_dir


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", required=True, type=Path)
    args = parser.parse_args()
    print(f"postprocessed={postprocess(args.case.resolve())}")


if __name__ == "__main__":
    main()
