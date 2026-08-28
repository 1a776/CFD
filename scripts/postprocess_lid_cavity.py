#!/usr/bin/env python3

"""Post-process one square lid-driven cavity case."""

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


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CASES_ROOT = PROJECT_ROOT / "cases" / "05_navier_stokes_equation"
DATA_ROOT = PROJECT_ROOT / "data" / "05_navier_stokes_equation" / "cases"
FIGURE_ROOT = PROJECT_ROOT / "figures" / "05_navier_stokes_equation" / "cases"


GHIA_U = {
    1000: [
        (1.0000, 1.00000), (0.9766, 0.65928), (0.9688, 0.57492),
        (0.9609, 0.51117), (0.9531, 0.46604), (0.8516, 0.33304),
        (0.7344, 0.18719), (0.6172, 0.05702), (0.5000, -0.06080),
        (0.4531, -0.10648), (0.2813, -0.27805), (0.1719, -0.38289),
        (0.1016, -0.29730), (0.0703, -0.22220), (0.0625, -0.20196),
        (0.0547, -0.18109), (0.0000, 0.00000),
    ],
    3200: [
        (1.0000, 1.00000), (0.9766, 0.53236), (0.9688, 0.48296),
        (0.9609, 0.46547), (0.9531, 0.46101), (0.8516, 0.34682),
        (0.7344, 0.19791), (0.6172, 0.07156), (0.5000, -0.04272),
        (0.4531, -0.086636), (0.2813, -0.24427), (0.1719, -0.34323),
        (0.1016, -0.41933), (0.0703, -0.37827), (0.0625, -0.35344),
        (0.0547, -0.32407), (0.0000, 0.00000),
    ],
    5000: [
        (1.0000, 1.00000), (0.9766, 0.48223), (0.9688, 0.46120),
        (0.9609, 0.45992), (0.9531, 0.46036), (0.8516, 0.33556),
        (0.7344, 0.20087), (0.6172, 0.08183), (0.5000, -0.03039),
        (0.4531, -0.07404), (0.2813, -0.22855), (0.1719, -0.33050),
        (0.1016, -0.40435), (0.0703, -0.43643), (0.0625, -0.42901),
        (0.0547, -0.41165), (0.0000, 0.00000),
    ],
}

GHIA_V = {
    1000: [
        (1.0000, 0.00000), (0.9688, -0.21388), (0.9609, -0.27669),
        (0.9531, -0.33714), (0.9453, -0.39188), (0.9063, -0.51500),
        (0.8594, -0.42665), (0.8047, -0.31966), (0.5000, 0.02526),
        (0.2344, 0.32235), (0.2266, 0.33075), (0.1563, 0.37095),
        (0.0938, 0.32627), (0.0781, 0.30353), (0.0703, 0.29012),
        (0.0625, 0.27485), (0.0000, 0.00000),
    ],
    3200: [
        (1.0000, 0.00000), (0.9688, -0.39017), (0.9609, -0.47425),
        (0.9531, -0.52357), (0.9453, -0.54053), (0.9063, -0.44307),
        (0.8594, -0.37401), (0.8047, -0.31184), (0.5000, 0.00999),
        (0.2344, 0.28188), (0.2266, 0.29030), (0.1563, 0.37119),
        (0.0938, 0.42768), (0.0781, 0.41906), (0.0703, 0.40917),
        (0.0625, 0.39560), (0.0000, 0.00000),
    ],
    5000: [
        (1.0000, 0.00000), (0.9688, -0.49774), (0.9609, -0.55069),
        (0.9531, -0.55408), (0.9453, -0.52876), (0.9063, -0.41442),
        (0.8594, -0.36214), (0.8047, -0.30018), (0.5000, 0.00945),
        (0.2344, 0.27280), (0.2266, 0.28066), (0.1563, 0.35368),
        (0.0938, 0.42951), (0.0781, 0.43648), (0.0703, 0.43329),
        (0.0625, 0.42447), (0.0000, 0.00000),
    ],
}


def _metadata(case: Path) -> dict:
    return json.loads((case / "metadata.json").read_text(encoding="utf-8"))


def _final_time(case: Path) -> str:
    times = []
    for path in case.iterdir():
        if not path.is_dir() or path.name in {"0", "0.orig", "constant", "system"}:
            continue
        try:
            value = float(path.name)
        except ValueError:
            continue
        if (path / "U").exists():
            times.append((value, path.name))
    if not times:
        raise RuntimeError(f"No numeric time directory with U found in {case}")
    return max(times)[1]


def _read_internal_field(path: Path, pattern: str) -> np.ndarray:
    text = path.read_text(encoding="utf-8")
    marker = re.search(pattern, text)
    if not marker:
        raise RuntimeError(f"Cannot read internalField from {path}")
    count = int(marker.group(1))
    body = text[marker.end():].split("\n)\n;\n\nboundaryField", 1)[0]
    return count, body


def _read_vector_field(path: Path) -> np.ndarray:
    count, body = _read_internal_field(
        path, r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\("
    )
    values = [
        (float(x), float(y), float(z))
        for x, y, z in re.findall(
            r"\(([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\)", body
        )
    ]
    if len(values) != count:
        raise RuntimeError(f"Expected {count} vectors in {path}, found {len(values)}")
    return np.asarray(values, dtype=float)


def _read_scalar_field(path: Path) -> np.ndarray:
    count, body = _read_internal_field(
        path, r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\("
    )
    values = [float(value) for value in re.findall(
        r"([-+0-9.eE]+)", body
    )]
    if len(values) != count:
        raise RuntimeError(f"Expected {count} scalars in {path}, found {len(values)}")
    return np.asarray(values, dtype=float)


def _load_case(case: Path) -> tuple[dict, str, np.ndarray, np.ndarray]:
    metadata = _metadata(case)
    time_name = _final_time(case)
    values = _read_vector_field(case / time_name / "U")
    centres_path = case / time_name / "C"
    if not centres_path.exists():
        centres_path = case / "constant" / "C"
    if not centres_path.exists():
        raise RuntimeError(f"Missing cell centres in {case / time_name / 'C'} and {centres_path}")
    centres = _read_vector_field(centres_path)
    if len(values) != len(centres):
        raise RuntimeError("U and C contain different numbers of cells")
    return metadata, time_name, values, centres[:, :2]


def _sample_centerlines(
    values: np.ndarray, centres: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    y = np.linspace(0.0, 1.0, 201)
    x = np.linspace(0.0, 1.0, 201)
    u_vertical = griddata(centres, values[:, 0], (np.full_like(y, 0.5), y), method="linear")
    v_horizontal = griddata(centres, values[:, 1], (x, np.full_like(x, 0.5)), method="linear")
    if np.any(~np.isfinite(u_vertical)):
        nearest = griddata(centres, values[:, 0], (np.full_like(y, 0.5), y), method="nearest")
        u_vertical = np.where(np.isfinite(u_vertical), u_vertical, nearest)
    if np.any(~np.isfinite(v_horizontal)):
        nearest = griddata(centres, values[:, 1], (x, np.full_like(x, 0.5)), method="nearest")
        v_horizontal = np.where(np.isfinite(v_horizontal), v_horizontal, nearest)

    u_vertical[0] = 0.0
    u_vertical[-1] = 1.0
    v_horizontal[0] = 0.0
    v_horizontal[-1] = 0.0
    return x, y, u_vertical, v_horizontal


def _write_csv(path: Path, header: list[str], rows: list[tuple[float, ...]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)


def _plot_field(
    values: np.ndarray, centres: np.ndarray, out: Path, title: str
) -> None:
    x = np.linspace(0.0, 1.0, 160)
    y = np.linspace(0.0, 1.0, 160)
    xx, yy = np.meshgrid(x, y)
    speed = np.hypot(values[:, 0], values[:, 1])
    speed_grid = griddata(centres, speed, (xx, yy), method="linear")
    ux_grid = griddata(centres, values[:, 0], (xx, yy), method="linear")
    uy_grid = griddata(centres, values[:, 1], (xx, yy), method="linear")
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    image = axes[0].pcolormesh(xx, yy, speed_grid, shading="auto", cmap="viridis")
    axes[0].set_title("|U|")
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("y")
    fig.colorbar(image, ax=axes[0])
    axes[1].streamplot(
        x, y, ux_grid, uy_grid, density=1.5, color=speed_grid, cmap="plasma"
    )
    axes[1].set_title("Streamlines")
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("y")
    fig.suptitle(title)
    fig.savefig(out, dpi=180)
    plt.close(fig)


def _plot_centerlines(
    x: np.ndarray,
    y: np.ndarray,
    u_vertical: np.ndarray,
    v_horizontal: np.ndarray,
    reynolds: int,
    out: Path,
) -> dict[str, float]:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    reference_u = np.asarray(GHIA_U.get(reynolds, []), dtype=float)
    reference_v = np.asarray(GHIA_V.get(reynolds, []), dtype=float)
    axes[0].plot(u_vertical, y, "o-", label="OpenFOAM")
    if len(reference_u):
        axes[0].plot(reference_u[:, 1], reference_u[:, 0], "k--", label="Ghia")
    axes[0].set_xlabel("u(0.5,y)")
    axes[0].set_ylabel("y")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    axes[1].plot(x, v_horizontal, "o-", label="OpenFOAM")
    if len(reference_v):
        axes[1].plot(reference_v[:, 0], reference_v[:, 1], "k--", label="Ghia")
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("v(x,0.5)")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    fig.suptitle(f"Re={reynolds} centerline comparison")
    fig.savefig(out, dpi=180)
    plt.close(fig)

    summary: dict[str, float] = {}
    if len(reference_u):
        numeric = np.interp(reference_u[:, 0], y, u_vertical)
        summary["u_max_abs_error_vs_ghia"] = float(np.max(np.abs(numeric - reference_u[:, 1])))
        summary["u_rmse_vs_ghia"] = float(np.sqrt(np.mean((numeric - reference_u[:, 1]) ** 2)))
    if len(reference_v):
        numeric = np.interp(reference_v[:, 0], x, v_horizontal)
        summary["v_max_abs_error_vs_ghia"] = float(np.max(np.abs(numeric - reference_v[:, 1])))
        summary["v_rmse_vs_ghia"] = float(np.sqrt(np.mean((numeric - reference_v[:, 1]) ** 2)))
    return summary


def postprocess(case: Path) -> Path:
    metadata, time_name, field, centres = _load_case(case)
    case_name = str(metadata.get("caseName", case.name))
    viscosity = float(metadata.get("viscosity", metadata.get("mu", 1.0)))
    reynolds = int(metadata.get("reynolds", round(1.0 / viscosity)))
    data_dir = DATA_ROOT / case_name
    figure_dir = FIGURE_ROOT / case_name
    data_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    x, y, u_vertical, v_horizontal = _sample_centerlines(field, centres)
    _write_csv(
        data_dir / "u_centerline.csv",
        ["y", "u"],
        [(float(yi), float(ui)) for yi, ui in zip(y, u_vertical)],
    )
    _write_csv(
        data_dir / "v_centerline.csv",
        ["x", "v"],
        [(float(xi), float(vi)) for xi, vi in zip(x, v_horizontal)],
    )
    _plot_field(
        field,
        centres,
        figure_dir / "field_and_streamlines.png",
        f"{case_name}, final time={time_name}",
    )
    comparison = _plot_centerlines(
        x,
        y,
        u_vertical,
        v_horizontal,
        reynolds,
        figure_dir / "centerline_comparison.png",
    )
    summary = {
        "caseName": case_name,
        "reynolds": reynolds,
        "viscosity": viscosity,
        "meshType": metadata.get("meshType", metadata.get("mesh", {}).get("type", "quad")),
        "cellsPerEdge": metadata.get("cellsPerEdge", metadata.get("resolution")),
        "cells": len(field),
        "finalTime": time_name,
        "uMin": float(np.min(field[:, 0])),
        "uMax": float(np.max(field[:, 0])),
        "vMin": float(np.min(field[:, 1])),
        "vMax": float(np.max(field[:, 1])),
        "centerline": comparison,
        "dataFiles": ["u_centerline.csv", "v_centerline.csv", "summary.json"],
        "figureFiles": ["field_and_streamlines.png", "centerline_comparison.png"],
    }
    (data_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    return data_dir


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", required=True, type=Path)
    args = parser.parse_args()
    output = postprocess(args.case.resolve())
    print(f"postprocessed={output}")


if __name__ == "__main__":
    main()
