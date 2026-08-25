#!/usr/bin/env python3

"""Collect, analyse, and plot all resolutions of one case family."""

from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
MPLCONFIGDIR = PROJECT_ROOT / "build" / "matplotlib"
MPLCONFIGDIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR))

from paths import (
    solver_analysis_dir,
    solver_analysis_figure_dir,
    solver_cases_dir,
    solver_data_dir,
    solver_figure_dir,
)


RAW_FIELDS = [
    "case",
    "caseName",
    "problem",
    "mesh",
    "resolution",
    "nominalH",
    "nCells",
    "finalTime",
    "finalTimeError",
    "normalizedL1",
    "normalizedL2",
    "normalizedLinf",
    "initialMass",
    "finalMass",
    "massChange",
    "normalizedMassError",
    "maxAbsResidualIntegral",
    "minCo",
    "maxCo",
    "targetCo",
    "initialAmplitude",
    "finalAmplitude",
    "minFinal",
    "maxFinal",
    "maxAbsFinal",
    "timeSteps",
    "meshOK",
    "boundedBelowInitial",
    "boundedAboveInitial",
    "solverEnded",
    "solverFatal",
]


def analysis_dirs(solver_family: str, case_name: str) -> tuple[Path, Path]:
    return (
        solver_analysis_dir(solver_family, case_name),
        solver_analysis_figure_dir(solver_family, case_name),
    )


def collect(solver_family: str, case_name: str, resolutions: list[int]) -> Path:
    data_dir, _ = analysis_dirs(solver_family, case_name)
    data_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    case_root = solver_cases_dir(solver_family) / case_name

    for resolution in sorted(resolutions):
        summary_path = (
            solver_data_dir(solver_family, case_name, resolution)
            / "summary.json"
        )
        if not summary_path.exists():
            raise RuntimeError(f"Missing per-case summary: {summary_path}")
        rows.append(json.loads(summary_path.read_text(encoding="utf-8")))

    raw_path = data_dir / "raw_results.csv"
    with raw_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=RAW_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in RAW_FIELDS})

    manifest = {
        "solverFamily": solver_family,
        "caseName": case_name,
        "resolutions": sorted(resolutions),
        "caseRoot": str(case_root),
        "caseDirectories": [str(case_root / f"N{n}") for n in sorted(resolutions)],
        "dataDirectories": [
            str(solver_data_dir(solver_family, case_name, n))
            for n in sorted(resolutions)
        ],
        "figureDirectories": [
            str(solver_figure_dir(solver_family, case_name, n))
            for n in sorted(resolutions)
        ],
        "rawResults": str(raw_path),
    }
    (data_dir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(f"rawResults={raw_path}")
    return raw_path


def _float(row: dict[str, str], name: str) -> float:
    return float(row[name])


def observed_order(previous_error: float, current_error: float, ratio: float) -> float:
    if previous_error <= 0.0 or current_error <= 0.0 or ratio <= 1.0:
        return float("nan")
    return math.log(previous_error / current_error) / math.log(ratio)


def analyse(solver_family: str, case_name: str) -> Path:
    data_dir, _ = analysis_dirs(solver_family, case_name)
    raw_path = data_dir / "raw_results.csv"
    with raw_path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    rows.sort(key=lambda row: int(row["resolution"]))
    if not rows:
        raise RuntimeError(f"No rows found in {raw_path}")

    summary_rows: list[dict[str, object]] = []
    previous: dict[str, str] | None = None
    for row in rows:
        resolution = int(row["resolution"])
        order_l1 = float("nan")
        order_l2 = float("nan")
        order_linf = float("nan")
        if previous is not None:
            ratio = resolution / int(previous["resolution"])
            order_l1 = observed_order(_float(previous, "normalizedL1"), _float(row, "normalizedL1"), ratio)
            order_l2 = observed_order(_float(previous, "normalizedL2"), _float(row, "normalizedL2"), ratio)
            order_linf = observed_order(_float(previous, "normalizedLinf"), _float(row, "normalizedLinf"), ratio)
        summary_rows.append(
            {
            "case": row.get("case", ""),
            "solverFamily": solver_family,
                "mesh": row.get("mesh", "quad"),
                "resolution": resolution,
                "nominalH": _float(row, "nominalH"),
                "nCells": row.get("nCells", ""),
                "normalizedL1": _float(row, "normalizedL1"),
                "normalizedL2": _float(row, "normalizedL2"),
                "normalizedLinf": _float(row, "normalizedLinf"),
                "observedOrderL1": order_l1,
                "observedOrderL2": order_l2,
                "observedOrderLinf": order_linf,
                "maxCo": _float(row, "maxCo"),
                "normalizedMassError": _float(row, "normalizedMassError"),
                "finalAmplitude": _float(row, "finalAmplitude"),
                "meshOK": row.get("meshOK", ""),
                "solverEnded": row.get("solverEnded", ""),
                "solverFatal": row.get("solverFatal", ""),
            }
        )
        previous = row

    mesh_type = str(summary_rows[0].get("mesh", "quad"))
    mesh_label = {
        "quad": "均匀结构化四边形",
        "tri": "三角形棱柱",
    }.get(mesh_type, mesh_type)
    summary_path = data_dir / "convergence_summary.csv"
    fields = list(summary_rows[0].keys())
    with summary_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in summary_rows:
            output = dict(row)
            for field in ("observedOrderL1", "observedOrderL2", "observedOrderLinf"):
                if isinstance(output[field], float) and math.isnan(output[field]):
                    output[field] = ""
            writer.writerow(output)

    valid_orders = [
        row["observedOrderL1"]
        for row in summary_rows
        if isinstance(row["observedOrderL1"], float)
        and not math.isnan(row["observedOrderL1"])
    ]
    first_error = float(summary_rows[0]["normalizedL1"])
    last_error = float(summary_rows[-1]["normalizedL1"])
    if valid_orders:
        order_text = ", ".join(f"{float(value):.6f}" for value in valid_orders)
        order_conclusion = f"L1 观察收敛阶依次为 `{order_text}`。"
    else:
        order_conclusion = "当前只有一个有效网格，暂时无法计算观察收敛阶。"

    report_lines = [
        f"# {case_name} 网格收敛性分析",
        "",
        f"本研究使用周期正弦波线性对流算例，比较不同 {mesh_label} 网格分辨率。",
        "所有网格使用相同的初始函数、速度、周期边界、CFL 和终止时间。",
        "",
        "$$p = \\frac{\\log(E_N/E_{2N})}{\\log(2)}$$",
        "",
        "其中 $E_N$ 是分辨率为 $N$ 时的误差，$p$ 是观察收敛阶。",
        "",
        "## 汇总表",
        "",
        "| N | cells | L1 | L2 | Linf | L1 order | maxCo | final amplitude |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in summary_rows:
        order = row["observedOrderL1"]
        order_text = "-" if isinstance(order, float) and math.isnan(order) else f"{float(order):.6f}"
        report_lines.append(
            f"| {row['resolution']} | {row['nCells']} | "
            f"{float(row['normalizedL1']):.8e} | {float(row['normalizedL2']):.8e} | "
            f"{float(row['normalizedLinf']):.8e} | {order_text} | "
            f"{float(row['maxCo']):.6f} | {float(row['finalAmplitude']):.8e} |"
        )
    report_lines.extend(
        [
            "",
            "## 结论",
            "",
            f"L1 误差从 `{first_error:.8e}` 变为 `{last_error:.8e}`。",
            order_conclusion,
            f"每个 N 的详细场数据、时间历史、日志和单案例图保存在 `{solver_family}/cases/<case>/Nxx`、`data/{solver_family}/cases/<case>/Nxx` 和 `figures/{solver_family}/cases/<case>/Nxx` 目录。",
            "",
        ]
    )
    report_path = data_dir / "analysis.md"
    report_path.write_text("\n".join(report_lines), encoding="utf-8")
    print(f"convergenceSummary={summary_path}")
    print(f"analysis={report_path}")
    return summary_path


def plot(solver_family: str, case_name: str) -> None:
    import matplotlib

    matplotlib.use("Agg")

    import matplotlib.pyplot as plt

    data_dir, figure_dir = analysis_dirs(solver_family, case_name)
    summary_path = data_dir / "convergence_summary.csv"
    figure_dir.mkdir(parents=True, exist_ok=True)
    with summary_path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    rows.sort(key=lambda row: int(row["resolution"]))
    if not rows:
        raise RuntimeError(f"No rows found in {summary_path}")

    n_values = [int(row["resolution"]) for row in rows]
    l1 = [float(row["normalizedL1"]) for row in rows]
    l2 = [float(row["normalizedL2"]) for row in rows]
    linf = [float(row["normalizedLinf"]) for row in rows]

    fig, axis = plt.subplots(figsize=(7.6, 5.2), constrained_layout=True)
    axis.loglog(n_values, l1, "o-", label="normalized L1")
    axis.loglog(n_values, l2, "s-", label="normalized L2")
    axis.loglog(n_values, linf, "^-", label="normalized Linf")
    reference = l1[-1] * n_values[-1]
    axis.loglog(n_values, [reference / n for n in n_values], "k--", label="slope 1 reference")
    axis.set_xlabel("N base intervals per direction")
    axis.set_ylabel("normalized error at t=1")
    axis.set_title(f"{case_name}: error convergence")
    axis.grid(True, which="both", alpha=0.25)
    axis.legend()
    fig.savefig(figure_dir / "convergence_errors.png", dpi=220)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(7.6, 5.2), constrained_layout=True)
    for field, label, marker in (
        ("observedOrderL1", "L1 order", "o"),
        ("observedOrderL2", "L2 order", "s"),
        ("observedOrderLinf", "Linf order", "^"),
    ):
        x = []
        y = []
        for row in rows:
            if row.get(field, ""):
                x.append(int(row["resolution"]))
                y.append(float(row[field]))
        if x:
            axis.plot(x, y, marker=marker, linestyle="-", label=label)
    axis.axhline(1.0, linestyle="--", color="black", label="first order")
    axis.set_xlabel("N base intervals per direction")
    axis.set_ylabel("observed order")
    axis.set_title(f"{case_name}: observed convergence order")
    axis.grid(True, alpha=0.25)
    axis.legend()
    fig.savefig(figure_dir / "convergence_order.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), constrained_layout=True)
    axes[0].plot(n_values, l1, "o-", label="normalized L1")
    axes[0].set_xlabel("N base intervals per direction")
    axes[0].set_ylabel("normalized L1 error")
    axes[0].grid(True, alpha=0.25)
    axes[0].legend()
    amplitudes = [float(row["finalAmplitude"]) for row in rows]
    axes[1].plot(n_values, amplitudes, "o-", label="final amplitude")
    axes[1].axhline(1.0, linestyle="--", color="black", label="exact amplitude")
    axes[1].set_xlabel("N base intervals per direction")
    axes[1].set_ylabel("amplitude")
    axes[1].grid(True, alpha=0.25)
    axes[1].legend()
    fig.suptitle(f"{case_name}: all-N comparison")
    fig.savefig(figure_dir / "all_N_comparison.png", dpi=220)
    plt.close(fig)

    print(f"analysisFigures={figure_dir}")
