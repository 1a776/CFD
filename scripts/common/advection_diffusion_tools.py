#!/usr/bin/env python3

"""Helpers for the periodic sine-wave advection-diffusion benchmark."""

from __future__ import annotations

import math
import re
from pathlib import Path


def exact_value(
    x: float,
    y: float,
    time_value: float,
    velocity: tuple[float, float, float],
    mu: float,
) -> float:
    """Return the exact solution of phi_t + div(U phi) - div(mu grad(phi)) = 0."""
    u, v, _ = velocity
    amplitude = math.exp(-8.0 * math.pi * math.pi * mu * time_value)
    phase = 2.0 * math.pi * (x + y - (u + v) * time_value)
    return amplitude * math.sin(phase)


def structured_values(
    nx: int,
    ny: int,
    domain: tuple[float, float, float, float],
    time_value: float,
    velocity: tuple[float, float, float],
    mu: float,
) -> list[float]:
    """Evaluate the exact solution at all structured cell centres."""
    xmin, xmax, ymin, ymax = domain
    dx = (xmax - xmin) / nx
    dy = (ymax - ymin) / ny
    return [
        exact_value(
            xmin + (i + 0.5) * dx,
            ymin + (j + 0.5) * dy,
            time_value,
            velocity,
            mu,
        )
        for j in range(ny)
        for i in range(nx)
    ]


def values_at_centres(
    centres: list[tuple[float, float, float]],
    time_value: float,
    velocity: tuple[float, float, float],
    mu: float,
) -> list[float]:
    """Evaluate the exact solution at arbitrary OpenFOAM cell centres."""
    return [
        exact_value(float(x), float(y), time_value, velocity, mu)
        for x, y, _ in centres
    ]


def rotating_peak_center(
    elapsed_time: float,
    initial_center: tuple[float, float] = (0.0, 0.5),
) -> tuple[float, float]:
    """Return the centre advected by U=(-y,x) after elapsed_time."""
    x0, y0 = initial_center
    return (
        x0 * math.cos(elapsed_time) - y0 * math.sin(elapsed_time),
        x0 * math.sin(elapsed_time) + y0 * math.cos(elapsed_time),
    )


def rotating_peak_exact_value(
    x: float,
    y: float,
    elapsed_time: float,
    epsilon: float,
    diffusion_start_time: float = math.pi / 2.0,
    initial_center: tuple[float, float] = (0.0, 0.5),
) -> float:
    """Exact rotating sharp-peak solution for the second advection-diffusion case."""
    physical_time = diffusion_start_time + elapsed_time
    xhat, yhat = rotating_peak_center(elapsed_time, initial_center)
    radius2 = (x - xhat) ** 2 + (y - yhat) ** 2
    return math.exp(-radius2 / (4.0 * epsilon * physical_time)) / (
        4.0 * math.pi * epsilon * physical_time
    )


def rotating_peak_values_at_centres(
    centres: list[tuple[float, float, float]],
    elapsed_time: float,
    epsilon: float,
    diffusion_start_time: float = math.pi / 2.0,
    initial_center: tuple[float, float] = (0.0, 0.5),
) -> list[float]:
    """Evaluate the rotating sharp-peak exact solution at OpenFOAM cell centres."""
    return [
        rotating_peak_exact_value(
            float(x),
            float(y),
            elapsed_time,
            epsilon,
            diffusion_start_time,
            initial_center,
        )
        for x, y, _ in centres
    ]


def rotating_peak_structured_values(
    nx: int,
    ny: int,
    domain: tuple[float, float, float, float],
    elapsed_time: float,
    epsilon: float,
    diffusion_start_time: float = math.pi / 2.0,
    initial_center: tuple[float, float] = (0.0, 0.5),
) -> list[float]:
    """Evaluate the rotating sharp-peak exact solution on structured cell centres."""
    xmin, xmax, ymin, ymax = domain
    dx = (xmax - xmin) / nx
    dy = (ymax - ymin) / ny
    return [
        rotating_peak_exact_value(
            xmin + (i + 0.5) * dx,
            ymin + (j + 0.5) * dy,
            elapsed_time,
            epsilon,
            diffusion_start_time,
            initial_center,
        )
        for j in range(ny)
        for i in range(nx)
    ]


def rotating_peak_boundary_field(
    patch_name: str,
    epsilon: float,
    diffusion_start_time: float,
    initial_center: tuple[float, float],
) -> str:
    """Return an analytic time-dependent Dirichlet boundary for the sharp peak."""
    safe_patch_name = "".join(
        character if character.isalnum() else "_"
        for character in patch_name
    )
    return f"""    {patch_name}
    {{
        // 原题边界条件：边界面上的 phi 也取旋转扩散尖峰解析解。
        // 数学公式：
        //   phi_b(x,y,tau) = 1/(4*pi*epsilon*(t0+tau))
        //                    * exp(-((x-xhat)^2+(y-yhat)^2)
        //                          /(4*epsilon*(t0+tau)))
        //   xhat = x0*cos(tau) - y0*sin(tau)
        //   yhat = x0*sin(tau) + y0*cos(tau)
        //
        // OpenFOAM 接口：
        //   type codedFixedValue 让边界值在每次 correctBoundaryConditions()
        //   时按当前 this->db().time().value() 重新计算。
        type            codedFixedValue;
        value           uniform 0;
        name            rotatingPeak_{safe_patch_name};
        code
        #{{
            const scalar epsilon = {epsilon:.16e};
            const scalar t0 = {diffusion_start_time:.16e};
            const scalar x0 = {initial_center[0]:.16e};
            const scalar y0 = {initial_center[1]:.16e};
            const scalar tau = this->db().time().value();
            const scalar physicalTime = t0 + tau;
            const scalar xhat = x0*cos(tau) - y0*sin(tau);
            const scalar yhat = x0*sin(tau) + y0*cos(tau);
            const vectorField& faceCentres = this->patch().Cf();
            scalarField values(faceCentres.size(), Zero);
            forAll(faceCentres, faceI)
            {{
                const scalar dx = faceCentres[faceI].x() - xhat;
                const scalar dy = faceCentres[faceI].y() - yhat;
                const scalar radius2 = dx*dx + dy*dy;
                values[faceI] =
                    exp(-radius2/(4.0*epsilon*physicalTime))
                   /(4.0*constant::mathematical::pi*epsilon*physicalTime);
            }}
            operator==(values);
        #}};
    }}"""


def rotation_velocity_values_at_centres(
    centres: list[tuple[float, float, float]],
) -> list[tuple[float, float, float]]:
    """Evaluate U=(-y,x,0) at OpenFOAM cell centres."""
    return [(-float(y), float(x), 0.0) for x, y, _ in centres]


def _field_header(object_name: str, class_name: str, dimensions: str) -> str:
    return f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /    O peration     |
    \\\\  /    A nd           |
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       {class_name};
    location    "0";
    object      {object_name};
}}

"""


def write_sine_initial_field(
    case: Path,
    nx: int,
    ny: int,
    domain: tuple[float, float, float, float],
    velocity: tuple[float, float, float],
    mu: float,
    field_name: str = "phi",
) -> Path:
    """Write the documented scalar initial condition into 0.orig."""
    values = structured_values(nx, ny, domain, 0.0, velocity, mu)
    body = "\n".join(f"    {value:.16e}" for value in values)
    output = case / "0.orig" / field_name
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        _field_header(field_name, "volScalarField", "[0 0 0 0 0 0 0]")
        + f"""// 第三题第一个算例的未知量：phi(x,y,t)。
//
// 对应控制方程：
//
//     d(phi)/dt + div(U phi) - div(mu grad(phi)) = 0
//
// 对应初始条件：
//
//     phi(x,y,0) = sin(2*pi*(x+y))
//
// internalField 中的 {len(values)} 个数是各有限体积单元中心的 phi(c,0)，
// 不是网格顶点值。scripts/prepare_case.py 会根据 system/blockMeshDict 的
// N x N 网格，把上式在每个单元中心逐个计算后写入这里。
dimensions      [0 0 0 0 0 0 0];

internalField   nonuniform List<scalar>
{len(values)}
(
{body}
);

boundaryField
{{
    // xMin/xMax 是一对 cyclic patch：
    // phi(0,y,t) = phi(1,y,t)。OpenFOAM 会配对两侧边界面，
    // fvc::div 和 fvc::laplacian 据此自动计算周期对流/扩散通量。
    xMin {{ type cyclic; }}
    xMax {{ type cyclic; }}

    // yMin/yMax 同样为周期边界：
    // phi(x,0,t) = phi(x,1,t)。
    yMin {{ type cyclic; }}
    yMax {{ type cyclic; }}

    // z 方向仅有一层单元。本二维问题不求解 z 方向变化。
    zMin {{ type empty; }}
    zMax {{ type empty; }}
}}
""",
        encoding="utf-8",
    )
    return output


def write_uniform_velocity_field(
    case: Path,
    velocity: tuple[float, float, float],
) -> Path:
    """Write the documented constant velocity field into 0.orig/U."""
    u, v, w = velocity
    output = case / "0.orig" / "U"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        _field_header("U", "volVectorField", "[0 1 -1 0 0 0 0]")
        + f"""// 对应第三题第一个算例的常速度：
//
//     U = ({u:g}, {v:g}, {w:g})
//
// 求解器用 fvc::flux(U) 计算每个面上的体积通量：
//
//     F_f = U_f dot S_f
//
// 然后用 fvc::div(faceFlux, phi) 构造对流散度项。
dimensions      [0 1 -1 0 0 0 0];

internalField   uniform ({u:.16g} {v:.16g} {w:.16g});

boundaryField
{{
    // 速度场必须与标量场使用相同的周期 patch 配对，
    // 才能使跨越 x/y 边界的面通量连续。
    xMin {{ type cyclic; }}
    xMax {{ type cyclic; }}
    yMin {{ type cyclic; }}
    yMax {{ type cyclic; }}

    // 二维 empty 面没有 z 向速度自由度。
    zMin {{ type empty; }}
    zMax {{ type empty; }}
}}
""",
        encoding="utf-8",
    )
    return output


def write_sine_initial_field_from_centres(
    case: Path,
    centres: list[tuple[float, float, float]],
    velocity: tuple[float, float, float],
    mu: float,
    field_name: str = "phi",
) -> Path:
    """Write the documented initial field for real triangular-cell centres."""
    values = values_at_centres(centres, 0.0, velocity, mu)
    body = "\n".join(f"    {value:.16e}" for value in values)
    output = case / "0.orig" / field_name
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        _field_header(field_name, "volScalarField", "[0 0 0 0 0 0 0]")
        + f"""// 第三题第一个三角形算例的未知量：phi(x,y,t)。
//
// 对应初始条件：
//
//     phi(x,y,0) = sin(2*pi*(x+y))
//
// 这是非结构化三角形网格。{len(values)} 个 internalField 数值按
// constant/C 中真实 OpenFOAM 单元中心逐个计算，不按数组 i,j 假定网格排序。
// Allrun 的顺序是：Gmsh 生成网格 -> gmshToFoam -> 写 cell centres ->
// scripts/prepare_case.py 读取 constant/C -> 写入本文件。
dimensions      [0 0 0 0 0 0 0];

internalField   nonuniform List<scalar>
{len(values)}
(
{body}
);

boundaryField
{{
    // xMin/xMax 和 yMin/yMax 是周期边界。cyclic patch 将两侧三角形
    // 边界面成对连接，使 fvc::div 和 fvc::laplacian 使用另一侧单元值。
    xMin {{ type cyclic; }}
    xMax {{ type cyclic; }}
    yMin {{ type cyclic; }}
    yMax {{ type cyclic; }}

    // zMin/zMax 是单层棱柱的二维 empty 面。
    zMin {{ type empty; }}
    zMax {{ type empty; }}
}}
""",
        encoding="utf-8",
    )
    return output


def write_rotating_peak_initial_field(
    case: Path,
    centres: list[tuple[float, float, float]],
    epsilon: float,
    field_name: str = "phi",
    diffusion_start_time: float = math.pi / 2.0,
    initial_center: tuple[float, float] = (0.0, 0.5),
    boundary_condition: str = "zeroDirichletApproximation",
) -> Path:
    """Write the second advection-diffusion case initial scalar field."""
    values = rotating_peak_values_at_centres(
        centres,
        0.0,
        epsilon,
        diffusion_start_time,
        initial_center,
    )
    body = "\n".join(f"    {value:.16e}" for value in values)
    output = case / "0.orig" / field_name
    output.parent.mkdir(parents=True, exist_ok=True)
    if boundary_condition == "analyticDirichlet":
        boundary_note = """// 原题说初始和边界条件由解析函数给出。
// 本案例使用 codedFixedValue 在 xMin/xMax/yMin/yMax 上按当前时间 tau
// 和边界面中心坐标计算解析 Dirichlet 值。"""
        boundary_entries = "\n".join(
            rotating_peak_boundary_field(
                patch_name,
                epsilon,
                diffusion_start_time,
                initial_center,
            )
            for patch_name in ("xMin", "xMax", "yMin", "yMax")
        )
    else:
        boundary_note = """// 原题说初始和边界条件由解析函数给出。这里先采用 fixedValue 0 的工程近似；
// 因为尖峰始终在区域内部旋转，边界解析值处于指数小量级，主要误差来自内部离散。"""
        boundary_entries = """    xMin { type fixedValue; value uniform 0; }
    xMax { type fixedValue; value uniform 0; }
    yMin { type fixedValue; value uniform 0; }
    yMax { type fixedValue; value uniform 0; }"""
    output.write_text(
        _field_header(field_name, "volScalarField", "[0 0 0 0 0 0 0]")
        + f"""// 第三题第二个算例的未知量：phi(x,y,tau)。
//
// OpenFOAM 计算时间记为 tau，范围是 0 <= tau <= 2*pi。
//
// 对应控制方程：
//
//     d(phi)/dt + div(U phi) - div(epsilon grad(phi)) = 0
//
// 对应旋转尖峰解析初值：
//
//     phi(x,y,0) = 1/(4*pi*epsilon*t0)
//                  * exp(-((x-x0)^2+(y-y0)^2)/(4*epsilon*t0))
//
// 其中：
//
//     epsilon = {epsilon:.16g}
//     t0 = {diffusion_start_time:.16g}
//     (x0,y0) = ({initial_center[0]:.16g},{initial_center[1]:.16g})
//
// internalField 的每个数都是单元中心值，不是网格顶点值。
// 对四边形网格，单元中心来自 blockMesh 的 N x N 规则网格；
// 对三角形网格，单元中心来自 gmshToFoam 后写出的 constant/C。
//
// 边界说明：
{boundary_note}
dimensions      [0 0 0 0 0 0 0];

internalField   nonuniform List<scalar>
{len(values)}
(
{body}
);

boundaryField
{{
{boundary_entries}
    zMin {{ type empty; }}
    zMax {{ type empty; }}
}}
""",
        encoding="utf-8",
    )
    return output


def write_rotating_velocity_field(
    case: Path,
    centres: list[tuple[float, float, float]],
) -> Path:
    """Write U=(-y,x,0) for the rotating sharp-peak case."""
    values = rotation_velocity_values_at_centres(centres)
    body = "\n".join(
        f"    ({u:.16e} {v:.16e} {w:.16e})" for u, v, w in values
    )
    output = case / "0.orig" / "U"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        _field_header("U", "volVectorField", "[0 1 -1 0 0 0 0]")
        + f"""// 第三题第二个算例的速度场：
//
//     U(x,y) = (-y, x, 0)
//
// 求解器用 fvc::flux(U) 计算面通量：
//
//     F_f = U_f dot S_f
//
// 然后用 fvc::div(faceFlux, phi) 构造对流项。
// internalField 是按单元中心坐标计算的非均匀速度。
dimensions      [0 1 -1 0 0 0 0];

internalField   nonuniform List<vector>
{len(values)}
(
{body}
);

boundaryField
{{
    // 速度边界用 zeroGradient 从内部外推到边界面；
    // 标量 phi 在边界上使用 fixedValue 0，所以边界标量通量主要由该 Dirichlet 值控制。
    xMin {{ type zeroGradient; }}
    xMax {{ type zeroGradient; }}
    yMin {{ type zeroGradient; }}
    yMax {{ type zeroGradient; }}
    zMin {{ type empty; }}
    zMax {{ type empty; }}
}}
""",
        encoding="utf-8",
    )
    return output


def parse_advection_diffusion_solver_log(path: Path) -> list[dict[str, float]]:
    """Extract completed time-step records from the third solver log."""
    if not path.exists():
        raise RuntimeError(f"Solver log not found: {path}")

    time_pattern = re.compile(
        r"Time = ([^\s]+)\s+step = (\d+)\s+deltaT = ([^\s]+)\s+"
        r"advectionDiffusionCo = ([^\s]+)"
    )
    rmin_pattern = re.compile(r"Rphi min\s*=\s*([^\s]+)")
    rmax_pattern = re.compile(r"Rphi max\s*=\s*([^\s]+)")
    phimin_pattern = re.compile(r"phi min\s*=\s*([^\s]+)")
    phimax_pattern = re.compile(r"phi max\s*=\s*([^\s]+)")
    mass_pattern = re.compile(r"mass\s*=\s*([^\s]+)")

    records: list[dict[str, float]] = []
    current: dict[str, float] | None = None
    for line in path.read_text(encoding="utf-8").splitlines():
        match = time_pattern.search(line)
        if match:
            current = {
                "time": float(match.group(1)),
                "step": float(match.group(2)),
                "deltaT": float(match.group(3)),
                "maxCo": float(match.group(4)),
            }
            continue
        if current is None:
            continue
        for pattern, name in (
            (rmin_pattern, "RphiMin"),
            (rmax_pattern, "RphiMax"),
            (phimin_pattern, "minT"),
            (phimax_pattern, "maxT"),
        ):
            match = pattern.search(line)
            if match:
                current[name] = float(match.group(1))
                if name == "maxT":
                    current["amplitude"] = 0.5 * (
                        current["maxT"] - current.get("minT", current["maxT"])
                    )
                break
        match = mass_pattern.search(line)
        if match:
            current["mass"] = float(match.group(1))
            records.append(current)
            current = None

    if not records:
        raise RuntimeError(f"No completed advection-diffusion records found in {path}")
    return records
