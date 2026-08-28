#!/usr/bin/env python3

"""Prepare and run OpenFOAM cases from JSON study configurations."""

from __future__ import annotations

import json
import os
import re
import shlex
import shutil
import subprocess
from pathlib import Path

from advection_rotation import (
    write_case_initial_field as write_rotation_initial_field,
    write_case_initial_field_from_centres as write_rotation_initial_field_from_centres,
    write_case_velocity_field as write_rotation_velocity_field,
    write_case_velocity_field_from_centres as write_rotation_velocity_field_from_centres,
)
from advection_sine import (
    write_case_initial_field as write_sine_initial_field,
    write_case_initial_field_from_centres,
)
from advection_diffusion_tools import (
    write_rotating_peak_initial_field as write_advection_diffusion_peak_initial_field,
    write_rotating_velocity_field as write_advection_diffusion_rotating_velocity_field,
    write_sine_initial_field as write_advection_diffusion_sine_initial_field,
    write_sine_initial_field_from_centres as write_advection_diffusion_sine_initial_field_from_centres,
    write_uniform_velocity_field as write_advection_diffusion_velocity_field,
)
from case_config import CaseConfig
from diffusion_tools import (
    write_discontinuous_initial_field,
    write_discontinuous_initial_field_from_centres,
    write_gaussian_initial_field,
    write_gaussian_initial_field_from_centres,
)
from poisson_tools import write_poisson_fields
from foam_fields import patch_uniform_vector_field, read_cell_geometry
from mesh_tools import mesh_resolution
from paths import PROJECT_ROOT


OUTER_PATCH_BOUNDARY = """boundary
(
    xMin
    {
        type patch;
        faces ((0 4 7 3));
    }
    xMax
    {
        type patch;
        faces ((1 2 6 5));
    }
    yMin
    {
        type patch;
        faces ((0 1 5 4));
    }
    yMax
    {
        type patch;
        faces ((3 7 6 2));
    }
    zMin
    {
        type empty;
        faces ((0 3 2 1));
    }
    zMax
    {
        type empty;
        faces ((4 5 6 7));
    }
);
"""

ADVECTION_DIFFUSION_PROBLEMS = {
    "sine_wave_advection_diffusion",
    "rotating_peak_advection_diffusion",
}

POISSON_PROBLEMS = {
    "poisson_manufactured",
}

LID_CAVITY_PROBLEMS = {
    "lid_driven_cavity",
    "tri_lid_driven_cavity",
}


def _structured_centres(
    nx: int,
    ny: int,
    domain: tuple[float, float, float, float],
    z: float,
) -> list[tuple[float, float, float]]:
    xmin, xmax, ymin, ymax = domain
    dx = (xmax - xmin) / nx
    dy = (ymax - ymin) / ny
    return [
        (xmin + (i + 0.5) * dx, ymin + (j + 0.5) * dy, z)
        for j in range(ny)
        for i in range(nx)
    ]


def _write_base_block_mesh(case: Path, config: CaseConfig, resolution: int) -> None:
    """Create the structured 2-D blockMesh dictionary without a template."""
    xmin, xmax, ymin, ymax = config.domain
    zmax = config.thickness
    boundary_type = "cyclic" if config.boundary_condition == "periodicXY" else "patch"
    neighbour = {
        "xMin": "neighbourPatch xMax;\n            transform translational;\n            separation (-1 0 0);",
        "xMax": "neighbourPatch xMin;\n            transform translational;\n            separation (1 0 0);",
        "yMin": "neighbourPatch yMax;\n            transform translational;\n            separation (0 -1 0);",
        "yMax": "neighbourPatch yMin;\n            transform translational;\n            separation (0 1 0);",
    }
    def patch(name: str, faces: str) -> str:
        if boundary_type == "cyclic":
            info = f"type cyclic;\n            {neighbour[name]}"
        else:
            info = "type patch;"
        return f"""    {name}
    {{
        {info}
        faces ({faces});
    }}"""

    path = case / "system" / "blockMeshDict"
    path.write_text(
        f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /    O  peration   |
    \\  /    A nd           |
     \\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       dictionary;
    location    "system";
    object      blockMeshDict;
}}

// 计算区域：Omega = [{xmin:g},{xmax:g}] x [{ymin:g},{ymax:g}]。
// 本文件由 JSON 直接生成，不依赖其它案例模板。
convertToMeters 1;

vertices
(
    ({xmin:g} {ymin:g} 0)
    ({xmax:g} {ymin:g} 0)
    ({xmax:g} {ymax:g} 0)
    ({xmin:g} {ymax:g} 0)
    ({xmin:g} {ymin:g} {zmax:g})
    ({xmax:g} {ymin:g} {zmax:g})
    ({xmax:g} {ymax:g} {zmax:g})
    ({xmin:g} {ymax:g} {zmax:g})
);

blocks
(
    // N={resolution}：每条边方向的单元数，二维单元数为 N*N。
    hex (0 1 2 3 4 5 6 7) ({resolution} {resolution} 1) simpleGrading (1 1 1)
);

edges
(
);

boundary
(
{patch("xMin", "(0 4 7 3)")}
{patch("xMax", "(1 2 6 5)")}
{patch("yMin", "(0 1 5 4)")}
{patch("yMax", "(3 7 6 2)")}
    zMin
    {{
        type empty;
        faces ((0 3 2 1));
    }}
    zMax
    {{
        type empty;
        faces ((4 5 6 7));
    }}
);
""",
        encoding="utf-8",
    )


def _write_base_system_files(case: Path, config: CaseConfig) -> None:
    """Create template-independent minimal system dictionaries."""
    system = case / "system"
    system.mkdir(parents=True, exist_ok=True)
    if config.problem in POISSON_PROBLEMS:
        return
    (system / "fvSchemes").write_text(
        """FoamFile
{
    format ascii;
    class dictionary;
    location "system";
    object fvSchemes;
}
ddtSchemes { default Euler; }
gradSchemes
{
    default Gauss linear;
}
divSchemes
{
    default none;
    div(phi,T) Gauss upwind;
    div(faceFlux,phi) Gauss upwind;
}
laplacianSchemes
{
    default none;
    laplacian(mu,phi) Gauss linear corrected;
}
interpolationSchemes
{
    default linear;
}
snGradSchemes
{
    default corrected;
}
""",
        encoding="utf-8",
    )
    (system / "controlDict").write_text(
        f"""FoamFile
{{
    format ascii;
    class dictionary;
    location "system";
    object controlDict;
}}
application {config.solver};
startFrom startTime;
startTime 0;
stopAt endTime;
endTime {format(config.end_time, ".17g")};
deltaT 1e-4;
writeControl runTime;
writeInterval 0.1;
purgeWrite 0;
writeFormat ascii;
writePrecision 12;
timePrecision 17;
writeCompression off;
runTimeModifiable no;
maxCo {config.max_co:g};
scalarField {config.scalar_field};
velocityField U;
""",
        encoding="utf-8",
    )
    (system / "fvSolution").write_text(
        """FoamFile
{
    format ascii;
    class dictionary;
    location "system";
    object fvSolution;
}
solvers
{
    // 显式求解器不组装线性方程组，此处无需配置线性求解器。
}
""",
        encoding="utf-8",
    )


def _write_base_velocity_field(case: Path, config: CaseConfig) -> None:
    """Create the basic advection velocity before patching its values."""
    if config.problem in LID_CAVITY_PROBLEMS:
        mesh_type = str(config.mesh_type).lower()
        empty_patches = (
            """    front
    {
        type empty;
    }
    back
    {
        type empty;
    }"""
            if mesh_type == "hybrid"
            else """    frontAndBack
    {
        type empty;
    }"""
        )
        if config.problem == "tri_lid_driven_cavity":
            wall_patches = """    leftWall
    {
        type noSlip;
    }
    rightWall
    {
        type noSlip;
    }
    movingTop
    {
        type fixedValue;
        value uniform (1 0 0);
    }"""
        else:
            wall_patches = """    leftWall
    {
        type noSlip;
    }
    rightWall
    {
        type noSlip;
    }
    bottomWall
    {
        type noSlip;
    }
    movingTop
    {
        type fixedValue;
        value uniform (1 0 0);
    }"""
        (case / "0.orig" / "U").write_text(
            f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /    O  peration     |
    \\\\  /    A nd           |
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       volVectorField;
    location    "0";
    object      U;
}}

// 第五题顶盖驱动腔流的速度初值。
//
//     U(x,y,0) = (0,0,0)
//
// movingTop 顶盖边界：
//
//     U = (1,0,0)
//
// 其余固壁边界：
//
//     U = (0,0,0)
dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField
{{
{wall_patches}
{empty_patches}
}}
""",
            encoding="utf-8",
        )
        return
    if config.problem not in {"sine_wave_advection", "solid_rotation_advection"}:
        return
    boundary = {
        name: ("cyclic" if config.boundary_condition == "periodicXY" else "zeroGradient")
        for name in ("xMin", "xMax", "yMin", "yMax")
    }
    boundary.update({"zMin": "empty", "zMax": "empty"})
    boundary_body = "\n".join(
        f"    {name}\n    {{\n        type {kind};\n    }}"
        for name, kind in boundary.items()
    )
    (case / "0.orig" / "U").write_text(
        f"""FoamFile
{{
    format ascii;
    class volVectorField;
    location "0";
    object U;
}}

dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField
{{
{boundary_body}
}}
""",
        encoding="utf-8",
    )


def _replace_or_append_dictionary_entry(path: Path, name: str, value: str) -> None:
    text = path.read_text(encoding="utf-8")
    pattern = re.compile(rf"^(\s*{re.escape(name)}\s+)([^;]+);", re.M)
    updated, count = pattern.subn(rf"\g<1>{value};", text, count=1)
    if count == 0:
        updated = text.rstrip() + f"\n{name:<16} {value};\n"
    path.write_text(updated, encoding="utf-8")


def _patch_fv_schemes(case: Path, config: CaseConfig) -> None:
    path = case / "system" / "fvSchemes"
    if config.problem in LID_CAVITY_PROBLEMS:
        is_projection = config.solver == "projectionFoamStudent"
        pressure_laplacian_entries = (
            f"    laplacian(dtCoeff,p) {config.laplacian_scheme};"
            if is_projection
            else (
                f"    laplacian((1|A(U)),p) {config.laplacian_scheme};\n"
                f"    laplacian(rAU,p) {config.laplacian_scheme};"
            )
        )
        path.write_text(
            f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /    O  peration   |
    \\\\  /    A nd           |
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSchemes;
}}

ddtSchemes
{{
    default Euler;
}}

gradSchemes
{{
    default Gauss linear;
}}

divSchemes
{{
    default none;
    div(phi,U) {config.div_scheme};
    div(phi,UStar) {config.div_scheme};
}}

laplacianSchemes
{{
    default none;
    laplacian(nu,U) {config.laplacian_scheme};
    laplacian(nu,UStar) {config.laplacian_scheme};
{pressure_laplacian_entries}
}}

interpolationSchemes
{{
    default linear;
}}

snGradSchemes
{{
    default corrected;
}}
""",
            encoding="utf-8",
        )
        return
    if config.problem in POISSON_PROBLEMS:
        path.write_text(
            f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /    O  peration   |
    \\\\  /    A nd           |
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSchemes;
}}

// 第四题 Poisson 方程是稳态问题，没有时间离散。
ddtSchemes
{{
    default         steadyState;
}}

gradSchemes
{{
    // corrected Laplacian 需要单元中心梯度。
    //     grad(phi)_c ≈ (1/V_c) sum_f phi_f S_cf
    default         Gauss linear;
}}

divSchemes
{{
    // 本题没有对流散度项。
    default         none;
}}

laplacianSchemes
{{
    // 对应求解器源码：
    //     fvm::laplacian(phi) == omega
    //
    // 对应有限体积形式：
    //     sum_f [ grad(phi)_f dot S_cf ] = V_c omega_c
    //
    // Gauss     -> 面通量求和；
    // linear    -> 面插值；
    // corrected -> 非正交修正。
    laplacian(phi)  {config.laplacian_scheme};
}}

interpolationSchemes
{{
    default         linear;
}}

snGradSchemes
{{
    // 控制边界面和非正交内部面的法向梯度。
    default         {config.sn_grad_scheme};
}}
""",
            encoding="utf-8",
        )
        return
    if config.problem in ADVECTION_DIFFUSION_PROBLEMS:
        case_label = (
            "第三题第一个周期正弦波案例"
            if config.problem == "sine_wave_advection_diffusion"
            else "第三题第二个旋转尖峰案例"
        )
        schemes_text = """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /    O peration     |
    \\\\  /    A nd           |
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSchemes;
}

ddtSchemes
{
    // 显式 Euler 更新写在求解器源码中：
    //     phi(c,n+1) = phi(c,n) + deltaT R(c,n)
    // 此处保留 Euler，作为该案例的时间离散记录。
    default         Euler;
}

gradSchemes
{
    // Gauss linear 的单元中心梯度会供 corrected 扩散修正使用。
    default         Gauss linear;
}

divSchemes
{
    default         none;

    // 对应源码：
    //     -fvc::div(faceFlux, phi, "div(faceFlux,phi)")
    //
    // 对应有限体积对流项：
    //     -1/V(c) sum_f [ F(cf) phi_f ]
    //     F(cf) = U_f dot S(cf)
    //
    // Gauss：用控制体表面通量和计算散度；
    // upwind：phi_f 取迎风侧单元值，是一阶有界格式。
    div(faceFlux,phi) Gauss upwind;
}

laplacianSchemes
{
    // 对应源码：
    //     fvc::laplacian(mu, phi)
    //
    // 对应有限体积扩散项：
    //     1/V(c) sum_f [ mu_f grad(phi)_f dot S(cf) ]
    //
    // linear：面上 mu/phi 的线性重构；
    // corrected：加入非正交网格的法向梯度修正。
    laplacian(mu,phi) Gauss linear corrected;
}

interpolationSchemes
{
    // 对流面值由 divSchemes 的 upwind 决定；
    // 该默认项服务于其它线性面插值。
    default         linear;
}

snGradSchemes
{
    // 与 laplacian(mu,phi) 的 corrected 选项配套，
    // 控制扩散通量中的法向梯度修正。
    default         corrected;
}
"""
        schemes_text = schemes_text.replace(
            "\nddtSchemes",
            "\n"
            f"// 当前文件对应：{case_label}。\n"
            "// 这里写“空间离散格式怎么选”，不写具体物理参数。\n"
            "// 具体入口关系：\n"
            "//     div(faceFlux,phi)  -> 对流项 div(U phi)，格式来自 JSON divScheme；\n"
            "//     laplacian(mu,phi) -> 扩散项 div(mu grad(phi))，格式来自 JSON laplacianScheme；\n"
            "//     snGradSchemes     -> 边界面和非正交面法向梯度修正，来自 JSON snGradScheme。\n\n"
            "ddtSchemes",
        )
        schemes_text = schemes_text.replace(
            "div(faceFlux,phi) Gauss upwind;",
            f"div(faceFlux,phi) {config.div_scheme};",
        )
        schemes_text = schemes_text.replace(
            "laplacian(mu,phi) Gauss linear corrected;",
            f"laplacian(mu,phi) {config.laplacian_scheme};",
        )
        schemes_text = schemes_text.replace(
            "default         corrected;\n}",
            f"default         {config.sn_grad_scheme};\n}}",
        )
        path.write_text(
            schemes_text,
            encoding="utf-8",
        )
        return
    text = path.read_text(encoding="utf-8")
    if config.problem.startswith("diffusion_"):
        updated, count = re.subn(
            r"^(\s*laplacian\(mu,phi\)\s+)([^;]+);",
            rf"\g<1>{config.laplacian_scheme};",
            text,
            count=1,
            flags=re.M,
        )
        if count != 1:
            raise RuntimeError(f"Cannot find laplacian(mu,phi) in {path}")
        updated, count = re.subn(
            r"(snGradSchemes\s*\{\s*.*?^\s*default\s+)([^;]+);",
            rf"\g<1>{config.sn_grad_scheme};",
            updated,
            count=1,
            flags=re.M | re.S,
        )
        if count != 1:
            raise RuntimeError(f"Cannot find snGradSchemes/default in {path}")
        path.write_text(updated, encoding="utf-8")
        return

    updated, count = re.subn(
        r"^(\s*div\(phi,T\)\s+)([^;]+);",
        rf"\g<1>{config.div_scheme};",
        text,
        count=1,
        flags=re.M,
    )
    if count != 1:
        raise RuntimeError(f"Cannot find div(phi,T) in {path}")

    if config.grad_t_scheme:
        grad_line = f"    grad(T)         {config.grad_t_scheme};"
        lines = updated.splitlines()
        grad_block_start: int | None = None
        grad_block_end: int | None = None
        existing_grad_line: int | None = None
        brace_depth = 0
        waiting_for_open_brace = False

        for line_number, line in enumerate(lines):
            if grad_block_start is None:
                if re.match(r"^\s*gradSchemes\s*\{\s*$", line):
                    grad_block_start = line_number
                    brace_depth = 1
                elif re.match(r"^\s*gradSchemes\s*$", line):
                    grad_block_start = line_number
                    waiting_for_open_brace = True
                continue

            if waiting_for_open_brace:
                if line.strip() == "{":
                    brace_depth = 1
                    waiting_for_open_brace = False
                continue

            if re.match(r"^\s*grad\(T\)\s+[^;]+;", line):
                existing_grad_line = line_number

            brace_depth += line.count("{") - line.count("}")
            if brace_depth == 0:
                grad_block_end = line_number
                break

        if grad_block_start is None or grad_block_end is None:
            raise RuntimeError(f"Cannot find complete gradSchemes block in {path}")

        if existing_grad_line is not None:
            lines[existing_grad_line] = grad_line
        else:
            # Insert immediately before the single closing brace of
            # gradSchemes. The previous regex kept that brace and added a
            # second one, which nested divSchemes incorrectly.
            lines.insert(grad_block_end, grad_line)

        updated = "\n".join(lines) + "\n"
    path.write_text(updated, encoding="utf-8")


def _patch_control_dict(case: Path, config: CaseConfig) -> None:
    path = case / "system" / "controlDict"
    if config.problem in POISSON_PROBLEMS:
        path.write_text(
            f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /    O  peration   |
    \\\\  /    A nd           |
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}}

// 第四题 Poisson 方程：
//     laplacian(phi) = omega
application     {config.solver};

// 这是稳态求解。endTime=0 不是物理终止时间，
// 而是保留 OpenFOAM 需要的时间目录管理接口。
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         0;
deltaT          1;

// 求解结束后写出 0/ 下的最终 phi。
writeControl    timeStep;
writeInterval   1;
purgeWrite      0;
writeFormat     ascii;
writePrecision  16;
timePrecision   12;
writeCompression off;
runTimeModifiable no;

// 求解器读取的字段名：
//     solutionField -> 0/phi，未知解；
//     sourceField   -> 0/omega，已知源项。
solutionField   {config.scalar_field};
sourceField     {config.source_field};

// corrected 非正交修正循环次数。
// 该值由 JSON 进入 cases，不是物理方程本身给定的参数。
nNonOrthogonalCorrectors {config.n_non_orthogonal_correctors};
""",
            encoding="utf-8",
        )
        return
    if config.problem in ADVECTION_DIFFUSION_PROBLEMS:
        case_time_description = (
            "第三题第一个案例在 t=1 比较数值解与解析解"
            if config.problem == "sine_wave_advection_diffusion"
            else "第三题第二个案例用 tau 从 0 旋转一圈到 2*pi，并与解析旋转尖峰比较"
        )
        path.write_text(
            f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /    O peration     |
    \\\\  /    A nd           |
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}}

// 应用程序入口。它必须与 build/03_advection_diffusion_equation/bin
// 中的可执行文件名一致。
application     {config.solver};

// 求解时间区间。{case_time_description}。
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {config.end_time:.17g};

// 初值 deltaT 是 OpenFOAM 标准字段；求解器每一步会根据
//
//     deltaT = alpha / max_c[(0.5 sum_f |F(cf)| + sum_f D(cf))/V(c)]
//
// 自动改写实际时间步。
deltaT          1e-4;

// 每 0.1 个物理时间写一次，便于观察衰减过程，同时避免显式细时间步
// 产生过多时间目录。求解器循环结束后还会强制写出 t=endTime。
writeControl    runTime;
writeInterval   0.1;
purgeWrite      0;

writeFormat     ascii;
writePrecision  12;
timePrecision   17;
writeCompression off;
runTimeModifiable no;

// 这两个字段决定求解器读取哪两个 0/ 文件：
// velocityField U   -> 0/U
// scalarField phi   -> 0/phi
velocityField   U;
scalarField     {config.scalar_field};

// alpha 是上式中的显式稳定安全系数。
// 它不是纯对流 CFL，而是对流和扩散共同构成的稳定监测量。
advectionDiffusionCo {config.advection_diffusion_co:.16g};

// 粗网格时额外限制最大时间步；细网格时通常由扩散稳定条件控制。
maxDeltaT       {config.max_delta_t:.16g};
""",
            encoding="utf-8",
        )
        return
    _replace_or_append_dictionary_entry(path, "application", config.solver)
    # Preserve the configured terminal time, especially values such as 2*pi.
    # The previous ":g" formatting rounded 2*pi to 6.28319 before OpenFOAM
    # read it, which shifted the requested final time.
    _replace_or_append_dictionary_entry(
        path,
        "endTime",
        format(config.end_time, ".17g"),
    )
    _replace_or_append_dictionary_entry(path, "maxCo", f"{config.max_co:g}")
    if config.problem in LID_CAVITY_PROBLEMS:
        _replace_or_append_dictionary_entry(path, "deltaT", f"{config.delta_t:.16g}")
        _replace_or_append_dictionary_entry(
            path,
            "adjustTimeStep",
            "true" if config.adjust_time_step else "false",
        )
        _replace_or_append_dictionary_entry(
            path,
            "writeInterval",
            f"{config.write_interval:.16g}",
        )
        if config.max_delta_t is not None:
            _replace_or_append_dictionary_entry(
                path,
                "maxDeltaT",
                f"{config.max_delta_t:g}",
            )
        _replace_or_append_dictionary_entry(
            path,
            "steadyStateControl",
            "true" if config.steady_state_control else "false",
        )
        _replace_or_append_dictionary_entry(
            path,
            "steadyVelocityTol",
            f"{config.steady_velocity_tol:.16g}",
        )
        _replace_or_append_dictionary_entry(
            path,
            "steadyMassTol",
            f"{config.steady_mass_tol:.16g}",
        )
        _replace_or_append_dictionary_entry(
            path,
            "minimumSteadySteps",
            str(config.minimum_steady_steps),
        )
        _replace_or_append_dictionary_entry(
            path,
            "requiredSteadySteps",
            str(config.required_steady_steps),
        )
    _replace_or_append_dictionary_entry(path, "scalarField", config.scalar_field)
    if config.problem.startswith("diffusion_"):
        _replace_or_append_dictionary_entry(path, "diffusionCo", f"{config.diffusion_co:g}")
        if config.max_delta_t is not None:
            _replace_or_append_dictionary_entry(path, "maxDeltaT", f"{config.max_delta_t:g}")
    # Keep the final time directory name consistent with the configured
    # terminal time, including values such as 2*pi.
    _replace_or_append_dictionary_entry(path, "timePrecision", "20")


def _write_transport_properties(case: Path, config: CaseConfig) -> Path:
    """Write constant/transportProperties for the diffusion solver."""
    path = case / "constant" / "transportProperties"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /    O peration     |
    \\\\  /    A nd           |
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       dictionary;
    location    "constant";
    object      transportProperties;
}}

// 对应求解器里的 dimensionedScalar mu。
// 数学上是扩散系数 μ，量纲是 L2/T。
mu              [0 2 -1 0 0 0 0] {config.diffusivity:.16g};
""",
        encoding="utf-8",
    )
    return path


def _patch_block_mesh_outer_boundaries(case: Path) -> None:
    """Use non-cyclic outer patches for compact support rotation profiles."""
    path = case / "system" / "blockMeshDict"
    text = path.read_text(encoding="utf-8")
    updated, count = re.subn(r"boundary\s*\(.*\)\s*;", OUTER_PATCH_BOUNDARY, text, count=1, flags=re.S)
    if count != 1:
        raise RuntimeError(f"Cannot replace boundary block in {path}")
    path.write_text(updated, encoding="utf-8")


def _patch_velocity_field(case: Path, config: CaseConfig) -> Path:
    """Apply the configured velocity to the generated case field."""
    if config.problem in LID_CAVITY_PROBLEMS:
        return case / "0.orig" / "U"
    if config.problem in POISSON_PROBLEMS:
        # Poisson 方程没有速度场 U；这里只是让统一的 case
        # preparation 流程跳过前三题的速度字段处理。
        return case / "0.orig" / config.scalar_field
    if config.problem.startswith("diffusion_"):
        return case / "0.orig" / config.scalar_field
    if config.problem == "sine_wave_advection_diffusion":
        return write_advection_diffusion_velocity_field(case, config.velocity)
    if config.problem == "rotating_peak_advection_diffusion":
        if config.mesh_type == "tri":
            if not (case / "constant" / "C").exists():
                return case / "0.orig" / "U"
            centres, _ = read_cell_geometry(case)
        else:
            nx, ny = mesh_resolution(case)
            centres = _structured_centres(nx, ny, config.domain, 0.5 * config.thickness)
        return write_advection_diffusion_rotating_velocity_field(case, centres)
    if config.problem == "solid_rotation_advection":
        if config.mesh_type == "tri":
            # The first preparation pass happens before Gmsh creates
            # constant/C.  The real cell-centred velocity is written later
            # by Allrun's --refresh-initial-only step.
            if not (case / "constant" / "C").exists():
                return case / "0.orig" / "U"
            centres, _ = read_cell_geometry(case)
            return write_rotation_velocity_field_from_centres(
                case,
                centres,
                config.velocity_model,
            )
        if config.mesh_type != "quad":
            raise NotImplementedError("Solid rotation velocity is currently implemented for quad mesh")
        nx, ny = mesh_resolution(case)
        return write_rotation_velocity_field(
            case,
            nx,
            ny,
            config.domain,
            config.velocity_model,
        )
    path = case / "0.orig" / "U"
    return patch_uniform_vector_field(path, config.velocity)


def _write_create_patch_dict(case: Path, config: CaseConfig) -> Path:
    """Create the final Gmsh patches for periodic or outer-boundary cases."""
    path = case / "system" / "createPatchDict"
    if config.boundary_condition == "periodicXY":
        horizontal_patches = """
    xMin
    {
        patchInfo
        {
            type cyclic;
            neighbourPatch xMax;
            transformType translational;
            separation (-1 0 0);
        }
        constructFrom patches;
        patches (xMinSource);
    }

    xMax
    {
        patchInfo
        {
            type cyclic;
            neighbourPatch xMin;
            transformType translational;
            separation (1 0 0);
        }
        constructFrom patches;
        patches (xMaxSource);
    }

    yMin
    {
        patchInfo
        {
            type cyclic;
            neighbourPatch yMax;
            transformType translational;
            separation (0 -1 0);
        }
        constructFrom patches;
        patches (yMinSource);
    }

    yMax
    {
        patchInfo
        {
            type cyclic;
            neighbourPatch yMin;
            transformType translational;
            separation (0 1 0);
        }
        constructFrom patches;
        patches (yMaxSource);
    }
"""
    else:
        horizontal_patches = """
    xMin
    {
        patchInfo { type patch; }
        constructFrom patches;
        patches (xMinSource);
    }

    xMax
    {
        patchInfo { type patch; }
        constructFrom patches;
        patches (xMaxSource);
    }

    yMin
    {
        patchInfo { type patch; }
        constructFrom patches;
        patches (yMinSource);
    }

    yMax
    {
        patchInfo { type patch; }
        constructFrom patches;
        patches (yMaxSource);
    }
"""
    path.write_text(
        f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /    O peration     | Mesh patch construction
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       dictionary;
    location    "system";
    object      createPatchDict;
}}

pointSync false;
writeCyclicMatch false;

patches
{{
    zMin
    {{
        patchInfo {{ type empty; }}
        constructFrom patches;
        patches (zMinSource);
    }}

    zMax
    {{
        patchInfo {{ type empty; }}
        constructFrom patches;
        patches (zMaxSource);
    }}
{horizontal_patches}
}}
""",
        encoding="utf-8",
    )
    return path


def _write_initial_field(case: Path, config: CaseConfig) -> Path:
    if config.problem in LID_CAVITY_PROBLEMS:
        mesh_type = str(config.mesh_type).lower()
        empty_patches = (
            """    front
    {
        type empty;
    }
    back
    {
        type empty;
    }"""
            if mesh_type == "hybrid"
            else """    frontAndBack
    {
        type empty;
    }"""
        )
        if config.problem == "tri_lid_driven_cavity":
            wall_patches = """    leftWall
    {
        type noSlip;
    }
    rightWall
    {
        type noSlip;
    }
    movingTop
    {
        type fixedValue;
        value uniform (1 0 0);
    }"""
            pressure_patches = """    leftWall
    {
        type zeroGradient;
    }
    rightWall
    {
        type zeroGradient;
    }
    movingTop
    {
        type zeroGradient;
    }"""
        else:
            wall_patches = """    leftWall
    {
        type noSlip;
    }
    rightWall
    {
        type noSlip;
    }
    bottomWall
    {
        type noSlip;
    }
    movingTop
    {
        type fixedValue;
        value uniform (1 0 0);
    }"""
            pressure_patches = """    leftWall
    {
        type zeroGradient;
    }
    rightWall
    {
        type zeroGradient;
    }
    bottomWall
    {
        type zeroGradient;
    }
    movingTop
    {
        type zeroGradient;
    }"""
        path = case / "0.orig" / "p"
        path.write_text(
            f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /    O  peration     |
    \\\\  /    A nd           |
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       volScalarField;
    location    "0";
    object      p;
}}

// 顶盖驱动腔流使用运动学压力 p = p_phys / rho。
//
// 初始压力：
//
//     p(x,y,0) = 0
//
// 固壁上采用零法向梯度作为压力边界初值。
dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{{
{pressure_patches}
{empty_patches}
}}
""",
            encoding="utf-8",
        )
        return path
    if config.problem in POISSON_PROBLEMS:
        if config.mesh_type == "tri":
            centres, _ = read_cell_geometry(case)
        else:
            nx, ny = mesh_resolution(case)
            centres = _structured_centres(
                nx,
                ny,
                config.domain,
                0.5 * config.thickness,
            )
        _, source_path = write_poisson_fields(
            case,
            centres,
            config.scalar_field,
            config.source_field,
        )
        return source_path
    if config.problem == "sine_wave_advection_diffusion":
        if config.mesh_type == "tri":
            centres, _ = read_cell_geometry(case)
            return write_advection_diffusion_sine_initial_field_from_centres(
                case,
                centres,
                config.velocity,
                config.diffusivity,
                config.scalar_field,
            )
        nx, ny = mesh_resolution(case)
        return write_advection_diffusion_sine_initial_field(
            case,
            nx,
            ny,
            config.domain,
            config.velocity,
            config.diffusivity,
            config.scalar_field,
        )
    if config.problem == "rotating_peak_advection_diffusion":
        initial_center_raw = config.initial_profile.get("center", [0.0, 0.5])
        if not isinstance(initial_center_raw, list | tuple) or len(initial_center_raw) != 2:
            raise RuntimeError("initialProfile.center must be [x0, y0]")
        initial_center = (float(initial_center_raw[0]), float(initial_center_raw[1]))
        diffusion_start_time = float(
            config.initial_profile.get("diffusionStartTime", 0.5 * 3.141592653589793)
        )
        if config.mesh_type == "tri":
            centres, _ = read_cell_geometry(case)
        else:
            nx, ny = mesh_resolution(case)
            centres = _structured_centres(nx, ny, config.domain, 0.5 * config.thickness)
        return write_advection_diffusion_peak_initial_field(
            case,
            centres,
            config.diffusivity,
            config.scalar_field,
            diffusion_start_time,
            initial_center,
            config.boundary_condition,
        )
    if config.problem == "diffusion_discontinuity":
        if config.mesh_type == "tri":
            centres, _ = read_cell_geometry(case)
            return write_discontinuous_initial_field_from_centres(
                case,
                centres,
                config.scalar_field,
            )
        nx, ny = mesh_resolution(case)
        return write_discontinuous_initial_field(
            case,
            nx,
            ny,
            config.domain,
            config.scalar_field,
        )
    if config.problem == "diffusion_gaussian":
        if config.mesh_type == "tri":
            centres, _ = read_cell_geometry(case)
            return write_gaussian_initial_field_from_centres(
                case,
                centres,
                config.scalar_field,
                config.diffusivity,
            )
        nx, ny = mesh_resolution(case)
        return write_gaussian_initial_field(
            case,
            nx,
            ny,
            config.domain,
            config.scalar_field,
            config.diffusivity,
        )
    if config.problem == "sine_wave_advection":
        if config.mesh_type == "quad":
            nx, ny = mesh_resolution(case)
            return write_sine_initial_field(case, nx, ny, config.velocity)
        if config.mesh_type == "tri":
            centres, _ = read_cell_geometry(case)
            return write_case_initial_field_from_centres(
                case,
                centres,
                config.velocity,
            )
    if config.problem == "solid_rotation_advection":
        if config.mesh_type == "tri":
            centres, _ = read_cell_geometry(case)
            return write_rotation_initial_field_from_centres(
                case,
                centres,
                config.initial_profile,
            )
        nx, ny = mesh_resolution(case)
        return write_rotation_initial_field(
            case,
            nx,
            ny,
            config.domain,
            config.initial_profile,
        )
    raise NotImplementedError(f"Unsupported problem for initial field: {config.problem}")


def _write_explicit_advection_diffusion_fv_solution(case: Path) -> None:
    """Document why the fully explicit solver has no linear-system controls."""
    (case / "system" / "fvSolution").write_text(
        """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /    O peration     |
    \\\\  /    A nd           |
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSolution;
}

solvers
{
    // 本求解器不调用 fvm::ddt、fvm::laplacian，也不组装 A*phi=b。
    //
    // 每一步直接计算：
    //     R = -div(U phi) + div(mu grad(phi))
    //     phi(new) = phi(old) + deltaT R
    //
    // 因而不存在需要在这里配置的线性方程组求解器或残差容限。
}
""",
        encoding="utf-8",
    )


def _write_lid_cavity_fv_solution(case: Path, config: CaseConfig) -> None:
    (case / "system" / "fvSolution").write_text(
        f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /    O  peration   |
    \\\\  /    A nd           |
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSolution;
}}

PISO
{{
    // 第五题 PISO 参数入口。
    nCorrectors 2;
    nNonOrthogonalCorrectors {config.n_non_orthogonal_correctors};
    pRefCell 0;
    pRefValue 0;
    momentumPredictor yes;
}}

solvers
{{
    U
    {{
        solver smoothSolver;
        smoother symGaussSeidel;
        tolerance {config.linear_tolerance:.16g};
        relTol 0;
    }}

    UStar
    {{
        $U;
    }}

    p
    {{
        solver GAMG;
        tolerance {config.linear_tolerance:.16g};
        relTol 0;
        smoother GaussSeidel;
    }}

    pFinal
    {{
        $p;
        relTol 0;
    }}
}}
""",
        encoding="utf-8",
    )


def _write_lid_cavity_physical_properties(case: Path, config: CaseConfig) -> None:
    (case / "constant" / "physicalProperties").write_text(
        f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /    O  peration   |
    \\\\  /    A nd           |
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       dictionary;
    location    "constant";
    object      physicalProperties;
}}

// 第五题顶盖驱动腔流的运动黏度。
//
// Reynolds 数定义为：
//
//     Re = U L / nu
//
// 这里取 U = 1, L = 1，因此：
//
//     Re = 1 / nu
nu {config.diffusivity:.16g};
""",
        encoding="utf-8",
    )


def _write_poisson_fv_solution(case: Path, config: CaseConfig) -> None:
    """Write linear-solver controls for the steady Poisson matrix."""
    (case / "system" / "fvSolution").write_text(
        f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /    O  peration   |
    \\\\  /    A nd           |
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSolution;
}}

solvers
{{
    // 对应数学方程组：
    //     A phi = b
    //
    // poissonFoamStudent 中的：
    //     phiEqn.solve()
    //
    // 会读取这里的线性求解器和容差。
    {config.scalar_field}
    {{
        solver          {config.linear_solver};
        tolerance       {config.linear_tolerance:.16g};
        relTol          0;
        smoother        GaussSeidel;
    }}
}}
""",
        encoding="utf-8",
    )


def _config_reference(config: CaseConfig) -> str:
    try:
        relative = config.path.relative_to(PROJECT_ROOT)
        return f"$projectRoot/{relative.as_posix()}"
    except ValueError:
        return config.path.as_posix()


def _write_case_scripts(case: Path, config: CaseConfig, resolution: int) -> None:
    config_reference = _config_reference(config)
    solver_path = (
        f"$projectRoot/build/{config.solver_family}/bin/{config.solver}"
    )
    hybrid_allrun_pre = None
    hybrid_create_patch = None
    if config.mesh_type == "hybrid":
        gmsh_python = config.gmsh_python or "/home/a776/vibeflow/python-env/bin/python"
        boundary_layer_thickness = (
            config.boundary_layer_thickness
            if config.boundary_layer_thickness is not None
            else 0.12
        )
        boundary_layer_count = (
            config.boundary_layer_count
            if config.boundary_layer_count is not None
            else max(2, round(resolution * boundary_layer_thickness))
        )
        hybrid_allrun_pre = f"""#!/bin/sh

set -eu

caseDir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
projectRoot=$(CDPATH= cd -- "$caseDir/../../../.." && pwd)
cd "$caseDir"

: "${{WM_PROJECT_DIR:?Please source /opt/openfoam14/etc/bashrc first}}"
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"

# 混合网格方腔流：
# 外圈为结构化四边形边界层，中心区域为三角形非结构网格。
# 该预处理与第一题混合网格一致，只是求解器换成了第五题的
# pisoFoamStudent。
gmshPython="${{VIBEFLOW_PYTHON:-{gmsh_python}}}"
if [ ! -x "$gmshPython" ]; then
    echo "Missing Gmsh Python interpreter: $gmshPython" >&2
    echo "Set VIBEFLOW_PYTHON to a Python environment containing gmsh." >&2
    exit 1
fi

sh "$caseDir/Allclean"
"""
        if config.problem == "tri_lid_driven_cavity":
            hybrid_allrun_pre = f"""#!/bin/sh

set -eu

caseDir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
projectRoot=$(CDPATH= cd -- "$caseDir/../../../.." && pwd)
cd "$caseDir"

: "${{WM_PROJECT_DIR:?Please source /opt/openfoam14/etc/bashrc first}}"
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"

# 第五题第二个算例：等边三角形顶盖驱动腔流。
# 外圈三条边生成边界层四边形带，中心保留三角形核心区。
gmshPython="${{VIBEFLOW_PYTHON:-{gmsh_python}}}"
if [ ! -x "$gmshPython" ]; then
    echo "Missing Gmsh Python interpreter: $gmshPython" >&2
    echo "Set VIBEFLOW_PYTHON to a Python environment containing gmsh." >&2
    exit 1
fi

sh "$caseDir/Allclean"
runApplication "$gmshPython" "$projectRoot/scripts/common/gmsh_triangular_cavity_hybrid_mesh.py" \
    --case "$caseDir" \
    --resolution {resolution} \
    --boundary-layer-thickness {boundary_layer_thickness:.16g} \
    --boundary-layer-count {boundary_layer_count} \
    --thickness {float(config.thickness):.16g}
"""
        else:
            hybrid_allrun_pre = f"""#!/bin/sh

set -eu

caseDir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
projectRoot=$(CDPATH= cd -- "$caseDir/../../../.." && pwd)
cd "$caseDir"

: "${{WM_PROJECT_DIR:?Please source /opt/openfoam14/etc/bashrc first}}"
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"

# 第五题第一个算例：方腔混合网格。
# 外圈为结构化四边形边界层，中心区域为三角形非结构网格。
gmshPython="${{VIBEFLOW_PYTHON:-{gmsh_python}}}"
if [ ! -x "$gmshPython" ]; then
    echo "Missing Gmsh Python interpreter: $gmshPython" >&2
    echo "Set VIBEFLOW_PYTHON to a Python environment containing gmsh." >&2
    exit 1
fi

sh "$caseDir/Allclean"
runApplication "$gmshPython" "$projectRoot/scripts/common/gmsh_hybrid_cavity_mesh.py" \
    --case "$caseDir" \
    --resolution {resolution} \
    --boundary-layer-thickness {boundary_layer_thickness:.16g} \
    --boundary-layer-count {boundary_layer_count} \
    --thickness {float(config.thickness):.16g}
"""
        hybrid_allrun_pre += f"""
runApplication -overwrite gmshToFoam "$caseDir/mesh/mesh.msh"
runApplication -overwrite createPatch
runApplication -overwrite checkMesh
runApplication -suffix centres foamPostProcess -constant -func writeCellCentres
runApplication -suffix volumes foamPostProcess -constant -func writeCellVolumes
python3 "$projectRoot/scripts/prepare_case.py" \
    --config "{config_reference}" \
    --N {resolution} \
    --refresh-initial-only
rm -rf "$caseDir/0"
cp -R "$caseDir/0.orig" "$caseDir/0"
"""
        hybrid_create_patch = """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /    O  peration     | Mesh patch construction
\\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       dictionary;
    location    "system";
    object      createPatchDict;
}

pointSync false;
writeCyclicMatch false;

patches
{
    front
    {
        patchInfo { type empty; }
        constructFrom patches;
        patches (frontSource);
    }
    back
    {
        patchInfo { type empty; }
        constructFrom patches;
        patches (backSource);
    }
}
"""
    poisson_run_comments = ""
    if config.problem in POISSON_PROBLEMS:
        poisson_run_comments = """# 第四题 Poisson 案例的执行顺序：
# 1. 清理上一次运行产生的网格、时间目录和日志；
# 2. 按 JSON 中的 N 生成四边形 blockMesh 或三角形 Gmsh 网格；
# 3. 生成与网格单元中心一致的 phi 初值和 omega 源项；
# 4. 运行 poissonFoamStudent，求解 fvm::laplacian(phi) == omega；
# 5. run_study.py 随后读取最终 phi，计算解析解误差并生成图表。
"""
    allclean = """#!/bin/sh

set -eu

caseDir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$caseDir"
rm -rf constant/polyMesh 0 postProcessing
for timeDir in ./[0-9]*; do
    [ "$timeDir" = "./0.orig" ] && continue
    [ -d "$timeDir" ] && rm -rf "$timeDir"
done
rm -f log.* run.batch.log
"""
    if config.mesh_type == "tri":
        allclean += """rm -rf mesh
rm -f constant/C constant/Cc* constant/Vc
"""
    if config.mesh_type == "quad":
        allrun_pre = None
        allrun = f"""#!/bin/sh

set -eu

caseDir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
    projectRoot=$(CDPATH= cd -- "$caseDir/../../../.." && pwd)
cd "$caseDir"

: "${{WM_PROJECT_DIR:?Please source /opt/openfoam14/etc/bashrc first}}"
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"

{poisson_run_comments}
solverPath="{solver_path}"
if [ ! -x "$solverPath" ]; then
    echo "Missing solver: $solverPath" >&2
    echo "Run: sh $projectRoot/scripts/build_student_solver.sh" >&2
    exit 1
fi

sh "$caseDir/Allclean"
python3 "$projectRoot/scripts/prepare_case.py" \\
    --config "{config_reference}" \\
    --N {resolution} \\
    --refresh-initial-only
cp -R "$caseDir/0.orig" "$caseDir/0"
runApplication -overwrite blockMesh
runApplication -overwrite checkMesh
runApplication -overwrite "$solverPath"
"""
    elif config.mesh_type == "tri":
        gmsh_python = config.gmsh_python or "/home/a776/vibeflow/python-env/bin/python"
        allrun_pre = f"""#!/bin/sh

set -eu

caseDir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
    projectRoot=$(CDPATH= cd -- "$caseDir/../../../.." && pwd)
cd "$caseDir"

: "${{WM_PROJECT_DIR:?Please source /opt/openfoam14/etc/bashrc first}}"
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"

# 三角形 Poisson 案例的预处理入口：
# Gmsh 生成二维三角形并拉伸为 OpenFOAM 薄棱柱；
# gmshToFoam 导入网格，createPatch 恢复四条物理边界；
# C/Vc 写出真实单元中心和体积，随后按这些中心生成 omega。
gmshPython="${{VIBEFLOW_PYTHON:-{gmsh_python}}}"
if [ ! -x "$gmshPython" ]; then
    echo "Missing Gmsh Python interpreter: $gmshPython" >&2
    echo "Set VIBEFLOW_PYTHON to a Python environment containing gmsh." >&2
    exit 1
fi

sh "$caseDir/Allclean"
runApplication "$gmshPython" "$projectRoot/scripts/common/gmsh_tri_mesh.py" \\
    --case "$caseDir" \\
    --resolution {resolution} \\
    --thickness {config.thickness:g} \\
    --xmin {config.domain[0]:g} \\
    --xmax {config.domain[1]:g} \\
    --ymin {config.domain[2]:g} \\
    --ymax {config.domain[3]:g}
runApplication -overwrite gmshToFoam "$caseDir/mesh/mesh.msh"
runApplication -overwrite createPatch
runApplication -overwrite checkMesh
runApplication -suffix centres foamPostProcess -constant -func writeCellCentres
runApplication -suffix volumes foamPostProcess -constant -func writeCellVolumes
python3 "$projectRoot/scripts/prepare_case.py" \\
    --config "{config_reference}" \\
    --N {resolution} \\
    --refresh-initial-only
rm -rf "$caseDir/0"
cp -R "$caseDir/0.orig" "$caseDir/0"
"""
        allrun = f"""#!/bin/sh

set -eu

caseDir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
    projectRoot=$(CDPATH= cd -- "$caseDir/../../../.." && pwd)
cd "$caseDir"

: "${{WM_PROJECT_DIR:?Please source /opt/openfoam14/etc/bashrc first}}"
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"

{poisson_run_comments}
solverPath="{solver_path}"
if [ ! -x "$solverPath" ]; then
    echo "Missing solver: $solverPath" >&2
    echo "Run: sh $projectRoot/scripts/build_student_solver.sh" >&2
    exit 1
fi

sh "$caseDir/Allrun.pre"
runApplication -overwrite "$solverPath"
"""
    elif config.mesh_type == "hybrid":
        allrun_pre = hybrid_allrun_pre
        allrun = f"""#!/bin/sh

set -eu

caseDir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
projectRoot=$(CDPATH= cd -- "$caseDir/../../../.." && pwd)
cd "$caseDir"

: "${{WM_PROJECT_DIR:?Please source /opt/openfoam14/etc/bashrc first}}"
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"

solverPath="{solver_path}"
if [ ! -x "$solverPath" ]; then
    echo "Missing solver: $solverPath" >&2
    echo "Run: sh $projectRoot/scripts/build_student_solver.sh" >&2
    exit 1
fi

sh "$caseDir/Allrun.pre"
runApplication -overwrite "$solverPath"
"""
    else:
        raise NotImplementedError(f"Unsupported mesh type: {config.mesh_type}")
    (case / "Allclean").write_text(allclean, encoding="utf-8")
    (case / "Allrun").write_text(allrun, encoding="utf-8")
    (case / "Allclean").chmod(0o755)
    (case / "Allrun").chmod(0o755)
    if allrun_pre is not None:
        (case / "Allrun.pre").write_text(allrun_pre, encoding="utf-8")
        (case / "Allrun.pre").chmod(0o755)
    (case / "case.foam").touch()


def prepare_case(
    config: CaseConfig,
    resolution: int,
    overwrite: bool = False,
    refresh_initial_only: bool = False,
) -> Path:
    """Prepare one configured OpenFOAM case directory."""
    config.require_implemented()
    if config.mesh_type not in {"quad", "tri", "hybrid"}:
        raise NotImplementedError(f"Unsupported mesh type: {config.mesh_type}")

    target = config.case_dir(resolution)
    if refresh_initial_only:
        if not target.exists():
            raise RuntimeError(f"Case directory does not exist: {target}")
        if config.problem in {"solid_rotation_advection", "rotating_peak_advection_diffusion"}:
            _patch_velocity_field(target, config)
        output = _write_initial_field(target, config)
        print(f"case={target}")
        print(f"resolution={resolution}")
        print(f"initialField={output}")
        return target

    if target.exists() and any(target.iterdir()):
        if not overwrite:
            raise RuntimeError(
                f"Case already exists: {target}. Use --overwrite to rebuild it."
            )
        shutil.rmtree(target)
    target.mkdir(parents=True, exist_ok=True)
    (target / "0.orig").mkdir()
    (target / "constant").mkdir()
    (target / "system").mkdir()
    _write_base_system_files(target, config)

    if config.mesh_type == "quad":
        _write_base_block_mesh(target, config, resolution)
        if config.problem == "sine_wave_advection_diffusion":
            block_mesh_path = target / "system" / "blockMeshDict"
            block_mesh_text = block_mesh_path.read_text(encoding="utf-8")
            block_mesh_text = re.sub(
                r"// 六面体控制体：\n"
                r"\s*//   x 方向 \d+ 个；\n"
                r"\s*//   y 方向 \d+ 个；\n"
                r"\s*//   z 方向 1 个。\n"
                r"\s*//\n"
                r"\s*// 总单元数是 \d+\*\d+\*1 = \d+。",
                "// 六面体控制体：\n"
                f"    //   x 方向 {resolution} 个；\n"
                f"    //   y 方向 {resolution} 个；\n"
                "    //   z 方向 1 个。\n"
                "    //\n"
                f"    // 总单元数是 {resolution}*{resolution}*1 = {resolution * resolution}。",
                block_mesh_text,
                count=1,
            )
            block_mesh_path.write_text(block_mesh_text, encoding="utf-8")
        if config.boundary_condition != "periodicXY":
            _patch_block_mesh_outer_boundaries(target)
    elif config.mesh_type == "tri":
        _write_create_patch_dict(target, config)
        (target / "mesh").mkdir(parents=True, exist_ok=True)
    else:
        # 混合网格方腔流的网格与 patch 由 Allrun.pre / gmsh 生成，
        # 这里仅保留 system 目录和后续写初值所需的结构。
        (target / "mesh").mkdir(parents=True, exist_ok=True)
        (target / "system" / "createPatchDict").write_text(
            """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /    O  peration     | Mesh patch construction
\\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       dictionary;
    location    "system";
    object      createPatchDict;
}

pointSync false;
writeCyclicMatch false;

patches
{
    front
    {
        patchInfo { type empty; }
        constructFrom patches;
        patches (frontSource);
    }
    back
    {
        patchInfo { type empty; }
        constructFrom patches;
        patches (backSource);
    }
}
""",
            encoding="utf-8",
        )
    if config.problem in POISSON_PROBLEMS:
        transport_path = target / "constant" / "transportProperties"
        if transport_path.exists():
            transport_path.unlink()
    _patch_fv_schemes(target, config)
    _patch_control_dict(target, config)
    if config.problem in LID_CAVITY_PROBLEMS:
        _write_lid_cavity_physical_properties(target, config)
        _write_lid_cavity_fv_solution(target, config)
    if config.problem.startswith("diffusion_") or config.problem in ADVECTION_DIFFUSION_PROBLEMS:
        _write_transport_properties(target, config)
    _write_base_velocity_field(target, config)
    _patch_velocity_field(target, config)
    if config.problem in ADVECTION_DIFFUSION_PROBLEMS:
        stale_scalar = target / "0.orig" / "T"
        if stale_scalar.exists():
            stale_scalar.unlink()
        _write_explicit_advection_diffusion_fv_solution(target)
    if config.problem in POISSON_PROBLEMS:
        _write_poisson_fv_solution(target, config)
    initial_field: Path | None = None
    if config.mesh_type == "quad" or config.problem in LID_CAVITY_PROBLEMS:
        initial_field = _write_initial_field(target, config)
    _write_case_scripts(target, config, resolution)

    metadata = {
        "solverFamily": config.solver_family,
        "caseName": config.case_name,
        "config": str(config.path),
        "equation": config.equation,
        "problem": config.problem,
        "meshType": config.mesh_type,
        "meshBackend": config.mesh_backend,
        "schemeName": config.scheme_name,
        "divScheme": config.div_scheme,
        "resolution": resolution,
        "solver": config.solver,
        "velocity": list(config.velocity),
        "domain": list(config.domain),
        "velocityModel": config.velocity_model,
        "initialProfile": config.initial_profile,
        "boundaryCondition": config.boundary_condition,
        "postprocess": config.postprocess,
        "scalarField": config.scalar_field,
        "sourceField": config.source_field,
        "linearSolver": config.linear_solver,
        "linearTolerance": config.linear_tolerance,
        "nNonOrthogonalCorrectors": config.n_non_orthogonal_correctors,
        "mu": config.diffusivity,
        "diffusionCo": config.diffusion_co,
        "advectionDiffusionCo": config.advection_diffusion_co,
        "laplacianScheme": config.laplacian_scheme,
        "snGradScheme": config.sn_grad_scheme,
        "endTime": config.end_time,
        "deltaT": config.delta_t,
        "adjustTimeStep": config.adjust_time_step,
        "writeInterval": config.write_interval,
        "maxCo": config.max_co,
        "maxDeltaT": config.max_delta_t,
        "steadyStateControl": config.steady_state_control,
        "steadyVelocityTol": config.steady_velocity_tol,
        "steadyMassTol": config.steady_mass_tol,
        "minimumSteadySteps": config.minimum_steady_steps,
        "requiredSteadySteps": config.required_steady_steps,
        "initialField": str(initial_field) if initial_field else "generated-after-mesh",
    }
    if config.problem == "lid_driven_cavity":
        metadata.update(
            {
                "reynolds": round(1.0 / config.diffusivity),
                "viscosity": config.diffusivity,
                "cellsPerEdge": resolution,
                "boundaryLayerThickness": config.boundary_layer_thickness,
                "boundaryLayerCount": config.boundary_layer_count,
            }
        )
    if config.problem == "tri_lid_driven_cavity":
        metadata.update(
            {
                "reynolds": round(1.0 / config.diffusivity),
                "viscosity": config.diffusivity,
                "cellsPerEdge": resolution,
                "boundaryLayerThickness": config.boundary_layer_thickness,
                "boundaryLayerCount": config.boundary_layer_count,
            }
        )
    (target / "metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(f"prepared={target}")
    print(f"caseName={config.case_name}")
    print(f"resolution={resolution}")
    print(f"scheme={config.div_scheme or config.laplacian_scheme}")
    return target


def run_case(case: Path, bashrc: Path) -> None:
    """Run one prepared case through its Allrun script."""
    command = (
        f"source {shlex.quote(str(bashrc))} && "
        f"sh {shlex.quote(str(case / 'Allrun'))}"
    )
    log_path = case / "run.batch.log"
    with log_path.open("w", encoding="utf-8") as stream:
        subprocess.run(
            ["bash", "-lc", command],
            cwd=case,
            stdout=stream,
            stderr=subprocess.STDOUT,
            check=True,
        )


def postprocess_configured_case(config: CaseConfig, resolution: int) -> None:
    """Post-process one configured case after it has run."""
    if config.problem == "tri_lid_driven_cavity":
        from postprocess_triangular_cavity import postprocess

        postprocess(config.case_dir(resolution))
        return
    if config.problem == "lid_driven_cavity":
        from postprocess_lid_cavity import postprocess as postprocess_lid_cavity_case

        postprocess_lid_cavity_case(config.case_dir(resolution))
        return

    if config.problem not in {
        "sine_wave_advection",
        "solid_rotation_advection",
        "diffusion_discontinuity",
        "diffusion_gaussian",
        "sine_wave_advection_diffusion",
        "rotating_peak_advection_diffusion",
        "poisson_manufactured",
    }:
        raise NotImplementedError(f"Unsupported postprocess problem: {config.problem}")
    from postprocess_case import postprocess_case

    target_monitor = (
        config.advection_diffusion_co
        if config.problem in ADVECTION_DIFFUSION_PROBLEMS
        else config.diffusion_co
        if config.problem.startswith("diffusion_")
        else config.max_co
    )
    postprocess_case(config.case_dir(resolution), PROJECT_ROOT, target_monitor)


def run_study(
    config: CaseConfig,
    resolutions: list[int],
    overwrite: bool,
    prepare_only: bool,
    bashrc: Path,
) -> None:
    """Prepare and optionally run all requested resolutions."""
    config.require_implemented()
    print(f"caseName={config.case_name}")
    print(f"problem={config.problem}")
    print(f"meshType={config.mesh_type}")
    print(f"scheme={config.div_scheme or config.laplacian_scheme}")
    print(f"resolutions={','.join(str(value) for value in resolutions)}")

    for resolution in resolutions:
        case = prepare_case(config, resolution, overwrite=overwrite)
        if prepare_only:
            continue
        print(f"running={case}")
        run_case(case, bashrc)
        postprocess_configured_case(config, resolution)

    if not prepare_only and config.problem in {
        "sine_wave_advection",
        "diffusion_discontinuity",
        "sine_wave_advection_diffusion",
        "rotating_peak_advection_diffusion",
        "poisson_manufactured",
    }:
        from study_analysis import analyse, collect, plot

        collect(config.solver_family, config.case_name, resolutions)
        analyse(config.solver_family, config.case_name)
        plot(config.solver_family, config.case_name)
