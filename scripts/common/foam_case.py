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
from foam_fields import patch_uniform_vector_field, read_cell_geometry
from mesh_tools import mesh_resolution, patch_block_mesh_resolution
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


def _copy_constant_inputs(template: Path, target: Path) -> None:
    source = template / "constant"
    destination = target / "constant"
    destination.mkdir(parents=True, exist_ok=True)
    if not source.exists():
        return
    for child in source.iterdir():
        if child.name == "polyMesh":
            continue
        target_child = destination / child.name
        if child.is_dir():
            shutil.copytree(child, target_child)
        else:
            shutil.copy2(child, target_child)


def _replace_or_append_dictionary_entry(path: Path, name: str, value: str) -> None:
    text = path.read_text(encoding="utf-8")
    pattern = re.compile(rf"^(\s*{re.escape(name)}\s+)([^;]+);", re.M)
    updated, count = pattern.subn(rf"\g<1>{value};", text, count=1)
    if count == 0:
        updated = text.rstrip() + f"\n{name:<16} {value};\n"
    path.write_text(updated, encoding="utf-8")


def _patch_fv_schemes(case: Path, config: CaseConfig) -> None:
    path = case / "system" / "fvSchemes"
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
    _replace_or_append_dictionary_entry(path, "scalarField", config.scalar_field)
    if config.problem.startswith("diffusion_"):
        _replace_or_append_dictionary_entry(path, "diffusionCo", f"{config.diffusion_co:g}")
        if config.max_delta_t is not None:
            _replace_or_append_dictionary_entry(path, "maxDeltaT", f"{config.max_delta_t:g}")
    # Keep the final time directory name consistent with the configured
    # terminal time, including values such as 2*pi.
    _replace_or_append_dictionary_entry(path, "timePrecision", "17")


def _patch_block_mesh_domain(case: Path, config: CaseConfig) -> None:
    """Patch a single-block 2-D box to the configured rectangular domain."""
    path = case / "system" / "blockMeshDict"
    text = path.read_text(encoding="utf-8")
    xmin, xmax, ymin, ymax = config.domain
    zmax = config.thickness
    vertices = f"""vertices
(
    ({xmin:g} {ymin:g} 0)
    ({xmax:g} {ymin:g} 0)
    ({xmax:g} {ymax:g} 0)
    ({xmin:g} {ymax:g} 0)
    ({xmin:g} {ymin:g} {zmax:g})
    ({xmax:g} {ymin:g} {zmax:g})
    ({xmax:g} {ymax:g} {zmax:g})
    ({xmin:g} {ymax:g} {zmax:g})
);"""
    updated, count = re.subn(r"vertices\s*\(.*?\)\s*;", vertices, text, count=1, flags=re.S)
    if count != 1:
        raise RuntimeError(f"Cannot patch vertices in {path}")
    path.write_text(updated, encoding="utf-8")


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
    """Apply the configured constant velocity to the template field."""
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
    if config.mesh_type not in {"quad", "tri"}:
        raise NotImplementedError(f"Unsupported mesh type: {config.mesh_type}")

    target = config.case_dir(resolution)
    template = config.template_case
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

    if not template.exists():
        raise RuntimeError(f"Template case is missing: {template}")

    if target == template:
        if overwrite and (target / "Allclean").exists():
            subprocess.run(["sh", str(target / "Allclean")], check=True)
    else:
        if target.exists() and any(target.iterdir()):
            if not overwrite:
                raise RuntimeError(
                    f"Case already exists: {target}. Use --overwrite to rebuild it."
                )
            shutil.rmtree(target)
        target.mkdir(parents=True, exist_ok=True)
        shutil.copytree(template / "0.orig", target / "0.orig")
        shutil.copytree(template / "system", target / "system")
        _copy_constant_inputs(template, target)

    if config.mesh_type == "quad":
        _patch_block_mesh_domain(target, config)
        patch_block_mesh_resolution(target, resolution)
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
    else:
        block_mesh_dict = target / "system" / "blockMeshDict"
        if block_mesh_dict.exists():
            block_mesh_dict.unlink()
        _write_create_patch_dict(target, config)
        (target / "mesh").mkdir(parents=True, exist_ok=True)
    _patch_fv_schemes(target, config)
    _patch_control_dict(target, config)
    if config.problem.startswith("diffusion_") or config.problem in ADVECTION_DIFFUSION_PROBLEMS:
        _write_transport_properties(target, config)
    _patch_velocity_field(target, config)
    if config.problem in ADVECTION_DIFFUSION_PROBLEMS:
        stale_scalar = target / "0.orig" / "T"
        if stale_scalar.exists():
            stale_scalar.unlink()
        _write_explicit_advection_diffusion_fv_solution(target)
    initial_field: Path | None = None
    if config.mesh_type == "quad":
        initial_field = _write_initial_field(target, config)
    _write_case_scripts(target, config, resolution)

    metadata = {
        "solverFamily": config.solver_family,
        "templateSolverFamily": config.template_solver_family,
        "templateCaseName": config.template_case_name,
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
        "mu": config.diffusivity,
        "diffusionCo": config.diffusion_co,
        "advectionDiffusionCo": config.advection_diffusion_co,
        "laplacianScheme": config.laplacian_scheme,
        "snGradScheme": config.sn_grad_scheme,
        "endTime": config.end_time,
        "maxCo": config.max_co,
        "maxDeltaT": config.max_delta_t,
        "initialField": str(initial_field) if initial_field else "generated-after-mesh",
    }
    (target / "metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(f"prepared={target}")
    print(f"caseName={config.case_name}")
    print(f"resolution={resolution}")
    print(f"scheme={config.div_scheme}")
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
    if config.problem not in {
        "sine_wave_advection",
        "solid_rotation_advection",
        "diffusion_discontinuity",
        "diffusion_gaussian",
        "sine_wave_advection_diffusion",
        "rotating_peak_advection_diffusion",
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
    print(f"scheme={config.div_scheme}")
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
    }:
        from study_analysis import analyse, collect, plot

        collect(config.solver_family, config.case_name, resolutions)
        analyse(config.solver_family, config.case_name)
        plot(config.solver_family, config.case_name)
