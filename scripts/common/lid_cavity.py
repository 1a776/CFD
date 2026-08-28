#!/usr/bin/env python3

"""Prepare and run one fifth-problem square lid-driven cavity case."""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CASES_ROOT = PROJECT_ROOT / "cases" / "05_navier_stokes_equation"


def load_case_config(path: Path) -> dict:
    data = json.loads(path.read_text(encoding="utf-8"))
    required = {
        "caseName",
        "solverFamily",
        "solver",
        "reynolds",
        "viscosity",
        "meshLevel",
        "cellsPerEdge",
    }
    missing = sorted(required - data.keys())
    if missing:
        raise RuntimeError(f"Missing configuration fields: {', '.join(missing)}")
    if data["solverFamily"] != "05_navier_stokes_equation":
        raise RuntimeError("This runner only accepts the fifth solver family")
    if data["solver"] not in {"projectionFoamStudent", "pisoFoamStudent"}:
        raise RuntimeError(
            "This runner only accepts projectionFoamStudent or pisoFoamStudent"
        )
    return data


def case_dir(config: dict) -> Path:
    return CASES_ROOT / str(config["caseName"])


def _block_mesh(config: dict) -> str:
    n = int(config["cellsPerEdge"])
    return f"""/*--------------------------------*- C++ -*----------------------------------*\\
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
    object      blockMeshDict;
}}

// 第五题第一个算例：方形顶盖驱动腔流。
//
// 数学区域：
//
//     Ω = [0,1] × [0,1]
//
// 本配置中的 N 表示每条边方向的单元数：
//
//     N = {n}
//
// 该案例使用规则四边形基准网格。后续若替换为混合非结构网格，
// 只需替换本文件和网格生成步骤，U、p、physicalProperties 及求解器接口不变。
convertToMeters 1;

vertices
(
    (0 0 0)
    (1 0 0)
    (1 1 0)
    (0 1 0)
    (0 0 0.1)
    (1 0 0.1)
    (1 1 0.1)
    (0 1 0.1)
);

blocks
(
    // 二维计算通过 z 方向一个单元和 empty 边界实现。
    hex (0 1 2 3 4 5 6 7) ({n} {n} 1) simpleGrading (1 1 1)
);

edges
(
);

boundary
(
    leftWall
    {{
        type wall;
        faces ((0 4 7 3));
    }}
    rightWall
    {{
        type wall;
        faces ((1 2 6 5));
    }}
    bottomWall
    {{
        type wall;
        faces ((0 1 5 4));
    }}
    movingTop
    {{
        type wall;
        faces ((3 7 6 2));
    }}
    frontAndBack
    {{
        type empty;
        faces
        (
            (0 3 2 1)
            (4 5 6 7)
        );
    }}
);
"""


def _hybrid_allrun_pre(config: dict) -> str:
    resolution = int(config["cellsPerEdge"])
    thickness = float(config.get("thickness", 0.1))
    layer_thickness = float(config.get("boundaryLayerThickness", 0.12))
    layer_count = int(
        config.get(
            "boundaryLayerCount",
            max(2, round(resolution * layer_thickness)),
        )
    )
    gmsh_python = config.get(
        "gmshPython", "/home/a776/vibeflow/python-env/bin/python"
    )
    return f"""#!/bin/sh

set -eu

caseDir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
projectRoot=$(CDPATH= cd -- "$caseDir/../../.." && pwd)
cd "$caseDir"

: "${{WM_PROJECT_DIR:?Please source /opt/openfoam14/etc/bashrc first}}"
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"

# 第五题第一案例的 Figure 2 混合网格：
# 外圈为结构化四边形边界层，中心区域为三角形非结构网格。
# 网格参数来自当前案例 JSON，不改变求解器和物理边界接口。
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
    --boundary-layer-thickness {layer_thickness:.16g} \
    --boundary-layer-count {layer_count} \
    --thickness {thickness:.16g}
runApplication -overwrite gmshToFoam "$caseDir/mesh/mesh.msh"
runApplication -overwrite createPatch
runApplication -overwrite checkMesh
"""


def _hybrid_create_patch() -> str:
    return """/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /    O  peration   | Mesh patch construction
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


def _fields(config: dict) -> tuple[str, str, str]:
    viscosity = float(config["viscosity"])
    mesh_type = str(config.get("meshType", "quad")).lower()
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
    u_field = """/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /    O  peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       volVectorField;
    location    "0";
    object      U;
}

// 初始条件：
//
//     U(x,y,0) = (0,0,0)
//
// 顶盖边界：
//
//     U(x,1,t) = (1,0,0)
//
// 其它三面：
//
//     U = (0,0,0)
dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField
{
    leftWall
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
    }
{empty_patches}
}
"""
    p_field = """/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /    O  peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       volScalarField;
    location    "0";
    object      p;
}

// 本项目使用运动学压力：
//
//     p_k = p / ρ
//
// 初始压力：
//
//     p_k(x,y,0) = 0
//
// 固壁法向压力条件：
//
//     ∂p_k/∂n = 0
//
// 压力参考自由度在 system/fvSolution/PISO 中设置。
dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{
    leftWall
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
    }
{empty_patches}
}
"""
    physical = f"""/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /    O  peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       dictionary;
    location    "constant";
    object      physicalProperties;
}}

// Reynolds 数定义：
//
//     Re = U L / ν
//
// 本案例取 U=1、L=1，因此：
//
//     Re = 1 / ν
//
// 当前工况：
//
//     Re = {config["reynolds"]}
//     ν  = {viscosity:.16g}
nu {viscosity:.16g};
"""
    u_field = u_field.replace("{empty_patches}", empty_patches)
    p_field = p_field.replace("{empty_patches}", empty_patches)
    return u_field, p_field, physical


def _system(config: dict) -> tuple[str, str, str]:
    delta_t = float(config.get("deltaT", 0.001))
    end_time = float(config.get("endTime", 10.0))
    max_co = float(config.get("maxCo", 0.2))
    write_interval = float(config.get("writeInterval", 1.0))
    steady_state_control = str(config.get("steadyStateControl", "true")).lower()
    steady_velocity_tol = float(config.get("steadyVelocityTol", 1e-6))
    steady_mass_tol = float(config.get("steadyMassTol", 1e-8))
    minimum_steady_steps = int(config.get("minimumSteadySteps", 1000))
    required_steady_steps = int(config.get("requiredSteadySteps", 20))
    div_scheme = str(config.get("divScheme", "Gauss linearUpwind grad(U)"))
    laplacian_scheme = str(config.get("laplacianScheme", "Gauss linear corrected"))
    n_corr = int(config.get("nCorrectors", 2))
    n_non_orth = int(config.get("nNonOrthogonalCorrectors", 1))
    tolerance = float(config.get("linearTolerance", 1e-8))
    solver_name = str(config["solver"])
    is_projection = solver_name == "projectionFoamStudent"
    adjust_time_step = str(
        config.get("adjustTimeStep", not is_projection)
    ).lower()
    max_delta_t = float(config.get("maxDeltaT", 0.01))

    control = f"""/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /    O  peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}}

application {solver_name};
startFrom startTime;
startTime 0;
stopAt endTime;
endTime {end_time:.16g};
deltaT {delta_t:.16g};
// PISO 求解器可读取下面三个字段做自适应时间步：
//     Co_max = max_c [ Δt / (2 V_c) Σ_f |Φ_f| ]
//     Δt_new = min(1.2 Δt_old, maxCo / Co_max Δt_old, maxDeltaT)
// projectionFoamStudent 当前仍使用固定 deltaT。
adjustTimeStep {adjust_time_step};
maxCo {max_co:.16g};
maxDeltaT {max_delta_t:.16g};
writeControl runTime;
writeInterval {write_interval:.16g};
purgeWrite 0;
writeFormat ascii;
writePrecision 12;
timePrecision 12;
writeCompression off;
runTimeModifiable no;

// 稳态控制：
//     max |Uⁿ - Uⁿ⁻¹| <= steadyVelocityTol
//     max |div(Uⁿ)|   <= steadyMassTol
// 两个条件连续满足 requiredSteadySteps 次，且至少推进
// minimumSteadySteps 个时间步后，求解器才提前结束。
steadyStateControl {steady_state_control};
steadyVelocityTol {steady_velocity_tol:.16g};
steadyMassTol {steady_mass_tol:.16g};
minimumSteadySteps {minimum_steady_steps};
requiredSteadySteps {required_steady_steps};
"""
    pressure_laplacian_key = "laplacian(dtCoeff,p)" if is_projection else "laplacian((1|A(U)),p)"
    pressure_laplacian_comment = (
        "    // projectionFoamStudent 的压力方程使用 dtCoeff 这一临时系数.\n"
        if is_projection
        else "    // pisoFoamStudent 的压力方程使用 (1|A(U)) 这一动量对角系数.\n"
    )
    div_projection_block = (
        "    // 预测速度 UStar 使用同一个对流离散格式.\n"
        "    // 这是源码 fvm::div(phi, UStar) 的字段名匹配入口.\n"
        f"    div(phi,UStar) {div_scheme};\n"
        if is_projection
        else ""
    )
    laplacian_projection_block = (
        "    // 预测速度方程实际调用 fvm::laplacian(nu, UStar).\n"
        f"    laplacian(nu,UStar) {laplacian_scheme};\n"
        if is_projection
        else ""
    )
    ustar_block = (
        "    // projectionFoamStudent 在预测步骤中组装并求解 UStarEqn,\n"
        "    // 因此线性求解器名称必须与矩阵字段 UStar 对应.\n"
        "    UStar\n"
        "    {\n"
        "        solver smoothSolver;\n"
        "        smoother symGaussSeidel;\n"
        f"        tolerance {tolerance:.16g};\n"
        "        relTol 0;\n"
        "    }\n"
        if is_projection
        else ""
    )
    pressure_laplacian_block = (
        pressure_laplacian_comment
        + f"    {pressure_laplacian_key} {laplacian_scheme};\n"
    )

    schemes = f"""/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /    O  peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSchemes;
}}

// 时间离散：
//
//     ∂U/∂t ≈ (Uⁿ⁺¹ - Uⁿ) / Δt
ddtSchemes
{{
    default Euler;
}}

// 单元中心梯度和非正交面法向梯度。
gradSchemes
{{
    default Gauss linear;
}}

divSchemes
{{
    default none;

    //     ∇·(U U) = 1/V_c Σ_f Φ_f U_f
    //     Φ_f = U_f · S_f
    div(phi,U) {div_scheme};
{div_projection_block}
}}

laplacianSchemes
{{
    default none;

    //     ∇·(ν∇U) = 1/V_c Σ_f (ν∇U)_f · S_cf
    laplacian(nu,U) {laplacian_scheme};
{laplacian_projection_block}

    // 压力方程对应的拉普拉斯项：
    //
    //     projectionFoamStudent:  ∇·(dtCoeff ∇p) = ∇·Φ*
    //     pisoFoamStudent:        ∇·((1|A(U)) ∇p) = ∇·Φ*
{pressure_laplacian_block}
}}

interpolationSchemes
{{
    default linear;
}}

snGradSchemes
{{
    default corrected;
}}
"""
    solution = f"""/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /    O  peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
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
    // 第五题的 PISO 参数入口。
    nCorrectors {n_corr};
    nNonOrthogonalCorrectors {n_non_orth};
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
        tolerance {tolerance:.16g};
        relTol 0;
    }}

{ustar_block}

    p
    {{
        solver GAMG;
        tolerance {tolerance:.16g};
        relTol 0;
        smoother GaussSeidel;
    }}

    // 末次压力校正使用 pFinal；PISO 与 projection 都保留这一入口。
    pFinal
    {{
        $p;
        relTol 0;
    }}
}}
"""
    return control, schemes, solution


def prepare(config_path: Path, overwrite: bool = False) -> Path:
    config = load_case_config(config_path)
    target = case_dir(config)
    if target.exists() and any(target.iterdir()):
        if not overwrite:
            raise RuntimeError(f"Case exists: {target}; use --overwrite")
        shutil.rmtree(target)

    (target / "0.orig").mkdir(parents=True)
    (target / "constant").mkdir()
    (target / "system").mkdir()
    mesh_type = str(config.get("meshType", "quad")).lower()
    if mesh_type == "hybrid":
        (target / "Allrun.pre").write_text(
            _hybrid_allrun_pre(config), encoding="utf-8"
        )
        (target / "Allrun.pre").chmod(0o755)
        (target / "system" / "createPatchDict").write_text(
            _hybrid_create_patch(), encoding="utf-8"
        )
    else:
        (target / "system" / "blockMeshDict").write_text(
            _block_mesh(config), encoding="utf-8"
        )
    u_field, p_field, physical_text = _fields(config)
    (target / "0.orig" / "U").write_text(u_field, encoding="utf-8")
    (target / "0.orig" / "p").write_text(p_field, encoding="utf-8")
    (target / "constant" / "physicalProperties").write_text(physical_text, encoding="utf-8")
    control, schemes, solution = _system(config)
    (target / "system" / "controlDict").write_text(control, encoding="utf-8")
    (target / "system" / "fvSchemes").write_text(schemes, encoding="utf-8")
    (target / "system" / "fvSolution").write_text(solution, encoding="utf-8")
    metadata = dict(config)
    metadata["casePath"] = str(target)
    metadata["mesh"] = {
        "type": mesh_type,
        "cellsPerEdge": config["cellsPerEdge"],
        "boundaryLayerThickness": config.get("boundaryLayerThickness"),
        "boundaryLayerCount": config.get("boundaryLayerCount"),
    }
    (target / "metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    (target / "case.foam").touch()
    _write_scripts(target, config)
    return target


def _write_scripts(target: Path, config: dict) -> None:
    solver = PROJECT_ROOT / "build" / config["solverFamily"] / "bin" / config["solver"]
    mesh_type = str(config.get("meshType", "quad")).lower()
    if mesh_type == "hybrid":
        mesh_commands = 'sh "$caseDir/Allrun.pre"'
    else:
        mesh_commands = """cp -R "$caseDir/0.orig" "$caseDir/0"
runApplication -overwrite blockMesh
runApplication -overwrite checkMesh"""
    allrun = f"""#!/bin/sh
set -eu
caseDir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$caseDir"
: "${{WM_PROJECT_DIR:?Please source /opt/openfoam14/etc/bashrc first}}"
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"
solverPath="{solver}"
if [ ! -x "$solverPath" ]; then
    echo "Missing solver: $solverPath" >&2
    echo "Run: sh {PROJECT_ROOT}/scripts/build_student_solver.sh" >&2
    exit 1
fi
sh "$caseDir/Allclean"
{mesh_commands}
rm -rf "$caseDir/0"
cp -R "$caseDir/0.orig" "$caseDir/0"
runApplication -overwrite "$solverPath"
runApplication -overwrite foamPostProcess -latestTime -func writeCellCentres
runApplication -overwrite foamPostProcess -latestTime -func writeCellVolumes
"""
    allclean = """#!/bin/sh
set -eu
caseDir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$caseDir"
rm -rf constant/polyMesh postProcessing
for timeDir in ./[0-9]*; do
    [ "$timeDir" = "./0.orig" ] && continue
    [ -d "$timeDir" ] && rm -rf "$timeDir"
done
# 保留 run.batch.log；该文件由外层 Python 在调用 Allrun 前打开。
rm -f log.*
"""
    (target / "Allrun").write_text(allrun, encoding="utf-8")
    (target / "Allclean").write_text(allclean, encoding="utf-8")
    (target / "Allrun").chmod(0o755)
    (target / "Allclean").chmod(0o755)


def run(config_path: Path, overwrite: bool = False) -> Path:
    config = load_case_config(config_path)
    target = prepare(config_path, overwrite=overwrite)
    log = target / "run.batch.log"
    with log.open("w", encoding="utf-8") as stream:
        subprocess.run(
            ["bash", "-lc", "source /opt/openfoam14/etc/bashrc && sh ./Allrun"],
            cwd=target,
            stdout=stream,
            stderr=subprocess.STDOUT,
            check=True,
        )
    subprocess.run(
        [
            "python3",
            str(PROJECT_ROOT / "scripts" / "postprocess_lid_cavity.py"),
            "--case",
            str(target),
        ],
        cwd=PROJECT_ROOT,
        check=True,
    )
    return target
