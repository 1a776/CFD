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
from case_config import CaseConfig
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
    text = path.read_text(encoding="utf-8")
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
    # Keep the final time directory name consistent with the configured
    # terminal time, including values such as 2*pi.
    _replace_or_append_dictionary_entry(path, "timePrecision", "17")


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
    if config.boundary_condition == "zeroScalarAtOuterBoundary":
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
    else:
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
        allrun = f"""#!/bin/sh

set -eu

caseDir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
    projectRoot=$(CDPATH= cd -- "$caseDir/../../../.." && pwd)
cd "$caseDir"

: "${{WM_PROJECT_DIR:?Please source /opt/openfoam14/etc/bashrc first}}"
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"

solverPath="{solver_path}"
gmshPython="${{VIBEFLOW_PYTHON:-{gmsh_python}}}"
if [ ! -x "$solverPath" ]; then
    echo "Missing solver: $solverPath" >&2
    echo "Run: sh $projectRoot/scripts/build_student_solver.sh" >&2
    exit 1
fi
if [ ! -x "$gmshPython" ]; then
    echo "Missing Gmsh Python interpreter: $gmshPython" >&2
    echo "Set VIBEFLOW_PYTHON to a Python environment containing gmsh." >&2
    exit 1
fi

sh "$caseDir/Allclean"
runApplication "$gmshPython" "$projectRoot/scripts/common/gmsh_tri_mesh.py" \\
    --case "$caseDir" \\
    --resolution {resolution} \\
    --thickness {config.thickness:g}
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
runApplication -overwrite "$solverPath"
"""
    else:
        raise NotImplementedError(f"Unsupported mesh type: {config.mesh_type}")
    (case / "Allclean").write_text(allclean, encoding="utf-8")
    (case / "Allrun").write_text(allrun, encoding="utf-8")
    (case / "Allclean").chmod(0o755)
    (case / "Allrun").chmod(0o755)
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
        if config.problem == "solid_rotation_advection":
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
        patch_block_mesh_resolution(target, resolution)
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
    _patch_velocity_field(target, config)
    initial_field: Path | None = None
    if config.mesh_type == "quad":
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
        "endTime": config.end_time,
        "maxCo": config.max_co,
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
    if config.problem not in {"sine_wave_advection", "solid_rotation_advection"}:
        raise NotImplementedError(f"Unsupported postprocess problem: {config.problem}")
    from postprocess_case import postprocess_case

    postprocess_case(config.case_dir(resolution), PROJECT_ROOT, config.max_co)


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

    if not prepare_only and config.problem == "sine_wave_advection":
        from study_analysis import analyse, collect, plot

        collect(config.solver_family, config.case_name, resolutions)
        analyse(config.solver_family, config.case_name)
        plot(config.solver_family, config.case_name)
