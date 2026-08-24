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

from advection_sine import write_case_initial_field as write_sine_initial_field
from case_config import CaseConfig
from mesh_tools import mesh_resolution, patch_block_mesh_resolution
from paths import PROJECT_ROOT


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
        if re.search(r"^\s*grad\(T\)\s+[^;]+;", updated, re.M):
            updated = re.sub(
                r"^(\s*)grad\(T\)\s+[^;]+;",
                grad_line,
                updated,
                count=1,
                flags=re.M,
            )
        else:
            def add_grad_t(match: re.Match[str]) -> str:
                return match.group(1).rstrip() + f"\n{grad_line}\n}}\n"

            updated, grad_count = re.subn(
                r"(gradSchemes\s*\{.*?^\})",
                add_grad_t,
                updated,
                count=1,
                flags=re.S | re.M,
            )
            if grad_count != 1:
                raise RuntimeError(f"Cannot find gradSchemes in {path}")
    path.write_text(updated, encoding="utf-8")


def _patch_control_dict(case: Path, config: CaseConfig) -> None:
    path = case / "system" / "controlDict"
    _replace_or_append_dictionary_entry(path, "application", config.solver)
    _replace_or_append_dictionary_entry(path, "endTime", f"{config.end_time:g}")
    _replace_or_append_dictionary_entry(path, "maxCo", f"{config.max_co:g}")


def _write_initial_field(case: Path, config: CaseConfig) -> Path:
    nx, ny = mesh_resolution(case)
    if config.problem == "sine_wave_advection":
        return write_sine_initial_field(case, nx, ny)
    raise NotImplementedError(f"Unsupported problem for initial field: {config.problem}")


def _config_reference(config: CaseConfig) -> str:
    try:
        relative = config.path.relative_to(PROJECT_ROOT)
        return f"$projectRoot/{relative.as_posix()}"
    except ValueError:
        return config.path.as_posix()


def _write_case_scripts(case: Path, config: CaseConfig, resolution: int) -> None:
    config_reference = _config_reference(config)
    solver_path = f"$projectRoot/build/bin/{config.solver}"
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
    allrun = f"""#!/bin/sh

set -eu

caseDir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
projectRoot=$(CDPATH= cd -- "$caseDir/../../.." && pwd)
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
    if config.mesh_type != "quad":
        raise NotImplementedError(f"Only quad mesh is implemented now: {config.mesh_type}")

    target = config.case_dir(resolution)
    template = config.template_case
    if refresh_initial_only:
        if not target.exists():
            raise RuntimeError(f"Case directory does not exist: {target}")
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

    patch_block_mesh_resolution(target, resolution)
    _patch_fv_schemes(target, config)
    _patch_control_dict(target, config)
    initial_field = _write_initial_field(target, config)
    _write_case_scripts(target, config, resolution)

    metadata = {
        "caseName": config.case_name,
        "config": str(config.path),
        "problem": config.problem,
        "meshType": config.mesh_type,
        "schemeName": config.scheme_name,
        "divScheme": config.div_scheme,
        "resolution": resolution,
        "solver": config.solver,
        "velocity": list(config.velocity),
        "endTime": config.end_time,
        "maxCo": config.max_co,
        "initialField": str(initial_field),
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
    if config.problem != "sine_wave_advection":
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

    if not prepare_only:
        from study_analysis import analyse, collect, plot

        collect(config.case_name, resolutions)
        analyse(config.case_name)
        plot(config.case_name)
