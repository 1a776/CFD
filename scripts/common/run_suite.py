#!/usr/bin/env python3

"""Generate and run N10/N20/N40/N80 cases for one case family."""

from __future__ import annotations

import argparse
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

from advection_tools import write_initial_field
from study_analysis import analyse, collect, plot


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RESOLUTIONS = (10, 20, 40, 80)
BLOCK_PATTERN = re.compile(
    r"(hex\s+\(0\s+1\s+2\s+3\s+4\s+5\s+6\s+7\)\s+)"
    r"\(\s*\d+\s+\d+\s+1\s*\)"
)


def parse_resolutions(value: str) -> list[int]:
    values = sorted({int(item.strip()) for item in value.split(",") if item.strip()})
    if not values or any(value <= 0 for value in values):
        raise argparse.ArgumentTypeError("resolutions must be positive integers")
    return values


def patch_block_mesh(case: Path, resolution: int) -> None:
    path = case / "system" / "blockMeshDict"
    text = path.read_text(encoding="utf-8")
    updated, count = BLOCK_PATTERN.subn(rf"\g<1>({resolution} {resolution} 1)", text, count=1)
    if count != 1:
        raise RuntimeError(f"Cannot patch blockMeshDict for {case}")
    path.write_text(updated, encoding="utf-8")


def copy_constant_inputs(template: Path, target: Path) -> None:
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


def write_case_scripts(case: Path, case_name: str) -> None:
    case_dir = case.as_posix()
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

solverPath="$projectRoot/build/bin/explicitAdvectionFoamStudent"
if [ ! -x "$solverPath" ]; then
    echo "Missing solver: $solverPath" >&2
    echo "Run: sh $projectRoot/scripts/build_student_solver.sh" >&2
    exit 1
fi

sh "$caseDir/Allclean"
python3 "$projectRoot/scripts/{case_name}/create_initial_fields.py" --case-dir "$caseDir"
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
    case_name: str,
    template: Path,
    resolution: int,
    overwrite: bool,
) -> Path:
    case_root = PROJECT_ROOT / "cases" / case_name
    target = case_root / f"N{resolution}"
    template = template.resolve()

    if target == template:
        if not target.exists():
            raise RuntimeError(f"Template case does not exist: {target}")
        if overwrite:
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
        copy_constant_inputs(template, target)

    patch_block_mesh(target, resolution)
    nx = resolution
    ny = resolution
    write_initial_field(target, nx, ny)
    write_case_scripts(target, case_name)
    (target / "metadata.json").write_text(
        "{\n"
        f'  "caseName": "{case_name}",\n'
        f'  "resolution": {resolution},\n'
        '  "initialFunction": "sin(2*pi*(x+y))",\n'
        '  "velocity": [1.0, 1.0],\n'
        '  "endTime": 1.0,\n'
        '  "maxCo": 0.2\n'
        "}\n",
        encoding="utf-8",
    )
    return target


def run_case(case: Path, bashrc: Path) -> None:
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


def run_suite(
    case_name: str,
    resolutions: list[int],
    overwrite: bool,
    prepare_only: bool,
    bashrc: Path,
) -> None:
    case_root = PROJECT_ROOT / "cases" / case_name
    template = case_root / "N20"
    if not template.exists():
        raise RuntimeError(f"The N20 template case is missing: {template}")

    print(f"caseFamily={case_name}")
    print(f"resolutions={','.join(str(value) for value in resolutions)}")
    for resolution in resolutions:
        case = prepare_case(case_name, template, resolution, overwrite)
        print(f"prepared={case}")
        if prepare_only:
            continue
        print(f"running={case}")
        run_case(case, bashrc)
        postprocess = PROJECT_ROOT / "scripts" / case_name / "plot_results.py"
        subprocess.run(
            [sys.executable, str(postprocess), "--case-dir", str(case)],
            cwd=PROJECT_ROOT,
            check=True,
        )

    if not prepare_only:
        collect(case_name, resolutions)
        analyse(case_name)
        plot(case_name)


def main(default_case_name: str) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--resolutions",
        type=parse_resolutions,
        default=list(DEFAULT_RESOLUTIONS),
        help="comma-separated N values, default: 10,20,40,80",
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument(
        "--bashrc",
        type=Path,
        default=Path(os.environ.get("OPENFOAM_BASHRC", "/opt/openfoam14/etc/bashrc")),
    )
    args = parser.parse_args()
    run_suite(
        default_case_name,
        args.resolutions,
        args.overwrite,
        args.prepare_only,
        args.bashrc.resolve(),
    )


if __name__ == "__main__":
    main("01_sine_wave_quad")
