#!/usr/bin/env python3

"""Print final errors for one 02 case resolution."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path


COMMON_DIR = Path(__file__).resolve().parents[1] / "common"
sys.path.insert(0, str(COMMON_DIR))

from advection_tools import (
    exact_values,
    latest_time,
    mesh_resolution,
    normalized_errors,
    read_scalar_field,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]
THICKNESS = 0.1


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-dir", type=Path, default=None)
    args = parser.parse_args()
    default_case = PROJECT_ROOT / "cases" / "02_sine_wave_quad_linearUpwind" / "N20"
    case = (
        args.case_dir.resolve()
        if args.case_dir is not None
        else Path(os.environ.get("VIBEFLOW_CASE_DIR", default_case))
        .expanduser()
        .resolve()
    )

    nx, ny = mesh_resolution(case)
    final_time, final_dir = latest_time(case)
    numerical = read_scalar_field(final_dir / "T")
    exact = exact_values(nx, ny, final_time)
    l1, l2, linf = normalized_errors(
        numerical, exact, THICKNESS / (nx * ny)
    )

    print(f"case={case}")
    print(f"finalTime={final_time:.16g}")
    print(f"normalizedL1={l1:.16e}")
    print(f"normalizedL2={l2:.16e}")
    print(f"normalizedLinf={linf:.16e}")
    print(f"minT={min(numerical):.16e}")
    print(f"maxT={max(numerical):.16e}")


if __name__ == "__main__":
    main()
