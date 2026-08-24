#!/usr/bin/env python3

"""Generate the sine-wave initial field for one 02 case resolution."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path


COMMON_DIR = Path(__file__).resolve().parents[1] / "common"
sys.path.insert(0, str(COMMON_DIR))

from advection_tools import mesh_resolution, write_initial_field


PROJECT_ROOT = Path(__file__).resolve().parents[2]


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
    output = write_initial_field(case, nx, ny)
    print(f"case={case}")
    print(f"resolution={nx}x{ny}")
    print(f"initialField={output}")


if __name__ == "__main__":
    main()
