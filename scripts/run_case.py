#!/usr/bin/env python3

"""Prepare, run, and optionally post-process one configured OpenFOAM case."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path


COMMON_DIR = Path(__file__).resolve().parent / "common"
sys.path.insert(0, str(COMMON_DIR))

from case_config import load_config
from foam_case import postprocess_configured_case, prepare_case, run_case


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="JSON config path or config stem")
    parser.add_argument("--N", type=int, required=True, help="cells per edge")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--no-prepare", action="store_true")
    parser.add_argument("--no-postprocess", action="store_true")
    parser.add_argument(
        "--bashrc",
        type=Path,
        default=Path(os.environ.get("OPENFOAM_BASHRC", "/opt/openfoam14/etc/bashrc")),
    )
    args = parser.parse_args()

    config = load_config(args.config)
    case = config.case_dir(args.N)
    if not args.no_prepare:
        case = prepare_case(config, args.N, overwrite=args.overwrite)

    run_case(case, args.bashrc.resolve())
    if not args.no_postprocess:
        postprocess_configured_case(config, args.N)


if __name__ == "__main__":
    main()
