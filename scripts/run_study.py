#!/usr/bin/env python3

"""Run a configured multi-resolution study."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path


COMMON_DIR = Path(__file__).resolve().parent / "common"
sys.path.insert(0, str(COMMON_DIR))

from case_config import load_config, parse_resolutions
from foam_case import run_study


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="JSON config path or config stem")
    parser.add_argument(
        "--resolutions",
        type=parse_resolutions,
        default=None,
        help="comma-separated N values; defaults to the config resolutions",
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument(
        "--bashrc",
        type=Path,
        default=Path(os.environ.get("OPENFOAM_BASHRC", "/opt/openfoam14/etc/bashrc")),
    )
    args = parser.parse_args()

    config = load_config(args.config)
    resolutions = args.resolutions if args.resolutions is not None else list(config.resolutions)
    run_study(
        config,
        resolutions=resolutions,
        overwrite=args.overwrite,
        prepare_only=args.prepare_only,
        bashrc=args.bashrc.resolve(),
    )


if __name__ == "__main__":
    main()
