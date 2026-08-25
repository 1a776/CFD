#!/usr/bin/env python3

"""Collect per-resolution summary.json files for one configured study."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


COMMON_DIR = Path(__file__).resolve().parent / "common"
sys.path.insert(0, str(COMMON_DIR))

from case_config import load_config, parse_resolutions


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="JSON config path or config stem")
    parser.add_argument(
        "--resolutions",
        type=parse_resolutions,
        default=None,
        help="comma-separated N values; defaults to the config resolutions",
    )
    args = parser.parse_args()

    config = load_config(args.config)
    resolutions = args.resolutions if args.resolutions is not None else list(config.resolutions)
    from study_analysis import collect

    collect(config.solver_family, config.case_name, resolutions)


if __name__ == "__main__":
    main()
