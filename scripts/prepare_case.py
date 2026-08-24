#!/usr/bin/env python3

"""Prepare one OpenFOAM case from a JSON configuration and one N value."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


COMMON_DIR = Path(__file__).resolve().parent / "common"
sys.path.insert(0, str(COMMON_DIR))

from case_config import load_config
from foam_case import prepare_case


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="JSON config path or config stem")
    parser.add_argument("--N", type=int, required=True, help="cells per edge")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--refresh-initial-only",
        action="store_true",
        help="only regenerate 0.orig/T for an already prepared case",
    )
    args = parser.parse_args()

    config = load_config(args.config)
    prepare_case(
        config,
        args.N,
        overwrite=args.overwrite,
        refresh_initial_only=args.refresh_initial_only,
    )


if __name__ == "__main__":
    main()
