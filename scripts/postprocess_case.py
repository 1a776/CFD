#!/usr/bin/env python3

"""Post-process one configured case resolution without running OpenFOAM."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


COMMON_DIR = Path(__file__).resolve().parent / "common"
sys.path.insert(0, str(COMMON_DIR))

from case_config import load_config
from foam_case import postprocess_configured_case


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="JSON config path or config stem")
    parser.add_argument("--N", type=int, required=True, help="cells per edge")
    args = parser.parse_args()

    config = load_config(args.config)
    postprocess_configured_case(config, args.N)


if __name__ == "__main__":
    main()
