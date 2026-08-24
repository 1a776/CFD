#!/usr/bin/env python3

"""Compute convergence tables for one configured study."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


COMMON_DIR = Path(__file__).resolve().parent / "common"
sys.path.insert(0, str(COMMON_DIR))

from case_config import load_config


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="JSON config path or config stem")
    args = parser.parse_args()

    config = load_config(args.config)
    from study_analysis import analyse

    analyse(config.case_name)


if __name__ == "__main__":
    main()
