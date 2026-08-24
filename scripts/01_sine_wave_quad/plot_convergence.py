#!/usr/bin/env python3

"""Rebuild the 01 case convergence tables and figures."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


COMMON_DIR = Path(__file__).resolve().parents[1] / "common"
sys.path.insert(0, str(COMMON_DIR))

from study_analysis import analyse, plot


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-name", default="01_sine_wave_quad")
    args = parser.parse_args()
    analyse(args.case_name)
    plot(args.case_name)


if __name__ == "__main__":
    main()
