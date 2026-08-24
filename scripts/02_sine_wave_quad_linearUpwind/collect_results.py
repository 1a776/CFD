#!/usr/bin/env python3

"""Collect per-N JSON summaries for the 02 case family."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


COMMON_DIR = Path(__file__).resolve().parents[1] / "common"
sys.path.insert(0, str(COMMON_DIR))

from run_suite import parse_resolutions
from study_analysis import collect


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--resolutions", type=parse_resolutions, default=[10, 20, 40, 80])
    args = parser.parse_args()
    collect("02_sine_wave_quad_linearUpwind", args.resolutions)


if __name__ == "__main__":
    main()
