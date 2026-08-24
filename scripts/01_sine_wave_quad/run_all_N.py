#!/usr/bin/env python3

"""Run N10/N20/N40/N80 for the first-order upwind case."""

from __future__ import annotations

import sys
from pathlib import Path


COMMON_DIR = Path(__file__).resolve().parents[1] / "common"
sys.path.insert(0, str(COMMON_DIR))

from run_suite import main


if __name__ == "__main__":
    main("01_sine_wave_quad")
