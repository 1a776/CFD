#!/usr/bin/env python3

"""Create data and figures for one 01 case resolution."""

from __future__ import annotations

import sys
from pathlib import Path


COMMON_DIR = Path(__file__).resolve().parents[1] / "common"
sys.path.insert(0, str(COMMON_DIR))

from postprocess_case import main


if __name__ == "__main__":
    main("01_sine_wave_quad")
