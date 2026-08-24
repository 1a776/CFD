#!/usr/bin/env python3

"""Shared project paths for the student CFD scripts."""

from __future__ import annotations

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
CONFIG_DIR = SCRIPTS_DIR / "configs"
CASES_DIR = PROJECT_ROOT / "cases"
DATA_DIR = PROJECT_ROOT / "data"
FIGURES_DIR = PROJECT_ROOT / "figures"
BUILD_DIR = PROJECT_ROOT / "build"


def project_path(*parts: str) -> Path:
    """Build an absolute path under the project root."""
    return PROJECT_ROOT.joinpath(*parts)
