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


def solver_cases_dir(solver_family: str) -> Path:
    return CASES_DIR / solver_family


def solver_case_dir(solver_family: str, case_name: str, resolution: int) -> Path:
    return solver_cases_dir(solver_family) / case_name / f"N{resolution}"


def solver_data_dir(
    solver_family: str, case_name: str, resolution: int | None = None
) -> Path:
    root = DATA_DIR / solver_family / "cases" / case_name
    return root if resolution is None else root / f"N{resolution}"


def solver_analysis_dir(solver_family: str, case_name: str) -> Path:
    return DATA_DIR / solver_family / "analysis" / case_name


def solver_figure_dir(
    solver_family: str, case_name: str, resolution: int | None = None
) -> Path:
    root = FIGURES_DIR / solver_family / "cases" / case_name
    return root if resolution is None else root / f"N{resolution}"


def solver_analysis_figure_dir(solver_family: str, case_name: str) -> Path:
    return FIGURES_DIR / solver_family / "analysis" / case_name


def solver_build_dir(solver_family: str) -> Path:
    return BUILD_DIR / solver_family


def project_path(*parts: str) -> Path:
    """Build an absolute path under the project root."""
    return PROJECT_ROOT.joinpath(*parts)
