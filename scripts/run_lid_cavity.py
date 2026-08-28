#!/usr/bin/env python3

"""Prepare and run exactly one fifth-problem lid-driven cavity case."""

from __future__ import annotations

import argparse
from pathlib import Path

from common.lid_cavity import run


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    target = run(args.config.resolve(), overwrite=args.overwrite)
    print(f"completed={target}")


if __name__ == "__main__":
    main()
