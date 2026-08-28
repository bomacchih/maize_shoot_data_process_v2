#!/usr/bin/env python3
"""Check the active Python environment against requirements-python.txt."""

from __future__ import annotations

import importlib.metadata
import re
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
REQUIREMENTS = PROJECT_ROOT / "data" / "requirements-python.txt"


def parse_requirement(line: str) -> tuple[str, str | None]:
    line = line.split("#", 1)[0].strip()
    if not line:
        return "", None
    match = re.fullmatch(
        r"([A-Za-z0-9_.-]+)(?:\[.*?\])?\s*(?:==\s*([^\s;]+))?(?:\s*;.*)?",
        line,
    )
    if not match:
        raise ValueError(f"Unsupported requirement syntax: {line}")
    return match.group(1), match.group(2)


def main() -> int:
    failures: list[str] = []
    warnings: list[str] = []
    passes: list[str] = []

    for raw_line in REQUIREMENTS.read_text(encoding="utf-8-sig").splitlines():
        name, pinned = parse_requirement(raw_line)
        if not name:
            continue
        try:
            observed = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            failures.append(f"{name} is not installed")
            continue
        if pinned is None:
            warnings.append(f"{name} {observed} is installed but not pinned")
        elif observed != pinned:
            failures.append(f"{name} is {observed}; required version is {pinned}")
        else:
            passes.append(f"{name} {observed}")

    print(f"[INFO] Python {sys.version.split()[0]}")
    for message in passes:
        print(f"[PASS] {message}")
    for message in warnings:
        print(f"[WARN] {message}")
    for message in failures:
        print(f"[FAIL] {message}")
    print(
        f"\nSummary: {len(passes)} pinned match(es), "
        f"{len(warnings)} unpinned installed package(s), "
        f"{len(failures)} missing/mismatched package(s)."
    )
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())

