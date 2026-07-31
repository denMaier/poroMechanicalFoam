#!/usr/bin/env python3
"""Measure the Cartesian checkerboard subspace on the 10-by-10 benchmark."""

from __future__ import annotations

import math
import pathlib
import re
import sys


def read_internal_field(path: pathlib.Path) -> list[float]:
    text = path.read_text()
    uniform = re.search(r"internalField\s+uniform\s+([^;]+);", text)
    if uniform:
        return [float(uniform.group(1))] * 100
    match = re.search(
        r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;",
        text,
        flags=re.DOTALL,
    )
    if not match:
        raise RuntimeError(f"cannot parse nonuniform scalar field {path}")
    count = int(match.group(1))
    values = [float(value) for value in match.group(2).split()]
    if len(values) != count:
        raise RuntimeError(f"expected {count} values in {path}, found {len(values)}")
    return values


def checkerboard_metric(values: list[float], nx: int = 10, ny: int = 10) -> float:
    if len(values) != nx * ny:
        raise RuntimeError(f"expected {nx*ny} cells, found {len(values)}")
    projections = [0.0, 0.0, 0.0]
    energy = 0.0
    for cell, value in enumerate(values):
        i = cell % nx
        j = cell // nx
        projections[0] += (-1.0 if i % 2 else 1.0) * value
        projections[1] += (-1.0 if j % 2 else 1.0) * value
        projections[2] += (-1.0 if (i + j) % 2 else 1.0) * value
        energy += value * value
    projected_energy = sum(projection * projection for projection in projections)
    return math.sqrt(projected_energy) / (
        math.sqrt(nx * ny) * math.sqrt(energy) + 1e-300
    )


def relative_l2(first: list[float], second: list[float]) -> float:
    numerator = math.sqrt(sum((a - b) ** 2 for a, b in zip(first, second)))
    denominator = math.sqrt(sum(a * a for a in first))
    return numerator / (denominator + 1e-300)


def main(paths: list[str]) -> None:
    if len(paths) != 3:
        raise SystemExit("usage: checkerboard_metric.py UNSTABILIZED EXPLICIT IMPLICIT")
    fields = [read_internal_field(pathlib.Path(path)) for path in paths]
    names = ("unstabilized", "explicit", "implicit")
    print("variant,checkerboardProjection,relativeL2ToImplicit")
    for name, field in zip(names, fields):
        print(
            f"{name},{checkerboard_metric(field):.12g},"
            f"{relative_l2(field, fields[2]):.12g}"
        )


if __name__ == "__main__":
    main(sys.argv[1:])
