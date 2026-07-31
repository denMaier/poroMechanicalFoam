#!/usr/bin/env python3
"""Plot the scale-2 stabilization-operator and solid-work pilot results."""

from __future__ import annotations

import csv
import pathlib
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from plot_matched_scale import rate_summary
from tau_fo_sweep import read_internal_field


DECOMPOSITION = re.compile(
    r"Pressure increment decomposition: iteration (\d+), "
    r"totalRms ([^,]+), checkerboardRms ([^,]+), "
    r"complementRms ([^,]+), checkerboardFraction ([^\s]+)"
)
SOLID_ITERATIONS = re.compile(
    r"^\s*(\d+)\s*: Converged - "
    r"(?:Relative residual|Absolute residual|Step norm relative) "
    r"tolerance met\.",
    re.MULTILINE,
)


def vector_field(path: pathlib.Path) -> np.ndarray:
    text = path.read_text()
    match = re.search(
        r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;",
        text,
        re.DOTALL,
    )
    if not match:
        raise RuntimeError(f"cannot parse vector field {path}")
    values = np.asarray(
        [
            [float(component) for component in vector.split()]
            for vector in re.findall(r"\(([^()]*)\)", match.group(2))
        ]
    )
    if len(values) != int(match.group(1)):
        raise RuntimeError(f"incorrect vector count in {path}")
    return values


def operator_row(directory: pathlib.Path, operator: str, factor: float) -> dict[str, float | str]:
    log = (directory / "tau_00_Fo_00" / "log.poroMechanicalFoam").read_text(
        errors="replace"
    )
    values = DECOMPOSITION.findall(log)
    complement = rate_summary(
        [(int(value[0]), float(value[3])) for value in values]
    )
    checkerboard = rate_summary(
        [(int(value[0]), float(value[2])) for value in values],
        noise_relative=1.0e-8,
    )
    summary = next(csv.DictReader((directory / "tauFoSweep.csv").open()))
    return {
        "kind": "operator",
        "operator": operator,
        "factor": factor,
        "complementRate": float(complement["rate"]),
        "complementLow": float(complement["low"]),
        "complementHigh": float(complement["high"]),
        "checkerboardRate": float(checkerboard["rate"]),
        "checkerboardToSource": float(summary["checkerboardToSource"]),
    }


def solid_row(
    directory: pathlib.Path,
    label: str,
    rtol: float,
    reference_pressure: np.ndarray,
    reference_displacement: np.ndarray,
) -> dict[str, float | str]:
    case = directory / "tau_00_Fo_00"
    log = (case / "log.poroMechanicalFoam").read_text(errors="replace")
    summary = next(csv.DictReader((directory / "tauFoSweep.csv").open()))
    pressure = read_internal_field(case / "10" / "poroFluid" / "p_rgh")
    displacement = vector_field(case / "10" / "solid" / "D")
    solid_iterations = [int(value) for value in SOLID_ITERATIONS.findall(log)]
    return {
        "kind": "solidWork",
        "label": label,
        "solidRtol": rtol,
        "outerIterations": int(summary["iterations"]),
        "solidCorrectors": sum(solid_iterations),
        "relativePressureDifference": np.linalg.norm(
            pressure - reference_pressure
        )
        / np.linalg.norm(reference_pressure),
        "relativeDisplacementDifference": np.linalg.norm(
            displacement - reference_displacement
        )
        / np.linalg.norm(reference_displacement),
    }


def main() -> None:
    root = pathlib.Path(__file__).resolve().parent
    output = root / "scale2DiagnosticPilots"
    output.mkdir(exist_ok=True)

    operator_rows = []
    for factor in (0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.5):
        token = f"{factor:g}".replace(".", "p")
        operator_rows.append(
            operator_row(
                root / f"scale2StabFactor{token}Pilot",
                "compact-minus-wide",
                factor,
            )
        )
    for factor in (0.08, 0.16, 0.32, 1.0):
        token = str(factor).replace(".", "p")
        operator_rows.append(
            operator_row(
                root / f"scale2PressureJumpFactor{token}Pilot",
                "pressure jump",
                factor,
            )
        )

    exact_directory = root / "scale2ConvergedSolidExact"
    exact_case = exact_directory / "tau_00_Fo_00"
    reference_pressure = read_internal_field(exact_case / "10" / "poroFluid" / "p_rgh")
    reference_displacement = vector_field(exact_case / "10" / "solid" / "D")
    solid_rows = [
        solid_row(
            exact_directory,
            r"$10^{-6}$",
            1.0e-6,
            reference_pressure,
            reference_displacement,
        ),
        solid_row(
            root / "scale2ConvergedSolidR1e3",
            r"$10^{-3}$",
            1.0e-3,
            reference_pressure,
            reference_displacement,
        ),
        solid_row(
            root / "scale2ConvergedSolidR1e2",
            r"$10^{-2}$",
            1.0e-2,
            reference_pressure,
            reference_displacement,
        ),
    ]

    fields = sorted({key for row in operator_rows + solid_rows for key in row})
    with (output / "pilotMetrics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(operator_rows + solid_rows)

    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.4), constrained_layout=True)
    rate_axis = axes[0]
    styles = {
        "compact-minus-wide": ("tab:blue", "o"),
        "pressure jump": ("tab:orange", "s"),
    }
    for operator, (color, marker) in styles.items():
        selected = sorted(
            (row for row in operator_rows if row["operator"] == operator),
            key=lambda row: float(row["factor"]),
        )
        factor = np.asarray([float(row["factor"]) for row in selected])
        complement = np.asarray([float(row["complementRate"]) for row in selected])
        checkerboard = np.asarray([float(row["checkerboardRate"]) for row in selected])
        rate_axis.plot(
            factor,
            complement,
            color=color,
            marker=marker,
            label=f"{operator}: complement",
        )
        rate_axis.plot(
            factor,
            checkerboard,
            color=color,
            marker=marker,
            linestyle="--",
            label=f"{operator}: checkerboard",
        )
    rate_axis.axhline(1.0, color="black", linewidth=1.2)
    rate_axis.axvline(
        9.0 * (11111.111111 / 19444.444444) / 32.0,
        color="0.4",
        linestyle=":",
        label=r"Aronson $c=1$ mapping",
    )
    rate_axis.set_xscale("log")
    rate_axis.set_ylim(0.6, 1.25)
    rate_axis.set_xlabel("stabilization factor")
    rate_axis.set_ylabel("measured outer rate")
    rate_axis.set_title(r"Hard point: $\tau=10^6$, $\mathrm{Fo}=10^{-5}$")
    rate_axis.grid(color="0.9", linewidth=0.5)
    rate_axis.legend(fontsize=8)

    work_axis = axes[1]
    positions = np.arange(len(solid_rows))
    correctors = np.asarray([float(row["solidCorrectors"]) for row in solid_rows])
    outer = np.asarray([float(row["outerIterations"]) for row in solid_rows])
    bars = work_axis.bar(positions, correctors, color="tab:blue", alpha=0.78)
    work_axis.set_xticks(positions, [str(row["label"]) for row in solid_rows])
    work_axis.set_xlabel("solid relative tolerance")
    work_axis.set_ylabel("total solid correctors", color="tab:blue")
    work_axis.tick_params(axis="y", labelcolor="tab:blue")
    work_axis.set_title(r"Converged point: $\tau=1$, $\mathrm{Fo}=0.1$")
    outer_axis = work_axis.twinx()
    outer_axis.plot(positions, outer, color="tab:red", marker="o", linewidth=2)
    outer_axis.set_ylabel("outer iterations", color="tab:red")
    outer_axis.tick_params(axis="y", labelcolor="tab:red")
    for bar, row in zip(bars, solid_rows):
        work_axis.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"{int(row['solidCorrectors'])}",
            ha="center",
            va="bottom",
            fontsize=9,
        )
    work_axis.text(
        0.02,
        0.97,
        "loose-solve field errors: $<7\\times10^{-11}$ relative",
        transform=work_axis.transAxes,
        va="top",
        fontsize=9,
    )
    work_axis.grid(axis="y", color="0.9", linewidth=0.5)

    fig.suptitle("Fixed-stress scale 2: stabilization dynamics and inner-solve work")
    path = output / "operatorAndSolidPilots.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    print(path)
    print(output / "pilotMetrics.csv")


if __name__ == "__main__":
    main()
