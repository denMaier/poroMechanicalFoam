#!/usr/bin/env python3
"""Generate modified Barry--Mercer tau--Fo surfaces."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import math
import os
import pathlib
import re
import shutil
import subprocess
from dataclasses import dataclass

import matplotlib

matplotlib.use("Agg")
from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np


DEFAULT_TAU = np.logspace(-2, 6, 9)
DEFAULT_FO = np.logspace(-6, 2, 9)
P_TOLERANCE = 1.0e-6
D_TOLERANCE = 1.0e-10


@dataclass(frozen=True)
class SweepPoint:
    tau_index: int
    fo_index: int
    tau: float
    fo: float
    storage: float
    conductivity: float
    case_name: str


def parse_values(text: str | None, default: np.ndarray) -> np.ndarray:
    if not text:
        return default
    values = [float(value) for value in text.replace(",", " ").split()]
    if not values or any(value <= 0 for value in values):
        raise ValueError("tau and Fo values must be positive")
    return np.asarray(values, dtype=float)


def run_checked(command: list[str], stdout=None) -> None:
    subprocess.run(command, check=True, stdout=stdout, stderr=subprocess.STDOUT)


def set_dictionary_entry(path: pathlib.Path, entry: str, value: str) -> None:
    run_checked(
        ["foamDictionary", str(path), "-entry", entry, "-set", value],
        stdout=subprocess.DEVNULL,
    )


def read_internal_field(path: pathlib.Path, count: int = 100) -> np.ndarray:
    text = path.read_text()
    uniform = re.search(r"internalField\s+uniform\s+([^;]+);", text)
    if uniform:
        return np.full(count, float(uniform.group(1)))
    match = re.search(
        r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;",
        text,
        flags=re.DOTALL,
    )
    if not match:
        raise RuntimeError(f"cannot parse pressure field {path}")
    values = np.asarray([float(value) for value in match.group(2).split()])
    if len(values) != int(match.group(1)):
        raise RuntimeError(f"incorrect value count in {path}")
    return values


def coupling_history(log_text: str) -> list[tuple[int, float, float]]:
    pattern = re.compile(
        r"\|\s*#\s+Iterations\s*\|\s*(\d+)\s*\|[^\n]*\n"
        r"\|\s*max\s+delta\(p_rgh\)\s*\|\s*([^\s|]+)\s*\|[^\n]*\n"
        r"\|\s*max\s+delta\(D\)\s*\|\s*([^\s|]+)",
    )
    return [
        (int(iteration), float(p_residual), float(d_residual))
        for iteration, p_residual, d_residual in pattern.findall(log_text)
    ]


def pressure_mode_history(log_text: str) -> list[tuple[int, float]]:
    pattern = re.compile(
        r"Pressure mode subspace diagnostic: iteration (\d+), "
        r"incrementNorm ([^\s]+)"
    )
    return [
        (int(iteration), float(increment_norm))
        for iteration, increment_norm in pattern.findall(log_text)
    ]


def scaled_update_history(
    history: list[tuple[int, float, float]],
) -> list[tuple[int, float]]:
    return [
        (
            iteration,
            max(p_residual / P_TOLERANCE, d_residual / D_TOLERANCE),
        )
        for iteration, p_residual, d_residual in history
        if iteration > 1
        and math.isfinite(p_residual)
        and math.isfinite(d_residual)
        and (p_residual > 0 or d_residual > 0)
    ]


def contraction_factor(history: list[tuple[int, float, float]]) -> float:
    usable = scaled_update_history(history)
    if len(usable) < 3:
        return math.nan

    # Fixed-stress updates can alternate between several spatial modes.  A
    # pointwise log fit aliases this beating when the spectral radius is near
    # one. Compare RMS envelopes over two adjacent windows instead.
    window = min(10, len(usable) // 2)
    earlier = usable[-2 * window : -window]
    later = usable[-window:]
    earlier_rms = math.sqrt(sum(value * value for _, value in earlier) / window)
    later_rms = math.sqrt(sum(value * value for _, value in later) / window)
    iteration_span = (
        sum(iteration for iteration, _ in later) / window
        - sum(iteration for iteration, _ in earlier) / window
    )
    return float((later_rms / earlier_rms) ** (1.0 / iteration_span))


def monotone_contraction(history: list[tuple[int, float, float]]) -> bool:
    """Return whether the scaled coupled update decreases late in the solve."""
    usable = scaled_update_history(history)
    if len(usable) < 3:
        return False

    late = usable[-min(20, len(usable)) :]
    return all(
        current < previous
        for (_, previous), (_, current) in zip(late, late[1:])
    )


def mode_contraction_factor(history: list[tuple[int, float]]) -> float:
    usable = [
        (iteration, increment_norm)
        for iteration, increment_norm in history
        if math.isfinite(increment_norm) and increment_norm > 0
    ]
    if len(usable) < 3:
        return math.nan

    # Once stabilization has removed the mode, its projected increment is
    # round-off around a constant final amplitude. Stop at that noise floor
    # instead of fitting the meaningless late sign changes.
    cutoff = usable[0][1] * 1.0e-8
    noise_start = next(
        (
            index
            for index, (_, value) in enumerate(usable[5:], start=5)
            if value < cutoff
        ),
        len(usable),
    )
    usable = usable[:noise_start]
    if len(usable) < 3:
        return math.nan

    # For an affine stationary iteration, projected increments obey
    # delta(a)_{k+1} = lambda*delta(a)_k. RMS blocks recover |lambda|
    # without being affected by a sign-alternating checkerboard mode.
    window = min(10, len(usable) // 2)
    earlier = usable[-2 * window : -window]
    later = usable[-window:]
    earlier_rms = math.sqrt(sum(value * value for _, value in earlier) / window)
    later_rms = math.sqrt(sum(value * value for _, value in later) / window)
    iteration_span = (
        sum(iteration for iteration, _ in later) / window
        - sum(iteration for iteration, _ in earlier) / window
    )
    return float((later_rms / earlier_rms) ** (1.0 / iteration_span))


def prepare_point(
    base: pathlib.Path,
    runs: pathlib.Path,
    point: SweepPoint,
    max_iterations: int,
    formulation: str,
    fixed_stress_scale: float,
    stab_factor: float,
    solid_rtol: float | None = None,
    solid_stol: float | None = None,
    solid_correctors: int | None = None,
    stabilization_type: str = "implicit",
    excluded_cell_zones: tuple[str, ...] = (),
) -> pathlib.Path:
    target = runs / point.case_name
    shutil.copytree(base, target)
    (target / f"{point.case_name}.foam").touch()

    coupling = target / "constant" / "poroCouplingProperties"
    hydraulic = target / "constant" / "poroFluid" / "poroHydraulicProperties"
    control = target / "system" / "controlDict"
    solid_properties = target / "constant" / "solid" / "solidProperties"

    set_dictionary_entry(control, "endTime", "10")
    stabilized = formulation == "implicit"
    coeffs = "poroSolidStabCoeffs" if stabilized else "poroSolidCoeffs"
    set_dictionary_entry(
        coupling,
        "poroSolidInterface",
        "poroSolidStab" if stabilized else "poroSolid",
    )
    set_dictionary_entry(coupling, "fixedStressStabilScale", f"{fixed_stress_scale:.16g}")
    if stabilized:
        set_dictionary_entry(
            coupling, f"{coeffs}/stabilizationType", stabilization_type
        )
        set_dictionary_entry(coupling, f"{coeffs}/stabFactor", f"{stab_factor:.16g}")
        set_dictionary_entry(
            coupling,
            f"{coeffs}/excludedCellZones",
            f"( {' '.join(excluded_cell_zones)} )",
        )
    set_dictionary_entry(
        coupling,
        f"{coeffs}/iterations",
        str(max_iterations),
    )
    # Record every outer iteration: the robust stability classification below
    # deliberately rejects oscillatory histories, including oscillations that
    # would be hidden by an even reporting interval.
    set_dictionary_entry(coupling, f"{coeffs}/infoFrequency", "1")
    set_dictionary_entry(
        hydraulic,
        "storageCoeffCoeffs/Ss",
        f"Ss [-1 1 2 0 0 0 0] {point.storage:.16g}",
    )
    set_dictionary_entry(
        hydraulic,
        "k",
        f"k [0 1 -1 0 0 0 0] {point.conductivity:.16g}",
    )
    solid_coeffs = "linearGeometryTotalDisplacementCoeffs"
    if solid_rtol is not None:
        set_dictionary_entry(
            solid_properties, f"{solid_coeffs}/rTol", f"{solid_rtol:.16g}"
        )
    if solid_stol is not None:
        set_dictionary_entry(
            solid_properties, f"{solid_coeffs}/sTol", f"{solid_stol:.16g}"
        )
    if solid_correctors is not None:
        set_dictionary_entry(
            solid_properties, f"{solid_coeffs}/nCorrectors", str(solid_correctors)
        )
    return target


def run_point(
    base: pathlib.Path,
    runs: pathlib.Path,
    point: SweepPoint,
    max_iterations: int,
    formulation: str,
    fixed_stress_scale: float,
    stab_factor: float,
    solid_rtol: float | None = None,
    solid_stol: float | None = None,
    solid_correctors: int | None = None,
    stabilization_type: str = "implicit",
    excluded_cell_zones: tuple[str, ...] = (),
) -> tuple[str, int]:
    target = prepare_point(
        base,
        runs,
        point,
        max_iterations,
        formulation,
        fixed_stress_scale,
        stab_factor,
        solid_rtol,
        solid_stol,
        solid_correctors,
        stabilization_type,
        excluded_cell_zones,
    )
    log_path = target / "log.poroMechanicalFoam"
    with log_path.open("w") as log:
        completed = subprocess.run(
            ["poroMechanicalFoam", "-case", str(target)],
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    return point.case_name, completed.returncode


def analyse_point(
    runs: pathlib.Path,
    point: SweepPoint,
    return_code: int,
    max_iterations: int,
    delta_t: float,
    source_period: float,
    source_volume: float,
) -> dict[str, float | int | str]:
    target = runs / point.case_name
    log_text = (target / "log.poroMechanicalFoam").read_text(errors="replace")
    history = coupling_history(log_text)
    mode_history = pressure_mode_history(log_text)
    iteration = history[-1][0] if history else max_iterations
    p_residual = history[-1][1] if history else math.inf
    d_residual = history[-1][2] if history else math.inf
    converged = (
        return_code == 0
        and iteration < max_iterations
        and p_residual <= P_TOLERANCE
        and d_residual <= D_TOLERANCE
    )

    pressure_path = target / "10" / "poroFluid" / "p_rgh"
    checkerboard = math.nan
    checkerboard_ratio = math.nan
    checkerboard_components = np.full(3, math.nan)
    if pressure_path.exists():
        pressure = read_internal_field(pressure_path)
        modes = np.asarray(
            [
                [1.0 if cell % 2 == 0 else -1.0 for cell in range(100)],
                [1.0 if (cell // 10) % 2 == 0 else -1.0 for cell in range(100)],
                [
                    1.0 if ((cell % 10 + cell // 10) % 2 == 0) else -1.0
                    for cell in range(100)
                ],
            ]
        )
        checkerboard_components = modes @ pressure / len(pressure)
        checkerboard = float(np.linalg.norm(checkerboard_components))
        source_rate = math.sin(2.0 * math.pi * delta_t / source_period)
        effective_storage = point.storage + 1.0 / 5555.555555555556
        source_pressure = abs(source_rate) * delta_t / (source_volume * effective_storage)
        checkerboard_ratio = checkerboard / source_pressure

    rho = contraction_factor(history)
    checkerboard_rho = mode_contraction_factor(mode_history)
    monotone = monotone_contraction(history)
    robustly_contracting = math.isfinite(rho) and rho < 1.0 and monotone
    return {
        "case": point.case_name,
        "tau": point.tau,
        "Fo": point.fo,
        "Ss": point.storage,
        "conductivity": point.conductivity,
        "returnCode": return_code,
        "iterations": iteration,
        "pResidual": p_residual,
        "DResidual": d_residual,
        "contractionFactor": rho,
        "checkerboardContractionFactor": checkerboard_rho,
        "envelopeContracting": int(math.isfinite(rho) and rho < 1.0),
        "monotone": int(monotone),
        "contracting": int(robustly_contracting),
        "converged": int(converged),
        "checkerboardAmplitude": checkerboard,
        "checkerboardXAmplitude": checkerboard_components[0],
        "checkerboardYAmplitude": checkerboard_components[1],
        "checkerboardXYAmplitude": checkerboard_components[2],
        "checkerboardToSource": checkerboard_ratio,
    }


def logarithmic_edges(values: np.ndarray) -> np.ndarray:
    if len(values) == 1:
        half_decade = math.sqrt(10.0)
        return np.asarray([values[0] / half_decade, values[0] * half_decade])
    edges = np.empty(len(values) + 1)
    edges[1:-1] = np.sqrt(values[:-1] * values[1:])
    edges[0] = values[0] / math.sqrt(values[1] / values[0])
    edges[-1] = values[-1] * math.sqrt(values[-1] / values[-2])
    return edges


def plot_surface(
    rows: list[dict[str, float | int | str]],
    tau_values: np.ndarray,
    fo_values: np.ndarray,
    output: pathlib.Path,
    formulation: str,
    fixed_stress_scale: float,
    stab_factor: float,
) -> None:
    shape = (len(tau_values), len(fo_values))
    checkerboard = np.full(shape, np.nan)
    contraction = np.full(shape, np.nan)
    checkerboard_contraction = np.full(shape, np.nan)
    tau_index = {value: index for index, value in enumerate(tau_values)}
    fo_index = {value: index for index, value in enumerate(fo_values)}

    for row in rows:
        i = tau_index[float(row["tau"])]
        j = fo_index[float(row["Fo"])]
        checkerboard[i, j] = float(row["checkerboardToSource"])
        contraction[i, j] = float(row["contractionFactor"])
        checkerboard_contraction[i, j] = float(
            row["checkerboardContractionFactor"]
        )

    x_edges, y_edges = np.meshgrid(
        logarithmic_edges(fo_values), logarithmic_edges(tau_values)
    )
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.2), constrained_layout=True)
    if formulation == "implicit":
        fig.suptitle(
            "Implicit compact-stencil stabilization "
            f"(fixed-stress scale {fixed_stress_scale:g}, factor {stab_factor:g})"
        )
    else:
        fig.suptitle(f"Unstabilized fixed stress (scale {fixed_stress_scale:g})")

    checker_log = np.log10(np.maximum(checkerboard, 1.0e-12))
    first = axes[0].pcolormesh(
        x_edges,
        y_edges,
        checker_log,
        shading="flat",
        cmap="viridis",
        vmin=-6,
        vmax=1,
    )
    axes[0].set_title("Checkerboard amplitude / source pressure")
    colorbar = fig.colorbar(first, ax=axes[0])
    colorbar.set_label(r"$\log_{10}(|p_{cb}|/p_{source})$")

    rate_norm = colors.TwoSlopeNorm(vmin=0.7, vcenter=1.0, vmax=1.3)
    second = axes[1].pcolormesh(
        x_edges,
        y_edges,
        contraction,
        shading="flat",
        cmap="coolwarm",
        norm=rate_norm,
    )
    axes[1].set_title("Measured fixed-stress contraction")
    colorbar = fig.colorbar(second, ax=axes[1])
    colorbar.set_label(r"$\rho$")

    third = axes[2].pcolormesh(
        x_edges,
        y_edges,
        checkerboard_contraction,
        shading="flat",
        cmap="coolwarm",
        norm=rate_norm,
    )
    axes[2].set_title("Checkerboard-mode contraction")
    colorbar = fig.colorbar(third, ax=axes[2])
    colorbar.set_label(r"$\rho_{cb}$")

    for axis in axes:
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlabel(r"Fourier number $\mathrm{Fo}$")
        axis.set_ylabel(r"Coupling strength $\tau=M_{Biot}/K_{vol}$")
        axis.grid(False)

    fig.savefig(output, dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tau-values", help="space- or comma-separated positive values")
    parser.add_argument("--fo-values", help="space- or comma-separated positive values")
    parser.add_argument("--jobs", type=int, default=min(4, os.cpu_count() or 1))
    parser.add_argument("--max-iterations", type=int, default=300)
    parser.add_argument(
        "--formulation",
        choices=("unstabilized", "implicit"),
        default="unstabilized",
    )
    parser.add_argument("--fixed-stress-scale", type=float)
    parser.add_argument("--stab-factor", type=float, default=1.0)
    parser.add_argument(
        "--stabilization-type",
        choices=("implicit", "explicit", "pressureJump"),
        default="implicit",
    )
    parser.add_argument("--solid-rtol", type=float)
    parser.add_argument("--solid-stol", type=float)
    parser.add_argument("--solid-correctors", type=int)
    parser.add_argument(
        "--excluded-cell-zones",
        help=(
            "space- or comma-separated cell zones whose touching faces "
            "are unstabilized"
        ),
    )
    parser.add_argument(
        "--output-directory",
        help="result-directory name below the benchmark directory",
    )
    parser.add_argument("--clean", action="store_true")
    parser.add_argument(
        "--reuse-existing",
        action="store_true",
        help="reanalyze existing default-grid cases without running the solver",
    )
    args = parser.parse_args()

    if args.fixed_stress_scale is None:
        fixed_stress_scale = 5.5 if args.formulation == "implicit" else 2.0
    else:
        fixed_stress_scale = args.fixed_stress_scale
    if fixed_stress_scale <= 0 or args.stab_factor <= 0:
        parser.error("fixed-stress scale and stabilization factor must be positive")
    if args.solid_rtol is not None and args.solid_rtol <= 0:
        parser.error("--solid-rtol must be positive")
    if args.solid_stol is not None and args.solid_stol <= 0:
        parser.error("--solid-stol must be positive")
    if args.solid_correctors is not None and args.solid_correctors <= 0:
        parser.error("--solid-correctors must be positive")
    excluded_cell_zones = tuple(
        value
        for value in (args.excluded_cell_zones or "").replace(",", " ").split()
        if value
    )

    case_dir = pathlib.Path(__file__).resolve().parent
    base = case_dir / "base"
    default_output = (
        "tauFoRunsImplicit" if args.formulation == "implicit" else "tauFoRuns"
    )
    output_name = args.output_directory or default_output
    output_path = pathlib.Path(output_name)
    if output_path.name != output_name or output_name in (".", ".."):
        parser.error("--output-directory must be a single directory name")
    runs = case_dir / output_name
    if args.clean and args.reuse_existing:
        parser.error("--clean and --reuse-existing cannot be combined")
    if args.clean and runs.exists():
        shutil.rmtree(runs)
    runs.mkdir(parents=True, exist_ok=True)

    if not args.reuse_existing:
        # The source and stabilization exclusion share a named cell zone.
        # Build it once on the shared solid mesh before copying case points.
        run_checked(
            ["blockMesh", "-case", str(base), "-region", "solid"],
            stdout=subprocess.DEVNULL,
        )
        run_checked(
            ["topoSet", "-case", str(base), "-region", "solid"],
            stdout=subprocess.DEVNULL,
        )
        run_checked(
            ["blockMesh", "-case", str(base), "-region", "poroFluid"],
            stdout=subprocess.DEVNULL,
        )

    tau_values = parse_values(args.tau_values, DEFAULT_TAU)
    fo_values = parse_values(args.fo_values, DEFAULT_FO)

    young = 10000.0
    poisson = 0.2
    bulk = young / (3.0 * (1.0 - 2.0 * poisson))
    delta_t = 10.0
    cell_size = 0.1
    source_volume = 0.001
    source_period = 200.0

    points: list[SweepPoint] = []
    for i, tau in enumerate(tau_values):
        storage = 1.0 / (tau * bulk)
        effective_storage = storage + 1.0 / bulk
        for j, fo in enumerate(fo_values):
            conductivity = fo * effective_storage * cell_size**2 / delta_t
            points.append(
                SweepPoint(
                    i,
                    j,
                    float(tau),
                    float(fo),
                    storage,
                    conductivity,
                    f"tau_{i:02d}_Fo_{j:02d}",
                )
            )

    return_codes: dict[str, int] = {}
    if args.reuse_existing:
        for point in points:
            log_path = runs / point.case_name / "log.poroMechanicalFoam"
            if not log_path.exists():
                parser.error(f"missing existing log: {log_path}")
            return_codes[point.case_name] = 0
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
            futures = {
                executor.submit(
                    run_point,
                    base,
                    runs,
                    point,
                    args.max_iterations,
                    args.formulation,
                    fixed_stress_scale,
                    args.stab_factor,
                    args.solid_rtol,
                    args.solid_stol,
                    args.solid_correctors,
                    args.stabilization_type,
                    excluded_cell_zones,
                ): point
                for point in points
            }
            for completed, future in enumerate(
                concurrent.futures.as_completed(futures), start=1
            ):
                case_name, return_code = future.result()
                return_codes[case_name] = return_code
                print(
                    f"[{completed}/{len(points)}] {case_name}: exit={return_code}",
                    flush=True,
                )

    rows = [
        analyse_point(
            runs,
            point,
            return_codes[point.case_name],
            args.max_iterations,
            delta_t,
            source_period,
            source_volume,
        )
        for point in points
    ]

    csv_path = runs / "tauFoSweep.csv"
    with csv_path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    plot_surface(
        rows,
        tau_values,
        fo_values,
        runs / "tauFoSurface.png",
        args.formulation,
        fixed_stress_scale,
        args.stab_factor,
    )
    print(csv_path)
    print(runs / "tauFoSurface.png")


if __name__ == "__main__":
    main()
