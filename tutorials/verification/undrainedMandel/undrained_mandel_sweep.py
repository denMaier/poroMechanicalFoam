#!/usr/bin/python3
"""Sweep the fully undrained first increment of a Mandel-type problem."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import math
import pathlib
import re
import shutil
import subprocess
from dataclasses import dataclass

import numpy as np


NX = 80
NZ = 16
YOUNG = 13_333_333.3333333
POISSON = 1.0 / 3.0
K_VOL = YOUNG / (3.0 * (1.0 - 2.0 * POISSON))
CELL_SIZE = 1.0 / NX
TOP_DISPLACEMENT = -1.0e-3
DEFAULT_TAU = (1.0, 3.0, 10.0, 30.0, 100.0)
DEFAULT_FO = (1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0)


@dataclass(frozen=True)
class Point:
    tau_index: int
    fo_index: int
    tau: float
    fo: float
    storage: float
    dt: float
    name: str


def values(text: str | None, default: tuple[float, ...]) -> list[float]:
    if text is None:
        return list(default)
    result = [float(item) for item in text.replace(",", " ").split()]
    if not result or any(item <= 0.0 for item in result):
        raise ValueError("sweep values must be strictly positive")
    return result


def run_checked(command: list[str], stdout=None) -> None:
    subprocess.run(command, check=True, stdout=stdout, stderr=subprocess.STDOUT)


def set_entry(path: pathlib.Path, entry: str, value: str) -> None:
    run_checked(
        ["foamDictionary", str(path), "-entry", entry, "-set", value],
        stdout=subprocess.DEVNULL,
    )


def mode_values(kind: str) -> np.ndarray:
    x = np.tile((-1.0) ** np.arange(NX), NZ)
    z = np.repeat((-1.0) ** np.arange(NZ), NX)
    return {"checkerboardX": x, "checkerboardZ": z, "checkerboardXZ": x * z}[kind]


def write_mode(path: pathlib.Path, name: str, data: np.ndarray) -> None:
    rows = [" ".join(f"{value:g}" for value in data[k : k + NX]) for k in range(0, len(data), NX)]
    path.write_text(
        "FoamFile\n"
        "{\n"
        "    version 2.0;\n"
        "    format ascii;\n"
        "    class volScalarField;\n"
        "    location \"constant/poroFluid\";\n"
        f"    object {name};\n"
        "}\n\n"
        "dimensions [0 0 0 0 0 0 0];\n\n"
        f"internalField nonuniform List<scalar>\n{len(data)}\n(\n"
        + "\n".join(rows)
        + "\n);\n\n"
        "boundaryField\n"
        "{\n"
        "    \"left|right|bottom|top\" { type zeroGradient; }\n"
        "    frontAndBack { type empty; }\n"
        "}\n"
    )


def build_base(base: pathlib.Path) -> None:
    run_checked(["blockMesh", "-case", str(base), "-region", "solid"], stdout=subprocess.DEVNULL)
    run_checked(["blockMesh", "-case", str(base), "-region", "poroFluid"], stdout=subprocess.DEVNULL)
    for name in ("checkerboardX", "checkerboardZ", "checkerboardXZ"):
        write_mode(base / "constant" / "poroFluid" / name, name, mode_values(name))


def prepare_case(
    base: pathlib.Path,
    runs: pathlib.Path,
    point: Point,
    formulation: str,
    max_iterations: int,
    fixed_stress_scale: float,
    stab_factor: float,
    load_variation: float,
) -> pathlib.Path:
    target = runs / point.name
    shutil.copytree(base, target)
    (target / f"{point.name}.foam").touch()

    control = target / "system" / "controlDict"
    coupling = target / "constant" / "poroCouplingProperties"
    hydraulic = target / "constant" / "poroFluid" / "poroHydraulicProperties"
    set_entry(control, "deltaT", f"{point.dt:.16g}")
    set_entry(control, "endTime", f"{point.dt:.16g}")
    set_entry(
        hydraulic,
        "storageCoeffCoeffs/Ss",
        f"Ss [-1 1 2 0 0 0 0] {point.storage:.16g}",
    )
    set_entry(coupling, "fixedStressStabilScale", f"{fixed_stress_scale:.16g}")

    stabilized = formulation != "unstabilized"
    coeffs = "poroSolidStabCoeffs" if stabilized else "poroSolidCoeffs"
    set_entry(coupling, "poroSolidInterface", "poroSolidStab" if stabilized else "poroSolid")
    set_entry(coupling, f"{coeffs}/iterations", str(max_iterations))
    set_entry(coupling, f"{coeffs}/infoFrequency", "1")
    if stabilized:
        set_entry(coupling, f"{coeffs}/stabilizationType", formulation)
        set_entry(coupling, f"{coeffs}/stabFactor", f"{stab_factor:.16g}")

    (target / "constant" / "topDisplacement.dat").write_text(
        "3\n"
        "(\n"
        "(0 (0 0 0))\n"
        f"({point.dt:.16g} (0 0 {TOP_DISPLACEMENT:.16g}))\n"
        f"({2.0 * point.dt:.16g} (0 0 {TOP_DISPLACEMENT:.16g}))\n"
        ")\n"
    )
    d_path = target / "0" / "solid" / "D"
    d_text = d_path.read_text()
    spatial_scale = [
        1.0 + load_variation * math.cos(math.pi * (i + 0.5) / NX)
        for i in range(NX)
    ]
    scale_entry = (
        "spatialScale nonuniform List<scalar>\n"
        f"        {NX}\n"
        "        (\n"
        + "\n".join(f"        {value:.16g}" for value in spatial_scale)
        + "\n        );"
    )
    d_text, count = re.subn(
        r"spatialScale\s+uniform\s+[^;]+;", scale_entry, d_text, count=1
    )
    if count != 1:
        raise RuntimeError(f"cannot replace top spatialScale in {d_path}")
    d_path.write_text(d_text)
    return target


def run_case(target: pathlib.Path) -> tuple[str, int]:
    log = target / "log.poroMechanicalFoam"
    with log.open("w") as output:
        result = subprocess.run(
            ["poroMechanicalFoam", "-case", str(target)],
            stdout=output,
            stderr=subprocess.STDOUT,
        )
    return target.name, result.returncode


def envelope_rate(history: list[tuple[int, float]]) -> float:
    usable = [(iteration, value) for iteration, value in history if iteration > 1 and value > 0 and math.isfinite(value)]
    if len(usable) < 6:
        return math.nan
    window = min(10, len(usable) // 2)
    first = usable[-2 * window : -window]
    second = usable[-window:]
    first_rms = math.sqrt(sum(value * value for _, value in first) / window)
    second_rms = math.sqrt(sum(value * value for _, value in second) / window)
    span = np.mean([iteration for iteration, _ in second]) - np.mean([iteration for iteration, _ in first])
    return float((second_rms / first_rms) ** (1.0 / span))


def checkerboard_rate(history: list[tuple[int, float]]) -> float:
    usable = [
        (iteration, value)
        for iteration, value in history
        if iteration > 1 and value > 0 and math.isfinite(value)
    ]
    if len(usable) < 6:
        return math.nan
    # Once the projected update has dropped eight orders below its initial
    # value, its late ratios describe sign-changing round-off rather than a
    # checkerboard eigenmode.  Fit the last envelopes before that floor.
    cutoff = usable[0][1] * 1.0e-8
    noise_start = next(
        (
            index
            for index, (_, value) in enumerate(usable[5:], start=5)
            if value < cutoff
        ),
        len(usable),
    )
    return envelope_rate(usable[:noise_start])


def parse_histories(text: str) -> dict[str, object]:
    updates = [
        (int(i), float(value))
        for i, value in re.findall(
            r"Pressure increment decomposition: iteration (\d+), totalRms ([^\s,]+)", text
        )
    ]
    checkerboard = [
        (int(i), float(value))
        for i, value in re.findall(
            r"Pressure mode subspace diagnostic: iteration (\d+), incrementNorm ([^\s]+)", text
        )
    ]
    equation = [
        (index + 1, float(value))
        for index, value in enumerate(
            re.findall(r"Pressure equation diagnostic: preMax [^,]+, preRms ([^,]+),", text)
        )
    ]
    coupling = [
        (int(i), float(p), float(d))
        for i, p, d in re.findall(
            r"\|\s*#\s+Iterations\s*\|\s*(\d+)\s*\|[^\n]*\n"
            r"\|\s*max\s+delta\(p_rgh\)\s*\|\s*([^\s|]+)\s*\|[^\n]*\n"
            r"\|\s*max\s+delta\(D\)\s*\|\s*([^\s|]+)",
            text,
        )
    ]
    return {
        "rho_update": envelope_rate(updates),
        "rho_pressure_equation": envelope_rate(equation),
        "rho_checkerboard": checkerboard_rate(checkerboard),
        "iterations": max((i for i, _, _ in coupling), default=0),
        "final_p_update": coupling[-1][1] if coupling else math.nan,
        "final_d_update": coupling[-1][2] if coupling else math.nan,
    }


def read_scalar_field(path: pathlib.Path) -> np.ndarray:
    text = path.read_text()
    uniform = re.search(r"internalField\s+uniform\s+([^;]+);", text)
    if uniform:
        return np.full(NX * NZ, float(uniform.group(1)))
    match = re.search(
        r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;",
        text,
        flags=re.DOTALL,
    )
    if not match:
        raise RuntimeError(f"cannot read internal field from {path}")
    data = np.asarray([float(value) for value in match.group(2).split()])
    if len(data) != int(match.group(1)):
        raise RuntimeError(f"wrong field size in {path}")
    return data


def final_pressure(case: pathlib.Path) -> np.ndarray:
    times = []
    for child in case.iterdir():
        try:
            time = float(child.name)
        except ValueError:
            continue
        path = child / "poroFluid" / "p_rgh"
        if time > 0.0 and path.exists():
            times.append((time, path))
    if not times:
        raise RuntimeError(f"no final pressure in {case}")
    return read_scalar_field(max(times)[1])


def pressure_metrics(p: np.ndarray) -> dict[str, float]:
    mean = float(np.mean(p))
    rms = float(np.sqrt(np.mean(p * p)))
    fluctuation_rms = float(np.sqrt(np.mean((p - mean) ** 2)))
    amplitudes = [float(np.mean(p * mode_values(name))) for name in ("checkerboardX", "checkerboardZ", "checkerboardXZ")]
    cb = float(np.linalg.norm(amplitudes))
    grid = p.reshape(NZ, NX)
    highpass_x = 0.25 * (grid[:, :-2] - 2.0 * grid[:, 1:-1] + grid[:, 2:])
    highpass_z = 0.25 * (grid[:-2, :] - 2.0 * grid[1:-1, :] + grid[2:, :])
    highpass_x_rms = float(np.sqrt(np.mean(highpass_x * highpass_x)))
    highpass_z_rms = float(np.sqrt(np.mean(highpass_z * highpass_z)))
    gridscale_rms = math.hypot(highpass_x_rms, highpass_z_rms)
    # The first high-pass row/column touches a boundary value.  Exclude two
    # cell layers when measuring bulk checkerboarding so the prescribed top
    # displacement response is not mislabeled as an interior pressure mode.
    interior_highpass_x = highpass_x[2:-2, 1:-1]
    interior_highpass_z = highpass_z[1:-1, 2:-2]
    interior_gridscale_rms = math.hypot(
        float(np.sqrt(np.mean(interior_highpass_x * interior_highpass_x))),
        float(np.sqrt(np.mean(interior_highpass_z * interior_highpass_z))),
    )
    scale = max(rms, np.finfo(float).tiny)
    fluctuation_scale = max(fluctuation_rms, np.finfo(float).tiny)
    return {
        "mean_pressure": mean,
        "pressure_rms": rms,
        "pressure_fluctuation_rms": fluctuation_rms,
        "checkerboard_amplitude": cb,
        "checkerboard_to_pressure_rms": cb / scale,
        "checkerboard_to_fluctuation_rms": cb / fluctuation_scale,
        "checkerboard_to_mean_pressure": cb / max(abs(mean), np.finfo(float).tiny),
        "checkerboard_x": amplitudes[0],
        "checkerboard_z": amplitudes[1],
        "checkerboard_xz": amplitudes[2],
        "highpass_x_rms": highpass_x_rms,
        "highpass_z_rms": highpass_z_rms,
        "gridscale_rms": gridscale_rms,
        "gridscale_to_pressure_rms": gridscale_rms / scale,
        "gridscale_to_fluctuation_rms": gridscale_rms / fluctuation_scale,
        "interior_gridscale_rms": interior_gridscale_rms,
        "interior_gridscale_to_pressure_rms": interior_gridscale_rms / scale,
        "interior_gridscale_to_fluctuation_rms":
            interior_gridscale_rms / fluctuation_scale,
    }


def analyse(
    case: pathlib.Path,
    point: Point,
    return_code: int,
    max_iterations: int,
    load_variation: float,
    formulation: str,
    fixed_stress_scale: float,
    stab_factor: float,
) -> dict[str, object]:
    log_text = (case / "log.poroMechanicalFoam").read_text(errors="replace")
    history = parse_histories(log_text)
    row: dict[str, object] = {
        "case": point.name,
        "tau": point.tau,
        "Fo": point.fo,
        "Ss": point.storage,
        "delta_t": point.dt,
        "load_variation": load_variation,
        "formulation": formulation,
        "fixed_stress_scale": fixed_stress_scale,
        "stab_factor": stab_factor if formulation != "unstabilized" else 0.0,
        "return_code": return_code,
        **history,
    }
    try:
        p = final_pressure(case)
        row.update(pressure_metrics(p))
    except (OSError, RuntimeError, ValueError):
        for name in (
            "mean_pressure", "pressure_rms", "pressure_fluctuation_rms",
            "checkerboard_amplitude", "checkerboard_to_pressure_rms",
            "checkerboard_to_fluctuation_rms",
            "checkerboard_to_mean_pressure", "checkerboard_x", "checkerboard_z",
            "checkerboard_xz", "highpass_x_rms", "highpass_z_rms",
            "gridscale_rms", "gridscale_to_pressure_rms",
            "gridscale_to_fluctuation_rms", "interior_gridscale_rms",
            "interior_gridscale_to_pressure_rms",
            "interior_gridscale_to_fluctuation_rms",
        ):
            row[name] = math.nan
    row["converged"] = bool(
        return_code == 0
        and int(row["iterations"]) < max_iterations
        and float(row["final_p_update"]) <= 1.0e-4
        and float(row["final_d_update"]) <= 1.0e-10
    )
    return row


def logarithmic_edges(data: list[float]) -> np.ndarray:
    values = np.asarray(data, dtype=float)
    if len(values) == 1:
        return np.asarray([values[0] / math.sqrt(10.0), values[0] * math.sqrt(10.0)])
    logs = np.log(values)
    interior = 0.5 * (logs[:-1] + logs[1:])
    return np.exp(np.concatenate(([logs[0] - (interior[0] - logs[0])], interior, [logs[-1] + (logs[-1] - interior[-1])])) )


def plot_summary(rows: list[dict[str, object]], tau_values: list[float], fo_values: list[float], path: pathlib.Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, Normalize

    metrics = (
        ("rho_update", r"$\rho_{\Delta p}$", "coolwarm", Normalize(0.7, 1.3)),
        ("rho_pressure_equation", r"$\rho_{r_p}$", "coolwarm", Normalize(0.7, 1.3)),
        ("rho_checkerboard", r"$\rho_{cb}$", "coolwarm", Normalize(0.7, 1.3)),
        ("interior_gridscale_to_pressure_rms", r"$\|H_hp\|_{int}/\|p\|_{RMS}$", "viridis", LogNorm(1e-12, 1.0)),
    )
    fig, axes = plt.subplots(1, 4, figsize=(15, 3.8), constrained_layout=True)
    for axis, (key, title, cmap, norm) in zip(axes, metrics):
        image = np.full((len(tau_values), len(fo_values)), np.nan)
        for row in rows:
            i = tau_values.index(float(row["tau"]))
            j = fo_values.index(float(row["Fo"]))
            value = float(row[key])
            if math.isfinite(value) and (not isinstance(norm, LogNorm) or value > 0):
                image[i, j] = value
        shown = axis.pcolormesh(
            logarithmic_edges(fo_values), logarithmic_edges(tau_values), image,
            shading="flat", cmap=cmap, norm=norm,
        )
        axis.set_title(title)
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlabel(r"Fourier number $Fo$")
        axis.set_ylabel(r"Coupling strength $\tau=M/K_{vol}$")
        fig.colorbar(shown, ax=axis, shrink=0.82)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_pressures(rows: list[dict[str, object]], runs: pathlib.Path, path: pathlib.Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plottable = []
    for row in rows:
        try:
            p = final_pressure(runs / str(row["case"])).reshape(NZ, NX)
        except (OSError, RuntimeError, ValueError):
            continue
        scale = max(abs(float(np.mean(p))), float(np.sqrt(np.mean(p * p))), np.finfo(float).tiny)
        plottable.append((row, (p - np.mean(p)) / scale))
    if not plottable:
        return
    ncols = min(3, len(plottable))
    nrows = math.ceil(len(plottable) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.3 * ncols, 2.5 * nrows), squeeze=False, constrained_layout=True)
    fields = [field for _, field in plottable]
    limit = max(np.quantile(np.abs(field), 0.995) for field in fields)
    limit = max(limit, 1e-12)
    for axis, (row, field) in zip(axes.flat, plottable):
        shown = axis.imshow(field, origin="lower", extent=(0, 1, 0, 0.2), aspect="auto", cmap="RdBu_r", vmin=-limit, vmax=limit)
        axis.set_title(
            rf"$\tau={float(row['tau']):.0e}$, $Fo={float(row['Fo']):.0e}$"
            + "\n"
            + rf"$\|H_hp\|_{{int}}/\|p\|={float(row['interior_gridscale_to_pressure_rms']):.2e}$"
        )
        axis.set_xlabel("x")
        axis.set_ylabel("z")
    for axis in axes.flat[len(plottable):]:
        axis.set_visible(False)
    fig.colorbar(shown, ax=axes, label=r"$(p-\bar p)/\|p\|_{RMS}$", shrink=0.8)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tau-values", help="values of tau=M/K_vol")
    parser.add_argument("--fo-values", help="values of the mesh Fourier number Fo")
    parser.add_argument(
        "--formulation",
        choices=("unstabilized", "implicit", "explicit", "pressureJump"),
        default="unstabilized",
    )
    parser.add_argument("--max-iterations", type=int, default=600)
    parser.add_argument("--fixed-stress-scale", type=float, default=1.5)
    parser.add_argument("--stab-factor", type=float, default=1.0)
    parser.add_argument("--conductivity", type=float, default=1.0e-12)
    parser.add_argument(
        "--load-variation",
        type=float,
        default=0.1,
        help="spatial multiplier amplitude in 1+epsilon*cos(pi*x/a)",
    )
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--output-directory")
    parser.add_argument("--clean", action="store_true")
    parser.add_argument(
        "--reuse-existing",
        action="store_true",
        help="reanalyse logs and fields in the output directory without rerunning",
    )
    parser.add_argument(
        "--resume-existing",
        action="store_true",
        help="reuse completed cases and run only missing tau--Fo points",
    )
    parser.add_argument(
        "--rerun-unresolved",
        action="store_true",
        help="reuse converged cases and rerun only absent or unresolved tau--Fo points",
    )
    args = parser.parse_args()
    if args.max_iterations <= 0 or args.jobs <= 0:
        parser.error("iteration and job counts must be positive")
    if args.fixed_stress_scale <= 0 or args.stab_factor <= 0 or args.conductivity <= 0:
        parser.error("scale and conductivity values must be positive")
    if not 0.0 <= args.load_variation < 1.0:
        parser.error("--load-variation must be in [0, 1)")

    tau_values = values(args.tau_values, DEFAULT_TAU)
    fo_values = values(args.fo_values, DEFAULT_FO)
    case_dir = pathlib.Path(__file__).resolve().parent
    base = case_dir / "base"
    output_name = args.output_directory or f"runs_{args.formulation}"
    if pathlib.Path(output_name).name != output_name:
        parser.error("output directory must be a single directory name")
    runs = case_dir / output_name
    if sum((args.clean, args.reuse_existing, args.resume_existing, args.rerun_unresolved)) > 1:
        parser.error(
            "--clean, --reuse-existing, --resume-existing and "
            "--rerun-unresolved are mutually exclusive"
        )
    if args.clean and runs.exists():
        shutil.rmtree(runs)
    if not (
        args.reuse_existing or args.resume_existing or args.rerun_unresolved
    ) and runs.exists() and any(runs.iterdir()):
        parser.error(f"output directory is not empty: {runs}; use --clean")
    runs.mkdir(parents=True, exist_ok=True)

    if not args.reuse_existing:
        set_entry(
            base / "constant" / "poroFluid" / "poroHydraulicProperties",
            "k",
            f"k [0 1 -1 0 0 0 0] {args.conductivity:.16g}",
        )
        build_base(base)

    points = []
    for i, tau in enumerate(tau_values):
        storage = 1.0 / (tau * K_VOL)
        effective_storage = storage + 1.0 / K_VOL
        for j, fo in enumerate(fo_values):
            dt = fo * effective_storage * CELL_SIZE**2 / args.conductivity
            points.append(Point(i, j, tau, fo, storage, dt, f"tau_{i:02d}_Fo_{j:02d}"))

    return_codes: dict[str, int] = {}
    if args.reuse_existing:
        for point in points:
            if not (runs / point.name / "log.poroMechanicalFoam").exists():
                parser.error(f"missing existing case {runs / point.name}")
            return_codes[point.name] = 0
    else:
        previous_convergence: dict[str, bool] = {}
        previous_summary = runs / "summary.csv"
        if args.rerun_unresolved and previous_summary.exists():
            with previous_summary.open(newline="") as source:
                previous_convergence = {
                    row["case"]: row["converged"] == "True"
                    for row in csv.DictReader(source)
                }
        cases = []
        for point in points:
            log = runs / point.name / "log.poroMechanicalFoam"
            reuse_converged = (
                args.rerun_unresolved
                and log.exists()
                and previous_convergence.get(point.name, False)
            )
            if (args.resume_existing and log.exists()) or reuse_converged:
                return_codes[point.name] = 0
            else:
                if args.rerun_unresolved and (runs / point.name).exists():
                    shutil.rmtree(runs / point.name)
                cases.append(
                    prepare_case(
                        base, runs, point, args.formulation, args.max_iterations,
                        args.fixed_stress_scale, args.stab_factor,
                        args.load_variation,
                    )
                )
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
            futures = {executor.submit(run_case, case): case for case in cases}
            for number, future in enumerate(concurrent.futures.as_completed(futures), 1):
                name, code = future.result()
                return_codes[name] = code
                print(f"[{number}/{len(cases)}] {name}: exit={code}", flush=True)

    rows = [
        analyse(
            runs / point.name,
            point,
            return_codes[point.name],
            args.max_iterations,
            args.load_variation,
            args.formulation,
            args.fixed_stress_scale,
            args.stab_factor,
        )
        for point in points
    ]
    csv_path = runs / "summary.csv"
    with csv_path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    plot_summary(rows, tau_values, fo_values, runs / "tauFoMap.png")
    plot_pressures(rows, runs, runs / "pressureFields.png")
    print(csv_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
