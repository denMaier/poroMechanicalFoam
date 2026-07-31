#!/usr/bin/env python3
"""Post-process the matched scale-2.5 sweeps without rerunning the solver."""

from __future__ import annotations

import argparse
import csv
import math
import pathlib
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np

from tau_fo_sweep import (
    D_TOLERANCE,
    P_TOLERANCE,
    coupling_history,
    pressure_mode_history,
    read_internal_field,
    scaled_update_history,
)


FORMULATIONS = ("Unstabilized", "Implicit stabilization")
DIRECTORIES = (
    "tauFoRunsUnstabilizedScale2p5",
    "tauFoRunsImplicitScale2p5",
)

RATE_MARGIN = 2.0e-3
CONTRACTING = 0
UNRESOLVED = 1
GROWING = 2


def rms_rate(history: list[tuple[int, float]]) -> float:
    """Return an RMS-envelope rate from a preselected iteration window."""
    if len(history) < 4:
        return math.nan
    window = min(10, len(history) // 2)
    earlier = history[-2 * window : -window]
    later = history[-window:]
    earlier_rms = math.sqrt(sum(value * value for _, value in earlier) / window)
    later_rms = math.sqrt(sum(value * value for _, value in later) / window)
    span = (
        sum(iteration for iteration, _ in later) / window
        - sum(iteration for iteration, _ in earlier) / window
    )
    if earlier_rms <= 0 or later_rms <= 0 or span <= 0:
        return math.nan
    return float((later_rms / earlier_rms) ** (1.0 / span))


def block_fit_rate(history: list[tuple[int, float]], block: int) -> float:
    """Fit a geometric rate to several RMS blocks, suppressing cycle phase."""
    if len(history) < 4 * block:
        return math.nan
    usable = history[-min(len(history), 8 * block) :]
    block_count = len(usable) // block
    usable = usable[-block_count * block :]
    centers = []
    rms_values = []
    for start in range(0, len(usable), block):
        segment = usable[start : start + block]
        centers.append(sum(iteration for iteration, _ in segment) / block)
        rms_values.append(
            math.sqrt(sum(value * value for _, value in segment) / block)
        )
    if min(rms_values) <= 0:
        return math.nan
    slope = np.polyfit(centers, np.log(rms_values), 1)[0]
    return float(math.exp(slope))


def rate_summary(
    history: list[tuple[int, float]], noise_relative: float | None = None
) -> dict[str, float | int]:
    """Return a multi-window rate interval and conservative classification."""
    usable = [
        (iteration, abs(value))
        for iteration, value in history
        if iteration > 1 and value > 0 and math.isfinite(value)
    ]
    floor_reached = False
    if noise_relative is not None and usable:
        cutoff = usable[0][1] * noise_relative
        noise_start = next(
            (
                index
                for index, (_, value) in enumerate(usable[5:], start=5)
                if value < cutoff
            ),
            len(usable),
        )
        floor_reached = noise_start < len(usable)
        usable = usable[:noise_start]

    estimates = []
    windows = (3, 5, 10, 20, 30, 40, 60) if noise_relative is not None else (10, 20, 30, 40, 60)
    for window in windows:
        if len(usable) >= 2 * window:
            estimate = rms_rate(usable[-2 * window :])
            if math.isfinite(estimate):
                estimates.append(estimate)
    for block in (20, 30):
        estimate = block_fit_rate(usable, block)
        if math.isfinite(estimate):
            estimates.append(estimate)

    # A strongly damped checkerboard can cross the relative diagnostic floor
    # before even the shortest pair of RMS windows is available.  That is
    # evidence of suppression, not an unresolved rate.  Retain an early
    # endpoint estimate so these cases are represented honestly in the maps.
    if floor_reached and len(usable) >= 2:
        span = usable[-1][0] - usable[0][0]
        if span > 0 and usable[0][1] > 0 and usable[-1][1] > 0:
            estimates.append((usable[-1][1] / usable[0][1]) ** (1.0 / span))

    if not estimates:
        return {
            "rate": math.nan,
            "low": math.nan,
            "high": math.nan,
            "state": UNRESOLVED,
            "floorReached": int(floor_reached),
        }
    low = min(estimates)
    high = max(estimates)
    central = float(np.median(estimates))
    if floor_reached or high < 1.0 - RATE_MARGIN:
        state = CONTRACTING
    elif low > 1.0 + RATE_MARGIN:
        state = GROWING
    else:
        state = UNRESOLVED
    return {
        "rate": central,
        "low": low,
        "high": high,
        "state": state,
        "floorReached": int(floor_reached),
    }


def finite_max(*values: float) -> float:
    """Return the maximum finite value, or NaN if none are finite."""
    finite = [value for value in values if math.isfinite(value)]
    return max(finite) if finite else math.nan


def component_update_histories(log_text: str) -> dict[str, list[tuple[int, float]]]:
    history = coupling_history(log_text)
    return {
        "pressure": [(iteration, pressure) for iteration, pressure, _ in history],
        "displacement": [
            (iteration, displacement) for iteration, _, displacement in history
        ],
    }


def mode_increment_histories(log_text: str) -> dict[str, list[tuple[int, float]]]:
    pattern = re.compile(
        r"Pressure mode diagnostic ([^:]+): iteration (\d+), "
        r"amplitude [^\s,]+, increment ([^\s,]+)"
    )
    histories: dict[str, list[tuple[int, float]]] = {}
    for name, iteration, increment in pattern.findall(log_text):
        try:
            value = float(increment)
        except ValueError:
            continue
        histories.setdefault(name, []).append((int(iteration), value))
    return histories


def pressure_equation_history(log_text: str) -> list[tuple[int, float]]:
    pattern = re.compile(
        r"Pressure equation diagnostic: preMax ([^\s,]+), "
        r"preRms ([^\s,]+), postMax ([^\s,]+), postRms ([^\s,]+), "
        r"linearInitial ([^\s,]+), linearFinal ([^\s,]+), "
        r"linearIterations (\d+)"
    )
    result = []
    for values in pattern.findall(log_text):
        pre_rms = float(values[1])
        linear_iterations = int(values[6])
        if linear_iterations > 0:
            result.append((len(result) + 1, pre_rms))
    return result


def combined_rate_summary(log_text: str) -> dict[str, object]:
    updates = component_update_histories(log_text)
    pressure = rate_summary(updates["pressure"])
    displacement = rate_summary(updates["displacement"])
    modes = {
        name: rate_summary(history, noise_relative=1.0e-8)
        for name, history in mode_increment_histories(log_text).items()
    }
    checkerboard_subspace = rate_summary(
        pressure_mode_history(log_text), noise_relative=1.0e-8
    )
    pressure_equation = rate_summary(pressure_equation_history(log_text))

    full_rate = finite_max(float(pressure["rate"]), float(displacement["rate"]))
    full_low = finite_max(float(pressure["low"]), float(displacement["low"]))
    full_high = finite_max(float(pressure["high"]), float(displacement["high"]))
    component_states = (int(pressure["state"]), int(displacement["state"]))
    if GROWING in component_states:
        full_state = GROWING
    elif component_states == (CONTRACTING, CONTRACTING):
        full_state = CONTRACTING
    else:
        full_state = UNRESOLVED

    checkerboard_rate = float(checkerboard_subspace["rate"])
    checkerboard_low = float(checkerboard_subspace["low"])
    checkerboard_high = float(checkerboard_subspace["high"])
    checkerboard_state = int(checkerboard_subspace["state"])

    return {
        "pressure": pressure,
        "displacement": displacement,
        "modes": modes,
        "checkerboardSubspace": checkerboard_subspace,
        "pressureEquation": pressure_equation,
        "fullRate": full_rate,
        "fullLow": full_low,
        "fullHigh": full_high,
        "fullState": full_state,
        "checkerboardRate": checkerboard_rate,
        "checkerboardLow": checkerboard_low,
        "checkerboardHigh": checkerboard_high,
        "checkerboardState": checkerboard_state,
    }


def matched_rates(log_text: str) -> tuple[float, float]:
    """Fit the coupled and checkerboard rates over the same iterations."""
    full = dict(scaled_update_history(coupling_history(log_text)))
    mode = dict(pressure_mode_history(log_text))
    common = sorted(
        iteration
        for iteration in full.keys() & mode.keys()
        if full[iteration] > 0
        and mode[iteration] > 0
        and math.isfinite(full[iteration])
        and math.isfinite(mode[iteration])
    )
    if len(common) < 4:
        return math.nan, math.nan

    # Do not fit the stabilized checkerboard after it has reached its
    # projection/noise floor. The identical cutoff is then applied to rho.
    initial_mode = mode[common[0]]
    cutoff = initial_mode * 1.0e-8
    noise_start = next(
        (
            index
            for index, iteration in enumerate(common[5:], start=5)
            if mode[iteration] < cutoff
        ),
        len(common),
    )
    common = common[:noise_start]
    return (
        rms_rate([(iteration, full[iteration]) for iteration in common]),
        rms_rate([(iteration, mode[iteration]) for iteration in common]),
    )


def mode_amplitude_history(log_text: str) -> dict[str, list[tuple[int, float]]]:
    pattern = re.compile(
        r"Pressure mode diagnostic ([^:]+): iteration (\d+), "
        r"amplitude ([^\s,]+)"
    )
    histories: dict[str, list[tuple[int, float]]] = {}
    for name, iteration, amplitude in pattern.findall(log_text):
        histories.setdefault(name, []).append((int(iteration), float(amplitude)))
    return histories


def checkerboard_metrics(pressure: np.ndarray) -> tuple[np.ndarray, float, float]:
    count = len(pressure)
    side = round(math.sqrt(count))
    if side * side != count:
        raise RuntimeError(f"expected a square Cartesian field, found {count} cells")
    cell = np.arange(count)
    i = cell % side
    j = cell // side
    modes = np.asarray(
        [(-1.0) ** i, (-1.0) ** j, (-1.0) ** (i + j)], dtype=float
    )
    fluctuation = pressure - np.mean(pressure)
    amplitudes = modes @ fluctuation / count
    reconstruction = amplitudes @ modes
    denominator = np.linalg.norm(fluctuation)
    energy_fraction = np.linalg.norm(reconstruction) / max(denominator, 1.0e-300)
    physical_range = np.percentile(pressure, 95) - np.percentile(pressure, 5)
    range_fraction = (
        2.0 * np.max(np.abs(reconstruction)) / max(physical_range, 1.0e-300)
    )
    return reconstruction, float(energy_fraction), float(range_fraction)


def load_sweep(
    directory: pathlib.Path,
    label: str,
    extended_directory: pathlib.Path | None = None,
) -> list[dict[str, object]]:
    with (directory / "tauFoSweep.csv").open() as stream:
        source_rows = list(csv.DictReader(stream))

    rows: list[dict[str, object]] = []
    for source in source_rows:
        case = str(source["case"])
        original_case_directory = directory / case
        extended_case_directory = (
            extended_directory / case if extended_directory is not None else None
        )
        case_directory = (
            extended_case_directory
            if extended_case_directory is not None
            and (extended_case_directory / "log.poroMechanicalFoam").exists()
            else original_case_directory
        )
        log_text = (case_directory / "log.poroMechanicalFoam").read_text(
            errors="replace"
        )
        robust = combined_rate_summary(log_text)
        rho_matched, rho_cb_matched = matched_rates(log_text)
        pressure = read_internal_field(case_directory / "10" / "poroFluid" / "p_rgh")
        reconstruction, energy_fraction, range_fraction = checkerboard_metrics(
            pressure
        )
        row: dict[str, object] = dict(source)
        coupling = coupling_history(log_text)
        if case_directory != original_case_directory and coupling:
            final_iteration, final_pressure_update, final_displacement_update = coupling[-1]
            row["iterations"] = final_iteration
            row["pResidual"] = final_pressure_update
            row["DResidual"] = final_displacement_update
            row["converged"] = int(
                final_pressure_update <= P_TOLERANCE
                and final_displacement_update <= D_TOLERANCE
            )

        full_state = int(robust["fullState"])
        checkerboard_state = int(robust["checkerboardState"])

        # Reaching the coupled stopping test is direct evidence that the full
        # update contracted, even when the history is too short (or contains
        # exact zeros) for the conservative multi-window estimator.  Preserve
        # the measured source rate for plotting in that case.
        converged = int(row["converged"])
        full_rate = float(robust["fullRate"])
        full_low = float(robust["fullLow"])
        full_high = float(robust["fullHigh"])
        if converged and full_state == UNRESOLVED:
            full_state = CONTRACTING
            fallback = float(source["contractionFactor"])
            if not math.isfinite(full_rate) and math.isfinite(fallback):
                full_rate = fallback
                full_low = fallback
                full_high = fallback

        if UNRESOLVED in (full_state, checkerboard_state):
            stability_regime = 4
        elif full_state == CONTRACTING and checkerboard_state == CONTRACTING:
            stability_regime = 0
        elif full_state == GROWING and checkerboard_state == CONTRACTING:
            stability_regime = 1
        elif full_state == GROWING and checkerboard_state == GROWING:
            stability_regime = 2
        else:
            stability_regime = 3

        row.update(
            {
                "formulation": label,
                "directory": directory,
                "logText": log_text,
                "pressure": pressure,
                "checkerboardField": reconstruction,
                "matchedContractionFactor": rho_matched,
                "matchedCheckerboardContractionFactor": rho_cb_matched,
                "checkerboardEnergyFraction": energy_fraction,
                "checkerboardRangeFraction": range_fraction,
                "pressureUpdateRate": robust["pressure"]["rate"],
                "pressureUpdateRateLow": robust["pressure"]["low"],
                "pressureUpdateRateHigh": robust["pressure"]["high"],
                "displacementUpdateRate": robust["displacement"]["rate"],
                "displacementUpdateRateLow": robust["displacement"]["low"],
                "displacementUpdateRateHigh": robust["displacement"]["high"],
                "pressureEquationRate": robust["pressureEquation"]["rate"],
                "pressureEquationRateLow": robust["pressureEquation"]["low"],
                "pressureEquationRateHigh": robust["pressureEquation"]["high"],
                "robustFullRate": full_rate,
                "robustFullRateLow": full_low,
                "robustFullRateHigh": full_high,
                "robustFullState": full_state,
                "robustCheckerboardRate": robust["checkerboardRate"],
                "robustCheckerboardRateLow": robust["checkerboardLow"],
                "robustCheckerboardRateHigh": robust["checkerboardHigh"],
                "robustCheckerboardState": checkerboard_state,
                "checkerboardDiagnosticFloorReached": robust[
                    "checkerboardSubspace"
                ]["floorReached"],
                "checkerboardRateDifference": (
                    float(robust["checkerboardRate"])
                    - full_rate
                ),
                "stabilityRegime": stability_regime,
                "extended": int(case_directory != original_case_directory),
            }
        )
        rows.append(row)
    return rows


def numeric(row: dict[str, object], key: str) -> float:
    return float(row[key])


def write_derived_csv(groups: list[list[dict[str, object]]], path: pathlib.Path) -> None:
    fields = [
        "formulation",
        "case",
        "tau",
        "Fo",
        "converged",
        "iterations",
        "contractionFactor",
        "checkerboardContractionFactor",
        "matchedContractionFactor",
        "matchedCheckerboardContractionFactor",
        "pressureUpdateRate",
        "pressureUpdateRateLow",
        "pressureUpdateRateHigh",
        "displacementUpdateRate",
        "displacementUpdateRateLow",
        "displacementUpdateRateHigh",
        "pressureEquationRate",
        "pressureEquationRateLow",
        "pressureEquationRateHigh",
        "robustFullRate",
        "robustFullRateLow",
        "robustFullRateHigh",
        "robustFullState",
        "robustCheckerboardRate",
        "robustCheckerboardRateLow",
        "robustCheckerboardRateHigh",
        "robustCheckerboardState",
        "checkerboardDiagnosticFloorReached",
        "checkerboardRateDifference",
        "stabilityRegime",
        "extended",
        "checkerboardToSource",
        "checkerboardEnergyFraction",
        "checkerboardRangeFraction",
    ]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for rows in groups:
            for row in rows:
                writer.writerow({field: row[field] for field in fields})


def plot_rate_classification(
    groups: list[list[dict[str, object]]], output: pathlib.Path
) -> None:
    all_rows = [row for rows in groups for row in rows]
    strengths = np.asarray(
        [max(numeric(row, "checkerboardEnergyFraction"), 1.0e-8) for row in all_rows]
    )
    norm = colors.LogNorm(vmin=max(1.0e-5, strengths.min()), vmax=strengths.max())
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.4), constrained_layout=True)
    mappable = None
    for axis, label, rows in zip(axes, FORMULATIONS, groups):
        for converged, marker, size in ((1, "o", 52), (0, "X", 70)):
            selected = [row for row in rows if int(row["converged"]) == converged]
            mappable = axis.scatter(
                [numeric(row, "robustFullRate") for row in selected],
                [numeric(row, "robustCheckerboardRate") for row in selected],
                c=[numeric(row, "checkerboardEnergyFraction") for row in selected],
                norm=norm,
                cmap="viridis",
                marker=marker,
                s=size,
                edgecolors="black",
                linewidths=0.45,
                label="converged" if converged else "iteration cap",
            )
        unresolved = [row for row in rows if int(row["stabilityRegime"]) == 4]
        if unresolved:
            x = np.asarray([numeric(row, "robustFullRate") for row in unresolved])
            y = np.asarray(
                [numeric(row, "robustCheckerboardRate") for row in unresolved]
            )
            axis.errorbar(
                x,
                y,
                xerr=np.asarray(
                    [
                        x
                        - np.asarray(
                            [numeric(row, "robustFullRateLow") for row in unresolved]
                        ),
                        np.asarray(
                            [numeric(row, "robustFullRateHigh") for row in unresolved]
                        )
                        - x,
                    ]
                ),
                yerr=np.asarray(
                    [
                        y
                        - np.asarray(
                            [
                                numeric(row, "robustCheckerboardRateLow")
                                for row in unresolved
                            ]
                        ),
                        np.asarray(
                            [
                                numeric(row, "robustCheckerboardRateHigh")
                                for row in unresolved
                            ]
                        )
                        - y,
                    ]
                ),
                linestyle="none",
                color="0.45",
                alpha=0.55,
                linewidth=0.7,
                zorder=0,
            )
        lower = 0.0
        upper = 1.3
        axis.plot([lower, upper], [lower, upper], color="0.45", linestyle=":")
        axis.axvline(1.0, color="black", linewidth=1.4)
        axis.axhline(1.0, color="magenta", linewidth=1.4)
        axis.set_xlim(0.6, upper)
        axis.set_ylim(lower, upper)
        axis.set_title(label)
        axis.set_xlabel(r"full coupled rate $\rho$")
        axis.set_ylabel(r"checkerboard rate $\rho_{cb}$")
        axis.grid(which="both", color="0.88", linewidth=0.5)
        axis.legend(loc="lower right")
    assert mappable is not None
    colorbar = fig.colorbar(mappable, ax=axes)
    colorbar.set_label(r"final checkerboard energy fraction $C_{cb}$")
    fig.suptitle("Matched fixed-stress scale 2.5: rate classification")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def grid(rows: list[dict[str, object]], key: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    taus = np.asarray(sorted({numeric(row, "tau") for row in rows}))
    fos = np.asarray(sorted({numeric(row, "Fo") for row in rows}))
    values = np.full((len(taus), len(fos)), np.nan)
    tau_index = {value: index for index, value in enumerate(taus)}
    fo_index = {value: index for index, value in enumerate(fos)}
    for row in rows:
        values[tau_index[numeric(row, "tau")], fo_index[numeric(row, "Fo")]] = numeric(
            row, key
        )
    return taus, fos, values


def logarithmic_edges(values: np.ndarray) -> np.ndarray:
    edges = np.empty(len(values) + 1)
    edges[1:-1] = np.sqrt(values[:-1] * values[1:])
    edges[0] = values[0] / math.sqrt(values[1] / values[0])
    edges[-1] = values[-1] * math.sqrt(values[-1] / values[-2])
    return edges


def plot_parameter_maps(groups: list[list[dict[str, object]]], output: pathlib.Path) -> None:
    fig, axes = plt.subplots(2, 4, figsize=(20.0, 9.2), constrained_layout=True)
    regime_cmap = colors.ListedColormap(
        ["#4daf4a", "#377eb8", "#e41a1c", "#ff7f00", "#969696"]
    )
    regime_norm = colors.BoundaryNorm(np.arange(-0.5, 5.5, 1.0), regime_cmap.N)
    source_strengths = [
        max(numeric(row, "checkerboardToSource"), 1.0e-8)
        for rows in groups
        for row in rows
    ]
    source_strength_norm = colors.LogNorm(
        vmin=max(1.0e-5, min(source_strengths)), vmax=max(source_strengths)
    )
    strength_values = [
        max(numeric(row, "checkerboardEnergyFraction"), 1.0e-8)
        for rows in groups
        for row in rows
    ]
    strength_norm = colors.LogNorm(
        vmin=max(1.0e-5, min(strength_values)), vmax=max(strength_values)
    )
    images = [None, None, None, None]
    iteration_max = max(
        numeric(row, "iterations") for rows in groups for row in rows
    )
    for row_index, (label, rows) in enumerate(zip(FORMULATIONS, groups)):
        tau, fo, _ = grid(rows, "checkerboardRateDifference")
        _, _, regime = grid(rows, "stabilityRegime")
        _, _, source_strength = grid(rows, "checkerboardToSource")
        _, _, strength = grid(rows, "checkerboardEnergyFraction")
        _, _, iterations = grid(rows, "iterations")
        _, _, converged = grid(rows, "converged")
        x_edges, y_edges = np.meshgrid(
            logarithmic_edges(fo), logarithmic_edges(tau)
        )
        images[0] = axes[row_index, 0].pcolormesh(
            x_edges,
            y_edges,
            regime,
            shading="flat",
            cmap=regime_cmap,
            norm=regime_norm,
        )
        images[1] = axes[row_index, 1].pcolormesh(
            x_edges,
            y_edges,
            source_strength,
            shading="flat",
            cmap="viridis",
            norm=source_strength_norm,
        )
        images[2] = axes[row_index, 2].pcolormesh(
            x_edges,
            y_edges,
            strength,
            shading="flat",
            cmap="viridis",
            norm=strength_norm,
        )
        images[3] = axes[row_index, 3].pcolormesh(
            x_edges,
            y_edges,
            iterations,
            shading="flat",
            cmap="magma",
            vmin=0,
            vmax=iteration_max,
        )
        for axis in axes[row_index]:
            failed_i, failed_j = np.where(converged < 0.5)
            axis.scatter(
                fo[failed_j], tau[failed_i], marker="x", color="white", s=27, linewidths=1.0
            )
            axis.set_xscale("log")
            axis.set_yscale("log")
            axis.set_xlabel(r"Fourier number $\mathrm{Fo}$")
            axis.set_ylabel(r"coupling strength $\tau$")
            axis.grid(False)
        axes[row_index, 0].set_ylabel(label + "\n" + r"coupling strength $\tau$")
    axes[0, 0].set_title("uncertainty-aware growth regime")
    axes[0, 1].set_title(r"amplitude $|p_{cb}|/p_{source}$")
    axes[0, 2].set_title(r"energy fraction $C_{cb}$")
    axes[0, 3].set_title("outer iterations (white x: capped)")
    labels = (
        "growth regime",
        r"$|p_{cb}|/p_{source}$",
        r"$C_{cb}$",
        "iterations",
    )
    for column, (image, label) in enumerate(zip(images, labels)):
        assert image is not None
        colorbar = fig.colorbar(image, ax=axes[:, column], label=label)
        if column == 0:
            colorbar.set_ticks([0, 1, 2, 3, 4])
            colorbar.set_ticklabels(
                [
                    "both contract",
                    "full grows; cb contracts",
                    "both grow",
                    "cb-only/inconsistent",
                    "unresolved",
                ]
            )
    fig.suptitle(
        "Matched fixed-stress scale 2.5: dynamics, final contamination, and convergence"
    )
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_residual_diagnostics(
    groups: list[list[dict[str, object]]], output: pathlib.Path
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.4), constrained_layout=True)
    mappable = None
    all_strengths = [
        max(numeric(row, "checkerboardToSource"), 1.0e-8)
        for rows in groups
        for row in rows
        if int(row["extended"])
    ]
    norm = colors.LogNorm(vmin=min(all_strengths), vmax=max(all_strengths))
    for axis, label, rows in zip(axes, FORMULATIONS, groups):
        selected = [
            row
            for row in rows
            if int(row["extended"])
            and math.isfinite(numeric(row, "pressureEquationRate"))
        ]
        x = np.asarray([numeric(row, "pressureUpdateRate") for row in selected])
        y = np.asarray([numeric(row, "pressureEquationRate") for row in selected])
        mappable = axis.scatter(
            x,
            y,
            c=[numeric(row, "checkerboardToSource") for row in selected],
            norm=norm,
            cmap="viridis",
            edgecolors="black",
            linewidths=0.45,
            s=58,
        )
        lower = min(0.9, float(x.min()), float(y.min()))
        upper = max(1.1, float(x.max()), float(y.max()))
        axis.plot([lower, upper], [lower, upper], color="0.4", linestyle=":")
        axis.axvline(1.0, color="black", linewidth=1.2)
        axis.axhline(1.0, color="black", linewidth=1.2)
        axis.set_xlim(lower, upper)
        axis.set_ylim(lower, upper)
        axis.set_aspect("equal", adjustable="box")
        axis.set_title(label)
        axis.set_xlabel(r"pressure-update rate $\rho_{\Delta p}$")
        axis.set_ylabel(r"pre-solve equation-defect rate $\rho_{r_p}$")
        axis.grid(color="0.9", linewidth=0.5)
    assert mappable is not None
    fig.colorbar(mappable, ax=axes, label=r"$|p_{cb}|/p_{source}$")
    fig.suptitle("Pressure update versus pressure-equation defect (extended cases)")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def select_representatives(rows: list[dict[str, object]]) -> list[str]:
    growing = [
        row
        for row in rows
        if int(row["robustCheckerboardState"]) == GROWING
    ]
    growing_case = max(
        growing, key=lambda row: numeric(row, "checkerboardEnergyFraction")
    )["case"]
    near_neutral = [
        row
        for row in rows
        if int(row["robustCheckerboardState"]) != GROWING
    ]
    neutral_case = min(
        near_neutral,
        key=lambda row: abs(numeric(row, "robustCheckerboardRate") - 1.0),
    )["case"]
    return [str(growing_case), str(neutral_case)]


def plot_histories(groups: list[list[dict[str, object]]], output: pathlib.Path) -> list[str]:
    cases = select_representatives(groups[0])
    lookups = [{str(row["case"]): row for row in rows} for rows in groups]
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 8.6), constrained_layout=True)
    line_colors = ("tab:red", "tab:blue")
    for column, case in enumerate(cases):
        for label, lookup, color in zip(FORMULATIONS, lookups, line_colors):
            row = lookup[case]
            log_text = str(row["logText"])
            full = scaled_update_history(coupling_history(log_text))
            mode = pressure_mode_history(log_text)
            mode_scale = next((value for _, value in mode if value > 0), 1.0)
            axes[0, column].semilogy(
                [iteration for iteration, _ in full],
                [value for _, value in full],
                color=color,
                label=label,
            )
            axes[1, column].semilogy(
                [iteration for iteration, value in mode if value > 0],
                [value / mode_scale for _, value in mode if value > 0],
                color=color,
                label=label,
            )
        row = lookups[0][case]
        axes[0, column].set_title(
            rf"$\tau={numeric(row, 'tau'):g}$, $\mathrm{{Fo}}={numeric(row, 'Fo'):g}$"
        )
        axes[0, column].axhline(1.0, color="0.35", linestyle=":", label="tolerance")
        axes[0, column].set_ylabel("scaled coupled update")
        axes[1, column].set_ylabel(r"$|\Delta a_{cb}|/|\Delta a_{cb}^{(1)}|$")
        axes[1, column].set_xlabel("outer iteration")
        for axis in axes[:, column]:
            axis.grid(which="both", color="0.88", linewidth=0.5)
            axis.legend()
    fig.suptitle("Representative matched-scale iteration histories")
    fig.savefig(output, dpi=220)
    plt.close(fig)
    return cases


def plot_fields(
    groups: list[list[dict[str, object]]], cases: list[str], output: pathlib.Path
) -> None:
    lookups = [{str(row["case"]): row for row in rows} for rows in groups]
    fig, axes = plt.subplots(2, 4, figsize=(15.6, 7.4), constrained_layout=True)
    titles = (
        "unstabilized pressure",
        "stabilized pressure",
        "unstabilized checkerboard",
        "stabilized checkerboard",
    )
    for axis, title in zip(axes[0], titles):
        axis.set_title(title)
    for row_index, case in enumerate(cases):
        unstabilized = lookups[0][case]
        stabilized = lookups[1][case]
        source_scale = max(numeric(unstabilized, "checkerboardAmplitude") / max(numeric(unstabilized, "checkerboardToSource"), 1.0e-300), 1.0e-300)
        fields = (
            np.asarray(unstabilized["pressure"]) - np.mean(unstabilized["pressure"]),
            np.asarray(stabilized["pressure"]) - np.mean(stabilized["pressure"]),
            np.asarray(unstabilized["checkerboardField"]),
            np.asarray(stabilized["checkerboardField"]),
        )
        normalized = [field.reshape((10, 10)) / source_scale for field in fields]
        limit = max(np.max(np.abs(field)) for field in normalized)
        for column, field in enumerate(normalized):
            image = axes[row_index, column].imshow(
                field,
                origin="lower",
                cmap="RdBu_r",
                vmin=-limit,
                vmax=limit,
                interpolation="nearest",
            )
            axes[row_index, column].set_xticks([])
            axes[row_index, column].set_yticks([])
        fig.colorbar(image, ax=axes[row_index, :], label=r"$p/p_{source}$")
        axes[row_index, 0].set_ylabel(
            rf"$\tau={numeric(unstabilized, 'tau'):g}$\n$\mathrm{{Fo}}={numeric(unstabilized, 'Fo'):g}$"
        )
    fig.suptitle("Final pressure and extracted checkerboard fields")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--unstabilized-directory", default=DIRECTORIES[0])
    parser.add_argument("--stabilized-directory", default=DIRECTORIES[1])
    parser.add_argument(
        "--extended-unstabilized-directory",
        default="matchedScale2p5Extended/unstabilized",
    )
    parser.add_argument(
        "--extended-stabilized-directory",
        default="matchedScale2p5Extended/implicit",
    )
    parser.add_argument("--output-directory", default="matchedScale2p5Plots")
    args = parser.parse_args()

    case_directory = pathlib.Path(__file__).resolve().parent
    directories = (
        case_directory / args.unstabilized_directory,
        case_directory / args.stabilized_directory,
    )
    extended_directories = (
        case_directory / args.extended_unstabilized_directory,
        case_directory / args.extended_stabilized_directory,
    )
    for directory in directories:
        if not directory.is_dir():
            parser.error(f"missing result directory: {directory}")
    output = case_directory / args.output_directory
    output.mkdir(parents=True, exist_ok=True)

    groups = [
        load_sweep(directory, label, extended)
        for directory, extended, label in zip(
            directories, extended_directories, FORMULATIONS
        )
    ]
    write_derived_csv(groups, output / "matchedScale2p5Metrics.csv")
    plot_rate_classification(groups, output / "rateClassification.png")
    plot_parameter_maps(groups, output / "parameterMaps.png")
    plot_residual_diagnostics(groups, output / "residualDiagnostics.png")
    representative_cases = plot_histories(groups, output / "iterationHistories.png")
    plot_fields(groups, representative_cases, output / "pressureFields.png")
    for path in sorted(output.iterdir()):
        print(path)


if __name__ == "__main__":
    main()
