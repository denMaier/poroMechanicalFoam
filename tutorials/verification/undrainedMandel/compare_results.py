#!/usr/bin/python3
"""Make matched tau--Fo plots for undrained Mandel sweeps."""

from __future__ import annotations

import argparse
import csv
import math
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import numpy as np

from undrained_mandel_sweep import NX, NZ, final_pressure, logarithmic_edges


def rows(path: pathlib.Path) -> dict[tuple[float, float], dict[str, str]]:
    with (path / "summary.csv").open() as handle:
        data = list(csv.DictReader(handle))
    return {(float(row["tau"]), float(row["Fo"])): row for row in data}


def dataset_label(data: dict[tuple[float, float], dict[str, str]]) -> str:
    row = next(iter(data.values()))
    formulation = row["formulation"]
    if formulation == "unstabilized":
        return "unstabilized"
    factor = float(row["stab_factor"])
    if formulation == "pressureJump":
        return f"paper pressure-jump, factor {factor:g}"
    if formulation == "explicit":
        return f"explicit compact-minus-wide, factor {factor:g}"
    return f"implicit compact-minus-wide, factor {factor:g}"


def surface(
    data: dict[tuple[float, float], dict[str, str]],
    tau_values: list[float],
    fo_values: list[float],
    metric: str,
) -> np.ndarray:
    result = np.full((len(tau_values), len(fo_values)), np.nan)
    for i, tau in enumerate(tau_values):
        for j, fo in enumerate(fo_values):
            value = float(data[(tau, fo)][metric])
            if math.isfinite(value):
                result[i, j] = value
    return result


def plot_surfaces(
    datasets: tuple[dict[tuple[float, float], dict[str, str]], ...],
    tau_values: list[float],
    fo_values: list[float],
    output: pathlib.Path,
) -> None:
    max_iterations = max(
        float(row["iterations"])
        for data in datasets
        for row in data.values()
    )
    metrics = (
        ("rho_update", r"$\rho_{\Delta p}$", "coolwarm", Normalize(0.7, 1.3)),
        (
            "rho_pressure_equation", r"$\rho_{r_p}$",
            "coolwarm", Normalize(0.7, 1.3),
        ),
        ("rho_checkerboard", r"$\rho_{cb}$", "coolwarm", Normalize(0.7, 1.3)),
        (
            "interior_gridscale_to_pressure_rms",
            r"$\|H_hp\|_{int}/\|p\|_{RMS}$",
            "viridis",
            LogNorm(1.0e-4, 1.0),
        ),
        ("iterations", "outer iterations", "magma", Normalize(0.0, max_iterations)),
    )
    fig, axes = plt.subplots(
        len(datasets), 5, figsize=(17.5, 3.5 * len(datasets)),
        squeeze=False, constrained_layout=True,
    )
    labels = tuple(dataset_label(data) for data in datasets)
    x_edges = logarithmic_edges(fo_values)
    y_edges = logarithmic_edges(tau_values)
    for row_index, (data, label) in enumerate(zip(datasets, labels)):
        for col_index, (metric, title, cmap, norm) in enumerate(metrics):
            axis = axes[row_index, col_index]
            image = surface(data, tau_values, fo_values, metric)
            shown = axis.pcolormesh(
                x_edges, y_edges, image, shading="flat", cmap=cmap, norm=norm
            )
            unresolved_x = []
            unresolved_y = []
            for tau in tau_values:
                for fo in fo_values:
                    if data[(tau, fo)]["converged"] != "True":
                        unresolved_x.append(fo)
                        unresolved_y.append(tau)
            if unresolved_x:
                axis.scatter(
                    unresolved_x, unresolved_y, marker="x", s=22,
                    linewidths=0.9, color="black",
                )
            axis.set_xscale("log")
            axis.set_yscale("log")
            axis.set_title(title)
            axis.set_xlabel(r"Fourier number $Fo$")
            if col_index == 0:
                axis.set_ylabel(
                    label + "\n" + r"Coupling strength $\tau=M/K_{vol}$"
                )
            fig.colorbar(shown, ax=axis, shrink=0.82)
    fig.savefig(output / "tauFoComparison.png", dpi=180)
    plt.close(fig)


def plot_selected_pressures(
    directories: tuple[pathlib.Path, ...],
    datasets: tuple[dict[tuple[float, float], dict[str, str]], ...],
    tau_values: list[float],
    fo_values: list[float],
    output: pathlib.Path,
) -> None:
    candidates = [
        (tau_values[0], fo_values[0]),
        (tau_values[len(tau_values) // 2], fo_values[0]),
        (tau_values[-1], fo_values[0]),
        (tau_values[-1], fo_values[-1]),
    ]
    selected = list(dict.fromkeys(candidates))
    normalized: list[list[np.ndarray]] = []
    for directory, data in zip(directories, datasets):
        fields = []
        for key in selected:
            p = final_pressure(directory / data[key]["case"]).reshape(NZ, NX)
            fields.append(
                (p - np.mean(p))
                / max(np.sqrt(np.mean(p * p)), np.finfo(float).tiny)
            )
        normalized.append(fields)
    limit = max(
        np.quantile(np.abs(field), 0.995)
        for fields in normalized
        for field in fields
    )
    fig, axes = plt.subplots(
        len(datasets), len(selected),
        figsize=(4.0 * len(selected), 2.5 * len(datasets)),
        squeeze=False, constrained_layout=True,
    )
    labels = tuple(dataset_label(data) for data in datasets)
    for row_index, (label, data, fields) in enumerate(
        zip(labels, datasets, normalized)
    ):
        for col_index, (key, field) in enumerate(zip(selected, fields)):
            axis = axes[row_index, col_index]
            shown = axis.imshow(
                field, origin="lower", extent=(0, 1, 0, 0.2), aspect="auto",
                cmap="RdBu_r", vmin=-limit, vmax=limit,
            )
            metric = float(data[key]["interior_gridscale_to_pressure_rms"])
            axis.set_title(
                rf"$\tau={key[0]:g}$, $Fo={key[1]:.0e}$"
                + "\n"
                + rf"$\|H_hp\|_{{int}}/\|p\|={metric:.3f}$"
            )
            axis.set_xlabel("x")
            if col_index == 0:
                axis.set_ylabel(f"{label}\nz")
    fig.colorbar(
        shown, ax=axes, label=r"$(p-\bar p)/\|p\|_{RMS}$", shrink=0.82
    )
    fig.savefig(output / "pressureComparison.png", dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("directories", nargs="+", type=pathlib.Path)
    parser.add_argument("--output-directory", type=pathlib.Path)
    args = parser.parse_args()
    if len(args.directories) < 2:
        parser.error("at least two sweep directories are required")
    directories = tuple(path.resolve() for path in args.directories)
    output = (
        args.output_directory or directories[0].parent / "tauFoComparison"
    ).resolve()
    output.mkdir(parents=True, exist_ok=True)
    datasets = tuple(rows(directory) for directory in directories)
    keys = set.intersection(*(set(data) for data in datasets))
    if not keys:
        parser.error("the summaries have no matched tau--Fo points")
    tau_values = sorted({key[0] for key in keys})
    fo_values = sorted({key[1] for key in keys})
    if len(keys) != len(tau_values) * len(fo_values):
        parser.error("matched points do not form a complete tau--Fo grid")
    plot_surfaces(datasets, tau_values, fo_values, output)
    plot_selected_pressures(directories, datasets, tau_values, fo_values, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
