#!/usr/bin/python3
"""Plot stabilization-factor sensitivity at one matched tau--Fo point."""

from __future__ import annotations

import argparse
import csv
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from undrained_mandel_sweep import NX, NZ, final_pressure


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("directories", nargs="+", type=pathlib.Path)
    parser.add_argument("--output-directory", type=pathlib.Path, required=True)
    args = parser.parse_args()
    output = args.output_directory.resolve()
    output.mkdir(parents=True, exist_ok=True)

    records = []
    for directory in args.directories:
        directory = directory.resolve()
        with (directory / "summary.csv").open() as handle:
            row = next(csv.DictReader(handle))
        records.append((float(row["stab_factor"]), directory, row))
    records.sort(key=lambda item: item[0])
    tau_fo = {(float(row["tau"]), float(row["Fo"])) for _, _, row in records}
    if len(tau_fo) != 1:
        parser.error("factor directories must describe one matched tau--Fo point")

    factors = [factor for factor, _, _ in records]
    gridscale = [
        float(row["interior_gridscale_to_pressure_rms"])
        for _, _, row in records
    ]
    iterations = [float(row["iterations"]) for _, _, row in records]
    rho_update = [float(row["rho_update"]) for _, _, row in records]
    rho_cb = [float(row["rho_checkerboard"]) for _, _, row in records]

    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.4), constrained_layout=True)
    axes[0].plot(factors, gridscale, "o-", color="#0072b2")
    axes[0].set_yscale("log")
    axes[0].set_ylabel(r"$\|H_hp\|_{int}/\|p\|_{RMS}$")
    axes[1].plot(factors, iterations, "o-", color="#d55e00")
    axes[1].set_ylabel("outer iterations")
    axes[2].plot(factors, rho_update, "o-", label=r"$\rho_{\Delta p}$")
    axes[2].plot(factors, rho_cb, "s--", label=r"$\rho_{cb}$")
    axes[2].set_ylim(0.7, 1.03)
    axes[2].set_ylabel("late contraction factor")
    axes[2].legend()
    for axis in axes:
        axis.set_xlabel("stabilization factor")
        axis.grid(alpha=0.25)
    tau, fo = next(iter(tau_fo))
    fig.suptitle(rf"Factor sensitivity at $\tau={tau:g}$, $Fo={fo:.0e}$")
    fig.savefig(output / "factorSensitivity.png", dpi=180)
    plt.close(fig)

    fields = []
    for _, directory, row in records:
        p = final_pressure(directory / row["case"]).reshape(NZ, NX)
        fields.append(
            (p - np.mean(p))
            / max(np.sqrt(np.mean(p * p)), np.finfo(float).tiny)
        )
    limit = max(np.quantile(np.abs(field), 0.995) for field in fields)
    ncols = 3
    nrows = int(np.ceil(len(fields) / ncols))
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(4.1 * ncols, 2.4 * nrows),
        squeeze=False, constrained_layout=True,
    )
    for axis, (factor, _, _), field, strength in zip(
        axes.flat, records, fields, gridscale
    ):
        shown = axis.imshow(
            field, origin="lower", extent=(0, 1, 0, 0.2), aspect="auto",
            cmap="RdBu_r", vmin=-limit, vmax=limit,
        )
        axis.set_title(
            f"factor {factor:g}\n"
            + rf"$\|H_hp\|_{{int}}/\|p\|={strength:.3f}$"
        )
        axis.set_xlabel("x")
        axis.set_ylabel("z")
    for axis in axes.flat[len(fields):]:
        axis.set_visible(False)
    fig.colorbar(
        shown, ax=axes, label=r"$(p-\bar p)/\|p\|_{RMS}$", shrink=0.82
    )
    fig.savefig(output / "factorPressureFields.png", dpi=180)
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
