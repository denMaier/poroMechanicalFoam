#!/usr/bin/env python3
"""Extend only unresolved matched-scale sweep points."""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import os
import pathlib

from plot_matched_scale import GROWING, combined_rate_summary
from tau_fo_sweep import SweepPoint, run_point


CONFIGURATIONS = (
    (
        "unstabilized",
        "tauFoRunsUnstabilizedScale2p5",
        "unstabilized",
    ),
    (
        "implicit",
        "tauFoRunsImplicitScale2p5",
        "implicit",
    ),
)


def selected_points(
    case_directory: pathlib.Path,
    source_name: str,
) -> list[tuple[SweepPoint, dict[str, str]]]:
    source_directory = case_directory / source_name
    with (source_directory / "tauFoSweep.csv").open() as stream:
        rows = list(csv.DictReader(stream))

    selected = []
    for row in rows:
        if int(row["converged"]):
            continue
        log_text = (
            source_directory / row["case"] / "log.poroMechanicalFoam"
        ).read_text(errors="replace")
        summary = combined_rate_summary(log_text)
        # The strongly checkerboard-growing region is already unambiguous at
        # 300 iterations. Spend the extension budget on all other capped
        # points, including apparently contracting points that merely need
        # more iterations to reach the absolute tolerances.
        if (
            int(summary["fullState"]) == GROWING
            and int(summary["checkerboardState"]) == GROWING
        ):
            continue
        case_tokens = row["case"].split("_")
        point = SweepPoint(
            tau_index=int(case_tokens[1]),
            fo_index=int(case_tokens[3]),
            tau=float(row["tau"]),
            fo=float(row["Fo"]),
            storage=float(row["Ss"]),
            conductivity=float(row["conductivity"]),
            case_name=row["case"],
        )
        metadata = dict(row)
        metadata["originalFullRate"] = str(summary["fullRate"])
        metadata["originalFullRateLow"] = str(summary["fullLow"])
        metadata["originalFullRateHigh"] = str(summary["fullHigh"])
        metadata["originalFullState"] = str(summary["fullState"])
        metadata["originalCheckerboardRate"] = str(summary["checkerboardRate"])
        metadata["originalCheckerboardState"] = str(summary["checkerboardState"])
        selected.append((point, metadata))
    return selected


def completed_log(path: pathlib.Path) -> bool:
    if not path.exists():
        return False
    return "\nEnd\n" in path.read_text(errors="replace")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-iterations", type=int, default=1200)
    parser.add_argument("--jobs", type=int, default=min(4, os.cpu_count() or 1))
    parser.add_argument("--output-directory", default="matchedScale2p5Extended")
    parser.add_argument(
        "--list-only",
        action="store_true",
        help="write the selection manifest without running cases",
    )
    args = parser.parse_args()
    if args.max_iterations <= 300:
        parser.error("--max-iterations must exceed the original cap of 300")
    if args.jobs <= 0:
        parser.error("--jobs must be positive")

    case_directory = pathlib.Path(__file__).resolve().parent
    base = case_directory / "base"
    output = case_directory / args.output_directory
    output.mkdir(parents=True, exist_ok=True)

    tasks = []
    manifest_rows = []
    for output_name, source_name, formulation in CONFIGURATIONS:
        runs = output / output_name
        runs.mkdir(parents=True, exist_ok=True)
        for point, metadata in selected_points(case_directory, source_name):
            metadata.update(
                {
                    "formulation": formulation,
                    "extendedMaxIterations": args.max_iterations,
                }
            )
            manifest_rows.append(metadata)
            log_path = runs / point.case_name / "log.poroMechanicalFoam"
            if not completed_log(log_path):
                tasks.append((runs, point, formulation))

    manifest_path = output / "extendedManifest.csv"
    with manifest_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(manifest_rows[0]))
        writer.writeheader()
        writer.writerows(manifest_rows)

    print(
        f"selected {len(manifest_rows)} points; "
        f"{len(tasks)} require execution; manifest: {manifest_path}",
        flush=True,
    )
    if args.list_only:
        return

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(
                run_point,
                base,
                runs,
                point,
                args.max_iterations,
                formulation,
                2.5,
                1.0,
            ): (runs.name, point.case_name)
            for runs, point, formulation in tasks
        }
        for completed, future in enumerate(
            concurrent.futures.as_completed(futures), start=1
        ):
            formulation, case_name = futures[future]
            _, return_code = future.result()
            print(
                f"[{completed}/{len(tasks)}] {formulation}/{case_name}: "
                f"exit={return_code}",
                flush=True,
            )


if __name__ == "__main__":
    main()
