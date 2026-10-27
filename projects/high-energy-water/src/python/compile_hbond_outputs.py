#!/usr/bin/env python3
"""Compile temperature-resolved hydrogen-bond output reproducibly.

The Fortran hydrogen-bond executable writes one ``dist.<temperature>`` file
per simulation temperature. This command validates those files, aggregates the
three reported hydrogen-bond channels, and estimates uncertainty by seeded
non-parametric bootstrap resampling.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np


EXPECTED_TEMPERATURES = (5, 15, 25, 35, 45, 55, 65, 75, 85, 95)


def bootstrap_error(
    data: np.ndarray,
    rng: np.random.Generator,
    iterations: int,
) -> float:
    """Return the bootstrap standard deviation of the sample mean."""
    if data.size == 0:
        raise ValueError("cannot bootstrap an empty series")
    indices = rng.integers(0, data.size, size=(iterations, data.size))
    return float(np.std(data[indices].mean(axis=1), ddof=1))


def read_distribution(path: Path) -> tuple[np.ndarray, ...]:
    """Read step, water-water, and the two water-host series from a file."""
    rows: list[tuple[int, float, float, float]] = []
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            stripped = line.strip()
            if not stripped or stripped.startswith("-"):
                continue
            parts = stripped.split()
            if len(parts) < 4:
                continue
            try:
                rows.append((int(parts[0]), *(float(value) for value in parts[1:4])))
            except ValueError:
                # The first non-data row is the Fortran header.
                continue
    if not rows:
        raise ValueError(f"no numeric hydrogen-bond rows found in {path}")
    columns = np.asarray(rows, dtype=float).T
    return columns[0].astype(int), columns[1], columns[2], columns[3]


def compile_outputs(
    input_dir: Path,
    output_dir: Path,
    seed: int = 20261027,
    iterations: int = 1_000,
) -> tuple[Path, Path]:
    """Compile expected temperature files and return the two output paths."""
    missing = [
        temp
        for temp in EXPECTED_TEMPERATURES
        if not (input_dir / f"dist.{temp}").is_file()
    ]
    if missing:
        raise FileNotFoundError(f"missing temperature files: {missing}")

    rng = np.random.default_rng(seed)
    step_set: set[int] = set()
    raw_data: dict[str, dict[int, float]] = {}
    average_row: dict[str, object] = {"Step": "AVG"}
    error_row: dict[str, object] = {"Step": "ERR"}
    summary_rows: list[list[float]] = []

    for temperature in EXPECTED_TEMPERATURES:
        steps, water_water, host_top, host_bottom = read_distribution(
            input_dir / f"dist.{temperature}"
        )
        step_set.update(int(step) for step in steps)
        series = {
            "HB_WW": water_water,
            "HB_WH_top": host_top,
            "HB_WH_bottom": host_bottom,
        }
        errors: dict[str, float] = {}
        for label, values in series.items():
            key = f"{temperature}_{label}"
            raw_data[key] = dict(zip(steps, values, strict=True))
            average_row[key] = float(np.mean(values))
            errors[label] = bootstrap_error(values, rng, iterations)
            error_row[key] = errors[label]

        summary_rows.append(
            [
                float(temperature),
                float(np.mean(water_water)), errors["HB_WW"],
                float(np.mean(host_top)), errors["HB_WH_top"],
                float(np.mean(host_bottom)), errors["HB_WH_bottom"],
            ]
        )

    wide_rows: list[dict[str, object]] = [average_row, error_row]
    for step in sorted(step_set):
        row: dict[str, object] = {"Step": step}
        for key, values in raw_data.items():
            row[key] = values.get(step, "")
        wide_rows.append(row)

    output_dir.mkdir(parents=True, exist_ok=True)
    wide_path = output_dir / "hbond_summary_sideways.csv"
    compact_path = output_dir / "hbond_summary_compact.csv"
    with wide_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(average_row))
        writer.writeheader()
        writer.writerows(wide_rows)

    compact_header = [
        "Temp_C",
        "HB_WW", "Bootstrap_SE_WW",
        "HB_WH_top", "Bootstrap_SE_WH_top",
        "HB_WH_bottom", "Bootstrap_SE_WH_bottom",
    ]
    with compact_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(compact_header)
        writer.writerows(summary_rows)
    return wide_path, compact_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=Path.cwd())
    parser.add_argument("--output-dir", type=Path, default=Path.cwd())
    parser.add_argument("--seed", type=int, default=20261027)
    parser.add_argument("--iterations", type=int, default=1_000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = compile_outputs(args.input_dir, args.output_dir, args.seed, args.iterations)
    print("Wrote:")
    for output in outputs:
        print(f"  {output}")


if __name__ == "__main__":
    main()
