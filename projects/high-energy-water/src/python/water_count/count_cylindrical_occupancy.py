#!/usr/bin/env python3
"""Count pocket waters and obtain block-averaged uncertainty across temperatures.

Run this script from a new ``water`` analysis directory. By default it expects
temperature trajectories in ../T5, ../T15, ..., ../T95 and the original Fortran
geometry file at ../top-bind. Paths are configurable from the command line.

Geometry follows hanal.f:
  1. Read two host atom groups and the cylindrical bin limits from top-bind.
  2. Unwrap those host atoms across periodic boundaries and center their union.
  3. Define the host axis from the two group centers.
  4. Count water oxygens satisfying z_lower < z < z_upper and r < radius.

MDAnalysis supplies coordinates in angstrom, so the top-bind geometry and all
reported distances are in angstrom.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import minimize_vectors


DEFAULT_TEMPERATURES = tuple(range(5, 100, 10))
DR_ANGSTROM = 0.25
GRID_HALF_LENGTH_ANGSTROM = 13.5


@dataclass(frozen=True)
class PocketGeometry:
    group1_zero_based: np.ndarray
    group2_zero_based: np.ndarray
    z_lower_angstrom: float
    z_upper_angstrom: float
    radius_angstrom: float


def parse_top_bind(path: Path) -> PocketGeometry:
    """Parse the fields consumed by the nbt=1 branch of hanal.f."""
    tokens: list[int] = []
    for raw_line in path.read_text().splitlines():
        # Fortran list-directed reads consume the leading numeric field(s) and
        # ignore the descriptive text that follows. Do not extract digits from
        # comments such as "(1 = top/bottom, 2 = cross)".
        for field in raw_line.split():
            try:
                tokens.append(int(field))
            except ValueError:
                break

    position = 0

    def take_one(label: str) -> int:
        nonlocal position
        if position >= len(tokens):
            raise ValueError(f"{path}: missing {label}")
        value = tokens[position]
        position += 1
        return value

    nbt = take_one("nbt")
    if nbt < 1:
        raise ValueError(f"{path}: nbt={nbt}; the Fortran cylinder requires two groups")

    groups: list[np.ndarray] = []
    for group_number in (1, 2):
        count = take_one(f"number of atoms in group {group_number}")
        indices = [take_one(f"atom in group {group_number}") for _ in range(count)]
        if not indices:
            raise ValueError(f"{path}: group {group_number} is empty")
        groups.append(np.asarray(indices, dtype=int) - 1)  # Fortran -> Python indexing

    _icf = take_one("icf")  # stochastic centering trials; unnecessary with explicit PBC unwrap
    izl = take_one("izl")
    izh = take_one("izh")
    irm = take_one("irm")

    return PocketGeometry(
        group1_zero_based=groups[0],
        group2_zero_based=groups[1],
        z_lower_angstrom=DR_ANGSTROM * (izl - 1) - GRID_HALF_LENGTH_ANGSTROM,
        z_upper_angstrom=DR_ANGSTROM * izh - GRID_HALF_LENGTH_ANGSTROM,
        radius_angstrom=DR_ANGSTROM * irm,
    )


def unwrap_relative(positions: np.ndarray, box: np.ndarray) -> np.ndarray:
    """Place a compact host selection in one periodic image around atom zero."""
    reference = positions[0]
    return reference + minimize_vectors(positions - reference, box=box)


def pocket_count(
    universe: mda.Universe,
    oxygen_indices: np.ndarray,
    geometry: PocketGeometry,
) -> int:
    """Count waters using the same cylinder coordinates and strict tests as hanal.f."""
    box = universe.dimensions
    if box is None or len(box) < 6 or np.any(np.asarray(box[:3]) <= 0):
        raise ValueError("Trajectory frame has no valid periodic box dimensions")

    orientation_indices = np.unique(
        np.concatenate((geometry.group1_zero_based, geometry.group2_zero_based))
    )
    host_positions = unwrap_relative(universe.atoms[orientation_indices].positions, box)
    index_to_local = {atom_index: local for local, atom_index in enumerate(orientation_indices)}
    group1_local = np.asarray([index_to_local[i] for i in geometry.group1_zero_based])
    group2_local = np.asarray([index_to_local[i] for i in geometry.group2_zero_based])

    host_center = host_positions.mean(axis=0)
    center1 = host_positions[group1_local].mean(axis=0)
    center2 = host_positions[group2_local].mean(axis=0)
    axis = center2 - center1
    axis_norm = np.linalg.norm(axis)
    if axis_norm == 0:
        raise ValueError("The two top-bind group centers coincide; host axis is undefined")
    axis /= axis_norm

    displacement = minimize_vectors(
        universe.atoms[oxygen_indices].positions - host_center,
        box=box,
    )
    axial = displacement @ axis
    radial_squared = np.einsum("ij,ij->i", displacement, displacement) - axial**2
    radial_squared = np.maximum(radial_squared, 0.0)  # suppress roundoff below zero

    inside = (
        (axial > geometry.z_lower_angstrom)
        & (axial < geometry.z_upper_angstrom)
        & (radial_squared < geometry.radius_angstrom**2)
    )
    return int(np.count_nonzero(inside))


def make_blocks(values: np.ndarray, requested_blocks: int) -> tuple[np.ndarray, int, int]:
    """Return equal contiguous block means, frames/block, and discarded tail size."""
    if requested_blocks < 2:
        raise ValueError("At least two blocks are required to calculate a block SEM")
    if len(values) < requested_blocks:
        raise ValueError(
            f"Only {len(values)} frames are available for {requested_blocks} blocks"
        )
    frames_per_block = len(values) // requested_blocks
    used = frames_per_block * requested_blocks
    block_means = values[:used].reshape(requested_blocks, frames_per_block).mean(axis=1)
    return block_means, frames_per_block, len(values) - used


def write_rows(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_temperatures(text: str) -> tuple[int, ...]:
    values = tuple(int(item.strip()) for item in text.split(",") if item.strip())
    if not values:
        raise argparse.ArgumentTypeError("temperature list is empty")
    return values


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Count unbiased pocket waters and calculate a block-averaged SEM."
    )
    parser.add_argument(
        "--temperatures",
        type=parse_temperatures,
        default=DEFAULT_TEMPERATURES,
        help="comma-separated Celsius temperatures (default: 5,15,...,95)",
    )
    parser.add_argument(
        "--directory-template",
        default="../{temperature}",
        help="temperature directory template (default: ../{temperature})",
    )
    parser.add_argument("--topology", default="prd.tpr", help="topology name in each directory")
    parser.add_argument("--trajectory", default="prd.xtc", help="trajectory name in each directory")
    parser.add_argument(
        "--top-bind",
        type=Path,
        default=Path("../top-bind"),
        help="original hanal.f top-bind file (default: ../top-bind)",
    )
    parser.add_argument("--host-atoms", type=int, default=126)
    parser.add_argument("--waters", type=int, default=2000)
    parser.add_argument("--sites-per-water", type=int, default=4)
    parser.add_argument(
        "--oxygen-offset",
        type=int,
        default=0,
        help="oxygen offset within each four-site water, zero based (default: 0)",
    )
    parser.add_argument(
        "--blocks",
        type=int,
        default=20,
        help="number of contiguous equal blocks per trajectory (default: 20)",
    )
    parser.add_argument(
        "--stride", type=int, default=1, help="analyze every STRIDE-th stored frame"
    )
    parser.add_argument(
        "--output-directory", type=Path, default=Path("water_count_results")
    )
    return parser.parse_args()


def main() -> int:
    args = arguments()
    geometry = parse_top_bind(args.top_bind.resolve())
    output_directory = args.output_directory.resolve()
    output_directory.mkdir(parents=True, exist_ok=True)

    first_oxygen = args.host_atoms + args.oxygen_offset
    oxygen_indices = first_oxygen + np.arange(args.waters) * args.sites_per_water
    summary_rows: list[dict] = []

    print(
        "Geometry from top-bind: "
        f"r < {geometry.radius_angstrom:.3f} A, "
        f"{geometry.z_lower_angstrom:.3f} < z < {geometry.z_upper_angstrom:.3f} A, "
        f"volume = {math.pi * geometry.radius_angstrom**2 * (geometry.z_upper_angstrom - geometry.z_lower_angstrom):.6f} A^3"
    )

    for temperature in args.temperatures:
        directory = Path(
            args.directory_template.format(temperature=temperature, temp=temperature)
        ).resolve()
        topology = directory / args.topology
        trajectory = directory / args.trajectory
        if not topology.is_file() or not trajectory.is_file():
            raise FileNotFoundError(
                f"T={temperature} C: expected {topology} and {trajectory}"
            )

        print(f"T={temperature} C: reading {trajectory}", flush=True)
        universe = mda.Universe(str(topology), str(trajectory))
        if oxygen_indices[-1] >= len(universe.atoms):
            raise IndexError(
                f"T={temperature} C has {len(universe.atoms)} atoms, but the expected last "
                f"water oxygen has zero-based index {oxygen_indices[-1]}"
            )
        if max(np.max(geometry.group1_zero_based), np.max(geometry.group2_zero_based)) >= args.host_atoms:
            raise IndexError("top-bind orientation group contains an atom beyond --host-atoms")

        counts: list[int] = []
        frame_rows: list[dict] = []
        for ts in universe.trajectory[:: args.stride]:
            count = pocket_count(universe, oxygen_indices, geometry)
            counts.append(count)
            frame_rows.append(
                {
                    "temperature_C": temperature,
                    "frame": ts.frame,
                    "time_ps": f"{ts.time:.8g}",
                    "water_count": count,
                }
            )

        count_array = np.asarray(counts, dtype=float)
        block_means, frames_per_block, discarded = make_blocks(count_array, args.blocks)
        mean_count = float(count_array.mean())
        instantaneous_sd = float(count_array.std(ddof=1))
        block_sd = float(block_means.std(ddof=1))
        block_sem = block_sd / math.sqrt(len(block_means))

        tag = f"T{temperature}"
        write_rows(
            output_directory / f"{tag}_frame_counts.tsv",
            ["temperature_C", "frame", "time_ps", "water_count"],
            frame_rows,
        )
        write_rows(
            output_directory / f"{tag}_block_means.tsv",
            ["temperature_C", "block", "first_analyzed_frame", "last_analyzed_frame", "mean_N"],
            [
                {
                    "temperature_C": temperature,
                    "block": block + 1,
                    "first_analyzed_frame": block * frames_per_block,
                    "last_analyzed_frame": (block + 1) * frames_per_block - 1,
                    "mean_N": f"{value:.10g}",
                }
                for block, value in enumerate(block_means)
            ],
        )

        unique, frequencies = np.unique(count_array.astype(int), return_counts=True)
        write_rows(
            output_directory / f"{tag}_Pn.tsv",
            ["N", "count", "P_N"],
            [
                {"N": int(n), "count": int(frequency), "P_N": f"{frequency / len(count_array):.10g}"}
                for n, frequency in zip(unique, frequencies)
            ],
        )

        summary_rows.append(
            {
                "temperature_C": temperature,
                "mean_N": f"{mean_count:.10g}",
                "instantaneous_SD_N": f"{instantaneous_sd:.10g}",
                "block_SEM_N": f"{block_sem:.10g}",
                "block_SD_N": f"{block_sd:.10g}",
                "analyzed_frames": len(count_array),
                "stride": args.stride,
                "blocks": len(block_means),
                "frames_per_block": frames_per_block,
                "discarded_tail_frames": discarded,
            }
        )
        print(
            f"T={temperature} C: <N>={mean_count:.6f}, "
            f"block SEM={block_sem:.6f} ({len(block_means)} blocks)",
            flush=True,
        )

    fields = [
        "temperature_C",
        "mean_N",
        "instantaneous_SD_N",
        "block_SEM_N",
        "block_SD_N",
        "analyzed_frames",
        "stride",
        "blocks",
        "frames_per_block",
        "discarded_tail_frames",
    ]
    summary_path = output_directory / "water_count_summary_by_temperature.tsv"
    write_rows(summary_path, fields, summary_rows)
    print(f"Complete. Summary: {summary_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise
