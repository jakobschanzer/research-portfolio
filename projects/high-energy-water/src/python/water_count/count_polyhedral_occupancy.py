#!/usr/bin/env python3
"""Framewise water occupancy and block SEM for the six-plane host cavity."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import minimize_vectors


TEMPERATURES = tuple(range(5, 100, 10))


@dataclass(frozen=True)
class VolumeDefinition:
    checks: np.ndarray          # (6,), zero-based host indices
    planes: np.ndarray          # (6, 3), zero-based host indices


def leading_integers(line: str) -> list[int]:
    values = []
    for field in line.split():
        try:
            values.append(int(field))
        except ValueError:
            break
    return values


def parse_top_bind(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read the two orientation groups used to center the host in janal.f."""
    tokens = []
    for line in path.read_text().splitlines():
        tokens.extend(leading_integers(line))
    pos = 0

    def take() -> int:
        nonlocal pos
        if pos >= len(tokens):
            raise ValueError(f"Incomplete top-bind file: {path}")
        value = tokens[pos]
        pos += 1
        return value

    if take() < 1:
        raise ValueError("top-bind does not define top/bottom orientation groups")
    groups = []
    for _ in range(2):
        count = take()
        groups.append(np.asarray([take() - 1 for _ in range(count)], dtype=int))
    return groups[0], groups[1]


def parse_vol(path: Path) -> VolumeDefinition:
    """Read [checks] and [planes] exactly as vol_init in vol.f does."""
    sections: dict[str, list[int]] = {}
    current = None
    for raw in path.read_text().splitlines():
        stripped = raw.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            current = stripped[1:-1].strip().lower()
            sections[current] = []
        elif current is not None:
            sections[current].extend(leading_integers(raw))
    checks = sections.get("checks", [])
    planes = sections.get("planes", [])
    if len(checks) != 6 or len(planes) != 18:
        raise ValueError(
            f"{path}: expected 6 check indices and 18 plane indices; "
            f"found {len(checks)} and {len(planes)}"
        )
    return VolumeDefinition(
        checks=np.asarray(checks, dtype=int) - 1,
        planes=(np.asarray(planes, dtype=int) - 1).reshape(6, 3),
    )


def unwrap_host(positions: np.ndarray, box: np.ndarray, reference_index: int) -> np.ndarray:
    reference = positions[reference_index]
    return reference + minimize_vectors(positions - reference, box=box)


def build_planes(
    centered_host: np.ndarray, definition: VolumeDefinition
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Reproduce vol_update_frame: cross(x2-x1, x3-x1), offset, check side."""
    normals = np.empty((6, 3), dtype=float)
    offsets = np.empty(6, dtype=float)
    check_signs = np.empty(6, dtype=float)
    for i, (a1, a2, a3) in enumerate(definition.planes):
        x1 = centered_host[a1]
        # Same cross-product order as vol.f: v2 × v1.
        normal = np.cross(centered_host[a2] - x1, centered_host[a3] - x1)
        magnitude = np.linalg.norm(normal)
        if magnitude == 0:
            raise ValueError(f"Plane {i + 1} is degenerate in this frame")
        normal /= magnitude
        offset = -float(normal @ x1)
        check_value = float(normal @ centered_host[definition.checks[i]] + offset)
        normals[i] = normal
        offsets[i] = offset
        check_signs[i] = 1.0 if check_value >= 0.0 else -1.0
    return normals, offsets, check_signs


def count_waters(
    universe: mda.Universe,
    host_count: int,
    oxygen_indices: np.ndarray,
    orientation_indices: np.ndarray,
    definition: VolumeDefinition,
) -> int:
    box = universe.dimensions
    if box is None or np.any(np.asarray(box[:3]) <= 0):
        raise ValueError("Frame has no valid periodic box")

    raw_host = universe.atoms[:host_count].positions
    host = unwrap_host(raw_host, box, int(orientation_indices[0]))
    host_center = host[orientation_indices].mean(axis=0)
    centered_host = host - host_center
    normals, offsets, check_signs = build_planes(centered_host, definition)

    water_positions = minimize_vectors(
        universe.atoms[oxygen_indices].positions - host_center, box=box
    )
    signed_plane_values = water_positions @ normals.T + offsets
    # vol_inside rejects only s*scheck < 0, so points on a plane are included.
    inside = np.all(signed_plane_values * check_signs >= 0.0, axis=1)
    return int(np.count_nonzero(inside))


def block_statistics(values: np.ndarray, blocks: int):
    if blocks < 2 or len(values) < blocks:
        raise ValueError("Need at least two blocks and at least one frame per block")
    frames_per_block = len(values) // blocks
    used = frames_per_block * blocks
    means = values[:used].reshape(blocks, frames_per_block).mean(axis=1)
    sem = float(means.std(ddof=1) / math.sqrt(blocks))
    return means, sem, frames_per_block, len(values) - used


def write_tsv(path: Path, headers: list[str], rows) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(headers)
        writer.writerows(rows)


def parse_temperatures(text: str):
    return tuple(int(value.strip()) for value in text.split(",") if value.strip())


def options():
    parser = argparse.ArgumentParser(
        description="Count waters inside the six-plane cavity defined by vol.f/vol."
    )
    parser.add_argument("--temperatures", type=parse_temperatures, default=TEMPERATURES)
    parser.add_argument("--directory-template", default="../{temperature}")
    parser.add_argument("--topology", default="prd.tpr")
    parser.add_argument("--trajectory", default="prd.xtc")
    parser.add_argument("--top-bind", type=Path, default=Path("../top-bind"))
    parser.add_argument("--vol", type=Path, default=Path("../vol"))
    parser.add_argument("--host-atoms", type=int, default=184)
    parser.add_argument("--waters", type=int, default=2000)
    parser.add_argument("--sites-per-water", type=int, default=4)
    parser.add_argument("--counterions", type=int, default=8)
    parser.add_argument("--oxygen-offset", type=int, default=0)
    parser.add_argument("--blocks", type=int, default=20)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--output-directory", type=Path, default=Path("water_count_results"))
    return parser.parse_args()


def main() -> int:
    args = options()
    group1, group2 = parse_top_bind(args.top_bind.resolve())
    orientation = np.unique(np.concatenate((group1, group2)))
    volume = parse_vol(args.vol.resolve())
    output = args.output_directory.resolve()
    output.mkdir(parents=True, exist_ok=True)

    oxygen_indices = (
        args.host_atoms
        + args.oxygen_offset
        + np.arange(args.waters) * args.sites_per_water
    )
    expected_atoms = (
        args.host_atoms + args.waters * args.sites_per_water + args.counterions
    )
    referenced = np.concatenate((orientation, volume.checks, volume.planes.ravel()))
    if np.any(referenced < 0) or np.any(referenced >= args.host_atoms):
        raise IndexError("top-bind or vol references an atom outside the host")

    print("Volume: intersection of 6 per-frame host planes")
    print("Check atoms (1-based):", " ".join(str(i + 1) for i in volume.checks))
    summary = []
    for temperature in args.temperatures:
        directory = Path(
            args.directory_template.format(temperature=temperature, temp=temperature)
        ).resolve()
        topology, trajectory = directory / args.topology, directory / args.trajectory
        if not topology.is_file() or not trajectory.is_file():
            raise FileNotFoundError(f"T={temperature}: expected {topology} and {trajectory}")

        print(f"T={temperature} C: reading {trajectory}", flush=True)
        universe = mda.Universe(str(topology), str(trajectory))
        if len(universe.atoms) != expected_atoms:
            raise ValueError(
                f"T={temperature}: found {len(universe.atoms)} atoms; expected "
                f"{args.host_atoms} + {args.waters}*{args.sites_per_water} + "
                f"{args.counterions} = {expected_atoms}"
            )

        counts, frame_rows = [], []
        for ts in universe.trajectory[:: args.stride]:
            value = count_waters(
                universe, args.host_atoms, oxygen_indices, orientation, volume
            )
            counts.append(value)
            frame_rows.append((temperature, ts.frame, f"{ts.time:.8g}", value))

        values = np.asarray(counts, dtype=float)
        block_means, sem, frames_per_block, discarded = block_statistics(values, args.blocks)
        mean, sd = float(values.mean()), float(values.std(ddof=1))
        tag = str(temperature)
        write_tsv(
            output / f"{tag}_frame_counts.tsv",
            ["temperature_C", "frame", "time_ps", "water_count"], frame_rows,
        )
        write_tsv(
            output / f"{tag}_block_means.tsv",
            ["temperature_C", "block", "first_analyzed_frame", "last_analyzed_frame", "mean_N"],
            [(temperature, i + 1, i * frames_per_block,
              (i + 1) * frames_per_block - 1, f"{x:.10g}")
             for i, x in enumerate(block_means)],
        )
        n_values, frequencies = np.unique(values.astype(int), return_counts=True)
        write_tsv(
            output / f"{tag}_Pn.tsv", ["N", "count", "P_N"],
            [(int(n), int(f), f"{f / len(values):.10g}")
             for n, f in zip(n_values, frequencies)],
        )
        summary.append((temperature, f"{mean:.10g}", f"{sd:.10g}", f"{sem:.10g}",
                        f"{block_means.std(ddof=1):.10g}", len(values), args.stride,
                        args.blocks, frames_per_block, discarded))
        print(f"T={temperature} C: <N>={mean:.6f}, block SEM={sem:.6f}", flush=True)

    write_tsv(
        output / "water_count_summary_by_temperature.tsv",
        ["temperature_C", "mean_N", "instantaneous_SD_N", "block_SEM_N",
         "block_SD_N", "analyzed_frames", "stride", "blocks",
         "frames_per_block", "discarded_tail_frames"], summary,
    )
    print("Complete:", output / "water_count_summary_by_temperature.tsv")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise
