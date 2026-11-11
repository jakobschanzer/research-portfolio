#!/usr/bin/env python3
"""Fail when the curated project contains unreviewed data or private paths."""

from __future__ import annotations

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
APPROVED_DATA_FILES = {
    "data/processed/pmf_thermodynamics.xlsx",
    "data/processed/direct_residual_decomposition.xlsx",
    "data/processed/spatial_maps.xlsx",
    "data/processed/occupancy_along_binding_coordinate.xlsx",
    "data/processed/mean_occupancy_vs_z.xlsx",
    "data/processed/cb7_occupancy_rerun.xlsx",
    "data/processed/hydrogen_bond_temperature.xlsx",
    "data/processed/unbiased_occupancy.xlsx",
    "data/derived/water_occupancy_by_host.xlsx",
    "data/derived/betacd_cluster_bootstrap_sem.tsv",
    "data/derived/drying_transition_analysis.xlsx",
    "data/derived/occupancy_moments.xlsx",
    "data/derived/pn_fattail_fits.xlsx",
}
REQUIRED_DATA_FILES = APPROVED_DATA_FILES - {
    "data/derived/drying_transition_analysis.xlsx",
    "data/derived/occupancy_moments.xlsx",
    "data/derived/pn_fattail_fits.xlsx",
}
TEXT_SUFFIXES = {".f", ".jl", ".md", ".py", ".sh", ".toml", ".txt", ".yml"}
PRIVATE_MARKERS = ("/Users/", "Box-Box", "Miscellaneous disorganization")


def main() -> None:
    failures: list[str] = []
    data_files = {
        path.relative_to(PROJECT_ROOT).as_posix()
        for path in (PROJECT_ROOT / "data").rglob("*")
        if path.is_file()
    }
    unexpected = sorted(data_files - APPROVED_DATA_FILES)
    missing = sorted(REQUIRED_DATA_FILES - data_files)
    if unexpected:
        failures.append(f"unreviewed data files: {unexpected}")
    if missing:
        failures.append(f"approved data files missing: {missing}")

    raw_dir = PROJECT_ROOT / "data" / "raw"
    if raw_dir.exists() and any(path.is_file() for path in raw_dir.rglob("*")):
        failures.append("raw simulation data must not be committed")

    for path in PROJECT_ROOT.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        if path.resolve() == Path(__file__).resolve():
            continue
        text = path.read_text(encoding="utf-8", errors="ignore")
        found = [marker for marker in PRIVATE_MARKERS if marker in text]
        if found:
            failures.append(f"private path marker in {path.relative_to(PROJECT_ROOT)}: {found}")

    if failures:
        raise SystemExit("Publication-boundary check failed:\n- " + "\n- ".join(failures))
    print(f"Publication boundary OK: {len(data_files)} reviewed data files; no private paths.")


if __name__ == "__main__":
    main()
