# Computational Molecular Science Portfolio

Research software and reproducible analysis by **Jakob L. Schanzer**.

This private portfolio highlights scientific computing across statistical mechanics, molecular simulation analysis, Fortran, Python, Julia, uncertainty quantification, and agent-assisted research engineering.

## Featured project: water in molecular cavities

[![Host systems: CB7, OA, and beta-cyclodextrin](projects/high-energy-water/figures/host_systems.png)](projects/high-energy-water/README.md)

### Does molecular confinement create “high-energy” water?

The project tests a common mechanistic claim in molecular recognition: that confined water is uniformly energetically unfavorable and that releasing it necessarily drives binding. It combines spatial thermodynamic maps, occupancy fluctuations, hydrogen-bond statistics, temperature response, and binding-coordinate analysis for three manuscript systems: **CB7, OA, and β-cyclodextrin**.

The central result is more nuanced. Populated hydration sites are, on average, enthalpically stabilized. Sharp collective drying occurs in rigid hydrophobic cavities, but drying alone is not evidence that expelling “high-energy water” drives association.

[Explore the project, figures, physics, and code →](projects/high-energy-water/README.md)

## Computational capabilities demonstrated

- Legacy and modernized **Fortran** analysis for cavity thermodynamics, hydrogen bonding, and host geometry
- Reproducible **Julia** figure generation and quantitative transition analysis
- Tested **Python** post-processing with deterministic bootstrap uncertainty estimates
- Spatial density, chemical-potential, enthalpy, and hydrogen-bond maps
- Statistical-mechanical interpretation connecting microscopic hydration to binding thermodynamics
- Geometry-aware analysis through a shared JANAL core and host-specific volume backends
- CI compilation/tests plus a constrained **Codex agent workflow** with explicit validation gates

## Portfolio roadmap

| Project | Focus | Status |
|---|---|---|
| [High-energy water](projects/high-energy-water/README.md) | Confined-water thermodynamics and collective drying | Featured; manuscript in preparation |
| Project 2 | Multiscale molecular modeling | Planned |
| Project 3 | Agentic scientific computing | Planned |

## Reproducibility

From the repository root:

```bash
make check
make test
```

Figure regeneration and data provenance are documented inside the featured project. Raw trajectories are intentionally excluded because of size and unpublished-work controls; compact processed inputs required for the displayed analyses are included.

> **Confidential research note:** this repository is private and contains results from a manuscript in preparation. Do not redistribute without the author’s permission.
