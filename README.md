# Jakob L. Schanzer — Research Portfolio

Computational molecular science, research software, and reproducible analysis by **Jakob L. Schanzer**.

This private portfolio highlights work across statistical mechanics, molecular simulation, scientific programming, uncertainty quantification, and agent-assisted research engineering. All underlying scientific work and all final scientific decisions are mine; Codex was used under my supervision to turn that work into an auditable software portfolio.

## Computational profile

| Area | Demonstrated capabilities |
|---|---|
| **Computational methods** | Molecular dynamics, replica-exchange molecular dynamics (REMD), umbrella sampling, WHAM, MBAR, spatial thermodynamics, occupancy-fluctuation analysis, and uncertainty quantification |
| **Scientific software** | Fortran analysis programs and geometry-specific adaptations, Python post-processing and tests, Julia analysis and visualization, version control, and continuous integration |
| **HPC-scale workflows** | Parallel replica and umbrella-window designs, multi-temperature simulation campaigns, and large-scale trajectory analysis |

## Current research: water in molecular cavities

[![Host systems: CB7, OA, and beta-cyclodextrin](projects/high-energy-water/figures/figure1_host_systems.png)](projects/high-energy-water/README.md)

### Does molecular confinement create “high-energy” water?

The project tests a common mechanistic claim in molecular recognition: that confined water is uniformly energetically unfavorable and that releasing it necessarily drives binding. It combines spatial thermodynamic maps, occupancy fluctuations, hydrogen-bond statistics, temperature response, and binding-coordinate analysis for three manuscript systems: **CB7, OA, and β-cyclodextrin**.

The central result is more nuanced. Populated hydration sites are, on average, enthalpically stabilized. Sharp collective drying occurs in rigid hydrophobic cavities, but drying alone is not evidence that expelling “high-energy water” drives association.

[Explore the project, figures, physics, and code →](projects/high-energy-water/README.md)

## Agent-assisted scientific software engineering

This repository is also a documented case study in using a coding agent for scientific research engineering.

- **What Codex did:** generated the portfolio structure around my existing research materials; created and iteratively refined the checked-in Julia figure pipeline; implemented reproducible Python post-processing, deterministic tests, publication-boundary checks, and continuous integration; and drafted documentation for my review.
- **What remained human-controlled:** I designed and ran the simulations, developed or adapted the core scientific analyses, supplied and interpreted the data, defined the physical claims, selected what could be disclosed, and reviewed every result and repository change.
- **How the agent was constrained:** a machine-readable policy limits the allowed systems and data; Fortran compilation, deterministic tests, figure regeneration, and a confidentiality scan serve as validation gates; commits and pushes require human review.

The complete control loop, example task contract, and limitations are documented in the [agentic workflow](projects/high-energy-water/docs/AGENTIC_WORKFLOW.md) and [validation record](projects/high-energy-water/docs/VALIDATION.md).

## Scientific software

| Code | Analysis performed |
|---|---|
| [`janal.f`](projects/high-energy-water/src/fortran/host_analysis/janal.f) | Processes independent unbiased trajectories and energy inputs to calculate cavity-water density, energy, chemical potential/test-particle insertion estimates, occupancy fluctuations, and axial–radial thermodynamic fields. |
| [`vol.f`](projects/high-energy-water/src/fortran/host_analysis/vol.f) | Supplies the host-following cavity-membership test used by JANAL. The checked-in six-plane implementation documents one geometry-specific production backend; other host geometries used their own `vol` input and adapted backend. |
| [`hb4.f`](projects/high-energy-water/src/fortran/hydrogen_bonds/hb4.f) | Counts water–water and water–host hydrogen bonds frame by frame under periodic boundary conditions. Ensemble counts from the bulk, host–guest, isolated-host, and isolated-guest simulations supplied the four terms used to calculate the binding hydrogen-bond balance, δ. |
| [`ring_angle.f`](projects/high-energy-water/src/fortran/ring_geometry/ring_angle.f) | Quantifies β-cyclodextrin ring orientation and inter-monomer spacing from trajectory frames. This geometry analysis supported selection and validation of the open-cavity rim restraints used in the restrained β-CD calculations; the associated free-energy offset was then obtained with a separate MBAR thermodynamic cycle. |

## Computational capabilities demonstrated

- Legacy and modernized **Fortran** analysis of independent unbiased trajectories—JANAL for cavity thermodynamics plus hydrogen-bond and host-geometry tools
- Replica-exchange umbrella sampling, WHAM reconstruction, and constant-heat-capacity thermodynamic decomposition
- Reproducible **Julia** figure generation and quantitative transition analysis
- Tested **Python** post-processing with deterministic bootstrap uncertainty estimates
- Spatial density, chemical-potential, enthalpy, and hydrogen-bond maps
- Statistical-mechanical interpretation connecting microscopic hydration to binding thermodynamics
- Geometry-aware analysis through a shared JANAL core and host-specific volume backends
- CI compilation/tests plus a constrained **Codex agent workflow** with explicit validation gates

## Reproducibility

From the repository root:

```bash
make check
make test
```

Figure regeneration and data provenance are documented inside the research project. Raw trajectories are intentionally excluded because of size and unpublished-work controls; compact processed inputs required for the displayed analyses are included.

> **Confidential research note:** this repository is private and contains results from a manuscript under peer review. External access should be granted only after explicit approval from the manuscript lead/PI, and permission to view does not imply permission to redistribute.
