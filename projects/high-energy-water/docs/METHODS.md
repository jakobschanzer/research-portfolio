# Methods overview

## Scientific systems

The curated portfolio covers the three host families presented in the manuscript: cucurbit[7]uril (CB7), octa-acid (OA), and β-cyclodextrin (β-CD). The β-CD analysis includes unrestrained and restrained protocols where the distinction is scientifically relevant.

## JANAL thermodynamic analysis

The Fortran JANAL program processes aligned molecular-simulation frames and energy inputs to accumulate cavity-resolved observables. Its outputs include water density, energy, chemical potential/test-particle insertion estimates, occupancy distributions, fluctuations, and hydrogen-bond statistics on axial–radial grids.

The analysis core is shared across systems. Cavity membership is delegated to a geometry-specific volume backend. This separation matters: a curved or polyhedral host cavity cannot be represented faithfully by a single universal boundary definition.

## Spatial maps

The manuscript maps water concentration, average hydrogen-bond coordination, and local enthalpy as functions of axial position `z` and radial position `r`. Display interpolation smooths the rendered grids but does not add physical resolution. Hydrogen-bond maps are expressed relative to a bulk reference at 25 °C.

## Binding thermodynamics

The potential of mean force is reported along the host–guest separation coordinate and decomposed into direct host–guest and solvent-mediated contributions. Enthalpic and entropic terms are interpreted together; a favorable binding free energy does not imply that each microscopic contribution is favorable.

## Occupancy and drying

The conditional distribution `P(N|z)` gives the probability of observing `N` cavity waters at guest position `z`. Mean occupancy and fluctuations are computed from that distribution. A 25–75% transition width is obtained by linear interpolation of the fraction-empty curve between adjacent 0.4 Å bins. “Sharp” and “gradual” refer to equilibrium change along the binding coordinate, not a kinetic rate.

## Uncertainty

Trajectory-derived standard errors use block or cluster bootstrap estimates where available. The temperature-resolved hydrogen-bond compiler uses a fixed random seed for reproducibility. The drying-transition uncertainty propagation uses the available pointwise standard errors; because cross-bin covariance was unavailable, this approximation is stated explicitly in the derived workbook.

## Interpretation boundary

The calculations support three claims:

1. Populated hydration sites can be enthalpically stabilized relative to bulk.
2. Rigid hydrophobic cavities can display collective, two-state-like drying.
3. Drying is not, by itself, evidence that release of unfavorable water drives binding.

They do not establish a universal mechanism for every host–guest system.

