# Data dictionary and provenance

All checked-in workbooks are compact, processed inputs for the three disclosed manuscript systems. Raw trajectories and undisclosed-system data are excluded.

| File | Contents | Consumed by |
|---|---|---|
| `data/processed/pmf_thermodynamics.xlsx` | WHAM-reconstructed PMFs and constant-heat-capacity-derived enthalpy, entropy, and heat-capacity profiles from the REMD/umbrella dataset | `analysis/pmf_decomposition/` |
| `data/processed/direct_residual_decomposition.xlsx` | Direct and operational residual components along the umbrella coordinate | `analysis/pmf_decomposition/` |
| `data/processed/spatial_maps.xlsx` | Axial–radial density, hydrogen-bond, and enthalpy grids derived from the unbiased simulations | `analysis/spatial_maps/` |
| `data/processed/mean_occupancy_vs_z.xlsx` | Authoritative corrected `⟨N⟩(z)` and SEM tables for CB7, OA, and restrained β-CD; the β-CD sheet comes from the frame-matched coordinate/water-count rerun | `analysis/occupancy/` |
| `data/processed/occupancy_along_binding_coordinate.xlsx` | Consolidated corrected `P(N|z)` distributions from the 100 ns conventional umbrella windows without replica exchange | `analysis/occupancy/`, `analysis/occupancy_distributions/` |
| `data/processed/cb7_occupancy_rerun.xlsx` | Retained corrected CB7 source package and provenance record; consolidated current inputs are in the two workbooks above | provenance only |
| `data/processed/hydrogen_bond_temperature.xlsx` | Temperature-dependent hydrogen-bond summaries from the unbiased simulations | `analysis/temperature_response/` |
| `data/processed/unbiased_occupancy.xlsx` | Authoritative occupancy values from the independent unbiased simulations, including the revised unrestrained β-CD 25 °C run | `analysis/temperature_response/` |
| `data/derived/water_occupancy_by_host.xlsx` | Host-level occupancy verification and block SEMs; the superseded unrestrained β-CD 25 °C uncertainty is intentionally blank | `analysis/temperature_response/` |
| `data/derived/betacd_cluster_bootstrap_sem.tsv` | Retained β-CD cluster-bootstrap SEM table from the preceding analysis version | provenance only |
| `data/derived/pn_fattail_fits.xlsx` | Generated audit of wet-basin quadratic fits, low-occupancy enhancement, all-window CB7 and β-CD state classifications, and the selected β-CD panel sequence | `analysis/occupancy_distributions/` |

## Units

- Position and distance: Å
- Temperature: °C in plotted temperature series; 298 K for the primary thermodynamic comparison
- Free energy and enthalpy: kJ mol⁻¹
- Density/concentration: see workbook column labels for normalization used by the source analysis
- Occupancy and hydrogen bonds: counts per cavity or water, as labeled in each output

## Generated outputs

`data/derived/drying_transition_analysis.xlsx`, `data/derived/occupancy_moments.xlsx`, and `data/derived/pn_fattail_fits.xlsx` are generated outputs accepted by the publication-boundary check. Review them before committing if the numerical analysis changes.
