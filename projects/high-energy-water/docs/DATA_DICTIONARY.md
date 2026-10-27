# Data dictionary and provenance

All checked-in workbooks are compact, processed inputs for the three disclosed manuscript systems. Raw trajectories and undisclosed-system data are excluded.

| File | Contents | Consumed by |
|---|---|---|
| `data/processed/pmf_thermodynamics.xlsx` | PMF, enthalpy, and entropy profiles for CB7, OA, and β-CD | `analysis/pmf_decomposition/` |
| `data/processed/direct_residual_decomposition.xlsx` | Direct and solvent-mediated force/free-energy components | `analysis/pmf_decomposition/` |
| `data/processed/spatial_maps.xlsx` | Axial–radial density, hydrogen-bond, and enthalpy grids | `analysis/spatial_maps/` |
| `data/processed/occupancy_along_binding_coordinate.xlsx` | Conditional occupancy distributions along the binding coordinate | `analysis/occupancy/` |
| `data/processed/hydrogen_bond_temperature.xlsx` | Temperature-dependent hydrogen-bond summaries | `analysis/temperature_response/` |
| `data/processed/unbiased_occupancy.xlsx` | Unbiased reference occupancy values | supporting analysis |
| `data/derived/water_occupancy_by_host.xlsx` | Host-level temperature/occupancy summary | `analysis/temperature_response/` |
| `data/derived/*_cluster_bootstrap_sem.tsv` | Cluster-bootstrap occupancy standard errors | `analysis/occupancy/` |

## Units

- Position and distance: Å
- Temperature: °C in plotted temperature series; 298 K for the primary thermodynamic comparison
- Free energy and enthalpy: kJ mol⁻¹
- Density/concentration: see workbook column labels for normalization used by the source analysis
- Occupancy and hydrogen bonds: counts per cavity or water, as labeled in each output

## Generated outputs

`data/derived/drying_transition_analysis.xlsx` and `data/derived/occupancy_moments.xlsx` are optional generated outputs accepted by the publication-boundary check. Review them before committing if the numerical analysis changes.
