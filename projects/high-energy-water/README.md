# Molecular Confinement Does Not Create High-Energy Water

**A computational study of hydration thermodynamics and collective drying in CB7, OA, and β-cyclodextrin**  
Jakob L. Schanzer · manuscript in preparation

## The scientific question

A common model of molecular recognition proposes that confinement creates energetically unfavorable—or “high-energy”—water and that releasing this water drives guest binding. This work separates two ideas that are often conflated:

1. **Energetics:** are the waters that actually occupy a cavity unfavorable relative to bulk?
2. **Drying:** does an approaching guest induce a collective transition between wet and dry cavity states?

The simulations show that **collective drying is real for CB7 and OA, but it is not evidence that their resident cavity waters were high in energy**. Restrained β-CD provides the contrasting case: it loses water gradually rather than through a comparable two-state transition.

## The manuscript story, figure by figure

### Figure 1 — Three cavities with different architectures

![Representative structures of OA, restrained beta-cyclodextrin, and CB7](figures/figure1_host_systems.png)

OA, restrained β-CD, and CB7 differ in shape, portal geometry, flexibility, and symmetry. Those architectural differences let us test which hydration behavior is general and which depends on the physical character of the pocket.

### Figure 2 — Favorable binding does not require a favorable solvent contribution

![PMF thermodynamic and direct-residual decomposition](figures/figure2_pmf_decomposition.png)

I obtained the temperature-dependent PMFs using **replica-exchange molecular dynamics combined with umbrella sampling**. Harmonic umbrella windows were placed every 0.04 nm along the host–guest separation coordinate. Exchanges were attempted across 20 replicas spanning 5–60 °C, and the biased distributions were combined with WHAM to reconstruct ΔG<sub>PMF</sub>(*z*,*T*).

At each position *z*, I fit the temperature-dependent PMF to a constant-heat-capacity model referenced to 298.15 K:

```text
ΔG_PMF(z,T) = ΔH_ref(z) − TΔS_ref(z)
             + ΔCp(z)[(T − T_ref) − T ln(T/T_ref)]
```

This fit separates the PMF into position-dependent enthalpic, entropic, and heat-capacity contributions: ΔH<sub>ref</sub>(*z*), ΔS<sub>ref</sub>(*z*), and ΔCp(*z*). The top row shows ΔG(*z*), ΔH(*z*), and −TΔS(*z*) at the reference temperature; the bottom row compares the PMF with the direct host–guest interaction and the operational non-direct residual.

- CB7 binds almost entirely enthalpically, and its direct host–guest interaction nearly reproduces the PMF.
- OA is also enthalpy dominated, but its residual opposes binding near the minimum.
- Restrained β-CD is the key counterexample: substantial water displacement accompanies favorable binding without a favorable binding enthalpy. Its residual also opposes binding.

The PMF well depths are approximately −110.4, −56.5, and −31.9 kJ mol⁻¹ for CB7, OA, and restrained β-CD. These are well depths relative to the separated state, not standard-state binding free energies.

### Table 1 — Hydrogen-bond deficiency is not a complete energetic model

![Hydrogen-bond accounting and binding-enthalpy comparison](figures/table1_hydrogen_bond_summary.png)

The hydration partition quantity, 𝒵, and the binding hydrogen-bond balance, δ, quantify water structure, but neither reproduces the full binding enthalpy. Hydrogen-bond counts omit bond-strength heterogeneity, host–water and other water–water interactions, guest desolvation, reorganization, and direct host–guest interactions. A hydrogen-bond deficit therefore cannot by itself establish a population of energetically unfavorable water.

### Figure 3 — Occupied cavity water is enthalpically stabilized on average

![Mean cavity occupancy as a function of temperature](figures/figure3_temperature_response.png)

Mean cavity occupancy decreases from 5 to 95 °C for every host model: approximately 6.40→5.36 waters for CB7, 4.20→3.04 for OA, 4.87→4.31 for restrained β-CD, and 3.69→2.92 for unrestrained β-CD. After accounting for the comparatively small thermal-expansion contribution, the negative temperature slopes imply a negative mean enthalpy for transferring water from bulk into each cavity.

This is the first direct challenge to the high-energy-water picture: **the resident cavity population is lower—not higher—in mean enthalpy than bulk water**.

### Figure 4 — Water occupies stabilized interfacial sites

![Spatial hydrogen-bond, density, and local-enthalpy maps](figures/figure4_spatial_maps.png)

The axial–radial maps show hydrogen-bond deviation from bulk (top), relative density (middle), and fitted local enthalpy (bottom). Dense hydration layers follow the host surfaces and coincide with the lowest local enthalpies. Regions with unfavorable fitted enthalpy have negligible occupancy and therefore do not constitute a substantial resident population available for release.

Some dense, low-enthalpy sites remain deficient in water–water hydrogen bonds. Other interactions—consistent with host–water stabilization—compensate for the missing water–water coordination. Reduced hydrogen bonding is therefore not synonymous with high energy.

### Figure 5 — The guest triggers sharp drying in CB7 and OA

![Mean cavity occupancy along the host-guest binding coordinate](figures/figure5_mean_occupancy_vs_z.png)

Mean occupancy identifies where water leaves as cyclohexane enters. OA passes from 25% to 75% empty over approximately 0.61 Å, and CB7 over 0.92 Å. Restrained β-CD responds more gradually and asymmetrically: approximately 2.21 Å through its negative-*z* narrower primary portal and 1.19 Å through its positive-*z* wider secondary portal.

These widths quantify an equilibrium response to guest position, not a kinetic drying rate. More importantly, a mean curve alone cannot distinguish concerted switching from the sequential loss of independent waters. That distinction requires the complete occupancy distributions.

### Figure 6 — P(N|z) reveals the drying mechanism

![Distance-resolved occupancy and log-probability distributions](figures/figure6_occupancy_distributions.svg)

Each column pairs the conditional occupancy distribution, P(N|*z*), with ln P(N|*z*) as the guest approaches. The schematics use one coordinate convention throughout: **negative *z* is left and positive *z* is right**; red is the host and blue is cyclohexane.

- **CB7:** near the portal, wet and dry populations coexist. Around *z* = 4.0 Å, probability is distributed between states centered near two and four waters, with a less probable intermediate occupancy.
- **OA:** at *z* = 3.6 Å, a wet maximum near *N* = 4 coexists with a dry *N* = 0 population. By *z* = 3.2 Å the dry state is most probable while the wet basin remains populated.
- **Restrained β-CD:** the distributions shift progressively toward lower *N* from either portal. Persistent wet/dry coexistence is absent, indicating continuous dehydration rather than the same collective transition.

The dashed red curves extrapolate a quadratic fit to the wet-state basin. Blue points mark low-occupancy probabilities enhanced above that near-Gaussian expectation. These non-Gaussian low-*N* tails develop before the cavity becomes dry and are the fluctuation signature of proximity to a drying transition. They are pronounced for CB7 and OA and develop more weakly and gradually for β-CD.

<details>
<summary><strong>Open host-specific P(N|z) panels for closer inspection</strong></summary>

#### CB7

![CB7 occupancy distributions](figures/cb7_pn_lnpn_paired.png)

#### OA

![OA occupancy distributions](figures/oa_pn_lnpn_paired.png)

#### Restrained β-CD

![Restrained beta-cyclodextrin occupancy distributions](figures/betacd_pn_lnpn_paired.png)

</details>

## Physical conclusion

The combined evidence supports a specific mechanism rather than a slogan:

- Populated hydration sites are enthalpically stabilized on average.
- CB7 and OA nevertheless approach a wet–dry phase boundary and undergo collective, two-state-like drying as the guest crosses the portal.
- Restrained β-CD dehydrates more continuously, demonstrating that collective drying depends on host architecture and wettability.
- Drying and high-energy-water release are physically distinct claims. Observing the first does not establish the second.

## Computational workflow

```text
simulation trajectories / energies
        │
        ├── replica-exchange umbrella sampling
        │       └── WHAM ΔG_PMF(z,T)
        │               └── constant-Cp fit → ΔH(z), ΔS(z), ΔCp(z)
        ├── JANAL (Fortran) ── density, energy, chemical potential,
        │                       fluctuations, occupancy distributions
        ├── H-bond analysis ── water–water and water–host networks
        └── geometry tools ─── host orientation and cavity membership
                                      │
                                      ▼
                             processed workbooks
                                      │
                                      ▼
                     Julia analysis + uncertainty propagation
                                      │
                                      ▼
                      manuscript figures and quantitative checks
```

The production design separates a shared JANAL analysis core from a geometry-specific cavity module. Each host uses its own `vol` input and, where required by its geometry, an adapted `vol.f` backend. The checked-in backend is a six-plane polyhedral implementation; it demonstrates that interface rather than claiming that one volume definition fits every host.

The repository reproduces the displayed analysis from compact processed workbooks. The underlying replica-exchange umbrella trajectories and full WHAM production directories are excluded because of size and manuscript confidentiality.

## Repository tour

| Path | Purpose |
|---|---|
| `src/fortran/host_analysis/` | Shared JANAL core and geometry-specific volume interface |
| `src/fortran/hydrogen_bonds/` | Hydrogen-bond and bridging-water analysis |
| `src/fortran/ring_geometry/` | β-CD ring orientation and spacing analysis |
| `src/python/` | Occupancy counting and deterministic aggregation |
| `analysis/` | Reproducible Julia figure and transition analyses |
| `data/processed/` | Compact inputs for the displayed manuscript analyses |
| `data/derived/` | Bootstrap errors, fit audits, and derived summary tables |
| `figures/` | Numbered manuscript figures plus detailed supporting panels |
| `docs/` | Methods, data provenance, validation, and agentic workflow |

## Reproduce the checks and figures

Requirements: GNU Fortran, Python 3.10+, and Julia 1.10+.

```bash
make check
make test

julia --project=. -e 'using Pkg; Pkg.instantiate()'
make figures
```

See [Methods](docs/METHODS.md), [Data dictionary](docs/DATA_DICTIONARY.md), [Validation](docs/VALIDATION.md), [Volume backends](docs/VOLUME_BACKENDS.md), and [Agentic workflow](docs/AGENTIC_WORKFLOW.md) for assumptions and limitations.

## Scope and confidentiality

Only the systems disclosed in the manuscript—CB7, OA, and β-CD, including the restrained β-CD protocol—are represented here. Raw trajectories are not included. This repository is private because the manuscript is in preparation; permission to view does not imply permission to redistribute.
