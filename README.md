# Janal & Volume Modules

This repo currently includes:
- `janal.f` — HostAnalysis core for cavity water energetics and statistics
- `vol.f`   — Volume module (generalized across cavity geometries)

## Build
# (example using gfortran; adjust flags as needed)
gfortran -O3 -fimplicit-none -o janal janal.f vol.f

## Inputs (expected files)
- `drive`, `temperature`, `top`, `top-bind`, `top-hb`, `pot`, `pot-sol`
See comments in `janal.f` for formats.

## Outputs (selected)
- `misc.out`, `muzr.out`, `enzr.out`, `dist.out`, `fluc.out`, `data.out`
2D fields: `gg`, `gh`, `gh2`, `gn`, `gn2`, `vn1`, `vn2`, `hd` (z × r grids)

## Notes
- Implements test particle insertion (TPI) for LJ and WCA inside pockets
- Computes H-bond distributions inside/outside the pocket
- Generalized across different cavity volumes
