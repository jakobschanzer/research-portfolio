# Geometry-specific cavity volumes

## Design

JANAL performs the common statistical-mechanical analysis, while a separate volume module answers whether a water position lies inside the host cavity. Production runs therefore pair:

1. the shared `janal.f` core;
2. a host-geometry-specific `vol.f` implementation; and
3. that implementation’s `vol` input file containing the relevant atom selections or parameters.

The volume implementation was adapted when the cavity geometry changed. Results should not be interpreted as though one boundary definition was reused unchanged across every host.

## Checked-in implementation

`src/fortran/host_analysis/vol.f` implements a six-plane convex region. At each frame it constructs oriented planes from host atom triplets, uses check atoms to determine the interior side, and tests candidate points against all six boundaries under periodic minimum-image conventions.

`vol_value()` currently returns `-1`, preserving JANAL’s legacy reported binding-volume behavior. The membership test is active when a valid `vol` file is present; absent or malformed input disables membership rather than silently constructing a default cavity.

This backend documents the interface and one production geometry pattern. It is not labeled as a universal backend, and a host assignment is deliberately not inferred without the original run configuration.

## Required provenance for another backend

Any future backend added to this repository should include its matching `vol` template, host assignment, atom-index convention, coordinate-frame definition, periodic-boundary treatment, and a simple inside/outside validation fixture.

