# Agentic scientific-software workflow

This repository demonstrates a constrained use of Codex for research engineering. The agent may organize code, improve reproducibility, write tests, regenerate figures, and surface assumptions; scientific interpretation and publication decisions remain human-controlled.

## Division of work in this repository

The underlying research predates the portfolio engineering: Jakob L. Schanzer designed and ran the simulations, developed or adapted the core scientific analysis programs, supplied the processed scientific results, and set the repository's publication boundary. Scientific interpretation was developed by Jakob with his advisor, Hank Ashbaugh, and their experimental collaborator, Bruce Gibb.

Working from those materials and under Jakob's direction, Codex:

- generated the repository and project structure;
- created and iteratively refined the checked-in Julia analysis and figure-generation scripts;
- implemented reproducible Python aggregation, deterministic tests, and numerical audit outputs;
- added Fortran compilation checks, continuous integration, and the automated publication-boundary scan; and
- drafted and revised the README, methods, provenance, validation, and workflow documentation.

Codex did not choose or run the production simulations, generate replacement scientific data, decide which unpublished systems to disclose, or independently determine the scientific conclusions. Jakob inspected the source materials, directed corrections, evaluated the physical meaning of the outputs with his human collaborators, and approved the final repository changes. The repository therefore demonstrates supervised agentic development with explicit scientific ownership, not autonomous authorship of the research.

## Control loop

```text
human objective
      ↓
inspect code + provenance
      ↓
propose a bounded change
      ↓
edit within publication boundary
      ↓
compile → test → regenerate affected analysis
      ↓
report evidence, warnings, and unresolved assumptions
      ↓
human review before commit/push
```

The root `AGENTS.md` is the machine-readable policy. It limits the allowed systems, forbids visibility/access changes, requires explicit labeling of synthetic data, and makes compilation/tests part of the acceptance criteria.

## Example task contract

**Objective:** make the temperature-dependent hydrogen-bond aggregation reproducible.

**Agent actions:** identify the unseeded bootstrap, add a seed and command-line interface, validate expected temperature files, emit portable CSV outputs, and create a deterministic synthetic unit test.

**Validation gate:** two runs using the same seed must be byte-identical; missing temperatures must fail clearly; scientific source files are never modified by the test.

**Human decision:** confirm that the estimator and reported uncertainty match the manuscript’s intended statistical interpretation.

This structure is designed for auditable scientific work: the agent accelerates implementation and verification without becoming an untraceable source of scientific claims.
