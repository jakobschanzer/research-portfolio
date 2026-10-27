# Agentic scientific-software workflow

This repository demonstrates a constrained use of Codex for research engineering. The agent may organize code, improve reproducibility, write tests, regenerate figures, and surface assumptions; scientific interpretation and publication decisions remain human-controlled.

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

