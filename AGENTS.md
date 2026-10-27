# Repository instructions for coding agents

## Scope

This is a private research portfolio for an unpublished manuscript. Preserve scientific provenance, reproducibility, and confidentiality.

## Publication boundary

- The only systems approved for this repository are CB7, OA, β-cyclodextrin, and its restrained manuscript variant.
- Do not add, infer, enumerate, or describe any other host systems, raw trajectories, private paths, or unpublished results.
- Do not change repository visibility or grant access.
- Do not push changes without an explicit human review request.

## Scientific integrity

- Never fabricate simulation data or replace production results with synthetic values without prominent labeling.
- Distinguish measured/simulated values, derived quantities, and illustrative examples.
- Preserve units and state assumptions beside transformations.
- Treat figure regeneration and numerical checks as acceptance criteria.
- Do not overstate conclusions beyond the manuscript evidence.

## Engineering workflow

1. Read the project README and relevant methods/data documentation.
2. Make the smallest auditable change.
3. Run `make check` and `make test` from the repository root.
4. For analysis changes, regenerate the affected figure and compare the numerical/visual result.
5. Summarize changed files, validation performed, and any unresolved scientific assumptions.

