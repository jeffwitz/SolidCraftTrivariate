# Bounded planar ring lookahead prototype

## Pre-registration

- Starting revision: `72e8e19c1bd4b25b41610a11d930f38543e40d40` from the isolated `main` checkout.
- Scope: add an opt-in candidate selector for convex sketch ring partitions; no default planner or recipe-format change.
- Frozen inputs: `seeded_convex_14`, `seeded_convex_04`, and `seeded_convex_26` from `docs/results/sketch-planner-polygon-evaluation/cases.json` (generator seed `20261126`).
- Candidate generation: compare the existing bounded host-separator alternatives and ring dynamic program globally rather than selecting tied-count alternatives by the extruded 3D score. The lookahead variant adds support offsets 1/4 and 3/4 (in addition to 1/3, 1/2, and 2/3), a host-search budget of 128 recursive visits and at most eight host plans. It also tries at most 128 exact shared-junction slides per recipe: for a non-input/non-guide graph vertex lying strictly between collinear neighbors, move it toward either neighbor by 1/8, 1/4, or 1/2 of that segment. Slides are rejected before construction if they leave any source segment containing the original vertex.
- Primary objective: minimize material block count. Secondary diagnostic: minimize the maximum 2D Jacobian condition number, sampled on each cell's closed 3×3 parameter grid. Deterministic first-candidate tie break within relative `1e-12`.
- Required gates: `partition_from_graph` exact rational area/volume and complete incidence checks; positive determinant lower bound on each emitted block; no input refusals. No candidate may replace a lower-block-count candidate.
- Baseline: default `partition` at the same starting revision. Report baseline and prototype block count, method, planar score, certificate lower bound, and elapsed time. The planar score is a diagnostic, not a continuous certificate.
- Exit: prototype API remains explicitly opt-in; report all three inputs, including no-change or worse outcomes, and run targeted Rust checks plus `make check` if the prototype is retained.

## Results

Pending implementation and paired measurements.
