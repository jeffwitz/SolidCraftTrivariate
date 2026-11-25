# Unequal-width filleted L regression

- Starting revision: `8eb6ea032f5834eaf01c206e9e3d8495cca6d2f6`.
- Production files: `src/orthogonal_sketch.rs`, `src/orthogonal_sketch_json.rs`, `src/sketch_fillet_chart.rs`, and `src/steinmetz_partition.rs` (the latter initializes the optional ledger field for other templates).
- Preconditions: one simple six-vertex orthogonal L, one re-entrant right-angle fillet, integer millimeter coordinates, positive arm clearances, no bores.
- Independent oracle: polygon area by rectangle inclusion-exclusion is `100*12 + 20*80 - 20*12 = 2560 mm^2`. The fillet adds `r^2*(1-pi/4)` per unit height, so expected volume is `30912 - 48*pi mm^3`.
- Required checks: positive certified determinant on every closed chart; conforming complete-face incidence; builder/audit replay; volume coefficients equal `30912` and `-48*pi`; the reason is replay-checked and can be omitted by historical models; unequal widths are accepted without changing historical equal-width outputs.
- Refusals: nonpositive arm clearance, invalid re-entrant vertex, and a radius not strictly smaller than either rectangular patch extent.
- Exit criteria: regression passes through the public interchange builder and independent audit, with no change to prior audit records.
