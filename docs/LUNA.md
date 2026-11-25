# Luna roadmap

The campaign develops the kernel and its audits in Rust. Increment L1 is implemented in `src/bezier.rs`, with native tests in `tests/bezier.rs`. Run `cargo test --locked --offline` before the cross-checks in `make check`. The Python description below records the historical starting point.

## L1 — Evaluate a native rational map

**Starting point.** Run `make check`, record `git rev-parse HEAD`, and inspect `build/affine_pair.json`. The initial proposal was to modify `trivariate/spline.py`, with tests and recipes under `tests/`. Preserve the API and affine-case results.

**Deliverable.** A tensor-product rational Bézier map with explicit degrees, supporting `(1,1,1)` and `(2,1,1)`, and evaluating `X` and `DX` on the closed cube. Fix the Jacobian-column convention, coefficient ordering, units, and face parameters before implementation. Distinguish invalid input from a family that is not yet supported.

**Independent recipe.** Extrude a quarter annulus with radii 1 and 2 mm through 3 mm along z. For each radius `r`, the arc's Euclidean controls are `(r,0)`, `(r,r)`, `(0,r)`, with weights `1`, `√2/2`, `1`. Choose axis order for a positive determinant, for example angle increasing in `u`, radius decreasing in `v`, and height increasing in `w`. Record how the irrational weight is represented.

**Acceptance criteria to freeze for this increment.** Check the affine case at vertices, faces, and interior points. For the curved case, check endpoints, angular midpoint, radius, and z on a grid that includes boundaries (absolute error at most `10⁻¹² mm` for these dimensions in double precision). Compare `DX` columns with independent derivatives and specify method and tolerance before the run. Test inconsistent dimensions and zero weights. Report pointwise tests separately from any proof over the continuous domain.

## Follow-on demonstrator

1. **Two curved blocks.** Split the quarter annulus along z, share the canonical face, verify permutations and reversals, save/reload, edit height, and test atomic rejection. Check exact subdivision and affine-field reproduction. Qualify weight and determinant bounds, injectivity, and partition before claiming J1 for this family.
2. **Bracket.** Build the legs, radius-6 fillet, and two diameter-8 holes from templates. Compare against the FieldCraft fixture (30 cells, 45 interfaces) and analytical volume `34080 − 616π mm³`. Justify extra subdivisions. Volume alone does not prove shape equality.
3. **J2.** Test H¹ conformity, rigid modes, affine field, equilibrium, and convergence under refinement; then compare with the finite-element reference.
4. **J3.** Verify dependency-based regeneration and sensitivities at a fixed regime. Report topology or selection events.
5. **J4.** Freeze corpus and thresholds; publish coverage including rejections, runtime, and memory at comparable precision.

Each delivery records commands, revision, parameters, results, limitations, and the next open criterion. The scientific milestones were open at the start.
