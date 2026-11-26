# Fixed-partition bore edge treatments

## Qualification scope

This study starts at `34b9720f059d9e93bb093edae557e77031f25582` and edits the published FieldCraft compact bracket at `docs/site/data/fieldcraft_x64_y64/fieldcraft_x64_y64.model.json`. The source carries the 60 × 50 × 40 mm L body, the R6 inside elbow, and two R4 bores. The unmodified model is first passed through `interchange::audit`; the edit then consumes those audited top-level block maps.

The two bore chart groups are `wing_east_` (axis Y, center [40,4,20] mm) and `wing_north_` (axis X, center [4,32,20] mm). Each treatment is independently parameterized in quarter-millimetre units. The accepted scope is a 0.25 mm end treatment on both mouths of each bore. Other widths, profiles, intersections, topology changes and imported geometries are outside this result.

## Fixed partition and checks

The ten logical chart IDs remain fixed. Each bore chart has three Bézier spans at its original axial stations `[0, 0.25, 7.75, 8]` mm. The two elbow charts receive the matching parameter splits `[0, 1, 31, 32]/32`. This gives 30 span elements while preserving the original `(u,v,w)` chart orientation. The exact whole-face audit finds all 53 span incidences, including all 11 original logical interfaces and 20 within-chart span joins.

Every sharp and chamfer span receives the outward-rounded closed-cube Jacobian certificate. The global lower bound for each state is serialized in `build/bracket-edge-edit-sharp.json` and `build/bracket-edge-edit-chamfer.json`. Input widths must satisfy `width = axial_length/32`, leave a positive middle span, and remain smaller than the bore radius. Source charts must match the expected ten-map inventory, exact circular wall identity, and four angular charts per bore. Rejection leaves the caller's input maps unchanged.

The chamfer profile changes the bore radius linearly from 4.00 mm to 4.25 mm over 0.25 mm at each end, then returns to the 4 mm cylinder. An independent exact-rational section oracle checks the entire circle, both ramp profiles, unchanged exterior restrictions and the added removal volume `49π/48 mm³` for the two bores together (about 3.20704250054 mm³). This is an exact geometry construction; the finite-radius predicates and floating-point certificate bounds retain their stated numerical method.

## Round fillet and singular boundary

The round branch uses a 0.25 mm exact rational quarter-circle profile with homogeneous weights `[1,1,2]` at the low end and `[2,1,1]` at the high end. It preserves all 30 spans and the 53 exact incidences. `certify_cap_contact` proves strict positive Jacobian away from the named bore-wall/cap contact rail and proves the determinant vanishes identically on that rail. The global determinant lower bound is therefore exactly zero; this branch is geometric-only and does not meet the closed-cube strict-regularity invariant. It must not be accepted as an analysis-ready block complex under the ordinary positive-Jacobian contract.

The singular set is limited to the exact bore-wall/cap contact rails in the 16 round end spans. It is not smoothed or reported as a positive certificate. The sharp and chamfer states remain strictly regular on their closed spans.

## Reproduction

From the repository root:

```sh
cargo test --offline --test bracket_edge_edit
cargo run --offline --example bracket_edge_edit -- sharp
cargo run --offline --example bracket_edge_edit -- chamfer
cargo run --offline --example bracket_edge_edit -- fillet
```

The example writes `build/bracket-edge-edit-{sharp,chamfer,fillet}.json`, including the source path, edge parameters, logical IDs, span intervals, degrees, homogeneous controls, incidence count, Jacobian scope and lower bound. The Rust API is `bracket_edge_edit::edit_compact_bracket(parent, bores, treatment)`. Use it only after auditing the parent model.
