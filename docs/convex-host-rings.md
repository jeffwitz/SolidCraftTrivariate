# Convex hosts and quadrilateral rings

Starting revision: 5252a4e. The automatic polygonal sketch planner accepts a convex stock polygon, disjoint convex cut polygons and their axial intervals. The trial consumes the recipe contours without construction guides. User guides, forced triangulation and Boolean nesting retain their existing paths.

The trial splits hosts, constructs rational quadrilateral rings and fills blind pocket floors. Existing candidates compete by certified block count and sampled condition on ties. Acceptance requires exact material area in every layer, conforming incidences and strictly positive closed-domain Jacobian bounds. Rebuild the frozen polygon corpus from its original recipes and compare against the previous kernel and the counts in the supplied Python study. Include fallbacks and all refusals in the results. Replaying older serialized models is a separate regression check.

## Sketch and volume construction

A separator trial first considers stock diagonals and support directions from the convex hulls of two cutter groups. A chord must strictly separate the groups. The planner tries balanced splits and, for up to six holes, another ordering that rewards clearance. It propagates the endpoints of later chords to neighboring hosts. Both searches are bounded to 64 recursive split nodes each. The clearance score uses a compensated chord length with explicit fused multiply-add operations. On the seeded case 05, different system implementations of `hypot` previously chose different separators at the same count of 57 blocks. The publication check now compares every rational control coordinate of the chosen partition.

For each host, the cyclic dynamic program can join one host side to one hole side or two host sides to a hole vertex. It compares balanced mandatory vertices with optional rational splits. Optional positions include quarters and projected cutter vertices rounded to sixteenths along a host edge. These fractions locate construction points on the unmodified side. Both traversals reuse the chosen points on shared edges. The native incidence and area checks assess the resulting geometry in exact arithmetic. The optimizer samples Jacobians to rank equal-count candidates.

For a blind cut the floor uses the same contour vertices as the walls. The planner tries rotated quadrilateral ladders and can bisect a side of an odd contour. A triangle uses its three edge midpoints and centroid. The selected footprint is extruded through the start and end of each cut. Earlier saved partitions remain auditable through their original planner.

## Reproduce

```sh
cargo build --release --locked --offline
mkdir -p build/demo
recipe=tests/regressions/overlapping_projections_diamonds.recipe.json
target/release/solidcraft-trivariate plan-sketch "$recipe" > build/demo/diamonds-plan.json
target/release/solidcraft-trivariate build "$recipe" build/demo/diamonds.json
target/release/solidcraft-trivariate audit build/demo/diamonds.json
python3 scripts/evaluate_convex_host_rings.py --work build/convex-ring-study
```

The diamond sketch builds in 10 blocks with 12 interfaces and exact volume 113,088 mm³. The previously archived construction used 72 blocks. The evaluation rebuilds the 45 frozen polygon recipes with the current planner.

## Measured results

The release executable built offline from the clean checkout `72e8e19` (kernel source revision `083fdfe`) passed the independent section and volume oracle on all 45 frozen polygon recipes and the 12 additional seeded recipes. The latter use two or three convex cutters and include blind intervals. The tests are reproducible public benchmarks. They measure the best certified candidate selected by this implementation.

| Set | Valid | Previous kernel blocks | Prototype blocks | Native automatic blocks |
|---|---:|---:|---:|---:|
| Frozen polygon corpus | 45/45 | 5,201 | — | 1,208 |
| Seeded convex subset of that corpus | 30/30 | 4,767 | 914 | 910 |
| Additional convex sketches | 12/12 | — | 442 | 420 |

For the 30-case subset the native solver matches the prototype count on 26 cases, improves three by two blocks each and uses two additional blocks on case 11. For case 10, two twelve-sided holes require 24 blocks in this ring construction. The previous kernel used 114. The measured maximum 3D Jacobian condition number falls from 49,185 to 130.3.

The block count falls on 30 of the 45 frozen sketches. The [independent extrusion diagnostics](results/polygon-extrusion-quality/README.md) rank in-plane distortion separately from the ratio of in-plane stretch to layer thickness. On eight sketches the maximum sampled 3D condition number increases by more than one part per million. The largest ratio is 2.69 on seeded case 17. Positive closed-domain determinant bounds, sampled conditioning and optimality are separate claims. This bounded automatic search targets polygonal sketches; intersecting or curved NURBS cuts use their existing builders. [Full per-case measurements and hashes](results/convex-host-rings.json)
