# Prepartitioning a NURBS sketch

An operation on the crossed plate specifies a sketch frame, direction, extrusion length, and loops to subtract. The same curve format describes a rectangular contour and a freeform contour. Each curve provides its degree, knots, and homogeneous poles `[wx, wy, w]`. [Two-sketch recipe](../tests/references/g2_cross_sketches.recipe.json)

## Deciding where to place construction lines

1. Extract spans over the curve's active interval. Rust function `extract_active_bezier` evaluates the B-spline's polar form in homogeneous coordinates using checked rational arithmetic. It accepts degrees 1 to 3, integer knots on a common scale, positive weights, and extended periodic sequences as well as curves clamped at the ends.
2. Examine both sides of every knot and the closure. At an interior knot of multiplicity m, a degree-p basis ordinarily guarantees C^(p−m) for m ≤ p. Poles and weights may yield greater regularity. `plan_loop` reports the knot, its multiplicity, the outgoing span, and geometric continuity measured exactly from the jets.
3. Preserve junction provenance. A position discontinuity or zero speed needs explicit handling. A tangent change marks a corner. A shared tangent with a curvature jump marks a G1 junction. A knot inserted to refine a straight edge may retain G2 continuity.
4. Add events imposed by **the other extrusion**. The transverse rectangle walls intersect the periodic profile at parameters that may lie inside spans. Each geometric event must carry identities from both sketches. These events decide where to subdivide spans and extend construction lines as faces. A representation knot and an intersection between operations have distinct provenance.

## Implemented geometric criterion

For two regular jets of an oriented planar curve, the test starts by checking position equality. It then checks that velocities are collinear and point in the same direction. If the right velocity is λ times the left with λ > 0, the junction is G2 when the difference between the right acceleration and λ² times the left acceleration is parallel to the tangent. This test allows both a speed change and a different tangential acceleration. It uses rational operations and reports arithmetic overflow. [Rust implementation](../src/sketch_junction.rs)

Tests cover a cubic rectangle with four corners and an additional subdivision, a rational circle and ellipse with double knots, a nonuniform periodic cubic loop, and two spans with homographic parameter changes. The geometric diagnostic applies to sketch curves. Qualifying their extrusion and intersection also requires volumetric regularity, interface, coverage, and fidelity checks. [Runnable tests](../tests/sketch_junction.rs)
