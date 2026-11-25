# Fixed-connectivity annulus

The experiment starts from revision `d25b6b6`. The annulus has radii 2 and 4 mm and height 8 mm. The operation treats all four circular edges together, with a 0.25 mm fillet radius or chamfer setback. Two experiments address two decisions the constructor can make from the sketch.

| Construction plan | Sharp | Chamfer | Fillet | Bézier audit cells per state |
|---|---|---|---|---:|
| Edit the gallery's four sectors | 4 volumes | 4 volumes | Strict-Jacobian obstruction with these incidences | 12 |
| Prepare a five-patch section and revolve it in two half-turns | 10 volumes | 10 volumes | 10 volumes | 128 |

## Editing the four original sectors

The prototype reads the four `full_annulus_v1` blocks and preserves their identities, orientations, angular poles, and weights up to a common homogeneous scale. It inserts two internal axial knots at 0.25 and 7.75 mm in each degree-(2,1,1) block. Both states have 24 poles per volume. At height z, radius varies between 2+d and 4−d. For the chamfer, d=max(0, 0.25−min(z,8−z)); the sharp state uses d=0. Editing moves two radial rows at each end. The ideal chamfered profile volume is 95.25π mm³.

All four sector interfaces match coefficient by coefficient. Bernstein certification on the twelve closed cells gives determinant lower bounds above 0.82079 mm³ in the sharp state and 0.69254 mm³ in the chamfered state, in local cell parameters. Radius, height, and angle each traverse their interval once, using the same dyadic circle approximation as the original constructor. [Constructor and tests](../examples/annulus_four_block_chamfer.rs) · [Sharp](results/annulus-fixed-topology/four-sharp.json) · [Chamfer](results/annulus-fixed-topology/four-chamfer.json)

## Fillet obstruction on these four volumes

Each original sector is the product of an angular quadrant and a square section. Its two angular end faces are internal interfaces. Two exterior physical faces meet at each of the four longitudinal edges from the section corners. After treating all four edges, the entire meridional boundary is smooth in tangent-plane sense, including at arc/line junctions.

In a tensor cell touching such an edge, take faces u=0 and v=0 and an interior point on the edge away from angular knots. The tangent derivatives of the two faces are respectively (Fv,Fw) and (Fu,Fw). If the physical boundary has one tangent plane, columns Fu, Fv, Fw all lie in that plane. The determinant therefore vanishes on the boundary. The same reasoning applies to one-sided derivatives of the adjacent cell, regardless of internal knot count. One such edge is enough to fail the project's criterion. The chamfer retains tangent breaks to which these parametric edges can terminate. [Rust witness for the rank defect and parameter changes](../tests/annulus_four_sector_obstruction.rs)

## Template prepared from the section and its revolution constructor

The constructor connects the midpoints of the four profile sides to an interior diamond centered at (3,4), with half-diagonals 0.5 and 2 mm. The four strips between the core and contour each contain a physical corner at the midpoint of their **map boundary**. The core and strips also exist in the sharp state. Each strip uses four tangential spans; two cover the editable region. There, depending on state, exterior poles describe two sharp-corner branches, two collinear segments, or two restrictions of a tangent rational circular arc. Setback and radius are 0.25 mm at all four corners.

The sweep uses two injective half-turns, each with two rational spans and a matched angular speed. Each volume meets the same neighbor on two distinct angular faces. All ten identities, incidences, degrees, knots, multiplicities, and weights are identical in all three states. Only wall pole coordinates change. The prototype JSON records all ten complete NURBS nets. [Rust constructor and audit](../examples/annulus_prepared_volumes.rs) · [Sharp](results/annulus-fixed-topology/ten-sharp.json) · [Chamfer](results/annulus-fixed-topology/ten-chamfer.json) · [Fillet](results/annulus-fixed-topology/ten-fillet.json)

For each state, all 128 extracted cells have determinant strictly greater than 0.0048828 mm³ in local coordinates. Homogeneous-coefficient checks find 352 cell interfaces and the expected 64 cell faces on the skin. In the section, the four strips cover disjoint meridional quadrants. In each quadrant, the strip boundary is simple (core side, two rails, outer contour). Positive determinant in every cell and oriented joins, together with this boundary of topological degree one, give a single cover. Revolution preserves it because r remains strictly positive and the two half-turns cover the circle once.

The oracle reads the JSON nets, extracts elements from knots, and integrates determinants with tensor Gauss quadrature of orders 12 and 16. It compares the sharp state against 96π. The treated states are compared with volumes measured previously from independent STEP re-imports. [Volume oracle](../tests/references/annulus_prepared_volumes.py)

| State | q16 volume (mm³) | Reference (mm³) |
|---|---:|---:|
| Sharp | 301.592894744620 | 301.592894744620 |
| Chamfer | 299.236700254428 | 299.236700254422 |
| Fillet | 300.581607414644 | 300.581607414643 |

Exact-fraction audit recovers all four profiles and their joins from the recorded data. For each fillet, the polynomial identity of a radius-1/4 mm circle holds over both complete spans. Tangents at both contacts and at the span junction have the expected direction. Among the 1,530 stored poles, 80 positions change from sharp to fillet and 120 from sharp to chamfer. [Profile check](../tests/references/annulus_prepared_profile.py)

## In what class are ten volumes minimal?

Consider a conforming quadrangulation of a disk section with smooth boundary, no double edges or identifications, and regular map corners. Let I be the number of interior vertices, B the number of boundary vertices, and q the number of quadrilaterals incident at a vertex. Euler's formula and side counting give F=I+B/2−1 and 4=Σᵢ(4−q)+Σᵦ(2−q). The image of a regular map corner has angle strictly less than π. Thus at least three maps must meet at an interior vertex and two at a smooth boundary vertex, so I≥4. A simple contour has an even number B≥4; hence **F≥5**. The diamond and its four strips attain this bound.

In a product construction of this section with angular intervals, one angular map on a closed interval would identify its two endpoints. Injectivity on the closed domain requires at least two intervals. The two half-turns attain this bound and give **ten volumes in this class**. Requiring four quadrants as distinct volumes gives twenty. The proof covers neither periodic maps with self-identification nor all 3D hexahedral decompositions not obtained by revolution.

## Replay and scope

```sh
cargo test --offline --test annulus_four_sector_obstruction
cargo test --offline --example annulus_four_block_chamfer
cargo test --offline --release --example annulus_prepared_volumes -- --nocapture
cargo run --offline --example annulus_prepared_volumes -- fillet 2 > build/annulus-ten-fillet.json
python3 tests/references/annulus_prepared_profile.py
python3 tests/references/annulus_prepared_volumes.py docs/results/annulus-fixed-topology/ten-{sharp,chamfer,fillet}.json
```

Qualification covers geometry and interfaces. The IGA space, degrees of freedom, and mechanical performance of the ten-volume template remain to be measured. The production path still needs to integrate the prototype's logical multi-span maps. The planner must count persistent volumes, spans, poles, and computation cost. The earlier 96-block single-span experiment used 96 geometric cells; this ten-volume template uses 128 cells in this implementation.
