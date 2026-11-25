# Two intersecting perpendicular bores

## Contract

In a box, two through-bores of common radius r have intersecting orthogonal axes parallel to two box directions. The canonical case uses X and Z. Their intersection is the local origin. Part data consists of stock bounds, axes, and radius. The constructor derives separator sheets from the equations of the two active surfaces. The partition uses a guard distance g strictly greater than r and strictly less than the distances to the four faces normal to X and Z. The two faces normal to Y must be farther than r.

## Construction rule

A point is material when both inequalities x²+y² ≥ r² and z²+y² ≥ r² hold. Subtracting them gives the planes x=z and x=−z. They indicate **which cylinder actually defines the wall** in each region. These planes cut both cylinders along branches of their intersection. The two tangent points (0,±r,0) are shared vertices.

In an octant, let m=min(|x|,|z|), d=max(|x|,|z|), e=|y|. In each dominance region, material is characterized by m²+e² ≥ r². For m≤g, the section is the rectangle [0,g] × [0,H_y] with a quarter-disk removed. Divide it into two curvilinear quadrilaterals. Their circular boundary is exactly

    C(q) = r (1−q², 2q) / (1+q²),   q ∈ [0,1/2] or [1/2,1].

The separator runs from (3r/5,4r/5) to (g,H_y). The outer boundary first follows the rectangle's vertical side, then its horizontal side. The homogeneous net joins each quadratic arc to that boundary by ruled interpolation. Extrude the section in the dominant-axis direction from d=m to the stock face d=H_d. Beyond the planes |x|=g and |z|=g, the remaining material is a box in each octant.

Each octant contains four curved blocks and one box, for 40 blocks total. Arcs and traces have integer homogeneous controls. The joins on diagonal planes carry exactly the same parameterization. At the poles, the extrusion-rule width is H_d, so an entire face is not collapsed to a point.

![Walls of both bores and separator in a section](figures/percages-secants.png)

Left: bore walls from the exported control nets, clipped to 7 mm from the intersection to show the junction. Right: the two regions to extrude in an octant of the reference case. The section axes use different scales.

## Why the maps cover the material

In the first sector, the outer-boundary point is B=(g,H_y t) for q=t/2; in the second, B=(g(1−t),H_y) for q=(1+t)/2. In both cases B−C has nonnegative components. Along the oriented arc, m decreases and e increases. The outer boundary has the same monotonicity. Homogeneous interpolation describes, after projection, segment C–B with a strictly increasing parameter in the radial direction. Differentiating in the angular direction, the derivative of the blending coefficient is parallel to the segment and vanishes from the determinant. The two remaining terms are the weighted tangents of C and B. Their determinant with B−C has a strict sign, including at the endpoints: the vertical edge uses g>r, and the horizontal edge uses H_y>r.

Each section therefore has a Jacobian of constant sign and a simple boundary traversed once. The topological degree on its interior is one in absolute value, so every interior point has one preimage. The two sections cover the rectangle outside the disk and share exactly the separator segment. For each section point, the extrusion rule is bijective on [m,H_d], with length at least H_d−g>0. The octant choice, then the dominant coordinate, separates block interiors. Points on equalities lie on the corresponding traces. This argument establishes coverage independently of numerical volume integration.

## Reproducible verification

From the repository root, run `cargo test --locked --offline --test steinmetz_core -- --nocapture`. The test reconstructs incidences, uses the outward-rounded Bernstein certificate on each closed cube, evaluates the walls, and integrates the map Jacobian. In the 60 × 40 × 30 mm box with r=3 mm, the volume common to the voids is 16 r³/3 = 144 mm³. At local height y, the common section of the two cylinders is a square of area 4(r²−y²), giving this volume by integration over [−r,r]. Expected material volume is 72,000 − 810π + 144 mm³.

Different radii or skew axes change the equations of the separator sheets. They require a separate qualification domain.

### Replay from the two features

```sh
mkdir -p build/steinmetz-demo
cargo run --locked --offline -- build tests/references/steinmetz_box.recipe.json build/steinmetz-demo/model.json
cargo run --locked --offline -- audit build/steinmetz-demo/model.json
cargo run --locked --offline -- export-step build/steinmetz-demo/model.json build/steinmetz-demo/perpendicular-bores.step
cargo run --locked --offline -- export-geometry-vtu build/steinmetz-demo/model.json 4 build/steinmetz-demo/blocks.vtu
cargo run --locked --offline -- export-analysis build/steinmetz-demo/model.json 4 build/steinmetz-demo/analysis.json
```

For this recipe, the persisted-map audit finds 68 interfaces, 104 exterior faces, and a shared Jacobian lower bound of 14.74559999999808 mm³. Independent map integration gives 69,599.309950592913 mm³. Each bore wall of length L has exposed area 2πrL − 8r²; the test checks these separately to 10⁻⁸ mm². Open CASCADE re-imports one valid closed solid of volume 69,599.30995051493 mm³, centered at (30,20,15) mm.

The exported IGA space has 554 scalar functions and 1,662 displacement components. On the exported model at Gauss order 4, an independent oracle checks 612 face-coefficient matches and exactly reproduces the affine field's parametric values and gradients on the rational test grid. The relative difference between assembled elastic energy and constant-strain energy is at most 2.06 × 10⁻¹⁵ for the stiffness diagnostic assembled at order 2. The export's parameter 4 sets volume quadrature. This test covers assembly and interfaces; it does not quantify stress concentration under load.
