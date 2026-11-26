# Conforming hexahedral reference for the crossed-tool union

The same ten-chart Rust union supplies the geometry. The reference calculation uses displacement-based finite elements with shared global nodes and full Gauss integration (2³ for HEX8, 3³ for twenty-node serendipity bricks). The elastic constants, left end clamp and right end traction of 1 MPa in the global x direction match the IGA study. The mean loaded-end displacement is the external work divided by the applied force of 40.22222222 N.

## Mesh construction

Sample each rational chart at increasing resolution in ordinate and width. At the two poles the first row is a triangular fan. Splitting each triangle into three quadrilaterals and each adjacent quadrilateral into four creates a conforming planar mesh. Extrude these quads in height, propagating the two footprint planes into the central profile in the analysis mesh. The element corners are eight distinct nodes even at the poles. The four internal footprints have exactly two incident element faces each at every mesh facet. The planar mesh approximates the spline boundary by chords. At equal mesh level HEX8 and HEX20 use identical faceted geometry; the quadratic elements have midside nodes on those same edges. The material volume converges to the native 6799.6047258023 mm³.

At every element, positive corner values of the affine planar Jacobian and positive corner values of the bilinear height certify a positive volume Jacobian across the cube for these extrusions. The solver checks Gauss-point determinants, manufactured affine strain energy, load, reaction and free-DOF residual. Mesh and displacement exports for ParaView use VTK type 12 or 25. The study independently assembles element stiffness and enforces strong nodal continuity. Only the material parameters and source geometry are shared with the IGA solver.

## Comparison

| Model | Perpendicular: mean end displacement (mm) | Inclined +3 mm: mean end displacement (mm) | Perpendicular: maximum force jump (N) | Inclined: maximum force jump (N) |
|---|---:|---:|---:|---:|
| IGA, ten charts, unrefined | 0.076785 | 0.080049 | 2.3124 | 2.1370 |
| IGA, ten charts, eight transverse intervals | 0.085631 | 0.166888 | 2.3368 | 2.8967 |
| HEX8, 8,184 elements | 0.088678 | 0.291584 | 1.2297 | 2.1265 |
| HEX8, 27,702 elements | 0.088958 | 0.300249 | 1.0352 | 1.7616 |
| HEX8, 65,760 elements | 0.089096 | 0.303561 | 0.9162 | 1.5334 |
| HEX20, 1,014 elements | 0.088881 | 0.306767 | 1.1706 | 1.7003 |
| HEX20, 8,184 elements | 0.089171 | 0.307772 | 0.8358 | 1.2167 |

The last two HEX20 meshes differ by 0.33% in displacement for either geometry. The transversely refined IGA solution gives 3.97% and 45.78% less displacement than the finer HEX20 solution, respectively. The coarse 22-chart matching-face control from the earlier study gives only 0.083168 mm in the inclined case. Simply extending the partial faces at the original approximation order did not resolve this discrepancy. The present comparison does not isolate the cause of the excess stiffness.

All stresses are recovered independently on the two incident elements. At integration points on the common chordal face, we compare the same physical tractions with opposite outward normals. Force jumps are integrated on the same eight original wing-face test functions per interface, pulled back to the corresponding subdivided face parameters. The reported maximum is the largest norm across the original 32 coefficient vectors. At the finer HEX20 level these jumps are still decreasing appreciably, so the comparisons of local forces are indicative. Using all 32 moment vectors and the finer HEX20 result as denominator, the relative Euclidean differences between the IGA and HEX20 forces on the wing/core sides are 19%/65% for the perpendicular tool and 40%/52% for the inclined tool. The discrete load and reaction balance does not quantify this local stress discrepancy.

## Reproduction

From the repository root, with NumPy and SciPy installed:

```sh
cargo build --locked --offline
study=docs/results/crossed-interface-iga
OPENBLAS_NUM_THREADS=1 python3 "$study/hex_reference.py" target/debug/solidcraft-trivariate build/s2-linear.json --levels 1 2 3 4 --rises 0 3
OPENBLAS_NUM_THREADS=1 python3 "$study/hex_reference.py" target/debug/solidcraft-trivariate build/s2-quadratic.json --quadratic --levels 1 2 --rises 0 3
OPENBLAS_NUM_THREADS=1 python3 "$study/hex_reference.py" target/debug/solidcraft-trivariate build/s2-vtk.json --levels 2 --rises -3 0 3 --vtk-directory build/s2-hex8
```
