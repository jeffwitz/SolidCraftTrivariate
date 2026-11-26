# Elasticity across partial interfaces in the S2 union

[Degree elevation, thickness refinement, and C¹ field study on the inclined union](bending-enrichment.md)

This independent NumPy/SciPy/SymPy study assembles the rational (3,1,1) basis on the ten chart control nets emitted by the native Rust CLI. Cases use transverse rises of −3, 0 and +3 mm across 100 mm. From the repository root run

```sh
cargo build --locked --offline
OPENBLAS_NUM_THREADS=1 python3 docs/results/crossed-interface-iga/run.py target/debug/solidcraft-trivariate docs/results/crossed-interface-iga/results.json
```

Seven full-face neighbours share trace coefficients. The two collapsed polar lines have single-valued displacement. Each of the four partial interfaces uses the physical height on the wing to evaluate the matching central trace at w=z/12. The common ordinate parameter is u. After multiplying by the square of the common cubic weight denominator, the trace difference has degree at most six in u and one in the wing height parameter. Equality at seven distinct u values and both wing-height endpoints therefore implies equality across that face in exact arithmetic. The script enforces these linear constraints by exact rational null-space reduction followed by floating-point orthonormalization and repeats the evaluations at 17 independent u locations. The four-node control underconstrains inclined interfaces: the largest held-out scalar jump for a unit coefficient vector in its admissible space is 0.0001121. With seven nodes this diagnostic stays below 8e−16.

The material is homogeneous and isotropic, E=1000 MPa and Poisson ratio 0.3. First, the free-body eigenspectrum is checked, using a rigid-mode threshold of 1e−10 times the largest eigenvalue. The affine uniaxial-stress test checks representability, energy and equilibrium against tractions integrated on the exposed boundary, subtracting the covered part of each central face. The second load case clamps the left end at x=0 and applies 1 MPa axial traction to the right end. Interface force is the sum of right-wing residual forces associated with virtual rigid translation. The reported tip observable is the traction-weighted mean x displacement of the loaded end.

| Rise (mm) | Independent scalar DOFs, ten charts | Mean end displacement, ten charts (mm) | Matching-face control, 22 charts (mm) | Relative difference |
|---:|---:|---:|---:|---:|
| −3 | 88 | 0.08004945 | 0.08316766 | −3.75% |
| 0 | 100 | 0.07678517 | 0.07937979 | −3.27% |
| +3 | 88 | 0.08004945 | 0.08316766 | −3.75% |

The matching control is built independently in the study by extending the two cutting planes through every central chart. It uses 172 independent scalar coefficients before clamping. Homogeneous control rows are checked for equality up to a constant scale at each declared interface. Both models have the same geometry, material and boundary loading. Both quadrature orders (8 and 14) agree on the end displacements to 2e−12 mm. The comparator is a second coarse approximation. This test establishes continuous displacement, force transfer and a measurable stiffness difference at these discretizations. Quantifying continuum displacement and local stress error requires analysis-space refinement. The Python assembler currently targets this union and its prescribed correspondences.

## Nodal force balance

Add `--nodal-forces nodal-forces.json` to the study command to write both force contributions and their signed sum in world XYZ (N). The test uses the cantilever load described above and all four wing/core interfaces.

For the assembled check, split the broken residual Ku−f into wing and core contributions. Select independent wing trace control coefficients by pivoted QR of the admissible clamped displacement space restricted to all four footprints. Lift a unit virtual displacement of each selected coefficient into that space. Its virtual work against each subdomain residual yields the conjugate forces. End coefficients shared by adjacent bands are represented once. Across 28 independent scalar coordinates for the perpendicular sweep and 16 for either inclined sweep, the largest vector imbalance over both quadratures is 2.29e−12 N.

For the stress check, evaluate strain from the physical gradient on each side separately and apply the isotropic constitutive law. Integrate each outward traction against the **same** eight rational basis functions of the wing face on its physical footprint. Each JSON row gives eight local vectors, indexed by 2*a+w with a=0..3 and w=0..1. These consistent face forces probe reconstructed stress and can jump even in a conforming displacement solution. They are face-local integrals, not point forces at vertices. With perpendicular tools the largest vector jump is 2.3124 N on the ten-chart model and 2.1778 N on the 22-chart model. For either tilt it is 2.1370 N and 2.3145 N respectively. The maxima occur on the loaded wing. Gauss orders 8 and 14 agree to eight decimal places in newtons. The affine stress control returns less than 6e−15 N jump. A continuum traction error estimate would require a refinement study.

## Refinement with fixed force probes

`refine_normal.py BINARY OUTPUT` inserts 1, 2, 4 and 8 analysis intervals along the transverse chart coordinate v. Add `--combined` to refine both the ordinate parameter u and v to 1, 2 and 4 intervals. Both studies keep the original ten geometric charts, the exact partial-face coupling, the load and the clamp. Children are used internally for basis evaluation and quadrature. The displacement space is continuous across new element faces and remains linear through the height w. On subdivided faces, the transpose of the rational Bézier restriction transfers the measured force moments back to the original eight test functions. Taking the maximum only after summing those contributions makes the comparison independent of the size of the test-function support.

The maximum force jump (N) in the perpendicular case changes from 2.3124 to 2.6111, 2.4170 and 2.3368 under transverse refinement. With either inclination it changes from 2.1370 to 2.4093, 2.4630 and 2.8967. Simultaneous u/v refinement to four intervals gives 2.4165 N and 2.4625 N respectively. At eight transverse intervals, mean loaded-end displacement is 0.08563 mm for perpendicular tools and 0.16689 mm for inclined tools. The maximum occurs on the loaded side at the first level and on the clamped side at the finest transverse level. The algebraic equilibrium residual stays below 4e−10. This experiment establishes neither convergence of stresses nor convergence of tip displacement in the full three-dimensional space.
