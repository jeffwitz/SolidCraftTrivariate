# Elasticity across partial interfaces in the S2 union

This independent NumPy/SciPy study assembles the rational (3,1,1) basis on the ten chart control nets emitted by the native Rust CLI. Cases use transverse rises of −3, 0 and +3 mm across 100 mm. From the repository root run

```sh
cargo build --locked --offline
OPENBLAS_NUM_THREADS=1 python3 docs/results/crossed-interface-iga/run.py target/debug/solidcraft-trivariate docs/results/crossed-interface-iga/results.json
```

Seven full-face neighbours share trace coefficients. The two collapsed polar lines have single-valued displacement. Each of the four partial interfaces uses the physical height on the wing to evaluate the matching central trace at w=z/12. The common ordinate parameter is u. After multiplying by the square of the common cubic weight denominator, the trace difference has degree at most six in u and one in the wing height parameter. Equality at seven distinct u values and both wing-height endpoints therefore implies equality across that face in exact arithmetic. The script enforces these linear constraints by numerical null-space reduction and repeats the evaluations at 17 independent u locations. The four-node control underconstrains inclined interfaces: the largest held-out scalar jump for a unit coefficient vector in its admissible space is 0.0001121. With seven nodes this diagnostic stays below 8e−16.

The material is homogeneous and isotropic, E=1000 MPa and Poisson ratio 0.3. First, the free-body eigenspectrum is checked, using a rigid-mode threshold of 1e−10 times the largest eigenvalue. The affine uniaxial-stress test checks representability, energy and equilibrium against tractions integrated on the exposed boundary, subtracting the covered part of each central face. The second load case clamps the left end at x=0 and applies 1 MPa axial traction to the right end. Interface force is the sum of right-wing residual forces associated with virtual rigid translation. The reported tip observable is the traction-weighted mean x displacement of the loaded end.

| Rise (mm) | Independent scalar DOFs, ten charts | Mean end displacement, ten charts (mm) | Matching-face control, 22 charts (mm) | Relative difference |
|---:|---:|---:|---:|---:|
| −3 | 88 | 0.08004945 | 0.08316766 | −3.75% |
| 0 | 100 | 0.07678517 | 0.07937979 | −3.27% |
| +3 | 88 | 0.08004945 | 0.08316766 | −3.75% |

The matching control is built independently in the study by extending the two cutting planes through every central chart. It uses 172 independent scalar coefficients before clamping. Homogeneous control rows are checked for equality up to a constant scale at each declared interface. Both models have the same geometry, material and boundary loading. Both quadrature orders (8 and 14) agree on the end displacements to 2e−12 mm. The comparator is a second coarse approximation. This test establishes continuous displacement, force transfer and a measurable stiffness difference at these discretizations. Quantifying continuum displacement and local stress error requires analysis-space refinement. The Python assembler currently targets this union and its prescribed correspondences.
