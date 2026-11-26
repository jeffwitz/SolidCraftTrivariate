# Tall-cylinder tensor covers and seam qualification

Baseline: main `279ec903`. The input is the published 15-cell through cylinder of radius 5 mm and height 36 mm. The two candidates preserve the square core and join opposite pairs of sectors or join all four sectors around a closed annulus. They use the original unit parameter intervals. The closed annulus identifies its two angular end traces.

## Acceptance criteria

* Rebuild the native source and both covers. Recover every homogeneous source coefficient up to one positive scale per span, all 15 source spans once per cover, with orientation preserved. Compare against independent rational basis evaluations.
* Preserve exterior faces, all core-to-annulus couplings and the explicitly identified self seam. Check sampled field traces and the assembled scalar H¹ diffusion solution against the 15-cell reference. Compare the linear system, solution and energy at relative tolerance 1e−8. Opening the seam must produce an observable change or a failed conformity check.
* Exercise both NURBS covers and the Bézier view on the public cylinder page. Scope is a C⁰ scalar H¹ discretization on this cylinder. The pre-existing geometry certificate concerns each one-sided span on its closed cube.

An additional three-dimensional elasticity test applies the exact torsion traction to the top face, clamps the bottom and leaves the cylindrical wall free. For shear modulus 1, radius 5 mm and twist rate 0.001 rad/mm, the expected end torque is G α π R⁴/2. The energy is half the torque times the total twist. Compare both integrated resultants and the discrete displacement with these analytical values. The displacement field is representable on this extruded geometry. Track quadrature sensitivity separately from the representation comparison.

## Geometry and topology

The three-volume cover has one square core and two adjoining half-rings. The two-volume cover has the same core plus one closed annulus. Its first and last angular control rows coincide. The published topology explicitly identifies the two end faces, with one pair of traces per axial layer. Exact rational extraction checks all 15 source cells and 60 evaluation witnesses per cover. Unit knot intervals preserve each source span's one-sided Jacobian and its native positive lower bound (approximately 0.21 mm³). The annulus is interpreted with its angular ends identified.

## Elasticity result

The three discretizations each have 204 vector coefficients, with 153 free after clamping the bottom. The solver evaluates rational basis gradients from each published lattice and assembles the full isotropic 3D stiffness. It reads the interpatch joins and the self seam from the cover. Under end torque, the analytical displacement is proportional to height and to distance from the cylinder axis, so it belongs to this space.

| Representation | Applied torque (N·mm) | Strain energy (N·mm) |
|---|---:|---:|
| 15 Bézier cells | 0.981747704252 | 0.017671458677 |
| 3 NURBS volumes | 0.981747704252 | 0.017671458677 |
| 2 NURBS volumes, seam identified | 0.981747704252 | 0.017671458677 |
| 2 volumes, seam coupling removed | 0.981747704252 | 0.020864050522 |

At order 16 the closed two-volume model recovers the analytical energy within 5.2 × 10⁻¹² relative. Raising quadrature to order 20 reduces that difference below 8 × 10⁻¹⁵ on the recorded machine. At order 16, omitting the three self-seam records retains all core interfaces but produces 216 vector coefficients and increases the energy by 18.07% under the same traction. The global residual and the support torque still balance. The analytical field detects the seam error even when global equilibrium passes.

[Elasticity results](../../site/box-cylinder/tall/cylinder-torsion.json) · [Scalar diffusion results](../../site/box-cylinder/tall/cylinder-iga.json)

Reproduce the geometry with `make check-published-geometry`. With NumPy and SciPy installed, `make check-cylinder-iga` runs both independent solvers and compares the numerical output with the published reports. The native generator lives in Rust. The analysis results here use the independent SciPy implementations. Fewer geometric volumes preserve the C⁰ field space and its integration spans in these tests.

## Scalar field and seam jump

An independent Poisson solve uses a spatially varying source, homogeneous values on the end disks and natural conditions on the wall. A separate numbering scheme matches coincident trace control locations. All three covers give 68 scalar basis functions. At the same five-point Gauss rule, the stiffness difference between the two-volume and 15-cell systems is below 2 × 10⁻¹⁵ in relative Frobenius norm. The scalar trace jump on the closed seam is zero at all 81 sampled points. Releasing the outer end row while preserving the core interfaces gives a maximum jump of 2.77 in the scalar field. The stiffness comparison measures equivalence of the representations at fixed quadrature; torsion above provides the independent analytical check.
