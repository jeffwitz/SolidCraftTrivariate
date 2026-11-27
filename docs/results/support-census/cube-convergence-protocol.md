# Uniform tricubic cube convergence: preregistered protocol

## Scope and frozen choices

Before numerical runs, this protocol fixes the domain, exact field, spaces,
quadrature and acceptance bands. The geometry is the single affine parent
cube `[0,1]^3 mm`, with no fitted or subdivided geometric maps. The material
is `E=2.6 MPa`, `nu=0.3` (`lambda=1.5 MPa`, `mu=1 MPa`). The exact
displacement is

```text
u(x,y,z) = (0, 0, sin(pi*z)) mm.
```

The body force is the analytic `-div(sigma(u))`, the exact Cauchy traction is
applied on every exterior face, and all displacement components are clamped
only on the bottom face `z=0`. The top is left natural; its exact traction is
nonzero. This avoids a top essential condition that would change the weak
problem. With the cube affine, its exact volume is 1 mm³ and every geometry
Jacobian determinant is exactly 1.

The primary field is tensor-product cubic B-spline, degree 3 and interior knot
multiplicity 1 (maximum continuity C2), on uniform tensor grids with `n=1,2,4`
cells per direction. The scalar coefficient counts are therefore
`(n+3)^3 = 64, 125, 343`; the grids have `1, 8, 64` elements. A full Bernstein
Q3 C0 control uses interior multiplicity 3 on the same three grids, giving
`(3*n+1)^3 = 64, 343, 2197` scalar coefficients. Both use the same exact
geometry, exact field weights and Gauss order 8. Gauss order 12 at the finest
level is the integration sensitivity check.

The implementation uses exact rational univariate extraction from
`parent_field::axis_rows`, `ExtractedElement::sample`, the variable-size
`assemble_sampled` operator, and the established `section_manufactured`
body/traction, constrained solve and error evaluator. It stores no global
dense stiffness matrix.

## Acceptance bands fixed before runs

- Each geometry sample has positive determinant and the integrated volume is
  within `2e-12 mm^3` of 1.0.
- Exact extraction rows preserve partition of unity; the affine geometry and
  affine displacement patch are reproduced to relative `2e-11` at interior,
  face and knot locations.
- The solve residual is below `1e-9`; force/work balance diagnostics are
  below `1e-8`.
- Both displacement and energy errors decrease strictly on `n=1,2,4`.
- Between `n=2` and `n=4`, the observed energy-error order is in `[2.4,3.6]`
  and displacement L2 order is in `[3.2,4.8]` for both spaces.
- At `n=4`, C2 has fewer scalar DOFs than C0 and its energy error is at most
  1.5 times the C0 error. The numerical error is not required to improve at
  equal element count because the two spaces are not the same dimension.
- At `n=4`, q8 versus q12 changes energy and L2 errors by at most `2e-8`
  relative and volume by at most `2e-12 mm^3`.

Every failed or indeterminate criterion remains in the output; no tolerance
is changed after observing results. This one analytic solution qualifies
only the affine cube, this material, these boundary data and these uniform
spaces. It does not establish curved-geometry convergence or stress accuracy.

## Frozen asymptotic supplement

If the first strict gate fails, retain its measurements unchanged and add an
`n=4` to `n=8` comparison at the same Gauss order 8, material, exact sine
field, geometry, boundary loading, and cubic spaces. Apply the same energy
order `[2.4,3.6]` and displacement L2 order `[3.2,4.8]` bands. Report the
supplement as a separate result; it does not replace or retroactively pass the
original `n=2` to `n=4` gate. Also report the q8/q12 `n=4` sensitivity already
specified above, the exact axial body resultant `7*pi N`, top traction
resultant `-3.5*pi N`, bottom reaction `-3.5*pi N`, and exact strain energy
`7*pi^2/8 N mm`. Give the sample covariance of `log(DOF)` and `log(error)` for
each space and error norm as a descriptive diagnostic only.
