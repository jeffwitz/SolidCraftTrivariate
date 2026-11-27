# Curved pocket mechanics protocol

## Frozen input and scope

- Frozen protocol revision: [00ca106](https://github.com/jeffwitz/SolidCraftTrivariate/commit/00ca106abe88ccc89b92cb27a97fe24b0100f581).
- Geometry: the published `plate_blind_pocket` recipe and model under
  `docs/site/data/plate_blind_pocket/`. The audited model has nine blocks,
  degree `(2,2,1)`, and a common positive weight on each control net.
- Existing bilinear benchmark and its public `solve(p,n,c1,order)` entry point
  remain unchanged in behavior. The curved study is an additional path.
- Study parameters: transverse field degree `p=3,4`, uniform XY spans
  `n=1,2,4,8`, both C0 and complete physical C1 across every shared chart edge.
  The axial field factor is exactly linear in physical `z`; no axial refinement.

## Fields, material, and loads

Primary manufactured displacement:

`u=(0,0,(z/20)*((x+2*y+10)/150)^5)` mm.

The affine control and reproduction check is
`u=(0,0,(z/20)*((x+2*y+10)/150))` mm. The p4 reproduction check also uses
`(z/20)*((x+2*y+10)/150)^2`; on a biquadratic chart this composes to degree at
most `(4,4,1)`.

Use the existing isotropic material convention: `E=2.6 MPa`, `nu=0.3`,
equivalently `lambda=1.5`, `mu=1` in solver units. Prescribe zero displacement
on the lowest exterior face. Apply the analytical `-div(sigma)` volume load
and `sigma*n` traction on all remaining exterior faces, including pocket walls
and floor. Include the constrained-face traction only in the independent
unrestricted virtual-work identity.

## C0/C1 constraints and evidence

C0 shares oriented edge value and tangent jets at every knot, plus the required
interior edge coefficient for quartic spans. C1 additionally equates a
physical gradient component on every open span. For the rational physical
gradient equations, clear the nonzero chart Jacobian denominators: with a
biquadratic map, the resulting univariate numerator degree is bounded by
`p+5`; therefore use at least `p+6` distinct points per span. Report this
bound, point count, exact constraint rank, and residual of the reconstructed
physical gradients at independent points including both span endpoints.

## Acceptance and independent oracle

- Reject any chart with unsupported degree, inconsistent homogeneous XY
  controls between z endpoints, nonconstant/zero weights, noninteger corner
  topology coordinates, or a nonpositive sampled Jacobian.
- Require the integrated volume to match the published model volume within
  `1e-10` relative and quadrature energy/load/work outputs to stabilize between
  the two highest available orders within `1e-7` relative.
- Require relative equilibrium residual below `1e-9` and unrestricted virtual
  work defect below `1e-7`.
- Compare geometry, exact field, energy, body resultant/moment, and separate
  top, wall, and floor traction integrals against the independent polynomial
  integration implemented by the coordinating agent. Do not treat quadrature
  self-convergence alone as an independent oracle.
- Preserve historical bilinear results and report DOFs, rank, quadrature order,
  error norms, residuals, wall/floor loads, and elapsed time for every run.

## Files owned by the mechanics implementation

`src/pocket_convergence.rs`, the compatible manufactured-profile extension in
`src/section_manufactured.rs`, `tests/pocket_convergence.rs`, and the matching
Rust module registration if needed. No commits or changes to geometry builders,
published model assets, or unrelated mechanics callers.
