# Blind-pocket lateral continuity sweep: preregistration

## Revision and scope

- Starting revision: `1c63c1280e1bef16fcf9a533b6fe70d6e0a8d158`.
- Operative scope: see [the pre-run amendment](pocket-lateral-c1-sweep-preregistration-amendment.md). It corrects the superseded Cartesian-grid description before any numerical sweep.
- Geometry source: `tests/references/blind_bilinear_pocket.recipe.json`, built by the repository CLI from the polygonal `sketch_layer_partition_v1` family. Its stock is 60 × 40 × 20 mm; the 20 × 10 mm pocket occupies `x=[20,40]`, `y=[15,25]` and removes the upper 6 mm. This is a new exact bilinear source model. The published `curved_sketch_layers_v1` model is not used: its maps are not bilinear, and its earlier mechanical measurements remain separate.
- New files: this protocol, the explicit source recipe, `tests/references/pocket_lateral_c1_sweep.py`, its focused tests, the result JSON, and the final report. No production Rust code or existing files are changed.
- Question: for polynomial degrees p=3 and 4, compare globally C0 and globally C1 lateral fields at n=1,2,4,8 subdivisions in each direction of each of the nine 20 × (15/10/15) mm parent rectangles.

## Mechanical model and reduction

Use a manufactured anti-plane displacement `u=(0,0,w(x,y))` with
`w=sin(pi*x/60) sin(pi*y/40)` mm. For homogeneous isotropic linear elasticity, the only nonzero stresses are `sigma_xz=mu*w_x` and `sigma_yz=mu*w_y`; set `mu=1 MPa`. This is an exact three-dimensional elasticity subproblem because `w_z=0` and the geometry is a straight extrusion. Integrating the 3D weak form in z gives the exact two-dimensional coefficient `t(x,y)`, equal to 20 mm outside the pocket footprint and 14 mm inside it. The four pocket walls have natural manufactured traction equal to the outward in-plane normal derivative times the removed 6 mm height. The external vertical stock faces are clamped (`w=0`); the top, floor, and bottom have zero anti-plane traction. The body force is `-mu*Delta(w)` over the solid.

The five distinct quadrilateral charts align exactly with the source partition. Each is split into n by n parametric subcells, and the varying bilinear Jacobian is retained in all integrals. C0 uses tensor-product open B-splines with interior knot multiplicity p. C1 uses the same knots with multiplicity p-1 within each chart and physical normal-gradient constraints on all eight oriented chart interfaces; the common trace supplies tangential-gradient continuity. These constraints include the footprint edges, slanted radial edges, and their vertices. The field is constant through physical thickness; the reduced coefficient is exactly 20 mm on the four outer charts and 14 mm on the footprint chart. The four exposed pocket walls carry the manufactured traction over the missing 6 mm. This is an exactly reduced anti-plane Ritz subproblem of 3D linear elasticity on the blind-pocket solid, not unrestricted 3D elasticity or the existing full 3D native solve.

## Measurements and controls

- Primary errors: relative 3D strain-energy norm and relative displacement L2 norm, with the exact thickness weight in both norms.
- Complexity: free scalar DOFs after homogeneous exterior Dirichlet conditions; vector-equivalent DOFs are three times this count. Report energy and L2 errors against DOFs as well as against h.
- Quadrature: tensor-product Gauss-Legendre integration, initial n=1 gate at orders 12 and 16. After correcting the two archived span-integration defects, p=4/C0/n=1 narrowly exceeded the 1e-9 error-observable threshold (relative energy-error change 1.92e-9). Therefore the full sweep uses independently assembled orders 16 and 24 and retains the same 1e-9 threshold for energy error, L2 error, and strain energy. This prospective increase is applied to all sweep cases.
- Solve: sparse direct solve; require normalized KKT/free residual below 1e-10. Require a positive smallest eigenvalue of the symmetric stiffness after essential conditions (sparse Lanczos check); for C1 also require the constraint residual below 1e-9.
- Independent checks: analytic MMS derivatives and source terms; independent high-order quadrature of exact norms; exact knot-multiplicity checks; numerical reconstruction of C1 jumps at all shared edges and vertices; compare order-12 and order-16 assemblies/results.
- The first n=1 solve in each p/continuity family must pass residual, definiteness, C1, and q16/q24 checks before the remaining sweep is accepted. Any failed case remains in the denominator and is reported as a failure.

## Thresholds and exit criteria

The sweep is complete only if all 16 (p, n, continuity) states solve, have residual below 1e-10, pass positive definiteness, and meet the q16/q24 observable tolerance. Report observed rates between consecutive n values, with no pre-imposed rate gate. Compare C0 and C1 at equal p,n and nearest equal DOF; do not claim C1 is better solely from asymptotic rate. The prior native inherited-space result on the curved nine-block model remains a separate historical measurement, not calibration data for this new bilinear model or reduced solver.

## Reproduction

`python3 tests/references/pocket_lateral_c1_sweep.py --output docs/results/pocket-lateral-c1-sweep.json`

The focused test command and environment are recorded in the final report. No tests are run before this protocol is written.
