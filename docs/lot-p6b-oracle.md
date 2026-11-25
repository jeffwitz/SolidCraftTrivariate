# P6b independent qualification protocol

This protocol freezes the independent checks for T-P6 and P7 before comparing
them with the Rust implementation. The reference scripts use only Python's
standard library and do not call the Rust constructors or copy their formulas
for the circle jets.

## Provenance and scope

- T-P6 source baseline: `ff2c36f`; implementation under review:
  `src/reparameterization.rs`.
- P7 source baseline: `7b61094` plus T-P6; implementation under review:
  `src/sketch_intersections.rs`.
- The checkout revision used to create a report is recorded in that report.
- This qualifies the two bounded analytic templates. It does not accept a
  sketch intersection as a volumetric feature or claim general CAD coverage.

## T-P6 oracle

For `a0+a1=0`, set `a=a0` and `a1=-a`. Direct polynomial integration gives
`J(a)=w0/3*(a-theta)^2 + w1/3*(-a+theta)^2 + 2*lambda*a^2`.
The unconstrained solution is `a*=rho*theta`, where
`rho=((w0+w1)/3)/((w0+w1)/3+2*lambda)`. The admissible interval is
`[-3/4,3/4]`. The reference independently clips the scalar minimizer, computes
the gradient, bound multiplier, tangent curvature, and evaluates the physical
map `x_i(u)=i+u+a_i*u*(1-u)`. It checks centered differences only in smooth
interior/strictly-active regimes. At activation, the domain map is continuous
and a unique two-sided derivative is intentionally not asserted. In the
32-epsilon uncertainty band the oracle requires the implementation to report
uncertainty and withhold sensitivity.

Frozen corpus: asymmetric interior weights, symmetric interior, positive and
negative strict activation, both exact activation thresholds, the adjacent
floating-point values on either side, regularization-only, flat nonunique,
zero weight on one side, invalid nonfinite/negative inputs, underflowed data
weight, and an overflowing finite objective. Acceptance targets are KKT
residual `<=1e-12`, map-jet centered-difference error `<=2e-8` for steps
`1e-4,5e-5,2.5e-5`, exact brick endpoints, minimum map speed `>=1/4`, and
exact regime/availability distinctions at thresholds.

As a supplementary geometry check after the frozen algebraic oracle, the P6
script also integrates `X^2+Y*Z` over a fixed non-affine sheared physical image
under the two reparameterized charts. Gauss orders 8 and 12 are compared at
`theta=0.2` and centered perturbations `±1e-3`; the exact physical sensitivity
is zero because the map remains a partition of the same image.

## P7 oracle

Inputs are two `[cx_um,cy_um,r_um]` integer triples. Classification is done
from exact integer `D²`, `(r1+r2)²-D²`, and `D²-(r1-r2)²`. For transverse
cases, the reference obtains the point from center-axis projection and a signed
perpendicular offset. It obtains each of the six point derivatives by solving
the 2-by-2 implicit system formed by the two circle constraint gradients,
not by differentiating the production closed form. The condition diagnostic
is independently calculated from the singular values of those rows.

Frozen corpus: all seven relation classes, external/internal tangency and
one-micrometre neighbors on either side, axis-aligned and oblique transverse
cases, parameter swaps, translation, quarter-turn rotation, reflection,
near-coincident maximum-radius circles at one-micrometre separation, center
and radius bounds, zero/negative radius, one-past-limit inputs, and `i64`
extremes. Ordinary-
scale transverse points must satisfy both circle constraints within
`2e-12 mm^2`; large-coordinate cases also report the residual relative to
squared radius to expose binary64 cancellation. Implicit jets must satisfy the
differentiated constraints within `2e-11 mm`.
One-micrometre centered differences must agree within `3e-6 mm/mm`; the
second-order Taylor ratio at steps 100 and 200 micrometres must be within 0.15
of four. Degenerate cases must have exact relation/cardinality and no jets.
