Geometry representation and certificates
========================================

Rational tensor-product maps
----------------------------

The native crate represents a volume block as a rational tensor-product
Bézier/NURBS map over the closed unit cube. Control points are homogeneous;
positive weights are required where the family contract says so. Physical
coordinates and dimensions are in millimeters. Parametric coordinates are
dimensionless. The public primitive API starts with
``RationalBezierMap::new`` and ``evaluate`` in the
`Bézier map implementation <https://github.com/jeffwitz/SolidCraftTrivariate/blob/main/src/bezier.rs>`_;
the crate documentation is the authoritative type-level API.

A model contains both geometry and topology. Geometry includes degrees, knots,
homogeneous controls and span domains. Topology includes block sides, oriented
face incidences, boundary status, provenance, and persistent identities where
the family provides them. Geometry equality alone does not define an
interface: restrictions must match under the recorded parameter orientation
and, for nested traces, under a certified refinement relation.

NURBS span extraction and exact simplification
-----------------------------------------------

The curve utilities in the
`NURBS curve module <https://github.com/jeffwitz/SolidCraftTrivariate/blob/main/src/nurbs_curve.rs>`_
extract active Bézier spans by exact knot insertion in integer/rational
arithmetic. The source degree and knot multiplicities are preserved as
provenance. A common integer scale represents rational knots; homogeneous
controls may be rescaled by a common nonzero factor without changing the
represented curve. Degree elevation and span subdivision are exact coefficient
transformations when they remain within checked arithmetic.

The crate also contains
`knot-reduction utilities <https://github.com/jeffwitz/SolidCraftTrivariate/blob/main/src/nurbs_knot_reduction.rs>`_.
Reduction is a separate operation from
exact extraction: a reduced representation must be accompanied by its
declared error contract and audit. A representation is not called exact merely
because sampled points agree. Unsupported degree, malformed/non-monotone knot
vectors, nonpositive weights where disallowed, inconsistent gauges, and
integer overflow are rejection cases. NURBS theory and the standard extraction
operations are covered in Piegl and Tiller, *The NURBS Book*, 2nd ed.,
Springer, DOI <https://doi.org/10.1007/978-3-642-59223-2>. This reference
describes the general representation; it does not certify any particular
SolidCraft model.

For a concrete bounded family, see ``tests/references/nurbs_strip.recipe.json``
and the
`NURBS strip reference <https://github.com/jeffwitz/SolidCraftTrivariate/blob/main/docs/nurbs-sketch.md>`_. The strip constructor requires compatible
monotone boundaries; it does not claim general reparameterization of arbitrary
NURBS profiles. For crossed profiles, see
``tests/references/nurbs_cross_oval.recipe.json`` and
``tests/references/nurbs_cross_bores.recipe.json``. The latter's exact event
contract is intentionally narrower than general algebraic surface
intersection.

Certificates and their domains
------------------------------

SolidCraft separates the following claims:

* **Representation fidelity:** the map equals its declared source geometry,
  or satisfies an explicitly bounded approximation contract.
* **Local regularity:** the denominator and Jacobian determinant stay strictly
  positive on a stated closed parametric domain.
* **Injectivity:** the map is one-to-one on that domain. A positive Jacobian
  alone is local and does not establish global injectivity.
* **Trace compatibility:** adjacent maps induce the same oriented geometry on
  their shared face, possibly after a supported exact subdivision.
* **Coverage:** the union of the active block images equals the target
  material, with disjoint interiors and the required boundary.
* **Analysis suitability:** the resulting discrete space is conforming and
  meets the stated integration or mechanics checks. This is not implied by
  geometric validity.

Every certificate must name its method, domain, assumptions, bounds, and
status. Bernstein-form interval bounds can prove polynomial positivity when
all coefficients on a subdomain are strictly positive. Subdivision may make
the bound decisive; if a fixed depth or arithmetic budget ends first, the
result is indeterminate rather than a proof. The underlying positivity methods
are established in the literature. Powers and Reznick bound the degree needed
in Pólya-type representations for polynomials positive on polyhedra.
Boudaoud, Caruso, and Roy study certificates in the Bernstein basis, including
degree elevation and subdivision as ways to obtain decisive coefficients.
These are related mathematical results, not specifications of the exact
interval tests implemented here. SolidCraft's certificates are the code and
recorded audit outputs for each family. See V. Powers and B. Reznick, “A new
bound for Pólya's theorem with applications to polynomials positive on
polyhedra,” *Journal of Pure and Applied Algebra* 164(1--2) (2001), 221--229,
`DOI 10.1016/S0022-4049(00)00155-9 <https://doi.org/10.1016/S0022-4049(00)00155-9>`_,
and F. Boudaoud, F. Caruso, and M.-F. Roy, “Certificates of Positivity in the
Bernstein Basis,” *Discrete & Computational Geometry* 39(4) (2008), 639--655,
`DOI 10.1007/s00454-007-9042-x <https://doi.org/10.1007/s00454-007-9042-x>`_.

The `spatial certificate implementation <https://github.com/jeffwitz/SolidCraftTrivariate/blob/main/src/spatial_certificate.rs>`_,
`planar checks <https://github.com/jeffwitz/SolidCraftTrivariate/blob/main/src/planar.rs>`_,
`annular audit <https://github.com/jeffwitz/SolidCraftTrivariate/blob/main/src/annular_audit.rs>`_,
`polygon-prism audit <https://github.com/jeffwitz/SolidCraftTrivariate/blob/main/src/polygon_prism_audit.rs>`_,
and `curved-sketch audit <https://github.com/jeffwitz/SolidCraftTrivariate/blob/main/src/curved_sketch_audit.rs>`_
implement these checks. A positive sampled
determinant is diagnostic only. ``tests/test_planar_certificate.py`` and
``tests/test_annular_certificate.py`` provide independent reference checks;
the Rust tests exercise production APIs.

Exact polygonal input contract
------------------------------

``sketch_layer_partition_v1`` uses integer millimeter vertices and straight
segments. It accepts simple contours subject to the stock, count, and cut
limits described in the project README. It does not treat a polygon as a
sampled curve. Curved contours use their own families, with explicit degrees,
knots, homogeneous controls, and parameter conventions. The polygonal source
recipe and selected model preserve loop and feature identifiers; the audit
does not silently repair invalid contours.

Local validity is not global coverage
--------------------------------------

For a block map ``F``, a strict lower bound on ``det(DF)`` over the closed cube
establishes orientation-preserving local regularity. It does not rule out
global self-overlap. Likewise, matching every declared interface does not
prove that blocks cover the intended solid. Coverage checks compare exact
section areas or family-specific analytic measures and incidence structure;
independent point classification and STEP readback are supplemental oracles.
An oracle's positive sample set is evidence for those samples, not a universal
proof.

Analysis export
---------------

``export-analysis`` creates a versioned conforming space package for the
supported geometry family and requested quadrature order. Geometry degree,
trace orientation, exact degree elevation, common DOF numbering, and volume
quadrature are recorded according to the package schema. Analysis suitability
is then checked independently through affine-field reproduction, interface
rank, equilibrium residual, and mechanics benchmarks where provided. Geometry
certificates do not automatically transfer to a different analysis basis or
quadrature rule. See the `NURBS space compiler <https://github.com/jeffwitz/SolidCraftTrivariate/blob/main/src/nurbs_space.rs>`_,
`tensor assembly <https://github.com/jeffwitz/SolidCraftTrivariate/blob/main/src/nurbs_tensor_assembly.rs>`_,
`analysis space <https://github.com/jeffwitz/SolidCraftTrivariate/blob/main/src/analysis_space.rs>`_,
and `solver package <https://github.com/jeffwitz/SolidCraftTrivariate/blob/main/src/solver_package.rs>`_.

Bernstein bounds and rational maps
----------------------------------

Johnen, Remacle and Geuzaine expand the polynomial Jacobian determinant in
a Bernstein basis and subdivide to tighten coefficient bounds. This gives
a direct reference for the closed-element validity test. For rational
SolidCraft charts the denominator also enters the test: a positive weight
bound permits the homogeneous numerator to determine the sign. A failed
strict-positivity proof at the subdivision budget is recorded separately
from a demonstrated inversion. A boundary-singular fillet must also be
identified separately from a regular chart.
`Johnen et al. (2013), Journal of Computational Physics <https://doi.org/10.1016/j.jcp.2012.08.051>`__.
