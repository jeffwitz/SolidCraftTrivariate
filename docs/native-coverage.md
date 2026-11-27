# Native and Python production coverage

This inventory distinguishes native construction from Python orchestration and
independent verification. Rust is the production kernel and the CLI delegates
to the same crate APIs. Python code under `tests/references` is retained as an
independent oracle; porting it does not make it a substitute for that oracle.

## Ported census operation

`scripts/census_transverse_transitions.py` classifies every model interface in
the support census. `scripts/check_polygon_join_mechanics.py` now obtains its
selection from the Rust CLI and compares the complete report against this
independent Fraction implementation before using the selected joins. The Rust crate now
exposes `transverse_transitions::classify`, and the CLI command
`classify-transverse-transitions MODEL OUTPUT` writes its per-interface JSON
records. The operation identifies axial/nonlateral links, unsupported frames,
and exact constant or affine first-jet transitions in either direction. A
successful first-jet identity is not a global reparameterization certificate.

The classifier reuses `constraint_system::Ratio` for exact rational arithmetic,
including checked cross-cancellation in products and checked `i128`
intermediates.
For a block whose homogeneous weights are one shared positive constant, it
normalizes the control net by that weight before forming jets; this is an exact
projective rescaling. Inputs that exceed checked arithmetic capacity currently
return an error and are not classified. The Python `Fraction` implementation
has unbounded integers and remains the independent reference for those cases.
Differential replay matched all records, statuses, and rational coefficients
for 70 models: the 68 support-census models (3,136 interfaces) and both named
models (50 interfaces). No arithmetic overflow occurred in this corpus.

For supported lateral polynomial sections of degree one through three extruded
linearly in z, with positive constant weights on each control net, each record also carries a separate `as_g1` result. It constructs exact
determinant-syzygy factors, reduces their polynomial gcd, and verifies the
split of total beta into two linear side factors. Admission checks coprime
alpha factors, positive alpha product, and nonvanishing section Jacobians on
the closed edge. The independent Python `Fraction` census agrees record for
record with the native CLI on all 68 polygon models. This is an interface-local
certificate; the constant/affine transition classification remains a separate
diagnostic. The [bilinear census](results/transverse-as-g1.md) and [curved pocket](results/curved-pocket-fields/README.md) give the admitted domains and mechanical results.

## Other Python invoked from project workflows

The Python commands in `Makefile` fall into these groups:

| Workflow | Role | Native relation |
| --- | --- | --- |
| `trivariate` package and its unit tests | Historical affine and annular reference workflow | The crate and `build`/`audit` CLI support the corresponding versioned families. Python provides comparison data. |
| `check_*`, `audit_*`, `census_*`, `measure_*` scripts | Corpus checks, independent calculations, diagnostics, and campaign measurements | These are validation/oracle programs, not replacement geometry builders. Several deliberately use exact `Fraction`, SciPy, or independent field assembly. |
| `build_box_cylinder_tall_patches.py` | Publishes native patch covers and checks the archived oracle | Rust chooses the cover via `plan-patches` or `cover-cylinder` and emits the homogeneous controls. Python verifies source-cell extraction and prepares the gallery files. |
| `prepare_manual.py`, gallery and plot scripts | Documentation, publication, and visualization support | Not production kernel functionality. |

## Other native production operations

The following production capabilities are implemented in Rust and exposed
through the CLI. The Python scripts and recipes around them prepare inputs,
replay examples, or serve as independent checks; they are not the production
constructors.

| Capability | Native entry point | Role / scope |
| --- | --- | --- |
| Crossed-tool bending solve | `crossed-bending` (`src/crossed_bending.rs`) | P(3,2,2) displacement field on the ten-chart union. With a +3 mm rise, four transverse spans, one thickness span and order 14, the native mean tip displacement is 0.301272 mm at 2,268 free vector DOFs. The separate published equal-DOF smooth-field comparison uses a different analysis space. |
| Fitted crossed-tool bending | `conforming-bending MODEL NV NW Q C0|C1 OUT` (`src/conforming_bending.rs`) | Exact recorded +3 mm S2 source only, fitted to the 22-chart conforming control. Q(3,2,2) with C0 coefficient joins on all chart faces; C1 removes uniform internal v/w cell knots only, so cross-chart traces remain C0. At q8, C0 2×2 gives 4,128 vector DOFs and 0.2972940001 mm; C1 2×2 gives 2,520 DOFs and 0.2896992614 mm. Equal-DOF C0 2×2 / C1 3×3 rows are 0.2972940001 / 0.2981775545 mm at 4,128 DOFs (88 / 198 cells). The separate C0 6×2 row reproduces the matched-IGA target 0.3055843318 mm at 11,112 DOFs (264 cells); this is distinct from the HEX20 reference. Geometry seam compatibility uses the recorded scaled floating tolerance; scalar C1 weights and knot-removal constraints use exact `Ratio` arithmetic. |
| Curved-cell patch planning | `plan-patches` (`src/gallery_patch_planner.rs`) | Plans candidate patch covers for the supported input family; bounded search and family admission still apply. |
| Axial composition | `compose-axial` (`src/composed_axial.rs`) | Composes aligned rectangular axial features with recorded source intervals and face provenance. |
| Cylinder covering | `cover-cylinder` (`src/cylinder_tensor_covers.rs`) | Assembles the supported square-core and angular-ring layouts, preserving source intervals. |
| Annulus boundary layouts | `annulus-edges` (`src/annulus_edges.rs`) | Produces the named linear, quadratic, and prepared layouts under the documented sharp/chamfer/fillet contracts. The original quadratic fillet has singular Jacobian rails; the prepared charts have their own certificates. |
| Transverse interface census | `classify-transverse-transitions` (`src/transverse_transitions.rs`) | Reports per-interface axial/frame status and supported exact constant/affine first-jet transitions. It does not certify a global reparameterization. |

## Deliberate research and generality limits

These are open generality or research problems, not identified Python-only
production operations awaiting a routine port:

* A complete arrangement algorithm for arbitrary algebraic/NURBS curves and
  unrestricted degeneracies is broader than the supported bounded sketch
  families and event contracts.
* Arbitrary partial, curved, or non-nested multipatch coupling is outside the
  exact whole-face and nested-trace patterns currently implemented.
* The one-level THB strip pilot is not a general residual estimator, automatic
  marking strategy, or proof of optimal adaptive convergence.
* Global feature identity matching across arbitrary topology changes and
  unrestricted corner-valence optimization are not implemented as general
  solvers.
* No general solver is claimed for topological hexahedral extension or for
  turning every admitted boundary into a rational tensor-product volume
  partition.
