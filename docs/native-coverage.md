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

## Other Python invoked from project workflows

The Python commands in `Makefile` fall into these groups:

| Workflow | Role | Native relation |
| --- | --- | --- |
| `trivariate` package and its unit tests | Historical affine and annular reference workflow | The crate and `build`/`audit` CLI support the corresponding versioned families. Keep Python as a reference. |
| `check_*`, `audit_*`, `census_*`, `measure_*` scripts | Corpus checks, independent calculations, diagnostics, and campaign measurements | These are validation/oracle programs, not replacement geometry builders. Several deliberately use exact `Fraction`, SciPy, or independent field assembly. |
| `build_box_cylinder_tall_patches.py` | Builds/checks patch input for the native patch assembler | The actual patch assembly is exposed in Rust and through `assemble-patches`; Python prepares and verifies the fixture. |
| `prepare_manual.py`, gallery and plot scripts | Documentation, publication, and visualization support | Not production kernel functionality. |

No other Python-only production operation was identified among the Python
commands used by `Makefile`. The remaining native coverage limit is arithmetic
range in `transverse_transitions`: checked `i128` is finite, unlike Python's
arbitrary-precision fractions. The measured census corpus had no overflow.
Other Python validators should stay independent unless a separate production
operation is identified and its contract is pre-registered.

## Other native production operations

The following production capabilities are implemented in Rust and exposed
through the CLI. The Python scripts and recipes around them prepare inputs,
replay examples, or serve as independent checks; they are not the production
constructors.

| Capability | Native entry point | Role / scope |
| --- | --- | --- |
| Curved-cell patch planning | `plan-patches` (`src/gallery_patch_planner.rs`) | Plans candidate patch covers for the supported input family; bounded search and family admission still apply. |
| Axial composition | `compose-axial` (`src/composed_axial.rs`) | Composes aligned rectangular axial features with recorded source intervals and face provenance. |
| Cylinder covering | `cover-cylinder` (`src/cylinder_tensor_covers.rs`) | Assembles the supported square-core and angular-ring layouts, preserving source intervals. |
| Annulus boundary layouts | `annulus-edges` (`src/annulus_edges.rs`) | Produces the named linear, quadratic, and prepared layouts under the documented sharp/chamfer/fillet contracts. The original quadratic fillet has singular Jacobian rails; the prepared charts have their own certificates. |
| Transverse interface census | `classify-transverse-transitions` (`src/transverse_transitions.rs`) | Reports per-interface axial/frame status and supported exact constant/affine first-jet transitions. It does not certify a global reparameterization. |

The census consumer `scripts/check_polygon_join_mechanics.py` now obtains the
selection from this Rust CLI and gates use on equality with its independent
Python `Fraction` report. The parity replay covered 70 models and 3,186
interfaces: 68 support-census models (3,136 interfaces) and two named models
(50 interfaces). All statuses and rational coefficients matched; no overflow
occurred. The exact `i128` boundary remains finite and is not inferred from
this corpus: out-of-range input is an error. No parity claim is made for
inputs that exceed that arithmetic range.

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

Each limitation should remain distinct from exact arithmetic overflow in the
transverse census: the former are scope/generalization boundaries, while the
latter is an implementation range limit with a currently passing corpus.
