Verification and evidence scope
===============================

Verification asks whether the implemented equations and algorithms are solved
as intended. It is distinct from validation against physical observations.
This project combines exact or interval-based geometric audits, independent
analytic oracles, manufactured-solution mechanics, mesh-convergence studies,
and comparisons with independent geometry or finite-element software. Each
method answers a different question; no single pass should be generalized
beyond its declared domain and discretization.

Geometric and algebraic checks
------------------------------

Native Rust constructors and audits live beside one another but are separate
implementations for each family (for example, ``src/curved_sketch.rs`` and
``src/curved_sketch_audit.rs``). Family audits check declared geometry
invariants such as coverage, exact volume, interfaces, and a named Jacobian
certificate. The Python references in ``tests/references`` provide independent
oracles for selected results. A sampled grid is a diagnostic only; a
continuous certificate must name its domain, method, and bounds. If a proof or
interval subdivision reaches its budget without a conclusion, the result is
indeterminate rather than accepted.

STEP exchange and readback
--------------------------

The exporter in ``src/step_boundary.rs`` serializes the audited boundary using
STEP AP214. ISO 10303-21 specifies the clear-text exchange encoding, while ISO
10303-42 defines generic geometric and topological representation resources;
neither standard makes a producer's export self-validating. The project gate
therefore distinguishes serialization from a readback in another CAD kernel.
Even agreement of volume, area, and moments is not a complete shape-identity
test: Kuzminykh and Hoffmann construct noncongruent solids sharing such
invariants. Readback checks should include closure/validity and classifications
at selected material, void, and boundary points, with reader tolerances stated.
See `ISO 10303-21:2016 <https://www.iso.org/standard/63141.html>`_,
`ISO 10303-42:2025 <https://www.iso.org/standard/91386.html>`_, and
A. Kuzminykh and C. Hoffmann, “On validating STEP product data exchange,”
*Computer-Aided Design* 40(2) (2008), 133--138,
`DOI 10.1016/j.cad.2007.09.004 <https://doi.org/10.1016/j.cad.2007.09.004>`_.

Preserve counterexamples and regressions under ``tests/counterexamples`` and
``tests/regressions``. An invalid input that is required to be rejected is
part of coverage, and rejection must leave the prior model or published
artifact unchanged. Rebuild tracked geometry from its source recipe when the
partition changes, then run::

   make check-published-geometry

Mechanics and the Method of Manufactured Solutions
---------------------------------------------------

The Method of Manufactured Solutions (MMS) chooses a smooth displacement
field, applies the continuum elasticity operator analytically to derive the
body force, and prescribes exact boundary data. The numerical solution can
then be compared with the chosen field in displacement and energy norms. MMS
is a code-verification technique: a successful rate study checks the coded
operator, loads, boundary conditions, quadrature, and solver together for
that case. It does not establish accuracy for arbitrary loads or validate a
physical material model.

The native T-junction example ``examples/manufactured_t.rs`` is a concrete
regression of a two-rate Q2 mesh. It uses a nonzero affine displacement plus
a smooth trigonometric perturbation on ``[0,2] × [0,2] × [0,1]``; the right
half is refined twice as finely in y, and the coarse interface is restricted
exactly onto two fine slave faces. The recipe and acceptance bands are recorded
in ``docs/campaign-run-ledger.md`` under M1. In particular, it requires an
affine patch test, a reduced residual check, observed L2 and H1 rates, and a
quadrature comparison. This is evidence for that mesh family and manufactured
field, not universal convergence.

For the S2 variable-size extracted-space roadmap study, use the regular affine
cube and its preregistered refinement levels and error norms in
``docs/recipe-driven-tsplines.md``. Do not substitute a singular geometry and
retain the smooth-case rate thresholds: the roadmap explicitly excludes that
inference.

Singular fillet caveat
----------------------

The four-sector annulus fillet has analytically identified Jacobian zeros on
four circles (at the two end planes and the stated radii). The proof of volume
and quadrant coverage does not make that map locally regular at those sets.
Consequently, results for the singular fillet are experimental evidence for
the declared parameterization and oracle, not qualification as a regular IGA
chart. See ``docs/annulus-original-edge-states.md`` and
``docs/results/inherited-continuity/README.md``. Do not use good volume
agreement to erase the singularity or imply mechanics qualification there.

Independent mechanics comparisons
---------------------------------

Comparisons among native trivariate, tetrahedral, or immersed spline solves
are meaningful only when the physical domain, boundary conditions, material,
loads, and observables are aligned. Report those definitions along with mesh
or basis size, quadrature, solver tolerance, residual, norm definition, and
resource limits. Different parameterizations can have different discrete
spaces even when their audited solid is the same. A field discrepancy is
therefore a discretization comparison, not necessarily a geometry defect.
The repository's campaign protocol and machine-readable measurements are
linked from ``docs/campaign-run-ledger.md``; preserve failed, timed-out, and
indeterminate trials in the denominator and state their limits.

As one current mechanics comparison, the crossed-tool bending study reports
0.305584 mm for the audited rational IGA model and 0.306089 mm for a matched
HEX20 mesh, both with 11,112 free vector degrees of freedom; their displacement
difference is 0.165% of the HEX20 value. A finer HEX20 reference with 125,004
free vector degrees of freedom gives 0.307772 mm. The study explicitly
qualifies global compliance for this geometry and loading, not local stress
accuracy or runtime. Reproduce the selected runs with
``make check-crossed-bending`` and see the
`measurement report <https://github.com/jeffwitz/SolidCraftTrivariate/blob/main/docs/results/crossed-interface-iga/bending-enrichment.md>`_.

Reproducible project gates from the repository root are::

   cargo test --locked --offline
   python3 -m unittest discover -s tests -v
   make check

``make check`` runs the Rust and Python suites, demonstrations, and published
geometry checks. Specialized mechanics studies have separate targets such as
``make check-local-thb`` and ``make check-cylinder-iga``; these should be
reported by name rather than implied by the general gate.

Primary verification reference
------------------------------

For the relation between a manufactured load and the exact field it
prescribes, see [Roache2009]_.

.. [Roache2009] P. J. Roache, “Fundamentals of Verification and Validation,”
   *Hermes*, 2009, Chapter 3: “Method of Manufactured Solutions.”
   `Author's manuscript <https://fenix.tecnico.ulisboa.pt/downloadFile/2815368242434470/FVV_Roache_2009.pdf>`_.
   MMS derives exact source terms from an assumed field and is used to check
   numerical implementation and observed order; it is not physical validation.
