Crossed-tool partial-interface mechanics
========================================

This reference describes the native elasticity check for the ten-chart S2
union of two crossing tools. Here ``S2`` names the crossed-tool geometry
experiment in ``docs/results/crossed-interface-iga``; it is unrelated to
roadmap stage labels elsewhere in this repository.

The native entry point is ``crossed_bending::analyze`` in
``src/crossed_bending.rs``. The filesystem command routes through the same
function:

.. code-block:: console

   cargo build --release --locked --offline
   mkdir -p build
   cli=target/release/solidcraft-trivariate
   "$cli" crossed-boolean tests/jobs/crossed_partial_rise3.recipe.json union build/inclined-cross.model.json
   "$cli" crossed-bending build/inclined-cross.model.json 1 1 8 build/crossed-bending-v1w1q8.json

The stored recipe is
``tests/jobs/crossed_partial_rise3.recipe.json``. It embeds the S2 sketch pair
and the positive 3 mm transverse rise. The API replays its embedded
``source_recipe`` through ``g2_cross_plate::boolean::report`` and
requires exact JSON equality with the saved model before accepting its charts.
Changing chart identifiers while leaving the source recipe untouched is
therefore rejected. The CLI publishes the analysis result atomically only
after a successful solve.

Model, field, and loading
-------------------------

The analysis raises the source geometry's degree-(3,1,1) homogeneous maps to a
Q(3,2,2) analysis basis. It subdivides the chart coordinates in ``v`` and
``w``; ``u`` remains unrefined. The implementation currently accepts
``1 <= refine_v <= 4``, ``1 <= refine_w <= 2`` and tensor Gauss orders 4
through 14. The chart geometry remains the replayed source geometry. Positive
homogeneous weights are required.

The constitutive model is small-strain, three-dimensional isotropic linear
elasticity with Young's modulus ``E = 1000 MPa`` and Poisson ratio ``nu = 0.3``.
The left end at ``x = 0`` is clamped. A uniform outward traction of
``1 MPa`` in the positive x direction is integrated over the right end; all
other boundaries are traction-free. The reported displacement is the
traction-weighted mean x displacement on that loaded end. Thus its work
conjugate is the reported total x force, ``40.2222222222 N`` for the checked
geometry.

Interior chart interfaces and within-chart subdivisions identify their full
Q(3,2,2) coefficient traces. The four wing-to-core interfaces are partial:
the wing and core parameterizations have different interface coverage. The
current route constrains these traces by overdetermined, preimage-aware
collocation and then computes a numerical nullspace of the resulting trace
matrix. It samples 21 points in ``u`` and 4 in each slave ``w`` cell. This is
not an exact mortar method, a proof of trace compatibility, or a general
nonmatching-interface algorithm. Its mapping assumptions are specific to the
recorded ten-chart S2 geometry, including the physical height relation used
to locate the master ``w`` cell. To resolve child-knot crossings, the native
path samples 65 points along each relevant edge, bisects sign-changing
brackets, then uses 21 points in each resulting ``u`` segment. At each ``u``
it splits the slave ``w`` interval at uniform child knots and at physical
``z = 12 k / refine_w`` crossings, and places four points on every resulting
subinterval.

Numerical parity
----------------

``partial_bending_study.py`` is the independent SciPy reference. It builds
its own tensor Bernstein basis, quadrature, global identification, collocated
trace matrix, singular-value nullspace and elasticity assembly. The checked
oracle data are in ``partial_bending_study.json``. Values below are for the
positive 3 mm transverse rise, ``refine_w = 1`` and quadrature order 8 unless
otherwise stated:

.. list-table:: Native parity against the independent S2 study
   :header-rows: 1

   * - ``refine_v``
     - Free vector DOFs
     - Independent tip mean (mm)
     - Native tip mean (mm)
     - Relative difference
   * - 1
     - 594
     - 0.2780814181868228
     - 0.2780814184139097
     - 8.15e-10
   * - 2
     - 1152
     - 0.2943430322223048
     - 0.2943430316017827
     - 1.74e-9
   * - 4
     - 2268
     - 0.3012717740825558
     - 0.3012717752323149
     - 3.81e-9

The two published ``refine_w = 2`` rows are also included in the independent
native CLI gate. Their partial trace matrices contain 840 rows and have rank
108, compared with 336 rows and rank 66 at ``refine_w = 1``:

.. list-table:: Positive-rise S2 with two w subdivisions
   :header-rows: 1

   * - ``refine_v``
     - Free vector DOFs
     - Independent tip mean (mm)
     - Native tip mean (mm)
     - Native trace rank
   * - 1
     - 996
     - 0.2784237912858544
     - 0.2784237914304856
     - 108
   * - 2
     - 1926
     - 0.2948563732415402
     - 0.2948563726691059
     - 108

The finest ``v=4,w=1`` row was also repeated with quadrature order 14: the
native tip mean is ``0.3012717756853457 mm`` versus ``0.3012717819154872 mm``
in the reference (relative difference ``2.07e-8``). The native q8-to-q14
change is ``4.53e-10 mm``. On the q14 run, the measured projected equilibrium
residual is ``1.06e-10``, affine patch energy error ``2.90e-13``, geometry
mismatch ``7.19e-14 mm``, integrated volume ``6799.604725802254 mm^3`` and
applied force ``40.22222222222224 N``. The matrix-free reduced stiffness
symmetry diagnostic is ``3.88e-17`` from six bilinear probes. At v=1,w=1,q8,
the dense reduced matrix's pre-symmetrization relative skew is measured
directly and is ``6.54e-16``; roundoff symmetrization is performed only after
this check. For the two v=1,w=2,q8 and v=2,w=2,q8 cases, measured symmetry
probes are ``4.86e-17`` and ``1.32e-16``, respectively; projected residuals
are ``1.03e-10`` and ``9.09e-11``.

The native test covers v=1,w=1,q8; v=4,w=1,q8/q14; and the two published
w=2 rows. It also checks a real CLI invocation and that malformed geometry
leaves an existing output file byte-for-byte intact. The dedicated CLI parity
script replays all 13 published rows, including zero and negative transverse
rise. Run the bounded native plus independent Python gate from the root:

.. code-block:: console

   cargo test --locked --offline --release --test crossed_bending
   make check-crossed-bending

Scope and failure interpretation
--------------------------------

The trace rank is determined numerically with relative cutoff ``1e-12``. The
native result reports that rank, the maximum sampled geometric gap, an affine
patch energy error, a projected equilibrium residual, and a stiffness
symmetry diagnostic. These are checks on this particular discretization and
run, not certificates of injectivity, global trace equality, mesh
convergence, or solution error. The Python study additionally reports
held-out trace samples, six rigid modes and an off-collocation displacement
jump; those independent diagnostics are not yet reproduced by the native
API. A matrix-free symmetry check uses six prescribed bilinear probes; it is
not an exhaustive matrix symmetry test.

The native route does not currently solve a manufactured-solution problem
(MMS) with a prescribed exact field and matching body force. Its affine energy
patch check and agreement with the independent S2 model are complementary
diagnostics, not an MMS substitute. They do not establish general-order
convergence for arbitrary loads or geometries.

Native parity is measured for rises ``-3, 0, +3 mm`` at ``(v,w)`` refinements
``(1,1)``, ``(2,1)``, and ``(4,1)``; positive-rise cases additionally cover
``(1,2)`` and ``(2,2)``. The finest supported inclined row for ``w=2`` in the
published study is v=2; no v=4,w=2 oracle row exists. Nonlinear material
response, contact, and arbitrary partial interfaces remain unqualified.
Singular fillets are outside this smooth crossed-tool mechanics case; their
stress convergence needs separate local refinement and an independent
reference.

Method context
--------------

The shared spline basis and geometry refinement follow the general
isogeometric-analysis formulation introduced by Hughes, Cottrell and Bazilevs
([CrossedHughes2005]_). Bézier extraction is the element-local representation used
by the native extracted-space assembly ([CrossedBorden2011]_). For nonmatching patch
interfaces, established isogeometric mortar methods formulate and analyze a
weak coupling with multiplier spaces; their stability depends on the choice
of trace and multiplier spaces ([CrossedBrivadis2015]_). The present S2 calculation
uses collocation for a fixed partial-interface geometry and does not claim
those mortar stability properties.

.. [CrossedHughes2005] T. J. R. Hughes, J. A. Cottrell and Y. Bazilevs,
   ``Isogeometric analysis: CAD, finite elements, NURBS, exact geometry and
   mesh refinement``, *Computer Methods in Applied Mechanics and Engineering*
   194 (2005), 4135--4195. DOI
   `10.1016/j.cma.2004.10.008 <https://doi.org/10.1016/j.cma.2004.10.008>`_.
.. [CrossedBorden2011] M. J. Borden, M. A. Scott, J. A. Evans and T. J. R. Hughes,
   ``Isogeometric finite element data structures based on Bézier extraction
   of NURBS``, *International Journal for Numerical Methods in Engineering*
   87 (2011), 15--47. DOI
   `10.1002/nme.2968 <https://doi.org/10.1002/nme.2968>`_.
.. [CrossedBrivadis2015] E. Brivadis, A. Buffa, B. Wohlmuth and L. Wunderlich,
   ``Isogeometric mortar methods``, *Computer Methods in Applied Mechanics and
   Engineering* 284 (2015), 292--319. DOI
   `10.1016/j.cma.2014.09.012 <https://doi.org/10.1016/j.cma.2014.09.012>`_.
