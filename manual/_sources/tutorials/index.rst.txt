Tutorials
=========

.. toctree::
   :maxdepth: 1

   solidcraft-plugin

Build and load a curved strip
-----------------------------

Start in a checkout of the repository with Cargo installed. This exercise
uses two rational sketch curves and a height of 5 mm. You will apply an axial
traction to one end and read back the support reaction.

.. code-block:: sh

   cargo build --release --locked --offline
   cli=target/release/solidcraft-trivariate
   out=build/manual-strip
   mkdir -p "$out"
   "$cli" build tests/references/nurbs_strip.recipe.json "$out/model.json"
   "$cli" audit "$out/model.json"
   "$cli" recover-spline-parent "$out/model.json" "$out/parent.json"

The source knot at 1 splits the interval from −2 to 4 into two cubic spans.
The model contains two Bézier cells. The recovered parent uses 20 homogeneous
controls. Its simple interior knot gives parametric C² continuity along the
sketch. The boundary names survive this change of representation.

.. code-block:: sh

   "$cli" refine-thb-strip "$out/model.json" tests/references/thb_interior.plan.json "$out/local.json"
   "$cli" solve-job "$out/local.json" tests/references/thb_end_traction.job.json "$out/result.json"
   python3 -c 'import json,sys; r=json.load(open(sys.argv[1])); print(r["free_vector_dofs"], r["force_n"], r["reaction_n"])' "$out/result.json"

The saved field has 544 scalar basis functions and 1,401 free displacement
coefficients after clamping. The resultant is about 9.19307 N along positive
z. The z reaction has the opposite sign. The job sets Young's modulus to
10 MPa and Poisson's ratio to zero, so the exact axial displacement for
this traction is z/10. The loaded end at z = 5 mm moves 0.5 mm.

To inspect the geometry independently of the solution, export its boundary.

.. code-block:: sh

   "$cli" export-step "$out/model.json" "$out/strip.step"

For the integration partition, basis weights and proofs behind this example,
see :doc:`../reference/spline-fields`.
