How-to guides
=============

Find a recipe your executable supports
--------------------------------------

Run ``solidcraft-trivariate capabilities`` and look up the recipe's
``family`` in :doc:`../reference/catalog`. Each example is taken from the
same source revision as the executable used to build this manual.

Build a bracket from its feature history
----------------------------------------

.. code-block:: sh

   cargo build --release --locked --offline
   cli=target/release/solidcraft-trivariate
   "$cli" import-solidcraft-sketch tests/references/solidcraft_sketch_bracket.timeline.json build/bracket.model.json
   "$cli" audit build/bracket.model.json
   "$cli" export-step build/bracket.model.json build/bracket.step

This importer retains the sketch, extrusion and bore intent. On the admitted
equal-leg history the compact construction gives ten cells. The historical
``import-solidcraft`` adapter consumes a different input format.
See :doc:`../reference/construction` for fallback conditions.

Test a continuity change on an existing part
--------------------------------------------

.. code-block:: sh

   cli=target/release/solidcraft-trivariate
   src=docs/results/s3-face-pocket.model.json
   "$cli" recover-pair-parents "$src" build/s3-pairs.json
   "$cli" export-inherited-space "$src" tests/references/matched_chart_jets.plan.json build/s3-field.json
   "$cli" audit-inherited-space build/s3-field.json
   "$cli" solve-job build/s3-field.json tests/references/s3_inherited_stress.job.json build/s3-solution.json

Compare the recovered pair count with the number of smooth field interfaces.
Then inspect ``energy_n_mm``, ``relative_residual`` and the force and reaction
vectors in the saved solution. To compare fields under a spatially varying
load, run ``make check-inherited-corpus``. Its manufactured solutions measure
error against the same physical displacement on each candidate space.

Rebuild the manual and API
--------------------------

.. code-block:: sh

   python3 -m venv .venv
   .venv/bin/pip install -r docs/rtd/requirements.txt
   make docs SPHINXBUILD=.venv/bin/sphinx-build
   python3 -m http.server --directory build/manual 8000

Open ``http://localhost:8000/``. The build generates the command and recipe
catalogue from the Rust executable and includes the crate's rustdoc. Read the
Docs uses the repository's ``.readthedocs.yaml`` to prepare the same inputs.
