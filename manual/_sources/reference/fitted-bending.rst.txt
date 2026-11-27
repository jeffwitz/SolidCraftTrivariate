Fitted crossed-tool bending
===========================

The fitted control uses the same physical union, clamp, traction and material as
:doc:`crossed-mechanics`. Along the inclined wing the core meets only part of
each wing-facing chart. Fitting the core to the two wing faces produces 22
charts with full interface traces. This gives a mechanical control for the
ten-chart partial-interface calculation. The source geometry for the measured
comparison is the recorded S2 union with a 3 mm rise over a 100 mm length.

Field construction
------------------

Each integration cell carries a rational displacement field of degree
(3,2,2). A full chart interface shares trace coefficients. Within each chart,
uniform subdivision in directions *v* and *w* supports two field choices:

* **C0** identifies trace values at each new subdivision.
* **C1** also identifies the parametric derivative across each new subdivision.
  If the two spans have the same parameter length and the same degree, the
  last interior, seam and first interior homogeneous controls satisfy
  :math:`H_L + H_R = 2 H_M`. The corresponding field relation is
  :math:`2 w_M d_M = w_L d_L + w_R d_R`, where :math:`w` is a geometry weight
  and :math:`d` is a displacement coefficient. The native signed lift applies
  these relations to the shared elasticity operator. The 22 original chart
  interfaces retain their shared value traces.

Reproduce the comparison
------------------------

.. code-block:: sh

   cargo build --release --locked --offline
   cli=target/release/solidcraft-trivariate
   mkdir -p build
   "$cli" crossed-boolean tests/jobs/crossed_partial_rise3.recipe.json union build/inclined-cross.json
   "$cli" conforming-bending build/inclined-cross.json 2 2 8 C0 build/fitted-c0.json
   "$cli" conforming-bending build/inclined-cross.json 3 3 8 C1 build/fitted-c1.json

At order 8 the independent solutions below give the mean axial displacement
of the loaded end. The table lists free vector displacement coefficients
after imposing the clamp and the field constraints.

.. list-table::
   :header-rows: 1

   * - Field
     - Subdivisions (*v*, *w*)
     - Integration cells
     - Free DOFs
     - Mean displacement (mm)
   * - C0
     - (2, 2)
     - 88
     - 4,128
     - 0.297294
   * - C1 within charts
     - (2, 2)
     - 88
     - 2,520
     - 0.289699
   * - C1 within charts
     - (3, 3)
     - 198
     - 4,128
     - 0.298178
   * - C0
     - (6, 2)
     - 264
     - 11,112
     - 0.305584

The first and third rows have the same field dimension. The third requires
110 additional integration cells. The last row can be compared with the
matched HEX20 result of 0.306089 mm and the finer HEX20 value of 0.307772 mm.
Agreement of this integral observable does not bound a local stress error.

The native regression suite is ``tests/conforming_bending.rs``. For the
independent reference see the
`continuity study <https://github.com/jeffwitz/SolidCraftTrivariate/blob/main/docs/results/crossed-interface-iga/smooth_bending.json>`_
and the
`matched IGA and HEX20 comparison <https://github.com/jeffwitz/SolidCraftTrivariate/blob/main/docs/results/crossed-interface-iga/bending-enrichment.md>`_.
For the relation between spline continuity, field dimension and refinement
see the isogeometric analysis formulation of
`Hughes, Cottrell and Bazilevs <https://doi.org/10.1016/j.cma.2004.10.008>`_.
