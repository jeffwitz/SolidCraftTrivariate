Explanation
===========

.. toctree::
   :maxdepth: 1

   plugin-boundaries

Why count cells, supports and unknowns separately?
--------------------------------------------------

Imagine drilling a circular hole through a plate. Several rational maps can
surround the hole, each with a nonsingular parameterization. Two neighboring
cells can share a complete face and still have different transverse
derivatives there. Their traces suffice to join continuous displacement
values. Matching physical gradients requires more information about both
maps.

A NURBS patch can group several compatible Bézier spans. With a knot of
multiplicity equal to its degree, the field can remain C⁰ across the knot.
This grouping reduces the number of design objects even when assembly still
integrates exactly the same spans. Recovering the original knot vector can
change the field dimension. Removing a redundant geometry knot is a separate
operation with a coefficient-level re-insertion test.

Consider the S3 plate. Of its 191 full faces, 62 pairs reconstruct a single
Bézier parent. A weaker first-derivative test admits 73 interfaces for the
selected C¹ field. The free field dimension drops from 2,679 to 1,710 on
the recorded test, while relative energy-norm error changes from 3.6405%
to 3.6417%. This measurement explains why the library exposes the choice
of field. A different load can reveal a much larger loss of accuracy.
`S3 measurements and counterexamples <https://github.com/jeffwitz/SolidCraftTrivariate/blob/main/docs/results/inherited-continuity/README.md>`_

Why retain a source face after refinement?
------------------------------------------

A traction is prescribed on a face of the designed part. Refining a field
may split its integration domain into many leaf faces. Without the mapping
back to the named face, a solver could integrate a different load after each
edit. The strip package saves the recipe-derived model and refinement plan.
Replaying the package reconstructs each exterior leaf's source selection
before either a clamp or a traction is applied.

The reference chapters describe the separate tests for
:doc:`geometric validity <../reference/geometry>`,
:doc:`basis nesting <../reference/spline-fields>` and
:doc:`mechanical error <../reference/verification>`.

Why a block count cannot select a partition by itself
-----------------------------------------------------

A separator between two polygonal holes determines the vertices available
to both surrounding quad rings. A wide separation margin can still place
one vertex almost on the extension of a hole edge. That cell becomes nearly
flat at its corner. For a planar Jacobian with singular values
:math:`s_1\geq s_2>0`, the condition number :math:`s_1/s_2` measures
directional distortion. Sampling it in 3D on a thin extrusion can chiefly
measure thickness. Report the planar value, corner angles and certified
volume determinant alongside the number of cells and the sample grid.

The native polygonal planner compares a bounded set of construction
candidates. A low count is an optimum only when a lower bound is reached
under the same interface and layering rules. Hinz and Buffa investigate
a complementary question: given the boundary contours and a parametric
multipatch domain, how can an elliptic PDE produce a valid interior map?
That approach provides a comparison for map quality once topology has
been chosen.
`Hinz and Buffa (2024) <https://doi.org/10.1007/s00366-024-01997-x>`__.
