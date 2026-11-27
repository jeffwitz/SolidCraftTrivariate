Spline fields, extraction, and local refinement
===============================================

Geometry splines and analysis splines have different jobs in this codebase.
The model's rational Bernstein control nets encode the geometry. Analysis
spaces attach field basis functions and global degrees of freedom to audited
blocks and interfaces. A model may therefore preserve exact geometry while an
analysis space is only a chosen approximation space over it.

Bézier extraction
-----------------

Bézier extraction rewrites a spline basis restricted to a knot cell in terms
of a local Bernstein basis and an extraction operator. It is an algebraic
change of local representation; it does not approximate the geometry or
create a new continuity guarantee at patch boundaries. The NURBS tensor
assembly code in ``src/nurbs_tensor_assembly.rs`` supports exact extraction
on its full-multiplicity cell grid. ``C0Tensor::extract`` returns the local
homogeneous net in patch orientation. The implementation checks the local
dimensions and orientation rather than relying on sampled values.

The original extraction data-structure paper is Scott et al. [Scott2011]_,
and the NURBS-specific finite-element formulation is Borden et al.
[Borden2011]_. The project uses extraction as a local assembly representation;
the cited papers do not certify this repository's model audits or numerical
results.

.. [Scott2011] M. A. Scott, M. J. Borden, C. V. Verhoosel, T. W. Sederberg,
   and T. J. R. Hughes, “Isogeometric finite element data structures based on
   Bézier extraction of T-splines,” ICES Report 10-45 (2010).
   `Author-hosted report <https://www.ices.utexas.edu/media/reports/2010/1045.pdf>`_.

.. [Borden2011] M. J. Borden, M. A. Scott, J. A. Evans, and T. J. R. Hughes,
   “Isogeometric finite element data structures based on Bézier extraction
   of NURBS,” *International Journal for Numerical Methods in Engineering*,
   87 (2011), 15–47. `DOI 10.1002/nme.2968
   <https://doi.org/10.1002/nme.2968>`_.

THB pilot and its limits
------------------------

``src/thb_strip.rs`` implements a deliberately narrow, one-level THB pilot
over a recovered cubic NURBS strip. Its public builders are
``build_thb_strip``, ``build_thb_strip_from_recipe``,
``build_thb_strip_from_recipe_cells``, and ``build_thb_strip_for_cells``.
The marked region must be a strict interior box; source degree, rational
weights, active-basis size, and refined-leaf count are checked. The module
constructs sparse local Bernstein extraction rows for the active truncated
basis and retains the prepared-space extraction on the same leaves for
comparison. Exact nesting and homogeneous geometry identity are checked
before a ``ThbStrip`` is returned.

This is a THB B-spline construction, not a T-spline construction, and it is
not a generic adaptive solver for every model family. In particular, THB
interior refinement does not resolve nonmatching or partial interface
footprints by itself. The current support-census report records that S2's
partial footprints require interface-specific handling. No result from this
pilot should be described as a general error estimator, an automatic marking
strategy, or a proof of optimal convergence.

Truncating coarse functions over refined regions preserves a partition
of unity in the THB construction [Giannelli2012]_.

.. [Giannelli2012] C. Giannelli, B. Jüttler, and H. Speleers,
   “THB-splines: the truncated basis for hierarchical splines,”
   *Computer Aided Geometric Design*, 29 (2012), 485–498.
   `DOI 10.1016/j.cagd.2012.03.025
   <https://doi.org/10.1016/j.cagd.2012.03.025>`_. The paper defines the
   truncated hierarchical basis and its properties; it is not a qualification
   of this implementation.

S2 context and references
-------------------------

S2 in ``docs/recipe-driven-tsplines.md`` asks for variable-size extracted
spaces and a non-affine manufactured elasticity solution on a regular affine
cube. The purpose of that manufactured case is to assess field approximation
independently of geometric approximation. A localized manufactured solution
or frozen geometric marking may seed an initial adaptive study; a new general
residual estimator is explicitly not required. Its convergence bands apply to
the prescribed smooth cube and must not be carried over to arbitrary local
refinement or singular physical solutions.

For adaptive THB analysis in a different, immersed setting, Divi et al.
[Divi2022]_ derive residual indicators with stabilization and boundary terms
for their model problems. That paper is useful context for what an estimator
requires, but neither its estimator nor its claims transfer to this project.

.. [Divi2022] S. C. Divi, P. H. van Zuijlen, T. Hoang, F. de Prenter,
   F. Auricchio, A. Reali, E. H. van Brummelen, and C. V. Verhoosel,
   “Residual-based error estimation and adaptivity for stabilized immersed
   isogeometric analysis using truncated hierarchical B-splines,”
   *Journal of Mechanics*, 38 (2022), 204–237.
   `DOI 10.1093/jom/ufac015 <https://doi.org/10.1093/jom/ufac015>`_.

Reproduce the local THB checks from the repository root with::

   make check-local-thb

The target builds the release CLI, runs ``thb_strip`` and
``fieldcraft_nurbs_h1`` Rust integration tests, runs the ignored cube
convergence replay, then checks independent THB mechanics, function values,
and tractions. Generated packages and measurements belong under ``build/``.
