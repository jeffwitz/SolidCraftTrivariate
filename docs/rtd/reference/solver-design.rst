Adaptive solves, pressure and design derivatives
=================================================

Marking and transfer
--------------------

An adaptive Q2 leaf carries three contributions to its squared indicator:
the volume residual, the traction jump on interior faces and the natural
boundary residual. The data structure also records whether an error upper
bound is available. Here it is false. The value depends on the current
mesh, loads and discrete solution.

For a fraction :math:`0<\theta\leq1`,
``AdaptiveResidual::bulk_mark`` sorts the squared indicators in descending
order and selects the shortest prefix containing at least that fraction of
their sum. Leaf identifiers order ties. A vanishing indicator selects no
leaves. Dörfler proved convergence for an adaptive Poisson algorithm.
The kernel uses its marking rule in elasticity. Error reduction on an
elastic part must be measured separately.
`Dörfler (1996) <https://doi.org/10.1137/0733054>`__.

The constrained multilevel solver assembles a coarse operator
:math:`P^T K P` using the nested Q2 transfer. It factors the coarsest
operator by Cholesky and adds an inverse reduced-diagonal correction at each
finer level. All levels use the same homogeneous mounting. An unsupported
mount or an inconsistent transfer fails during construction. For timings,
include preconditioner setup, iterations and assembly in addition to
the residual. For the linear iteration underlying the solver, see
`Hestenes and Stiefel (1952) <https://doi.org/10.6028/jres.049.044>`__.

Nearly incompressible materials
-------------------------------

As Poisson's ratio approaches one half, a displacement discretization can
overconstrain the volume strain. The optional ``projected_bulk`` operator
projects the Q2 displacement divergence onto an element-local Q1 pressure
space with eight coefficients. It condenses the pressure variables before
the global solve. Its recovery operator returns coefficients of mean stress
with tension positive, in MPa. Compare reaction equilibrium, the true 3D
residual, and a displacement or stress observable against a reference under
the same loading. Elguedj, Bazilevs, Calo and Hughes derive small- and
finite-strain projection methods for higher-order NURBS. The kernel's
qualified law in this module is small-strain isotropic elasticity.
`Author manuscript <https://www.ices.utexas.edu/media/reports/2007/0726.pdf>`__.

Design derivatives
------------------

``differentiate-job MODEL JOB DIRECTION OUTPUT`` holds the chart
connectivity and quadrature rule fixed while differentiating geometry,
material, loads and supported constraint equations. For a reduced
equilibrium :math:`K_r q=f_r`, the tangent satisfies
:math:`K_r q'=f'_r-K'_r q`. For one scalar observable the adjoint uses the
transpose of the same reduced operator. Moving prescribed values and trace
coefficients enter through the derivative of the lift. Record both
objective derivatives and replay the central difference inside the same
topological regime. A change of active features calls for a new regime.
Wall, Frenzel and Cyron describe shape sensitivities with NURBS design
variables in isogeometric structural optimisation.
`Wall et al. (2008) <https://doi.org/10.1016/j.cma.2008.01.025>`__.

Executable checks
-----------------

.. code-block:: sh

   cargo test --locked --offline --test multilevel
   cargo test --locked --offline --test adaptive_enrichment
   cargo test --locked --offline --test mechanical_sensitivity

The symbol-level contracts are in the crate modules
``adaptive_estimator``, ``multilevel``, ``projected_bulk``,
``elasticity_sensitivity`` and ``mechanical_sensitivity``.
