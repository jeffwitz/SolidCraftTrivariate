FieldCraft analysis exchange
============================

The desktop bridge connects a SolidCraft command plugin to the native Rust
mechanical solver and to FieldCraft's results viewer. The plugin's WebAssembly
component can call declared SolidCraft commands and exchange files in the
approved companion directory. The desktop launcher runs the native processes.

Domain and load
---------------

A solve takes an explicit trivariate **recipe** and a **mechanical job**. The
job selects support and traction surfaces using feature provenance and, where
specified, a plane. An imported STEP boundary alone has no trivariate basis.
A recipe solve concerns the named model even when a different CAD document
is active. After a host-only Boolean, fillet or chamfer, the displayed body can
differ from the last complete trivariate model.

The reference bracket has ten blocks and eleven internal interfaces. The
sharp stock is clamped on the outer profile at :math:`y=0` and loaded on the
outer profile at :math:`y=50` mm. A uniform traction of 6.25 MPa points in
positive :math:`x`. The linear elastic material has :math:`E=210000` MPa and
:math:`\\nu=0.3`. With the Q2 field, order-eight Gauss integration and no
refinement, the run has 171 scalar unknowns before essential constraints,
a 2000 N applied resultant, a support reaction of approximately −2000 N in
:math:`x`, and 135.44475 N mm of strain energy. These values identify this
discretization. See :doc:`solver-design` for the mechanical formulation.

Fields and units
----------------

.. list-table::
   :header-rows: 1
   :widths: 26 28 46

   * - Quantity
     - Native result or VTU
     - FieldCraft exchange
   * - Position and displacement
     - mm
     - Converted to metres in the archive.
   * - Cauchy stress and von Mises stress
     - MPa
     - Converted to Pa in the archive.
   * - Reaction and energy
     - N and N mm
     - Read global equilibrium and energy from the native result JSON.

Each sampled hexahedron is a display element within a volume chart. Samples
on opposite sides of an interface can have different one-sided stresses.
The displacement vector and stress fields are evaluated on the solved charts.
A stress maximum on this plotting grid is a sampled maximum.

Validation
----------

A replay should check three separate results: the solver's residual and
reaction balance, agreement of the sampled fields after a FieldCraft
read/export round trip, and a contour drawn by the released FieldCraft CLI.
Each catches a different failure, including wrong unit conversion or loss
of the solved result during archive loading. See the
`FieldCraft 0.8.1 results codec <https://github.com/tangency-labs/fieldcraft/blob/v0.8.1/crates/io/src/results.rs>`_
and the `SolidCraft plugin build <https://github.com/jeffwitz/SolidCraftTrivariate/actions/workflows/solidcraft-trivariate-plugin.yml>`_.

On the reference bracket at four display subdivisions per chart direction,
the released FieldCraft reader returned all 1,250 chart samples in source
order. Comparing its ASCII VTU with the native VTU gave maximum absolute
differences of :math:`1.4\times10^{-17}` mm for displacement and
:math:`2.9\times10^{-14}` MPa over all six stress components. The
round-trip check explicitly permutes the six symmetric stress entries into
VTK order (xx, yy, zz, xy, yz, xz). A scalar von Mises comparison alone
cannot catch every axis swap.
