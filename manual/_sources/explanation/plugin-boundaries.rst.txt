Which model does the solver see?
================================

A face extrusion on the bracket illustrates the importance of identifying
the operand. The Rust kernel can construct the spline extrusion tool.
SolidCraft can combine its STEP boundary with an edited bracket and display
the resulting feature. The tool still describes only the added or removed
volume. An elasticity solve on that tool would impose its supports and loads
on a different body from the one in the viewport.

The geometry path
-----------------

A native trivariate model contains rational maps, their interfaces and
boundary identifiers. The STEP boundary supplies SolidCraft's renderer with
faces it can tessellate. It does not carry the interior map network. A native
body edit can therefore change the CAD boundary without producing a new
trivariate analysis domain. Preserve the recipe, model and their lineage
alongside the SolidCraft design.

The analysis path
-----------------

A mechanical job ties its material, support and load to a particular volume
model and named boundary faces. The Rust solver assembles a displacement field
on those volume maps. A separate visualization samples this field for display.
Increasing the density of display samples smooths the picture. Changing the
analysis space or its quadrature changes the numerical experiment.

This distinction matters at a constrained edge. Two display elements can
share the same displacement sample while the one-sided stress differs.
Interpret displayed nodal values using the field definitions in
:doc:`../reference/fieldcraft-plugin`. To change the physics, edit and rerun
the mechanical job and retain the result and job together.
