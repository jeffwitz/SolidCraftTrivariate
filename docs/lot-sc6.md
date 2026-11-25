# SC6: import circular sketch bores from SolidCraft

## Starting revision and scope

- Starting repository revision: `e6b9979d59990426cafe58d8ece628f16f424a35`.
- Native target: G20 family `polygon_circular_bores_v1` from the frozen G20
  protocol in `docs/campagne-20261125.md` (constructor revision `f005c19`).
- Host evidence: SolidCraft CLI `solidcraft/66`, version `0.48.1`; `solidcraft-cli
  script FILE.sc --out HISTORY.json` extracts editable history without
  rebuilding the document.
- Scope: one closed straight-line outer profile with one to four exact circular
  sketch curves, extruded once as a positive, one-sided new body. This adapter
  imports the saved geometric snapshot; it does not preserve editability.

## Accepted host history

The document has exactly two unsuppressed `timeline.add` features in order:

1. A sketch on the origin XY plane. Its active geometry contains one simple
   closed contour of non-construction straight line curves and one to four
   non-construction `circle` curves. The contour graph is formed by the stored
   line endpoint indices. Each circle stores a point index `c` for its centre
   and an exact integral positive radius `r` in millimetres. The referenced
   centre point must exist and all source coordinates must be finite exact
   integral millimetres. Coincident constraints are accepted only when their
   endpoint coordinates agree exactly; other constraints are excluded.
2. An extrusion referencing that sketch. It selects exactly the outer contour
   through `profiles:{"type":"curves","loops":[[...outer line IDs...]]}`;
   the circles are unselected inner loops interpreted by the host as voids.
   The extrusion is `new_body`, positive, one-sided `total_distance`, with an
   exact positive integral-millimetre distance and no offset, taper, thin wall,
   or extra option.

The adapter sets G20's numerical guard parameter to `guard_mm = 2 * radius_mm`.
This is the smallest admissible guard in G20 and does not change the source
solid. The independent G20 constructor/audit must still verify that the
resulting guard diamonds lie strictly inside the polygon, are disjoint, respect
the 32-vertex budget, and fit the certified domain. A real host solid can
therefore be excluded when this required native partition cannot be certified.

## Refusals and publication

Refuse any other feature count/order, host format/version, plane, sketch curve
kind, construction curve, open/branched/self-intersecting contour, duplicate
identity, unresolved expression, fractional or out-of-domain coordinate,
non-integral radius, invalid circle centre reference, noncoincident constraint,
ambiguous/inconsistent extrusion profile selection, negative/symmetric/tapered
extrusion, offset, Boolean feature, or guard placement rejected by G20. The
adapter builds and audits a candidate in memory and returns only the audited
G20 model. CLI publication remains atomic: a refusal preserves an existing
output. Direct-command SC1 import is unchanged. Timeline dispatch must route a
two-feature sketch/extrude history containing a circle curve to this adapter;
polygon-only histories continue through SC5.

## Qualification thresholds and independent oracle

Freeze actual editable `.sc` documents and their extracted histories for one,
two, and four holes, plus an oblique/concave outer contour where admissible.
Generate the host STEP from the saved `.sc` and the native STEP from the
imported/audited model. Reopen both with FreeCAD and require one valid solid,
volume consistent with `(outer area - sum(pi*r^2))*height`, host/native volume
and area differences below `1e-8 mm^3` and `1e-8 mm^2`, centre-of-mass
difference below `1e-9 mm`, symmetric-difference volume below `1e-8 mm^3`, and
matching interior, exterior, and bore point classifications. The host solid
and native builder are not used as each other's geometric oracle.

Exercise mutations to contour closure, circle radius/centre, circle type,
construction state, constraint type, profile selection, extrusion direction,
offset/taper, feature reference, and feature count. Verify refusals preserve a
pre-existing output. Record valid host cases excluded by the native G20 guard
domain separately from malformed files; a refusal on a frozen required valid
case counts against coverage. Keep extracted histories, STEP files, oracle
inputs, and reports under `build/` or `tests/references/` as appropriate.

## Explicit exclusions

SC6 does not import sketch arcs, rounded profiles, multiple polygonal islands,
circle arrays made by features, side-drilled `solid.hole` operations, Boolean
feature timelines, non-XY/attached planes, parameter expressions, or
non-integral sketch dimensions. It does not infer material topology from
sampled points. Expanding this subset requires new real host documents and
separate differential evidence before broadening the adapter.
