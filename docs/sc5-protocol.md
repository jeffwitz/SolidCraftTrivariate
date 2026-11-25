# SC5 protocol: closed polygon sketch extrusion import

## Scope and host evidence

SC5 translates a canonical SolidCraft `solidcraft/66` history into either the
native `polygon_prism_v1` or `polygonal_perforated_prism_v1` family. The host
history is extracted with `solidcraft-cli script FILE.sc --out HISTORY.json`;
this command reads the saved editable document and does not rebuild its solid.
The host reports version `0.48.1` in the `about` record. The available command
registry documents `sketch.create`, `sketch.line`, `sketch.finish`, and
`solid.extrude`; the importer accepts only a smaller, frozen subset below.

## Accepted history

The history has exactly two active `timeline.add` commands, in dependency
order:

1. An unsuppressed `sketch` feature on the origin `XY` plane. Its sketch
   contains only non-construction `line` curves. Curves refer to point indices
   in the stored point array; their directed endpoints form one or more closed
   cycles. Point coordinates are finite, exact integral millimetres in the
   native G17/G18 domain. A simple G17 contour has 3..24 vertices. G18 accepts
   one to four strict hole cycles with 3..32 total vertices. Loops may be
   concave. The stored curve graph, not array order or
   geometric proximity, defines adjacency.
2. An unsuppressed `extrude` feature referencing that sketch, with all profiles,
   `new_body`, positive one-sided `total_distance`, literal positive integral
   millimetre distance, no offset, no taper, and no thin-wall mode. The output
   is an extrusion along +Z with no top shift. Width, height, orientation, and
   profile containment come from the recorded sketch and are revalidated by
   the native constructor and independent audit.

The outer loop is the unique cycle strictly containing all other cycles. Hole
cycles must be pairwise disjoint and strictly inside it. The importer assigns
stable profile labels `outer` and `hole0`...`hole3`; holes are ordered by their
smallest safe source vertex ID. Vertex IDs derive from the source sketch point IDs after
the adapter's documented safe-name conversion. Loop winding and cycle start
are normalized by the native constructors; no coordinate-based welding is
performed.

## Refusals

Refuse malformed or unsupported metadata, extra/missing/suppressed features,
unresolved or non-literal parameters, non-XY/offset/attached planes, arcs or
construction geometry, open or branched endpoint graphs, unsupported sketch
constraints, non-finite/fractional/out-of-range coordinates, duplicate source
IDs, self-intersection, loop contact, nested or exterior holes, collinear
adjacent vertices, ambiguous profile selection, and any extrusion option
outside the accepted subset. Reject before publishing a candidate. Call the
versioned native builder for the selected family and its independent audit;
only return the audited model. The caller's existing output remains untouched
on any refusal.

## Qualification plan

Before adoption, freeze real `.sc` sources and extracted histories for a
triangle, a concave simple polygon, and a closed profile containing one or more
holes. For each source, compare host-produced STEP to the native model's STEP
by reopened geometry, volume, and point classification; do not treat equal
volume alone as proof of equal shape. Preserve the source file, extracted
history, command-line invocations, host version, STEP comparison tolerance,
and per-case admission result. Exercise direct and timeline imports against
the same recipe, then mutate feature references, curve closure, point values,
extrusion flags, and loop relationships. Every refusal must preserve a
pre-existing output. Report unsupported but valid SolidCraft histories as
explicit exclusions rather than successful coverage.

## Version boundary

The initial adapter is deliberately specific to the `solidcraft/66` stored
history shape and host version captured above. It does not evaluate parameters,
constraints, sketch solvers, arbitrary planes, draft, symmetric/reverse
extrusion, face attachments, or boolean operations. Extending that boundary
requires a new protocol revision and differential host cases before code
changes.
