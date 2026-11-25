# SC9 protocol: native SolidCraft counterbores into G23

## Revision, route, and scope

Starting revision: `4d5ae4c156eeacc0706bbbb188f2001a87d056e7`, with G23's
`polygon_counterbores_v1` constructor. The adapter is isolated in
`src/solidcraft_counterbore.rs`. The bridge route accepts exactly a saved
SolidCraft 0.48.1 history with one initial `box` feature and one `hole`
feature whose hole type is `counterbore`; the equivalent direct command script
has exactly `solid.box` and `solid.hole`. SC8 simple blind-hole histories and
all previously qualified SC1–SC8 histories remain routed as before.

The SolidCraft 0.48.1 registry and real host confirm one `solid.hole` operation
with pilot `diameter`, `type: counterbore`, `cb_diameter` and `cb_depth`; pilot
depth defaults to through-all. The native host feature was used directly; no
coaxial two-hole approximation is involved. The saved `.sc` documents reopen
as a single body. The adapter is `src/solidcraft_counterbore.rs` and offers
`import_counterbore_timeline` and `import_counterbore_direct`; the bridge routes
`counterbore` separately from SC8's `simple` hole subset.

## Frozen witnesses

Use the SC8 40×40×12 mm rectangular box with an 8 mm through pilot and 12 mm
counterbore diameter at depth 4 mm. Save and freeze host histories and `.sc`
files for a centered top opening, a decentered top opening, and a bottom
opening only if the host exposes the reverse direction. Host faces `F1:+z`
and `F1:-z` determine the opening side; replay must confirm their direction and
shoulder positions. The exact G23 fields are pilot `r=4`, `R=6`, guard `G=12`,
depth 4, and opening `end` from the top or `start` from the bottom.

The analytic volume is
`19200 - pi*(4^2*12 + (6^2-4^2)*4) = 19200 - 272*pi mm^3`.
For a top entry the shoulder lies at `z=8` and the wider wall spans `[8,12]`;
for a bottom entry it lies at `z=4` and the wider wall spans `[0,4]`. The
shoulder area is `pi*(6^2-4^2)=20*pi mm^2`, with no bore end cap because the
pilot passes through.

## Qualification and refusals

Replay and reopen each real `.sc` file with SolidCraft, export STEP, import the
same history through G23, and compare host and native STEP in OCCT. Three pairs
passed: centered top, decentered top, and decentered bottom. Each is one valid
solid with the expected bounds and analytic volume `19200 - 272*pi mm^3`.
Symmetric Boolean differences were zero. Volume matched the analytic value to
about `4.5e-9 mm^3`; total area matched the independent decomposition
`5120 + 80*pi mm^2` within `2.3e-9 mm^2`. Centroid errors were below
`1.7e-12 mm`; twelve inside/outside points per pair passed. The independent
surface terms are pilot wall `64*pi`, counterbore wall `48*pi`, and shoulder
`20*pi mm^2`; the shoulder lies at z=8 from the top and z=4 from the bottom.
See [SC9 measurements](results/solidcraft-counterbore-20261125.json).

The adapter uses the host face/direction, center, pilot and counterbore
diameters, and counterbore depth. A host-accepted 13 mm counterbore diameter
was refused because G23 requires integer radii; CLI refusal preserved an
existing output sentinel byte-for-byte. Unsupported fields, finite pilot
depth, conical points, odd diameters, nonintegral placement, and inconsistent
face/direction are rejected before model publication.

Reject a host-accepted but native-out-of-contract counterbore (for example an
odd integer counterbore diameter, which gives a half-millimetre radius outside
G23's integer-radius lattice), a conical/drilled tip, a non-through pilot, or
inconsistent face/direction fields. Verify that a refused import leaves an
existing output model byte-for-byte unchanged. A host-rendered feature is
qualified only for these recorded histories and the exact G23 continuous audit
domain.
