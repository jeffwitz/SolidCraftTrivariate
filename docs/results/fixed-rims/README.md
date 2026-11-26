# Edge treatments on the original trivariate partitions

The S3 and FieldCraft experiments retain the volume identities of their
source models. An axial knot introduces additional Bézier spans inside a
chart. Surface counts in a shell export therefore differ from the number
of logical sides or volumes.

## S3 pocket stock

Both outside horizontal perimeters of the stock are treated at 0.25 mm.
The eight straight edges meet in mitres at the corners. The transverse cut,
NURBS through-hole and pocket remain present. Both exchange files were
generated from the modified trivariate maps and independently re-imported
as one valid closed solid.

- [Chamfer STEP](s3-stock-chamfer.step)
- [Fillet STEP](s3-stock-fillet.step)
- [Exact section and Jacobian checks](../s3-stock-rim.md)
- [STEP shell measurements and contact probes](s3-step-check.json)

A separate [offset-rail experiment](../s3-offset-rail.md) concerns the top
lip of the curved hole. Its error applies to a contact curve. The stock-rim
STEP files do not apply any treatment to that hole.

## Regularity at a smooth boundary

At a map edge where both incident faces are exterior and receive a common
tangent plane, both cross-edge derivatives and the edge tangent lie in that
plane. The Jacobian of that chart vanishes on the edge. The [actual-model
audit](../fixed-edge-rank.md) identifies the relevant incidence pattern in
both source models, and the [original annulus study](../../annulus-four-block-topology.md)
provides the prepared, closed-domain regular alternative for that geometry.

## FieldCraft compact bracket

The four circular mouths of the two transverse bores are treated. Each
state has ten original chart identities and eleven original neighbors. Two
axial knots in each wing are extended through the two elbow charts. The
three states share the same sequence of axial intervals [0,1,31,32]/32.
They have 30 Bézier spans and 53 matching span interfaces. The bend radius
and the remaining exterior edges are unchanged.

On the published x64/y64 bracket the bores have a 4 mm radius. Exact
identities on the sixteen angular patches recover the 0.25 mm chamfer and
the toroidal fillet sections, along with the contact circles at each mouth.
The removed volumes are 49π/48 = 3.20704250 mm³ and
101π/48 − 17π²/32 = 1.36720720 mm³, respectively. Every chamfer span has a
positive Jacobian on its closed cube. The fillet has fourteen such spans
and sixteen spans with a zero on just the cap contact rail.

- [Chamfer derivation and oracle](../../../tests/references/bracket_edge_edit_oracle.md)
- [Fillet derivation and oracle](../../../tests/references/bracket_fillet_rank_oracle.md)

Generate a logical model with its original IDs and normalized parent
parameter intervals, or export the same boundary in STEP:

```sh
cargo build --release --locked --offline
model=docs/site/data/fieldcraft_x64_y64/fieldcraft_x64_y64.model.json
target/release/solidcraft-trivariate edit-fieldcraft-bore-rims "$model" chamfer build/bracket-chamfer.json
target/release/solidcraft-trivariate edit-fieldcraft-bore-rims "$model" fillet build/bracket-fillet.step
```

The [five-variant replay](fieldcraft-series.json) uses all published
FieldCraft bore-radius pairs, including the 51/16 mm radius. For each of the
15 emitted states the independent oracle checks the normalized knot
intervals, integer Bernstein signs on every span and the integrated volume.
All five sharp-to-treated volume differences agree with the analytic
bore-mouth integrals within 1e−7 mm³. Run
`python3 tests/references/fieldcraft_rim_series.py --cli target/release/solidcraft-trivariate --work build/rim-series`.

On x64/y64, FreeCAD/OCCT imported one closed, valid solid in each state.
All 288 probes per state, just inside and outside the expected material
boundary at both mouths of both bores, agreed with the analytical reference.
[Chamfer STEP](fieldcraft-chamfer.step) · [Fillet STEP](fieldcraft-fillet.step) ·
[CAD measurements](fieldcraft-step-check.json)

The [S3 contact-rail controls](s3-contact-rails.json) cover only the named
internal lip span. They are separately downloadable for checking the
approximation bounds against the original sketch.
