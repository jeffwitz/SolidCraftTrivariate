# Independent comparison of sketch partitions

## Euler lower bounds for the blind pocket

The frozen case `plate_blind_pocket` is a 60 × 40 × 20 mm plate with an
interior 20 × 10 mm pocket from z=14 to z=20. Its exact target volume is
`48 000 − 20 × 10 × 6 = 46 800 mm³`.

For a conforming quadrilateral mesh of a planar region, let `B` and `I` be
its boundary and interior vertex counts, `E` its edge count, and `Q` its
quadrilateral count. Since each interior edge is incident to two quads and
each boundary edge to one, `4Q = 2E − B`. Combining this with Euler's formula
`V − E + Q = χ`, where `V = B + I`, gives

```text
Q = B/2 + I − χ.
```

The lower slab is a disk (`χ=1`). Its rectangular outer boundary requires at
least four vertices. Matching the four corners of the upper pocket boundary
requires at least four interior vertices in this lower section. Therefore
`Q ≥ 4/2 + 4 − 1 = 5`. These four locations are the reentrant pocket corners
when viewed on the upper boundary; treating them as vertices of the lower
section is essential to the bound.

The upper slab is an annulus (`χ=0`). Its outer rectangle and pocket contour
require at least eight boundary vertices, four on each contour, so `Q ≥ 8/2
= 4`. For the layer-aligned family tested here, every block is cut at the
feature event z=14 and uses a conforming all-quadrilateral section on each
slab. In that family, the lower bound is nine blocks, split 5 below and 4
above. It is not a global minimum over all possible 3D block complexes:
blocks that span z=14 outside the pocket footprint could reduce the total,
but would require partial/interface raccords that this family has not
qualified. The Euler count does not certify the geometric maps, their
positive Jacobians, or a 3D partition by itself.

## C-section control

The frozen `c_section_two_holes` section has one connected C-shaped material
component and two square through-holes. Its boundary has 16 required corner
vertices: eight on the outer C and four on each hole. Here `χ=1−h=−1`, with
`h=2`, so `Q = B/2 + I + h − 1 ≥ 8 + 2 − 1 = 9`. The exact C-section area is
`20×20 − 14×12 − 2×2×2 = 224 mm²`; at height 7 mm its exact volume is
`1 568 mm³`.

## Oracle and recorded baseline

`tests/references/pocket_partition_comparison.py` invokes the CLI on the
versioned recipe corpus, then independently reads each emitted model JSON.
It reconstructs the trilinear block corners using rational arithmetic,
checks exact planar edge incidences and Euler characteristics, checks section
containment and pairwise interior disjointness, and compares exact summed
section area and volume with shoelace references derived from the recipes.
The fixed cases include the blind pocket, a centered and a translated
rectangular pocket, the C-section with two holes, and an L-shaped bracket
with a hole and blind pocket.

Run from the repository root:

```sh
rtk python3 tests/references/pocket_partition_comparison.py \
  --cli target/release/solidcraft-trivariate \
  --work build/pocket-partition-comparison-final \
  --output build/pocket-partition-comparison-final/report.json
```

At the recorded baseline (`d25b6b6`, release CLI SHA-256
`5b13b9287a30cd69d2bfba16662fded47943cd1048d4bf7fc1fad36dbebce379`), all
five cases matched exact recipe volumes and section coverage. The blind
pocket used 17 blocks (9 lower, 8 upper), and the C-section used 35 blocks;
both therefore failed to attain the proposed nine-block lower bounds. The
remaining counts were 26 for the centered rectangle, 26 for the translated
rectangle, and 77 for the L-shaped bracket. This baseline does not qualify a
later candidate.

On the guided candidate release CLI, the same command passed all six cases. The
CLI hash used for this run is recorded in its generated JSON report; the
integration commit will rerun the command against its own release binary.
The blind pocket has 5 lower and 4 upper
cells, with exact section areas 2 400 and 2 200 mm² and volume 46 800 mm³.
The lower incidence counts are `(V,E,Q)=(8,12,5)` and `χ=1`; the annular
section counts are `(8,12,4)` and `χ=0`. It attains the conditional nine
block bound. The C-shaped section also attains its Euler bound: 9 cells,
`(V,E,Q)=(16,26,9)`, `χ=−1`, area 224 mm², and volume 1 568 mm³. The centered
and translated rectangles each use 14 blocks and have volume 9 504 mm³; the
L-shaped polygonal prism with a through-hole and blind pocket uses 26 blocks
and has volume 5 084 mm³. Each case has exact section area conservation,
zero pairwise interior-overlap witnesses, and zero section area gap.

The same 27-point-per-block Jacobian condition-number sampling used by the
partition-cost reference gives these historical-to-guided maxima: blind
pocket `3.33333 → 10.57210`, C-section `8.88748 → 121.90145`, and L-shaped
polygonal bracket `5.00000 → 62.98412`. These are sampled `κ₂(J)` diagnostics;
the increase indicates greater map distortion in the guided output and is
not a continuous bound. The translated C-hole variant moves the lower 2 × 2
hole +1 mm in x, away from the C opening. Five of the nine block IDs remain
common across the edit. The local edit check derives its expected unaffected
set geometrically: a block
must have no control point at any of the original lower-hole corners. Exactly
5 blocks meet that condition, and their IDs, degrees, and homogeneous controls
must remain identical. The interface network has 10 edges before and after;
the report counts exact interfaces between those unchanged blocks.

The [integration report](results/pocket-partition-comparison.json) records the CLI
hash, per-layer incidences, exact areas, volumes, sampled condition numbers,
and the translated-hole ID/control/interface comparison.

## Limits

The oracle covers polygonal straight-edged sections and trilinear blocks. It
checks non-overlap with a separating-axis test using exact rational
projections on all eight quadrilateral edge normals, after requiring strict
convexity; it also rejects planar T-junctions. For a translated top profile,
it verifies every upper corner's displacement against the recipe's prescribed
linear shift. Its coverage argument is exact in the section plane: cells lie within the
expected material, their interiors do not overlap, and their exact areas sum
to the recipe area. This is not a continuous 3D injectivity or Jacobian proof.
The Euler minima assume a conforming quadrilateral mesh, layer alignment at
feature events, and preservation of the specified contour corners;
nonconforming T-junctions, extra feature boundaries, cross-event blocks, or a
different topology require a new count. The L-shaped case
is a polygonal prism proxy, not a qualification of a curved mechanical
bracket.
