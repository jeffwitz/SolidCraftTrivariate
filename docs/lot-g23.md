# G23 counterbore constructor protocol

## Scope and revision

Starting revision: `d7d764affcb9db12c87e37e7e0b1ff026a336818` (G21 certified
baseline). Production scope is limited to `src/counterbores.rs`; constructor
tests are in `tests/counterbores.rs`. The candidate family is
`polygon_counterbores_v1`. JSON/CLI adaptation and the independent continuous
audit/OCCT oracle are separate responsibilities.

The input has one polygonal outer profile, integer height and terminal shear,
and one to four disjoint counterbores. Each bore has an ID, integer center,
pilot radius `r`, guard radius `G`, counterbore radius `R`, integer depth `d`,
and an opening (`start` or `end`). Require `1 <= r < R`, `G >= 2R`,
`0 < d < H`, all guards strictly contained and mutually disjoint, and the
existing G21 profile/shear bounds. The pilot hole passes through the entire
height; the larger counterbore removes material over depth `d` from the named
end. No cone, threads, nonintegral dimensions, or arbitrary bore directions
are included.

## Exact layered construction

Use G21's certified guard/profile complex and exact shear map. Create a global
z event at the counterbore shoulder (`d` for `start`, `H-d` for `end`) for each
bore. Split every G21 shell block at all events with homogeneous rational
restriction so neighboring traces remain whole-face compatible. The guard to
counterbore-radius ring (G to R) remains in every layer. On the layer where
only the pilot radius is removed, add eight concentric R-to-r material sectors;
there is no material core inside r. On the larger-counterbore layer, omit that
annular material. This leaves the circular radius-R wall exposed only through
the counterbore depth and the radius-r wall exposed only through the remainder
of the pilot bore. The R-to-r sectors' end face at the event plane is the flat
shoulder.

Stable IDs:

- outer G21 cells retain their root semantic ID and use `_z{z0}_{z1}` for the
  layer interval;
- counterbore ring:
  `bore_{id}_counter_sector_{east|north|west|south}_half_{0|1}_z{z0}_{z1}`;
- pilot-only annulus:
  `bore_{id}_pilot_sector_{east|north|west|south}_half_{0|1}_z{z0}_{z1}`.

The two circular families have degrees `[2,1,1]`, controls ordered with `w`
fastest, `u` increasing around the arc, and `v` moving from G to R or R to r.
Reuse G21's exact rational arc parameterization and odd-quadrant paired
reversal/reflection needed for trace identity. Restrict in `w` with exact
integer homogeneous numerators and one positive common factor per block.

Design-face semantics are `bore/{id}/counterbore_wall` on the R side only in
the wider removed layer, `bore/{id}/pilot_wall` on the r side only in the
pilot-only layer, and `bore/{id}/shoulder` on the annular end at the shoulder
plane. Externally visible shell caps keep `start` and `end`; all other boundary
faces must receive explicit stable classifications. Internal R interfaces
between the G-to-R guard ring and R-to-r annulus are not selectable.

The ideal volume is

`A_outer * H - sum(pi * (r^2 * H + (R^2-r^2) * d))`.

## Exit criteria and refusal coverage

The independent audit must reconstruct its expected layers from the recipe,
certify the planar maps continuously, compare circle/shoulder traces exactly,
and enumerate the full incidence and selection inventory. Verify positive
weights and Jacobians, exact G21 reuse after layer restriction, stable ordering,
input permutation, translation/shear, start/end openings, and 1..4 disjoint
holes. Reject overlapping guards, guard/profile contact, invalid radii/depth,
duplicate IDs, wrong traces/tags, and overflow. STEP/OCCT qualification compares
volume, center of mass, bounds and point classifications; a partially tested
case is not promoted to a general counterbore guarantee.
