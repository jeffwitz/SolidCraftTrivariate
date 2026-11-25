# G21 interval-bore constructor protocol

## Frozen scope

Base revision: `05f3b7b`. The family is `polygon_interval_bores_v1`; the frozen
input corpus is `tests/references/interval_bores_corpus.json` (30 recipes), and
the independent partition/geometry oracle is `tests/references/interval_bores.py`.
The builder is `src/interval_bores.rs`; its audit is kept independent in
`src/interval_bores_audit.rs`.

Accepted recipes use an exact polygonal outer profile, one to four disjoint
guarded circular bores, an integer height, an integer terminal shear, and an
integer removal interval `[lo, hi]` for each bore. The layer partition includes
all interval endpoints. The guard annulus is present at every height; the disk
is filled by eight circle-to-inner-diamond sectors and four Q1 core cells only
where the material remains. A through bore, either blind cap, a sealed cavity,
coincident layer events, concave profiles, shear, ordering, and the frozen
invalid recipes remain in scope. No claim extends beyond these admitted
profiles and the continuous certificates available to the audit.

## Orientation and trace regression

The first implementation exposed a ring chart trace defect on the smallest
frozen `through` recipe. A separate sealed-cavity witness exercises the central
Q1 cells, which are absent from a through-bore recipe. Both regression recipes
are preserved under `tests/regressions/`; the independent audit must certify a
strictly positive determinant lower bound (`> 0.0 mm^3`), exact incidence, and
the expected volume coefficient. The old behavior fails the positive-Jacobian
assertion on the cavity's core blocks and fails incidence on the ring-only
through case.

The central Q1 cells use proper counter-clockwise rotations in the sequence
east, north, west, south. The ring sectors preserve G20's exact conic trace
convention: its quarter-arc endpoint weights are asymmetric (2, 2, 4), so
north and south pair a reversed angular parameter with the matching axis
reflection. This keeps angle increasing and preserves projective radial-seam
traces. For the central Q1 cell,
the east control order is `[east, south-east, north-east, center]`, which is
`[P00, P01, P10, P11]` in the `(u*2+v)*2+w` layout; the other cells are its
proper rotations. The inner-diamond edge is represented as a Q2 degree
elevation of a straight Q1 edge at common weight, with controls
`[[2r,0,4], [r,r,4], [0,2r,4]]`. This makes its restricted trace agree exactly
with the Q1 core edge; a geometrically straight but nonlinearly parameterized
edge is not sufficient for incidence.

Before adoption, run the constructor tests, independent audit tests, the frozen
30-case corpus, and the STEP/OCCT integration oracle. Keep any constructor vs
oracle disagreement as a counterexample and report rejected mandatory cases
as coverage failures.
