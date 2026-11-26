# Independent bore-chamfer oracle

## Result

The frozen 10-chart FieldCraft bracket has two separate radius-4 mm through-bores, each spanning 8 mm axially. The generated fixed-partition edit adds a quarter-millimetre chamfer at both axial mouths of each bore. The independent oracle passes all four mouth sections.

There are four physical circular mouth edges in the target (two bores times two axial ends). The oracle checks all four, represented by 16 quarter-arc chart sections. The bore-wall test is an exact polynomial identity in the rational Bézier angle and axial parameters, after clearing the homogeneous denominator. It proves the radius changes linearly from 4.25 mm at each mouth to 4 mm at the end of the 0.25 mm chamfer support, and remains 4 mm in the middle span. The four source sector arcs per bore are checked to be four connected, consistently oriented quarter-circles, so the calculation covers each complete circle exactly once. The untouched radial surface of every bore chart and both elbow maps are independently compared with their exact source restrictions.

The cross-section integral for one mouth is

```text
∫₀ʷ π[(R + w − s)² − R²] ds = π(Rw² + w³/3).
```

For `R=4 mm`, `w=1/4 mm`, each mouth adds `49π/192 mm³`. Thus each bore's two mouths add `49π/96 mm³`, and both bores add exactly

```text
49π/48 mm³ = 3.2070425005395804 mm³.
```

The oracle also independently certifies the tensor Bernstein determinant of all 30 emitted Bezier elements, with no negative coefficient and no all-zero boundary stratum. The serialized construction reports 10 logical charts and 53 interfaces; the oracle checks the exact geometric restrictions, while the interface count is the constructor's declared inventory.

## Reproduction and provenance

From the repository root:

```sh
cargo run --offline --example bracket_edge_edit
python3 tests/references/bracket_edge_edit_oracle.py
```

The first command writes `build/bracket-edge-edit.json`. The Python oracle independently reads the source `docs/site/data/fieldcraft_x64_y64/fieldcraft_x64_y64.model.json` and that output, then writes `build/bracket-edge-edit-independent.json`. The source SHA-256 is `82fb0d3f7a0977a47c705886d13632900dc62429be572231fc4f152a38688e09`; the output SHA-256 audited here is `406dfdcfe8913bba52863e8a1b564dd43321505a517a1165f8b7a77493e8b7b5`.

The derivative certificate reuses the independent exact integer Bernstein implementation in `tests/references/s3_stock_rim.py`. Radius identities, seam restrictions, and the analytic volume calculation are implemented directly in `tests/references/bracket_edge_edit_oracle.py` with Python `Fraction` arithmetic. The fixed-radius result applies to this exact source model and these two 0.25 mm end chamfers; it does not qualify a round fillet or a different topology.
