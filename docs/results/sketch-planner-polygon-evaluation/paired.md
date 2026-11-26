# Paired polygonal sketch partition evaluation

These paired scores are historical observations at `4e8a394`. The [current kernel results](../../convex-host-rings.md) cover the same 45 inputs and 12 additional sketches.

Public synthetic corpus SHA-256: `f10964bb0417a702a394f66aa96e148b6ead67685a3a9f9f541cea8d6ab98e7b`.
Baseline executable SHA-256: `ea0bbf8609bb575cc5ac5f8dec6c8c1659a85050a95f7f0db86869d15160786c` (`baseline-e5bd5c6`).
Candidate executable SHA-256: `684c7c6a89423f798648800e28c53d1065c23b8aaa7e5407409714cea4a36f70` (`candidate-4e8a3944`).
The baseline source is `e5bd5c6`; the frozen candidate source is `4e8a3944cd5aa85f8ba0984b109b6649b2d54b6f`. The 45-case public synthetic bundle SHA-256 was recorded before candidate scoring.
All measurements are paired on the same 45 recipes: 5 fixed convex through-hole, 5 fixed convex blind-pocket, 5 fixed concave stock/cut, and 30 seeded convex cases split evenly between through-holes and blind pockets. Integer millimetres are within the input coordinate contract.
Each admitted model is independently checked with exact rational section containment, edge incidence, non-overlap, area conservation, and exact volume, plus the CLI closed-domain determinant certificate. Each build and audit operation had a 60 second timeout.

## Paired outcomes by family

| family | n | baseline_valid | candidate_valid | new | regressions | common_success | common_failure | median_block_delta | baseline_euler_gap | candidate_euler_gap |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| concave_stock_or_cut | 5 | 5 | 5 | 0 | 0 | 5 | 0 | 0 | 132 | 132 |
| convex_blind_pockets | 5 | 5 | 5 | 0 | 0 | 5 | 0 | -40 | 328 | 181 |
| convex_through_holes | 5 | 5 | 5 | 0 | 0 | 5 | 0 | -22 | 116 | 12 |
| seeded_convex_blind_pockets | 15 | 15 | 15 | 0 | 0 | 15 | 0 | 0 | 3575 | 3114 |
| seeded_convex_through_holes | 15 | 15 | 15 | 0 | 0 | 15 | 0 | 0 | 1441 | 1170 |

A negative block delta means the candidate used fewer blocks. Rejections and timeouts remain in the denominator. All 45 paired inputs passed geometry on both executables. Total blocks fell from 6,184 to 5,201 (−15.9%); there were no new admissions or regressions. Euler gap is observed blocks minus the bound using mandatory input contour vertices; the conditional bound including inserted host vertices and exact per-case gaps are in the JSON.

## Construction methods

Among the 30 seeded convex cases (stock polygons have 5–12 sides), the baseline used triangulation on all 30. The candidate selected radial construction on 9/10 one-cut cases; `seeded_convex_15` is the remaining one-cut case kept on triangulation. All 20 multiple-cut cases retained triangulation. The input-vertex Euler gap totals for the one-cut stratum are 974→242; the multiple-cut stratum is 4042→4042.

### concave_stock_or_cut

Baseline: `{'orthogonal_feature_footprint': 4, 'guided_quad_footprint': 1}`

Candidate: `{'orthogonal_feature_footprint': 4, 'guided_quad_footprint': 1}`

### convex_blind_pockets

Baseline: `{'triangulated_feature_footprint': 3, 'rectangular_hosts_y': 1, 'rectangular_hosts_x': 1}`

Candidate: `{'radial_convex_footprint': 3, 'rectangular_hosts_y': 1, 'rectangular_hosts_x': 1}`

### convex_through_holes

Baseline: `{'triangulated_feature_footprint': 3, 'rectangular_hosts_x': 2}`

Candidate: `{'radial_convex_footprint': 2, 'rectangular_hosts_x': 2, 'convex_host_separators': 1}`

### seeded_convex_blind_pockets

Baseline: `{'triangulated_feature_footprint': 15}`

Candidate: `{'triangulated_feature_footprint': 11, 'radial_convex_footprint': 4}`

### seeded_convex_through_holes

Baseline: `{'triangulated_feature_footprint': 15}`

Candidate: `{'radial_convex_footprint': 5, 'triangulated_feature_footprint': 10}`

## Sampled map condition

The condition metric is diagnostic: for every block, the spectral condition number of the Jacobian was evaluated on the closed 5×5×5 parameter grid. The reported statistics take each solid’s maximum over its blocks. This is not a continuous regularity certificate.

| executable | n | median | p95 | worst |
| --- | --- | --- | --- | --- |
| baseline-e5bd5c6 | 45 | 285.97290288907664 | 4534.385859639157 | 49185.00826535531 |
| candidate-4e8a3944 | 45 | 218.0454462587007 | 4534.385859639157 | 49185.00826535531 |

The candidate reduced blocks sharply on some cases while worsening sampled condition. The largest paired rise is `blind_triangle_shallow`: 105→56 blocks and κmax 53.94→908.70 (16.84×). Seeded cases `_03`, `_09`, `_21`, and `_12` rise 1.41–1.84×. These remain visible even though the aggregate median improves.
Per-case block counts, observed timing, construction method, Euler bounds, pass/refusal status, and per-block 5³ condition samples are in the adjacent JSON. Candidate timing overlapped a separate cargo test process at startup; observed wall totals (235.60 s baseline, 209.23 s candidate) are descriptive, not a controlled speed claim.

## Scope and limits

The corpus is deterministic public synthetic input, not the externally owned private holdout. It covers convex polygonal through-holes, convex blind pockets, and concave stock or cut profiles. Nonparallel and stacked 3D bores are assessed separately. The Euler bound applies only to connected conforming quadrilateral sections under the stated mandatory-corner assumptions. Jacobian samples are not used as a regularity certificate.
