# Polygonal sketch planner evaluation

This directory preserves the frozen inputs and the paired measurements at `4e8a394`. The [current kernel results](../../convex-host-rings.md) cover the same 45 inputs and 12 additional sketches.

This directory contains a paired evaluation of the baseline planner and the frozen candidate on 45 deterministic public synthetic recipes. It is separate from any externally owned private holdout.

## Frozen inputs and executables

- Baseline source: `e5bd5c6`; executable SHA-256: `ea0bbf8609bb575cc5ac5f8dec6c8c1659a85050a95f7f0db86869d15160786c`.
- Candidate source: `4e8a3944cd5aa85f8ba0984b109b6649b2d54b6f`; executable SHA-256: `684c7c6a89423f798648800e28c53d1065c23b8aaa7e5407409714cea4a36f70`.
- Frozen recipe corpus SHA-256: `f10964bb0417a702a394f66aa96e148b6ead67685a3a9f9f541cea8d6ab98e7b`.

The recipes use integer millimetres inside the kernel's coordinate range. The corpus contains five fixed convex through-hole cases, five convex blind-pocket cases, five concave stock or cut cases, and 30 seeded convex cases (15 through-hole and 15 blind-pocket). The seeded profiles use five to twelve sides and vary cut count, scale, placement, slant, and orientation. The generator seed is 20261126.

## Results

See [paired.md](paired.md) for the paired summary and [paired.json](paired.json) for per-case comparisons. `baseline.json` and `candidate.json` hold the per-executable observations, including build/audit result, raw construction method, exact volume, block count, conditional and input-vertex Euler bounds, and sampled condition numbers. `baseline-model-backward-replay.json` reports the separate audit of all 45 archived baseline models through the candidate executable; all 45 returned a valid audit. This replay does not change fresh-recipe geometry scores.

The candidate and baseline both passed the geometric gate on all 45 inputs, with no refusals or timeouts. The candidate used 5,201 blocks versus 6,184 for baseline, a reduction of 983 blocks (15.9%). Fifteen cases used fewer blocks and thirty retained their previous count. For the seeded profiles, all 30 baseline models used `triangulated_feature_footprint`. The candidate selected `radial_convex_footprint` on 9 of 10 one-cut cases; the remaining one-cut case (`seeded_convex_15`) and all 20 multiple-cut cases retained triangulation. The nine radial conversions reduced their summed input-vertex Euler gap from 771 to 39 blocks. All 20 multiple-cut cases retained their baseline Euler gap.

There is a map-quality tradeoff. The maximum sampled 2-norm Jacobian condition across each solid had median 285.97 at baseline and 218.05 for the candidate; nearest-rank p95 and worst remained 4,534.39 and 49,185.01. `blind_triangle_shallow` used 56 instead of 105 blocks, while its sampled maximum rose from 53.94 to 908.70. Other case-level increases are listed in `paired.json`.

The condition number is sampled on the same closed 5×5×5 parameter grid in every block. It is a diagnostic and does not replace the continuous closed-domain determinant certificate. The Euler relation applies to connected conforming quadrilateral sections. Input-only bounds use mandatory contour vertices; conditional bounds also count the inserted vertices present in the audited partition.

The baseline runner took 235.60 seconds and candidate runner 209.23 seconds. The candidate run overlapped a separate Cargo test process at startup, so these elapsed times are descriptive and do not establish a speedup.

## Reproduction

Run these commands from the repository root. The executable paths should point to release binaries built from the stated revisions. The runner locates `tests/references/pocket_partition_comparison.py` relative to this repository, imports its exact section predicates, and records each executable's SHA-256 before and after the run.

```sh
python3 docs/results/sketch-planner-polygon-evaluation/generate.py
python3 docs/results/sketch-planner-polygon-evaluation/run.py \
  --cli /path/to/e5bd5c6/target/release/solidcraft-trivariate \
  --tag baseline-e5bd5c6 \
  --work /tmp/sketch-planner-baseline \
  --output /tmp/sketch-planner-baseline.json \
  --timeout 60
python3 docs/results/sketch-planner-polygon-evaluation/run.py \
  --cli /path/to/4e8a3944/target/release/solidcraft-trivariate \
  --tag candidate-4e8a3944 \
  --work /tmp/sketch-planner-candidate \
  --output /tmp/sketch-planner-candidate.json \
  --timeout 60
python3 docs/results/sketch-planner-polygon-evaluation/replay.py \
  --cli /path/to/4e8a3944/target/release/solidcraft-trivariate \
  --baseline-models /tmp/sketch-planner-baseline \
  --baseline-result /tmp/sketch-planner-baseline.json \
  --output /tmp/sketch-planner-backward-replay.json
python3 docs/results/sketch-planner-polygon-evaluation/report.py \
  --baseline /tmp/sketch-planner-baseline.json \
  --candidate /tmp/sketch-planner-candidate.json \
  --corpus docs/results/sketch-planner-polygon-evaluation/cases.json \
  --json /tmp/sketch-planner-paired.json \
  --markdown /tmp/sketch-planner-paired.md
```

Each build and audit operation has a 60-second timeout. Rejections and timeouts remain in the per-family denominator. The generator emits the corpus hash; it must match the frozen hash above before comparing results.
