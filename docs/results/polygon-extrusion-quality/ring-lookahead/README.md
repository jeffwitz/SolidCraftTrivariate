# Paired evaluation: planar ring lookahead

Candidate source revision: `f406011839fdf6dd92f3d7857c53b2b35c90231e`; binary SHA-256: `e2ac9c06f1148c914f53c13ad05e4f71c6900847c758647bcd7e0470cbd55bfc`.
Baseline source revision: `72e8e19c1bd4b25b41610a11d930f38543e40d40`; binary SHA-256: `8f33eeba801e38e56fbe52ad6ccb96a68a28d2f20187a7d67df869af0c0f9ac2`.
Baseline ledger SHA-256: `bc532a128c4fda221df6f6b3e555bae1a391b997b88a16075911041ee7918a44`; evaluator SHA-256: `18387fe3ae81f7fd01bb0f143c8ce5d3feaccbdd1b9be1fc1d31bf004ff807a7`.

## Cohort and geometric gate

The cohort contains 57 fixed recipes: 45 historical cases and 12 supplements. The candidate recipe adds only the explicit `partition_strategy=ring_planar_lookahead` selector. Both sides use the frozen `run.py` build/audit and exact Fraction section checks.

Geometry passes: baseline 57/57; candidate 57/57.
Block counts: 57 equal, 0 decreased, 0 increased; total 1628 → 1628 (delta 0).
Block-map sets: 29 changed and 28 unchanged (comparison ignores IDs and partition-plan metadata; it hashes sorted degrees and homogeneous controls).
Common-success build+audit cost: median 0.9253 s baseline / 2.1751 s candidate; cohort sums 159.400 s / 398.124 s. These are descriptive single-run timings.

## Independent quality comparison

All ratios are adverse-direction ratios where values above 1 mean worse: candidate/baseline for condition numbers and thickness aspect, baseline/candidate for minimum corner sine. Metrics use the independent planar evaluator; extrema are sampled, not continuous bounds. A relative change must exceed 1e−10 to be counted as worsening.

| Metric | Paired n | Median adverse ratio | p95 | Maximum | Improved | Unchanged | Worsened |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Planar spectral κ | 57 | 0.992588 | 1 | 1 | 29 | 28 | 0 |
| Normal thickness aspect | 57 | 1 | 1.3131 | 1.51557 | 2 | 41 | 14 |
| Full 3D spectral κ | 57 | 1 | 1.31301 | 1.5157 | 9 | 36 | 12 |
| Minimum corner sine | 57 | 1 | 1.08468 | 1.69703 | 16 | 35 | 6 |

### Extreme witnesses

- `seeded_convex_16`: 118.061 → 178.929 (adverse ratio 1.51557); candidate witness block `quad_341d1_180d1_405d1_254d1_6767d14_5409d28_15079d7_7326d7_z0_15`, parameter `[1.0, 1.0, 0.25]`.
- `seeded_convex_16`: 118.062 → 178.947 (adverse ratio 1.5157); candidate witness block `quad_341d1_180d1_405d1_254d1_6767d14_5409d28_15079d7_7326d7_z0_15`, parameter `[1.0, 1.0, 0.25]`.
- `fresh_convex_04`: 0.277161 → 0.163321 (adverse ratio 1.69703); candidate witness block `quad_42d1_25d1_423d8_100d1_56d1_69d1_56d1_84d1_z0_10`, parameter `[0.0, 0.0, 0.0]`.
- `seeded_convex_14`: 229.388 → 31.0621 (adverse ratio 0.135413); candidate witness block `quad_-1638d1_-609d1_-2931446d1817_-1040473d1817_-838d1_-457d1_-758d1_-339d1_z0_13`, parameter `[0.0, 0.0, 0.0]`.


Geometry regression cases: none.

## Per-family counts

| Family | n | Geometry baseline/candidate | Equal blocks | Block delta |
| --- | ---: | ---: | ---: | ---: |
| concave_stock_or_cut | 5 | 5/5 | 5 | +0 |
| convex_blind_pockets | 5 | 5/5 | 5 | +0 |
| convex_through_holes | 5 | 5/5 | 5 | +0 |
| fresh_convex_host_rings | 12 | 12/12 | 12 | +0 |
| seeded_convex_blind_pockets | 15 | 15/15 | 15 | +0 |
| seeded_convex_through_holes | 15 | 15/15 | 15 | +0 |

## Cases with at least one quality regression

| Case | Block Δ | Planar κ ratio | Thickness ratio | 3D κ ratio | Corner-sine adverse ratio |
| --- | ---: | ---: | ---: | ---: | ---: |
| blind_triangle_shallow | 0 | 0.853369 | 1.42852 | 1.42852 | 1.55842 |
| blind_two_levels | 0 | 0.992588 | 1.00757 | 1.00757 | 1.04152 |
| seeded_convex_01 | 0 | 0.55229 | 1.16011 | 1.16011 | 0.950968 |
| seeded_convex_02 | 0 | 0.873659 | 1 | 1 | 1.08468 |
| seeded_convex_04 | 0 | 0.254986 | 1.25003 | 1.25003 | 0.431726 |
| seeded_convex_10 | 0 | 0.546792 | 1.06502 | 1.06524 | 0.407151 |
| seeded_convex_11 | 0 | 0.791967 | 1.3131 | 1.31301 | 0.799399 |
| seeded_convex_13 | 0 | 0.573651 | 1.05993 | 1.05996 | 1 |
| seeded_convex_16 | 0 | 0.703485 | 1.51557 | 1.5157 | 1.06309 |
| seeded_convex_19 | 0 | 0.582351 | 1.11777 | 1.11778 | 1 |
| seeded_convex_23 | 0 | 0.539208 | 1.02974 | 1.02975 | 1.08093 |
| seeded_convex_26 | 0 | 0.600135 | 1.11274 | 1.02357 | 0.810514 |
| seeded_convex_28 | 0 | 0.272 | 1.26517 | 1.26389 | 0.534669 |
| fresh_convex_01 | 0 | 0.810995 | 1.03751 | 0.810995 | 1 |
| fresh_convex_03 | 0 | 0.655395 | 1.00786 | 0.655395 | 0.969905 |
| fresh_convex_04 | 0 | 0.828387 | 1 | 0.828387 | 1.69703 |

## Reproduction

Run the comparison script from the checkout root with the frozen release binary, baseline review ledger, and baseline model directory:

```sh
python3 scripts/compare_ring_lookahead.py \
  --cli /path/to/solidcraft-trivariate \
  --baseline-ledger /path/to/clean-ledger.json \
  --baseline-models /path/to/replay \
  --work build/ring-lookahead-candidate \
  --output-dir docs/results/polygon-extrusion-quality/ring-lookahead
```

Machine-readable paired results are in `paired-results.json`; frozen supplement recipes are in `supplemental-cases.json`. The candidate command also writes its detailed raw run ledger under the selected work directory.
