# Planar and thickness quality of the 45 polygonal gallery models

This report records an independent Python evaluation of the 45 gallery models in the frozen polygonal planner corpus. The case IDs come from `docs/results/sketch-planner-polygon-evaluation/paired.json`; later gallery additions are excluded. The evaluated model JSON files are listed with SHA-256 digests in `results.json`.

## Method

The evaluator reads the homogeneous trilinear block controls directly. It rejects blocks whose degrees are not `[1, 1, 1]`, which do not have eight controls, or whose control weights are not one common positive value. It does not call the Rust builder, audit, or quality routines.

At every point of the closed 5×5×5 parameter grid, it computes:

- **Planar spectral condition number**: `sigma_max / sigma_min` of the 3×2 Jacobian formed by the in-plane parameter derivatives. The smaller singular value is evaluated as `area / sigma_max`, avoiding cancellation in the small eigenvalue.
- **Thickness-direction aspect**: `sigma_max / abs(dw · n)`, where `n` is the unit normal to the local planar section and `dw` is the third Jacobian column. This measures planar stretch against normal thickness and is insensitive to in-plane shear in the extrusion direction.
- **Maximum 3D spectral condition number**: largest divided by smallest singular value of the full 3×3 Jacobian. This remains a separate measure and includes skew.

The minimum corner sine is the minimum sine between the two in-plane tangent vectors at the eight parametric corners. The report also records the smallest sampled Jacobian determinant. Every case metric has a block and parameter witness. Ties follow manifest order and block order.

These are sampled diagnostics, not continuous quality bounds. A poor value identifies a witness to inspect; a good sampled value does not certify an interior minimum. Determinant certificates and topological validity are outside this report.

## Results

All 45 frozen case IDs were present in the gallery manifest and evaluated: 1,208 active blocks in total, with 5 to 74 blocks per case. For the first three metrics, the table gives median / nearest-rank p95 / maximum. Corner sine is adverse in the lower tail, so its row gives median / nearest-rank p05 / minimum.

| Metric | Median | p95 | Maximum and witness case |
| --- | ---: | ---: | --- |
| Planar spectral condition number | 27.172 | 123.942 | 229.388, `bench_poly_seeded_convex_14` |
| Thickness-direction aspect | 123.362 | 596.151 | 892.902, `bench_poly_seeded_convex_09` |
| Maximum 3D spectral condition number | 130.286 | 596.151 | 892.902, `bench_poly_seeded_convex_09` |
| Minimum corner sine | 0.17844 | 0.037475 (p05) | 0.010746 minimum, `bench_poly_seeded_convex_14` |

The planar-condition witness is block `quad_-158d1_161d1_-17482d179_11069d179_-78d1_279d1_28891d358_363551d716_z0_13` at parameter `[1, 1, 0]`. Its planar singular values are 312.788 mm and 1.36357 mm. The same case has the minimum corner sine, witnessed at that block's `[1, 1, 0]` corner.

The maximum thickness-direction aspect is in block `quad_-1218d1_-1018d1_-226d1_-322d1_222d1_-148d1_1022d1_-148d1_z0_3` of `bench_poly_seeded_convex_09` at `[0, 0, 0]`: planar `sigma_max` is 2,678.707 mm and normal thickness per parameter unit is 3.000 mm. The maximum 3D condition number occurs at the same location and is 892.902. This separation makes clear that the thin direction accounts for nearly all the 3D conditioning in this witness.

## Reproduction

From the repository root:

```sh
python3 scripts/evaluate_polygon_extrusion_quality.py \
  --manifest docs/site/experiments.json \
  --models docs/site/data \
  --baseline-cases docs/results/sketch-planner-polygon-evaluation/paired.json \
  --output docs/results/polygon-extrusion-quality/results.json
```

The JSON records the source revision, Python version, evaluator and manifest hashes, each model hash, per-case extrema, and witnesses. `summarize_model(model_dict)` exposes the same case-level summary for another gallery or benchmark integration.


## Geometry at the largest in-plane condition numbers

The following rows trace the two incident face edges at the worst sampled planar corner back to the source sketch. Vertex indices start at zero in the recipe. The angle is the interior angle of the material cell, measured between its two incident edge vectors. An acute angle and its supplement have the same sine.

| Case | In-plane κ | Interior angle | Sketch vertex at the corner | Other ends of the two incident edges |
| --- | ---: | ---: | --- | --- |
| seeded_convex_14 | 229.39 | 179.38° | cut_2[6] | cut_2[5] and an interior construction point near (80.701, 507.753) |
| seeded_convex_26 | 126.54 | 177.07° | cut_2[10] | cut_2[0] and stock[9] |
| seeded_convex_04 | 123.94 | 177.94° | cut_2[4] | cut_2[3] and stock[5] |

Case 14 is a candidate for changing a separator or the points introduced on it. In cases 26 and 04 the cell connects directly to an original stock vertex. An alternative partition could change which hole and stock vertices share a cell. The witness alone does not identify the construction decision that first made this pairing unavoidable.
