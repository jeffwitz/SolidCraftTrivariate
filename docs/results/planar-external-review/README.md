# External polygon partition review

## Scope and reproducibility

This is an external set of 11 polygonal recipes supplied in `heldout_11_recipes.zip`. It is separate from the 45 published polygon cases and the 12 additional cases. The prototype used mixed_sizes during development; the other ten recipes were withheld. The exact recipe bytes are preserved under [`inputs/`](inputs/); [`inputs.sha256.json`](inputs.sha256.json) records the archive SHA-256 and every recipe SHA-256. The ZIP itself and supplied model file are not copied into this result directory.

The supplied archive SHA-256 is `308f53a93925dbe2aa86759390c6a21943813a92f36889a9221b104adc6b4f89`. The standalone `mixed_sizes.recipe.json` supplied alongside it matches the archive entry byte for byte.

Both CLIs ran sequentially against the same 11 source recipes. The baseline used the clean release CLI from revision `72e8e19c1bd4b25b41610a11d930f38543e40d40` (SHA-256 `8f33eeba801e38e56fbe52ad6ccb96a68a28d2f20187a7d67df869af0c0f9ac2`). The frozen prototype used revision `f406011839fdf6dd92f3d7857c53b2b35c90231e` (SHA-256 `e2ac9c06f1148c914f53c13ad05e4f71c6900847c758647bcd7e0470cbd55bfc`) with `partition_strategy: ring_planar_lookahead` inserted into a temporary build copy of each recipe. The versioned input files were not edited. [The per-case JSON report](paired_evaluation.json) includes exact input and build-recipe hashes, model and STEP hashes, native build/audit reports, per-layer coverage, independent topology/provenance results, all sampled metrics, wall times, and failure details.

The independent geometry audit uses exact `Fraction` polygon moments and volume, Shapely unions and pairwise intersections on every axial layer, a separately constructed FreeCAD/OCCT prism-minus-prisms reference, and imported STEP validity, solid count, volume, centroid, and deterministic point classification. Planar quality metrics were computed by `scripts/evaluate_polygon_extrusion_quality.py` (SHA-256 `018a75fcf2c379012e25afee2678d858105cfe29b845d3812fb630fe0fddec14`): the maximum planar spectral condition number and thickness aspect use a 5 × 5 × 5 parameter grid per trilinear block; minimum corner sine uses the eight corners. These samples are diagnostics and do not certify a continuous condition bound.

Recorded environment: FreeCAD 26.3.0, OCCT 7.8.1, Python 3.11.14, NumPy 1.26.4, and Shapely 2.1.2. Per-case timings below are sequential wall-clock seconds. Build, explicit native audit, STEP export, and independent audit/metric costs are reported separately in the JSON.

`run_phase.py` is the versioned FreeCADCmd runner. Set `SCT_REVIEW_PHASE`, `SCT_REVIEW_CLI`, `SCT_REVIEW_CHECKOUT`, `SCT_REVIEW_WORK`, `SCT_REVIEW_OUTPUT`, and `SCT_REVIEW_METRICS_SCRIPT`, then run `freecadcmd docs/results/planar-external-review/run_phase.py`. For the prototype phase, also set `SCT_REVIEW_STRATEGY=ring_planar_lookahead`; for baseline parity, set `SCT_REVIEW_SUPPLIED_MODEL` to the provided mixed-size model. Generated models, STEP files, and intermediate phase reports belong under `build/` and are not part of this versioned input/result set.

## Results

Both strategies built and passed the native and independent checks for all 11 recipes. They produced 219 blocks in total each, with identical exact material volumes for every recipe. Cross-section symmetric difference and pairwise overlap were zero in all layers; the largest floating-point area-sum discrepancy was 1.82 × 10⁻¹² mm². Across each phase, 3,706 independent classification probes had zero mismatches. The largest CAD/STEP volume discrepancy was 6.67 × 10⁻¹⁶ relative, and the largest centroid discrepancy was 1.02 × 10⁻¹² mm. No cases failed or were refused.

| Recipe | Blocks, both | Max planar κ, baseline → lookahead | Min corner sine, baseline → lookahead | Build + audit seconds, baseline → lookahead |
|---|---:|---:|---:|---:|
| L_stock_two_holes_one_arm | 17 | 3.7500 → 3.7500 | 1.0000 → 1.0000 | 0.0244 → 0.5231 |
| blind_two_depths_plus_hole | 48 | 11.3389 → 11.3389 | 0.3069 → 0.3069 | 1.1125 → 2.3063 |
| grid_8_square_holes | 32 | 7.7205 → 6.6569 | 0.3162 → 0.3681 | 1.1944 → 3.5920 |
| hole_1mm_from_edge | 4 | 109.0008 → 109.0008 | 0.0333 → 0.0333 | 0.0361 → 0.0806 |
| mixed_sizes | 14 | 558.9996 → 45.4834 | 0.0037 → 0.0656 | 1.4240 → 4.6310 |
| pockets_top_and_bottom | 28 | 7.0923 → 6.8094 | 0.5145 → 0.5300 | 0.1525 → 0.7529 |
| skewed_hexagon_div40 | 12 | 14.2666 → 14.2666 | 0.2418 → 0.2418 | 0.9959 → 2.0167 |
| thin_web_1mm | 8 | 133.3342 → 133.3342 | 0.0250 → 0.0250 | 0.0949 → 0.1452 |
| three_holes_triangle_layout | 16 | 14.8024 → 12.8764 | 0.1750 → 0.2119 | 2.6212 → 4.7646 |
| two_16gons | 32 | 8.7352 → 8.7352 | 0.3764 → 0.3764 | 3.4571 → 7.2961 |
| two_triangles | 8 | 27.2133 → 22.8291 | 0.0948 → 0.0906 | 2.3321 → 4.3582 |

The lookahead reduced the sampled maximum planar condition number on five cases and left six unchanged; it worsened none on the eleven external cases. It left block counts unchanged. The minimum corner sine improved on four cases, was unchanged on six, and decreased slightly on `two_triangles`; that case's sampled planar condition number still improved. Full metric witnesses are retained per case in the JSON.

Aggregate costs were:

| Phase | Build | Explicit audit | STEP export | Independent audit + metrics |
|---|---:|---:|---:|---:|
| Baseline | 9.108 s | 4.337 s | 4.442 s | 91.729 s |
| Ring lookahead | 19.518 s | 10.949 s | 11.020 s | 147.745 s |

Build plus explicit audit took 13.445 s baseline and 30.467 s with lookahead (2.27×). These are one-pass totals, with no warmup or repeated timing runs.

## Supplied mixed-size model

The supplied `mixed_sizes.model.json` has SHA-256 `949649c718a79485a1e1147737e088253d46bea07b355db727a7a4556731397c`. Building its matching recipe with the clean baseline CLI produced a byte-identical model with the same SHA-256. The lookahead output has a different model hash because it selected a different partition, and it independently passed exact-volume, coverage, topology, provenance, STEP, and classification checks for the same physical material region. The supplied model itself is not duplicated here.

## Failures

There were zero build refusals and zero native or independent audit failures across both phases. The failure list is empty in `paired_evaluation.json`.
