# Geometry recipe replay inventory

This report inventories standalone geometry recipes and recipe cases embedded in JSON corpora under `tests/references/` and `tests/regressions/`. It records builds and audits with the current CLI. Corpus fields named `expected`, `accepted` or `admitted` set each case's expected build result. A refused build is also checked to preserve a pre-existing output file. The run writes generated files under `/tmp`; it does not publish or replace model assets.

## Reproduce

From the repository root, with the locked offline release binary already built:

```sh
python3 tests/references/replay_geometry_recipes.py \
  --cli target/release/solidcraft-trivariate \
  --work /tmp/solidcraft-recipe-replay \
  --output docs/results/reference-recipe-replay.json
```

The script hashes the CLI before and after replay, records the checkout revision and all case outcomes, and returns nonzero if an expected build or audit fails. The report JSON is the per-recipe evidence and includes the intentional refusal rows.

## Result

At revision `3e48ae3d5761f2c838ab30a86435c605463fa6b7`, the release CLI replayed **234 cases**: 167 expected admissions and 67 corpus-declared refusals. All 234 matched their expected build results; every admitted model passed `audit`, and each refusal preserved its target file. The CLI SHA-256 and per-case records are in [reference-recipe-replay.json](reference-recipe-replay.json).

The replay manifest has SHA-256 `ee5867b6424ddbe95349d68f51d2feca0c432a7d5010d565409f73dc89dbab21`; the runner SHA-256 is `58c25670391d86dfd99ed4a5aa604c637f3271d0b52e077706936b780cf1f74e`. These bind the result to the exact inputs and script even though the checkout has an uncommitted review overlay.

The S3 principal and open-side recipes each embed a parent model. Their embedded parents compare exactly, as parsed JSON values, with the model rebuilt from `tests/references/g2_cross_sketches.recipe.json`. Both S3 models also build and audit successfully.

## Stale and non-build inputs

The standalone `tests/references/g2_cross_plate.recipe.json` used a retired top-level rectangle/profile schema and failed with `unexpected object fields`. It is now expressed with the current sketch recipe schema by reusing the already-existing `g2_cross_sketches.recipe.json` geometry. The stock, profile poles, through directions and rational rectangle coordinates are unchanged; the canonical recipe builds and audits. The full independent junction script was rerun against this path and the frozen CLI: it verified four seam generators, four seam arcs, C0 IGA with q16 residual `2.87 × 10⁻¹⁵`, and all five production variants.

The renamed `tests/references/g2_plate_jets.input.json` was independently reanalyzed by `g2_plate_jets.py`; its curve-jet result remains `verified` with the same source-byte hash. It is kept as an analysis input rather than presented as a production solid recipe.

Three files carry family labels but are fixtures for other commands, not build recipes, so they are listed separately rather than counted as expected refusals:

- `tests/references/g2_plate_jets.input.json` is the historical S1 curve-analysis input. Its `g2_plate_v1` geometry builder is not registered in the current production CLI; the curve jet script can analyze its profile without building a solid.
- `tests/regressions/pocket_wall_internal_face.json` is an audit face-selection assertion.
- `tests/regressions/rational_sector_endpoint_y_control.json` is a homogeneous-control mutation witness for the rational-sector audit.

The replay inventory excludes only these exact paths and fails if their family labels change. Their purpose and exclusion are serialized in the report JSON.

## Current FieldCraft STEP qualification

The five independent SolidCraft source command files are preserved under `tests/references/fieldcraft_scripts/`. They were extracted from `docs/results/fieldcraft-replay/fieldcraft-native.zip` at revision `3e48ae3d5761f2c838ab30a86435c605463fa6b7`; their bytes match the recorded source hashes:

| Input | SHA-256 |
|---|---|
| `x40-y80.script.json` | `830ee65bd8ad0ba9231e7522d403815f36d4edc30ddb3853902484827e54e924` |
| `x50-y74.script.json` | `445614af81985a823ef4ee1019409bad2ccce0d36602ab33ece7153417fa6af0` |
| `x51-y74.script.json` | `1351082f4b871ad3260419b7c641d4358ee3c941a784001336bb26f7cd0ca03c` |
| `x64-y64.script.json` | `b9cfc6a449ca4ec54022b9cf7b8ef46adc2017548b93e9441cb10bdb0565ee7f` |
| `x80-y40.script.json` | `89d24566cd6ae4dd2e8298ec5c17109be9c30b50ea810c35cf32bdf02a572e3a` |

The oracle defaults to these scripts and derives bore-probe centers and radii from them, independently of the current production recipes. The command sequence encoded in each source is `sketch.create`, `sketch.line`, `sketch.finish`, `solid.extrude`, `solid.fillet`, and two `solid.hole` operations. This is the recorded SolidCraft design command sequence. The archive and host timing records do not retain a shell invocation for exporting the independent host STEP, so no export command is inferred here.

The separate [current FieldCraft oracle report](fieldcraft-compact-current-oracle.json) rebuilds and audits the five published recipes under `docs/site/data/fieldcraft_*`, checks exact parsed-model parity with the published files, then compares each current STEP against the independent host STEP under `<host-dir>/<case>/host.step`. All five passed. Each model has 10 blocks, 11 interfaces, 38 exterior sides and exact homogeneous trace equality. The two independent classification grids test 7,680 points across the five brackets with zero mismatches. Another 80 bore-radius probes and 40 tangency-plane probes all pass; each bounding-box difference is zero. The maximum adaptive OCCT mass-property differences are `2.19 × 10⁻¹¹ mm³` in volume and `1.64 × 10⁻⁷ mm²` in area. The report records the FreeCAD version, CLI hash and per-input hashes, including the source script used for each probe set.

From the repository root in a checkout containing the independent STEP inputs at `build/solidcraft-bridge/`:

```sh
freecadcmd tests/references/fieldcraft_compact_model_oracle.FCMacro
```

The macro accepts `FIELDCRAFT_*` environment overrides for the CLI, work directory, source-script directory, host STEP directory, OCCT volume helper, selected cases, and output. The default source-script directory is the tracked `tests/references/fieldcraft_scripts/`; `--source-dir` can select another directory explicitly. The oracle also rebuilds each current recipe with `target/release/solidcraft-trivariate` and checks it against the published model before running STEP comparisons.

## Limits

This replay establishes CLI build and model audit outcomes for the listed recipes. It does not rerun external CAD, STEP, IGA, or mechanics oracles; those remain covered by their dedicated scripts and measured reports. The corpus-declared refusal outcome is not an independent proof that each rejection is scientifically necessary.
