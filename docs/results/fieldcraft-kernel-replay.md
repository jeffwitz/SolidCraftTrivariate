# FieldCraft bracket replay: independent STEP qualification

## Scope and provenance

Five reconstructed FieldCraft bracket variants were compared with their matching SolidCraft host STEP outputs. The upstream recipe and solver expectations are pinned to FieldCraft revision `35e834652d8548384c2da8ef1b28d588c9af2ffd` in `references.lock`. The original part is a 60 × 50 × 40 mm L bracket with 8 mm legs, an inside R6 fillet, and two orthogonal through bores. Variant labels identify the paired host script and replay recipe; radii below are stated explicitly to avoid relying on those labels.

This report describes the archived grid. The current constructor selects the [compact partition](fieldcraft-compact.md). The archived candidates are all admitted 48-block models with 74 shared interfaces and 140 exterior sides. The corresponding host B-Rep has 11 faces. Both representations import as one valid solid with one shell. The different face and edge counts reflect the trivariate block boundary partition; they are not expected to match the merged host B-Rep topology.

## Independent checks

`tests/references/fieldcraft_kernel_replay.py` reads each host and candidate STEP with FreeCAD/OCCT. It exports each imported shape to a BREP in a case-specific directory and invokes the separately compiled `occ-volume` helper, which computes adaptive OCCT Gauss–Kronrod volume and surface area at tolerances 1e−9 and 1e−12. It also checks the source cylinder radii, exact analytic volume, bounding box, 384 uniform interior points, and 16 radial points per pair of bores. The radial points lie at r−0.02 mm (expected void) and r+0.02 mm (expected material), at four cardinal directions for each bore. The script returns a failing status if a gate fails.

For these orthogonal through holes, the independent volume formula is

`V = 34080 − π(360 + 8(rx² + ry²)) mm³`,

where `rx` is the radius of the bore along X and `ry` the radius along Y. Each source host volume, replay STEP volume, and analytic value agrees within the recorded integration error and the 1e−8 relative gate. The candidate-to-host volume discrepancy is at most 4.8e−9 mm³ (1.5e−13 relative); the surface-area discrepancy is at most 6.51e−7 mm² (6.2e−11 relative). All five have identical bounding boxes, zero mismatches on the 384-point interior grid, and zero failures among the 16 targeted bore probes.

| Case | Bore radii X / Y (mm) | Adaptive volume (mm³) | Candidate − host volume (mm³) | Candidate − host area (mm²) | Grid / radial probes |
|---|---:|---:|---:|---:|---:|
| x40-y80 | 5 / 2.5 | 32163.62848131003 | −2.07e−10 | −6.51e−7 | 384/0; 16/0 |
| x50-y74 | 4.625 / 3.125 | 32165.98467580065 | +2.44e−10 | −6.51e−7 | 384/0; 16/0 |
| x51-y74 | 4.625 / 3.1875 | 32156.06902398781 | +2.91e−10 | −6.51e−7 | 384/0; 16/0 |
| x64-y64 | 4 / 4 | 32144.77892538703 | −1.64e−9 | −6.50e−7 | 384/0; 16/0 |
| x80-y40 | 2.5 / 5 | 32163.62848130545 | −4.78e−9 | −6.49e−7 | 384/0; 16/0 |

The swapped x40-y80 and x80-y40 cases have equal volumes; the bore-directed probes separately confirm that each radius is on the intended axis. Candidate STEP SHA-256 values are listed in each JSON result under `build/fieldcraft-oracle/`; the same records contain host STEP, source script, recipe, model, and audit hashes.

## Mechanical reference

FieldCraft provides the `filleted-bracket-fillet-stress` job specification and expected values, but not its original finite-element mesh/result artifact in this checkout. The same support, 2 kN +X load, steel material, and R6 bracket are replayed on the baseline x64-y64 trivariate model by `tests/references/mechanical_jobs/fieldcraft_bracket_hex_pa5.json`. The Q2 solve and independent SciPy reassembly are recorded in `docs/results/fieldcraft-bracket-pa5-mechanics.md`: force reaction is −2000 N, relative force balance 1.10e−11, moment balance 5.81e−12, and SciPy displacement coefficients agree to 5.33e−13 relative.

This verifies the transferred load/support job and equilibrium. It does not compare the FieldCraft tet10 peak von Mises stress (about 213 MPa on the fine mesh) against the Q2 result: the available Q2 record contains block-mean stresses rather than the same nodal fillet-face peak, and the original tet mesh and nodal field are unavailable.

The default `FreeCAD.Shape.Volume` value for the imported candidate STEP is not used for acceptance. The piecewise rational B-spline boundary requires adaptive mass integration; the script preserves the default value only as a diagnostic in its raw shape summary. All reported acceptance volumes come from `occ-volume` on the exported BREP.

## Reproduction

Run from the repository root after building the release CLI and with the FreeCAD Python environment and `occ-volume` binary available. The independent host export for each case must already be present under `build/solidcraft-bridge/<case>/host.step`. The supported case IDs are `x40-y80 x50-y74 x51-y74 x64-y64 x80-y40`:

```sh
: "${FC_ENV:?Set FC_ENV to the FreeCAD Python environment}"
: "${FC_BUILD:?Set FC_BUILD to the FreeCAD build containing lib and Mod/Part}"
case_id=x64-y64
host="build/solidcraft-bridge/$case_id"
archive=build/fieldcraft-archive/fieldcraft-native
unzip -oq docs/results/fieldcraft-replay/fieldcraft-native.zip -d build/fieldcraft-archive
variant="$archive/variants/$case_id"
LD_LIBRARY_PATH="$FC_BUILD/lib:$FC_ENV/lib" \
PYTHONPATH="$FC_BUILD/lib:$FC_BUILD/Mod/Part" \
  "$FC_ENV/bin/python" tests/references/fieldcraft_kernel_replay.py \
  "$host/host.step" "$variant/$case_id.step" \
  --output "build/fieldcraft-oracle/$case_id.json" \
  --occ-volume build/solidcraft-bridge/occ-volume --fast \
  --source-script "$archive/replay/inputs/$case_id.script.json" \
  --candidate-recipe "$variant/$case_id.recipe.json"
```

Set `case_id` to the case being checked. BREP intermediates are written under `build/fieldcraft-oracle/C/` so concurrent runs cannot overwrite one another. Full run records and hashes are regenerated in `build/fieldcraft-oracle/`; that directory is build output, not a committed reference.
