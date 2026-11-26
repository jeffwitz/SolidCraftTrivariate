# Independent STEP readback for orthogonal sketch cases

From the repository root, set `OCCT_PREFIX` to the installation containing the OCCT headers and libraries used by FreeCAD. Build the independent mass-properties helper and run the four cases.

```sh
mkdir -p build/solidcraft-bridge
c++ -std=c++17 -I"$OCCT_PREFIX/include/opencascade" tests/references/occ_volume.cpp \
    -L"$OCCT_PREFIX/lib" -Wl,-rpath,"$OCCT_PREFIX/lib" \
    -lTKTopAlgo -lTKBRep -lTKG3d -lTKMath -lTKernel -o build/solidcraft-bridge/occ-volume
cargo build --locked --offline --bin solidcraft-trivariate
target/debug/solidcraft-trivariate build tests/references/unequal_fillet_l.recipe.json build/unequal-fillet-l.model.json
target/debug/solidcraft-trivariate export-step build/unequal-fillet-l.model.json build/unequal-fillet-l.step
target/debug/solidcraft-trivariate build tests/references/slender_plate_host.recipe.json build/slender-plate-host.model.json
target/debug/solidcraft-trivariate export-step build/slender-plate-host.model.json build/slender-plate-host.step
target/debug/solidcraft-trivariate build tests/references/unequal_fillet_cross_bores.recipe.json build/unequal-fillet-cross-bores.model.json
target/debug/solidcraft-trivariate export-step build/unequal-fillet-cross-bores.model.json build/unequal-fillet-cross-bores.step
target/debug/solidcraft-trivariate build tests/regressions/fieldcraft_three_bores.recipe.json build/fieldcraft-three-bores.json
target/debug/solidcraft-trivariate export-step build/fieldcraft-three-bores.json build/fieldcraft-three-bores.step
OCC_VOLUME_BIN=build/solidcraft-bridge/occ-volume freecadcmd --console tests/references/audit_orthogonal_steps.FCMacro
```

The macro reads AP214 STEP with FreeCAD 26.3.0 / OCCT and checks validity, one solid, one closed solid, one shell, face count, and volume. It exports the imported B-Reps and invokes the separate `occ-volume` helper for adaptive Gauss–Kronrod mass integration at relative tolerances 1e−9 and 1e−12. The default `Shape.Volume` is retained as a diagnostic, not as the acceptance volume. Analytic targets are derived from the integer sketch dimensions and circular feature radii.

| Case | STEP solids / shells / faces | STEP volume (mm³) | Analytic volume (mm³) | Relative error | Result |
| --- | ---: | ---: | ---: | ---: | --- |
| Unequal filleted L, no bores | 1 / 1 / 18 | 30761.203552627685 | `30912 − 48π` = 30761.20355262769 | <1e−15 | Valid, closed, volume agrees |
| Slender plate with one through bore | 1 / 1 / 24 | 14060.707993412296 | `14400 − 108π` = 14060.707993412303 | 5.2e−16 | Valid, closed, volume agrees |
| Unequal filleted L with transverse bores in both arms | 1 / 1 / 38 | 30359.079692968182 | `30912 − 176π` = 30359.079692968196 | 4.8e−16 | Valid, closed, volume agrees |
| FieldCraft bracket with an extra R2 bore | 1 / 1 / 52 | 32044.24796047382 | `34080 − 648π` = 32044.247960473815 | 1.14e−16 | Valid, closed, volume agrees |

On the ten-block unequal bracket the default FreeCAD property gives 30359.079691341878 mm³. The earlier 48-block export had 140 boundary faces and gave 30358.52082001161 mm³ with the default property. Both representations pass adaptive integration against the same analytic target. Full measurements are written by the macro to `build/orthogonal-step-freecad-audit.json`.
