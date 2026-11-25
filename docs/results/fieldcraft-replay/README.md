# Historical 48-block FieldCraft bracket models

[Current 10-block partition and comparison](../fieldcraft-compact.md)

[Download the five archived models, recipes, STEP files, and audits](./fieldcraft-native.zip)

The original bracket has 60 and 50 mm arms, an 8 mm wall, a 40 mm extrusion, and an inside fillet of radius 6 mm. The archive contains five 48-block models: the original x64-y64 case and four campaign variants that change only the bore diameters. The first bore is centered at x=40 mm and drilled along Y. The second is centered at y=32 mm and drilled along X. Both are centered at z=20 mm.

| Case | Bore diameter along Y (mm) | Bore diameter along X (mm) | Analytical volume (mm³) |
|---|---:|---:|---:|
| x64-y64, original | 8 | 8 | 32,144.778925 |
| x40-y80 | 5 | 10 | 32,163.628481 |
| x80-y40 | 10 | 5 | 32,163.628481 |
| x50-y74 | 6.25 | 9.25 | 32,165.984676 |
| x51-y74 | 6.375 | 9.25 | 32,156.069024 |

The historical Rust constructor takes the L outline, fillet, and the centers, axes, and radii of both bores. Fractional coordinates use a sixteenth-millimeter grid followed by exact rational scaling. The two regions around the bores share the same Z levels. Each archived model has 48 blocks and 74 conforming interfaces. These counts describe the historical partition and do not establish minimality. The current sketch-based construction and its 10-block results are documented separately above.

FreeCAD/OCCT re-imported each STEP as a valid one-shell solid. The five STEP files agree with their corresponding SolidCraft reference files at all 1,920 grid points and 80 probes placed ±0.02 mm from the bore walls. Adaptive mass-property integration found a maximum volume difference below `5 × 10⁻⁹ mm³` for each corresponding pair.

After extracting the archive, the CLI can rebuild each solid with `build`, check it with `audit`, and export it with `export-step`. The archive README describes the input and results for the mechanical check on the original part.

[Independent check of the five exports](../fieldcraft-kernel-replay.md) · [FieldCraft model source](https://github.com/tangency-labs/fieldcraft/tree/35e834652d8548384c2da8ef1b28d588c9af2ffd/tests/corpus/filleted-bracket) · [Bracket mechanical replay](../fieldcraft-bracket-pa5-mechanics.md)
