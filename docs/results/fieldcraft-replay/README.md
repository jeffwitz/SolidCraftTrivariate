# FieldCraft sketch models

The current replay bundle contains five two-bore brackets and a four-sector control plate. The CLI builds the brackets from the L-shaped sketch, inside fillet and bore axes. Each has ten blocks and eleven interfaces. The ZIP includes the source recipe, trivariate model, native audit and STEP export for all six solids.

[Download the current bundle](fieldcraft-compact.zip) · [Inspect the original bracket](https://jeffwitz.github.io/SolidCraftTrivariate/?model=fieldcraft_x64_y64)

| Case | Bore diameter along Y (mm) | Bore diameter along X (mm) | Analytical volume (mm³) |
|---|---:|---:|---:|
| x64-y64, original | 8 | 8 | 32,144.778925 |
| x40-y80 | 5 | 10 | 32,163.628481 |
| x80-y40 | 10 | 5 | 32,163.628481 |
| x50-y74 | 6.25 | 9.25 | 32,165.984676 |
| x51-y74 | 6.375 | 9.25 | 32,156.069024 |

To rebuild one bracket directly from its tracked recipe:

```sh
cargo build --release --locked --offline
mkdir -p build/fieldcraft-replay
cli=target/release/solidcraft-trivariate
"$cli" build docs/site/data/fieldcraft_x64_y64/fieldcraft_x64_y64.recipe.json build/fieldcraft-replay/bracket.json
"$cli" audit build/fieldcraft-replay/bracket.json
"$cli" export-step build/fieldcraft-replay/bracket.json build/fieldcraft-replay/bracket.step
```

[Volume and STEP qualification](../fieldcraft-compact.md) · [Mechanical analysis](../fieldcraft-bracket-pa5-mechanics.md)
