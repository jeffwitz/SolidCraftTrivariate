# Crossed sketch tool Booleans

The first milestone uses the archived S2 crossed-tool recipe. The acceptance target is a disjoint collection of rational volume charts for intersection, both differences and union. A is the axial cubic profile tool and B is the transverse rectangular tool. Compare integrated chart volumes to the independent SolidCore measurements within 1e-6 cubic millimetres. The construction requires paired monotone profile arcs with identical rational ordinate and weight controls up to scale. Unresolved pairings and unsupported directions must return an error before writing output. Test the four material classes in physical space separately from the volume check.

## Run

The public Rust API constructs one `g2_cross_plate::boolean::CrossedTools` arrangement with `from_recipe`. Call `charts` on it to select any of the four operations without repeating curve extraction or intersection isolation. The CLI accepts both S2 sketch encodings.

```sh
cargo run --locked --offline -- crossed-boolean tests/references/g2_cross_sketches.recipe.json intersection build/crossed-intersection.json
cargo test --locked --offline --test crossed_tool_boolean
python3 tests/references/crossed_tool_oracle.py target/debug/solidcraft-trivariate
```

To incline the transverse sweep, wrap the same sketch recipe with the fields `family: "inclined_crossed_tools_v1"`, `sketch_pair` and integer `transverse_rise_mm`. The displacement from start to finish is (stock width, 0, rise) in millimetres. The vertical extrusion stays fixed. Both transverse caps must remain strictly between the stock caps over the full width.

## Construction and qualification

The two rectangular ordinates split the periodic cubic at exact rational parameters. Pair the left and right arcs over each ordinate interval, retaining their two directed source-parameter intervals. Ascending projected ordinate controls with positive weights give a strictly increasing ordinate in the open span. Common homogeneous weight and ordinate rows make the interpolation linear in x. The control differences certify positive interior width. Selecting the occupied columns and merging consecutive retained height intervals gives degree (3,1,1) trivariate charts in a positive integer homogeneous gauge. The tilted transverse height planes use the exact expression z = z_start + rise * x / width. The complete profile column survives unchanged in a union.

| Operation | Charts | Volume in cubic millimetres |
|---|---:|---:|
| intersection | 2 | 712.3493494135 |
| a-minus-b | 8 | 2777.3825035800 |
| b-minus-a | 4 | 3309.8728728087 |
| union | 10 | 6799.6047258023 |

The independent oracle reads the serialized controls and evaluates the source sketch with SciPy de Boor evaluation. Each of 12 combinations (four operations at rises −3, 0 and +3 mm) passes 5,145 physical probes, including offsets of 0.0001 mm about curved and inclined faces. Each probe requires exactly one covering chart in material and zero in void. Gauss orders 12 and 20 reproduce the volume targets within 1e−6 cubic millimetres. The archived SolidCore run supplies the perpendicular reference. For the tilted cases, translating each transverse z interval by a function of x preserves its height and its intersection volume with the full-height axial tool. Point tests distinguish the three geometries.

The output schema `crossed_tool_charts_v1` stores the source recipe, crossings in the base sketch, paired source intervals and volume controls. The chart count is a construction result. Some interfaces are partial faces after height merging, and the two polar edges of the profile can collapse. This first schema has no assembled interface constraint space, certified mechanical result or STEP export. Its admission domain requires exactly two resolved crossings per transverse ordinate and matching left/right ordinate parametrizations. Free-form surface intersections, tangencies, coincident faces, full containment and chained cuts remain open.
