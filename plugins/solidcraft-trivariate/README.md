# SolidCraft Trivariate Plugin

This is a separately installable SolidCraft command plugin. It snapshots the active
SolidCraft timeline for the native `solidcraft-trivariate` companion, displays a
companion-produced native sketch preview on explicit request, and inserts the companion's
STEP result through SolidCraft's own STEP importer.

## Commands

| Entry | Action |
| --- | --- |
| `snapshot-timeline` | Calls `timeline.export_script`, refuses a script with omitted design content, and writes `timeline.json` to the approved `companion` folder. |
| `extrude-active-sketch` | Builds the selected pre-extrusion sketch as a new native-kernel body using `height_mm` (default 12, allowed 1–64 mm), audits it, and inserts its STEP payload as a native import feature. |
| `extrude-active-timeline` | Uses the existing sketch-plus-extrusion timeline importer. This is separate from sketch-only extrusion so it cannot duplicate an existing body by default. |
| `build-recipe` | Reads a JSON family recipe from the approved companion folder, builds and audits it with the shared Rust constructor, exports STEP, and inserts it as a native import feature. For the qualified `fieldcraft_x64_y64` recipe, `rim` can be `sharp` (default), `chamfer`, or `fillet`. |
| `show-sketch-stage` | Reads `partition.json`, verifies the complete source script against the open document, adds the companion's ordinary sketch feature with `timeline.add`, and shows it with `browser.visibility`. |
| `hide-sketch-stage` | Rechecks the current timeline and preview feature, then hides that sketch with `browser.visibility`. |
| `insert-result` | Reads `result.step` from the approved `companion` folder and adds it through SolidCraft's native timeline import feature, then sets the active viewport to native `shaded` style. |

There is no automatic preview creation or visibility change when installing, enabling, or
opening a design. Before the first explicit Show Sketch Stage action, the plugin has added
no preview feature. The preview descriptor is restricted to one plane and construction
line curves. SolidCraft renders and tessellates the native sketch and imported STEP with
its own renderer and kernel; this plugin creates no mesh, tessellation, or custom viewport
content. Inserting the result changes the active view style to shaded, a global display
setting; users can restore another style in SolidCraft's View menu.

The companion descriptor contract is `solidcraft_sketch_preview_v1`. It contains
`source_feature_id`, the exact `source_feature` and `source_script`, a `feature` sketch
with name `SCT Partition <source_feature_id>` and its assigned feature id, and optional
`internal_edges`. A preview installed by this plugin has its id recorded back into the
descriptor. Before every show or hide, the plugin verifies that the open timeline is the
recorded source timeline after removing only that exact preview `timeline.add` step. A
changed document or changed preview is refused. Companion files are read and written only
inside the folder the user approves as `companion` during installation.

The imported STEP text is carried in the native `timeline.add` feature, so host display and
tessellation use SolidCraft's reader and renderer without an external path or generated mesh.
`timeline.export_script` is SolidCraft's actual replayable command history. SolidCraft
0.48.1 refuses designs containing components or assemblies, and the plugin refuses any
export whose `omitted` list is nonempty. Visibility changes are saved with the design and
are undoable, per the host's `browser.visibility` contract.

## Build and check

Build a core Wasm module and wrap it as a component. SolidCraft does not link WASI
interfaces, so a `wasm32-wasip2` component is refused by the host:

```sh
rustup target add wasm32-unknown-unknown
cd plugins/solidcraft-trivariate
cargo build --locked --release --target wasm32-unknown-unknown
wasm-tools component new target/wasm32-unknown-unknown/release/solidcraft_trivariate_plugin.wasm -o solidcraft-trivariate.wasm
solidcraft-cli plugin check .
solidcraft-cli plugin test . --require-cases
```

`extrude-active-sketch` reads the specified sketch before any host extrusion and accepts
the native line-and-arc profile subset supported by `solidcraft_polygon::extrude_source_sketch`.
`extrude-active-timeline` is the separate path for an already extruded sketch. The plugin
does not make a separate mesh; SolidCraft imports the STEP text as a native timeline
feature and renders the resulting body itself.

The component imports the `craft:plugin@0.2.0` `commands` and `files` interfaces from
SolidCraft's plugin contract. The workflow at
`.github/workflows/solidcraft-trivariate-plugin.yml` builds it, checks it with the
SolidCraft 0.48.1 CLI, tests the entries, and packages `plugin.json` with the component.

## Upstream contract provenance

The minimal WIT contract in `wit/plugin.wit` is adapted from SolidCraft's
`crates/craft-plugin/wit/plugin.wit` and `docs/plugins.md` at upstream revision
`765f26e0d1fc7c66a0d03ea2dae28ee3b65270b8` (SolidCraft 0.48.1). Upstream is MIT/Apache-2.0
licensed. This plugin uses the documented WIT imports and host commands; it contains no
copied SolidCraft implementation code. Keep the upstream origin and revision with any
redistribution of the WIT-derived interface.
