# SolidCraft plugin icons

The plugin package supplies original vector illustrations for its commands. In SolidCraft
0.48.1, `x-solidcraft-icon` names a relative `.svg` path inside the plugin package; it does
not accept a built-in icon identifier. The host checks the SVG against its supported subset
when loading the package. An absent or rejected image falls back to the generic plugin icon,
and the command remains available.

| Plugin command | Package icon | Meaning |
| --- | --- | --- |
| `extrude-active-sketch` | `icons/extrude.svg` | Sketch extruded into a solid |
| `extrude` | `icons/extrude.svg` | Sketch extrusion |
| `drill` | `icons/drill.svg` | Drill bit and work surface |
| `boolean` | `icons/boolean.svg` | Overlapping operands |
| `chamfer` | `icons/chamfer.svg` | Clipped corner profile |
| `fillet` | `icons/fillet.svg` | Rounded corner profile |
| `fieldcraft-iga-solve-display` | `icons/compute.svg` | Analysed patch and solved field |

The shared extrusion drawing reflects one operation across two plugin entries. Commands
without a declared image use SolidCraft's generic plugin glyph. A missing or refused custom
image also falls back to that glyph, though the release checker requires every declared image
to be present. `scripts/check_solidcraft_plugin_icons.py` checks this exact mapping, the
package-local paths, and the SVG constraints before packaging. The plugin workflow also runs the host's
`solidcraft-cli plugin check` and captures the installed command group and dialogs in its
desktop screenshot artifact.

## Native icon identifiers and permissions

SolidCraft 0.48.1 has native operation icon masters named `extrude`, `hole`, `combine`,
`chamfer`, and `fillet`. Their source files are the version-pinned
[extrude](https://github.com/tangency-labs/solidcraft/blob/765f26e0d1fc7c66a0d03ea2dae28ee3b65270b8/crates/look/icons/design/svg/model/extrude.svg),
[hole](https://github.com/tangency-labs/solidcraft/blob/765f26e0d1fc7c66a0d03ea2dae28ee3b65270b8/crates/look/icons/design/svg/model/hole.svg),
[combine](https://github.com/tangency-labs/solidcraft/blob/765f26e0d1fc7c66a0d03ea2dae28ee3b65270b8/crates/look/icons/design/svg/model/combine.svg),
[chamfer](https://github.com/tangency-labs/solidcraft/blob/765f26e0d1fc7c66a0d03ea2dae28ee3b65270b8/crates/look/icons/design/svg/model/chamfer.svg),
and [fillet](https://github.com/tangency-labs/solidcraft/blob/765f26e0d1fc7c66a0d03ea2dae28ee3b65270b8/crates/look/icons/design/svg/model/fillet.svg)
masters at revision `765f26e0d1fc7c66a0d03ea2dae28ee3b65270b8`. They identify the host's
operation names, but the plugin API does not bind to those IDs.

The native icon masters are Brand Material under SolidCraft's
[version-pinned `LICENSE-BRAND`](https://github.com/tangency-labs/solidcraft/blob/765f26e0d1fc7c66a0d03ea2dae28ee3b65270b8/LICENSE-BRAND).
Its sections 3(a) and 3(b) require written maintainer permission before using them in another
program or distributing copies outside SolidCraft. This plugin therefore does not copy or
modify that artwork. Its six SVGs are separate plugin illustrations using `currentColor`
and the plugin icon contract's 64 × 64 vector canvas.

Plugin commands display in their own group. The 0.48.1 ribbon implementation explicitly
keeps plugin buttons out of the native-command pinning flow
([`toolbar.rs`](https://github.com/tangency-labs/solidcraft/blob/765f26e0d1fc7c66a0d03ea2dae28ee3b65270b8/crates/ui-egui/src/toolbar.rs#L2129)).
The documentation's ribbon pin feature applies to built-in commands; it is not a supported
way to move or rebind a plugin operation to a stock Sketch button. The plugin exposes its own
command buttons and uses the command image assigned by its manifest.
