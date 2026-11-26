# Gallery asset regeneration

Run from the repository root after building the current release CLI:

```sh
cargo build --release --locked --offline
python3 scripts/gallery_catalog.py
```

The script reads the tracked 68-entry `experiments.json`, rebuilds each model
from its linked recipe, audits it with the Rust CLI, exports a STEP file, and
refreshes the audit and cache-versioned links. Polygon gallery recipes use thicker extrusions for viewing. Each original height and pocket interval was multiplied by the smallest integer of at least two that brings the stock height to at least one sixth of its longest XY span. The frozen benchmark retains its original dimensions. Displayed thickness, volume and quality diagnostics are measured on the downloadable gallery geometry. For polygonal extrusions the audit contains the in-plane condition number, minimum corner sine and ratio of in-plane stretch to normal layer thickness. The WebGL panel displays each diagnostic separately alongside the 3D condition number. The downloadable audit gives their block IDs and parametric witnesses. It preserves prior independent readbacks as
historical records when they belong to an older STEP export. Repeated runs
reuse existing STEP bytes when the only difference is the export timestamp.

The browser check can then be run with:

```sh
node tests/references/gallery_browser.mjs build/gallery-qa
```
