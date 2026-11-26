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

The [NURBS comparison](compare.html) shows both audited layouts for each of
the 46 polygon benchmarks. To regenerate the compact display geometry, run
`python3 scripts/check_sketch_nurbs.py --cli target/release/solidcraft-trivariate`
followed by `python3 scripts/gallery_nurbs_comparison.py`. Integer homogeneous
control points, span incidence and patch membership are retained in the viewer data.


### Reading a partition

The **NURBS volumes** view assigns one colour to each volume. Use **Bézier cells** to see its constituent polynomial or rational spans. A single-span volume has the same extent in both views. With a cut enabled, the colours on the section follow the same assignment as the outer surface.

For each of the 46 polygon benchmarks, the catalogue compares the two qualified candidates by volume count, then cell count. An exact tie keeps the existing cell layout. Both candidates remain inspectable in the synchronized comparison. For example, seeded 14 has 18 volumes over 30 cells on its original layout and 14 volumes over 34 cells with reselected separators. These counts describe the candidates evaluated by the kernel. They do not establish a global minimum.


For a curved model, or a selected polygon layout consisting entirely of single-span volumes, the gallery also tests joining adjacent compatible cells in pairs. Where available, download the **NURBS representation** alongside the native model. The representation supplies the homogeneous control lattice, knot vectors, and source-cell parameter maps. The replay checks every elevated control coefficient and evaluates rational witnesses on both sides of a seam. Internal knots may retain C0 continuity. Connections to neighbouring volumes are recorded on the original cell faces, which can be subfaces of a combined volume. This geometry check does not qualify a new conforming IGA analysis space.
