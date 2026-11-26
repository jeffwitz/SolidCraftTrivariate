# Fixed-partition bore-edge edit preregistration

## Starting point and scope

- Starting revision: `34b9720f059d9e93bb093edae557e77031f25582`.
- Target: the actual compact FieldCraft bracket in `docs/site/data/fieldcraft_x64_y64/fieldcraft_x64_y64.model.json`, with the R6 elbow, two R4 transverse bores, 60 × 50 × 40 mm stock, in millimetres. The unequal-leg fixture in `tests/references/unequal_fillet_cross_bores.recipe.json` is secondary only.
- Owned files: `src/bracket_edge_edit.rs`, `examples/bracket_edge_edit.rs`, focused Rust tests and this report.
- The existing ten logical blocks and their adjacency are fixed. Each operation must preserve this partition and every unaffected control net. No remeshing, block insertion, trimming, or external CAD Boolean is allowed.
- Implement chamfers first on both axial ends of both bores. Then make an explicit rational round-fillet attempt on the same four circular edges; singular/unsupported domains must be reported explicitly.

## Preconditions and method

Use the repository's exact rational bore charts and the two axial endpoint traces. The requested edge treatment is a piecewise rational control-net edit over the existing chart spans. Preserve exact source geometry away from the treatment support; interfaces are derived from common edited traces, never matched by sampling alone.

Before accepting an edit, check (1) positive Jacobian on every closed span using the native spatial certificate, (2) exact homogeneous equality on every declared interface after consistent orientation/normalization, (3) unchanged block count and topology, (4) independent section/volume and radius checks, and (5) unchanged control data outside the declared support. A sampled geometry check is diagnostic only.

## Frozen thresholds and rejection cases

- Jacobian: certificate status must be `Valid` with a strictly positive lower bound on every span; `Unknown` or `Invalid` is rejection.
- Interfaces: exact integer homogeneous equality after common scaling; no tolerance-based acceptance.
- Analytic checks: edge radius and axial endpoint location within `1e-10 mm`; volume against an independently integrated exact/rational section formula within relative `2e-10`.
- Reject zero/negative chamfer width, width reaching the bore radius, overlapping treatment supports, nonintegral/unrepresentable chart subdivisions, any changed nonincident trace, and any Jacobian/interface indeterminate.
- Fillet is admitted only where a nondegenerate whole-block neighborhood exists on the fixed partition. If a circular edge intersects a singular chart edge/corner or requires a changed partition, record that domain as unsupported; do not silently regularize or approximate it.

## Independent oracle and exit criteria

The oracle will independently evaluate the source and edited rational maps using exact fractions at endpoints and interior parameters, derive section radii and volume independently, and enumerate all ten-block incidences from face control nets. Exit requires a reproducible example, focused tests for positive and rejected cases, certificate and interface summaries, a clear fillet-domain result, and a report of limitations. No claim is made beyond the tested fixed partition and parameter range.

## Correction before implementation

The first preregistration draft named the unequal-leg regression fixture. After initial implementation work had begun, but before qualifying or externally measuring that fixture, the parent clarified that the target is the published FieldCraft R6/R4 compact bracket. Its primary fixture is now authoritative. The unequal-leg recipe remains a secondary optional case. No acceptance threshold was changed.
