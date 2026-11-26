# Blind evaluation of sketch partitioning

This is the baseline comparison frozen at `4e8a394`. The [current kernel results](convex-host-rings.md) cover the same 45 inputs and 12 additional sketches.

The reference executable is built from `e5bd5c6`. The evaluation owner keeps the generator, seeds, cases and measurements outside the implementation workspace until the candidate revision is frozen. The owner can publish a SHA-256 commitment to a private bundle containing the generator version, seeds and scoring script before development. Specify the family definitions, units, tolerances, time budget and acceptance thresholds in advance. Disclose the bundle after scoring so that others can reproduce the comparison. Cases disclosed during development become regressions for subsequent runs.

## Frozen comparison

The independent evaluator retained the 45 recipes until the kernel was frozen and then published the generator, inputs and paired scores. Both executables passed all 45 geometry gates. The new planner reduced the block count on 15 cases and kept it on 30. Among the 30 randomly generated convex stocks it improved nine of ten single-footprint cases and used the historical triangulation on all twenty multi-footprint cases. Five solids had a larger sampled maximum condition number. [Polygon evaluation and data](results/sketch-planner-polygon-evaluation/README.md) · [Ten public 3D controls](results/propagated-bores.md)

The release checks passed 464 Rust tests, the separately invoked release-mode bracket stiffness test, 58 Python tests, the demonstration commands and strict Clippy on all targets. The candidate also replayed the 45 archived polygon models, the six baseline bracket controls that previously built successfully and the eleven models in the published gallery.

## Paired comparison

Run both executables on exactly the same recipes under the same limits. Report counts and paired changes separately for convex polygonal through-holes, blind pockets, concave stock or cut profiles, and nonparallel or stacked bores. Within each family, vary feature count, vertex count, clearance, aspect ratio, relative placement and orientation. Include rotations or rational-coordinate equivalents supported by the input format. Mark inputs outside its coordinate contract separately. A rejected valid input, timeout or fallback remains in its family's denominator.

Report three distinct outcomes for each admissible input: whether a model passed independent geometric checks, which construction produced it, and its cost and quality. List admitted compact models, valid fallback models, rejected cases and timeouts. Any invalid model returned as valid fails the geometric gate. Record regressions on inputs the reference handled.

- Check material coverage, absence of overlap, boundary provenance and complete interface matching independently of candidate selection. Require the kernel's closed-domain regularity certificate. Exact layer incidence and area checks apply to straight prismatic cases. Use independently defined material probes and volume checks for pierced 3D cases.
- Record block count, the applicable lower bound and their difference. For a conforming connected quadrilateral section use `Q = B/2 + I − χ`. Report a bound based on mandatory input vertices separately from one conditional on inserted host vertices. Include event layers and required pocket-floor traces when applicable. Leave the bound unspecified for a 3D complex without a justified counting argument.
- Sample the spectral condition number on the same closed 5 × 5 × 5 parameter grid in every block of both outputs. Summarize the per-solid maxima by median, 95th percentile and worst value. Compare common successes pairwise and report newly admitted or lost cases separately. Publish block-count improvements alongside any deterioration in map quality.

## Development controls

Use public synthetic inputs independent of the holdout. For polygonal sections, search separating lines from convex features and test their intersections with the current host. A separating line can meet adjacent host edges. Verify the resulting cell topology, convexity and full traces before admitting a quadrilateral ring. On blind pockets the floor must also admit a conforming cap, including boundary parity. Treat concave stock and holes as a separate extension.

For stacked bores, propagate a requested cut through mapped face coordinates. Carry axis permutations, reversals, knot parameters and feature provenance. Split a rational Bézier span in homogeneous coordinates. A trace entering a neighbor must be an admitted isoparametric cut there before exact tensor subdivision can be used. Close the cut set across all incident faces, deduplicate requests, enforce a declared block budget and certify both children and their interfaces. Record the first incompatible trace or budget refusal. Compare the final count with other certified plans before selection.


## Run on independent private sketches

The kernel evaluated here is frozen at `4e8a3944cd5aa85f8ba0984b109b6649b2d54b6f`. The next commit supplies the three bracket regression recipes and has identical Rust sources. Build the reference and candidate into separate directories. This leaves the active checkout and the private corpus intact.

```sh
git worktree add --detach ../sct-baseline e5bd5c6
git worktree add --detach ../sct-candidate 15b2626
cargo build --release --locked --offline --manifest-path ../sct-baseline/Cargo.toml
cargo build --release --locked --offline --manifest-path ../sct-candidate/Cargo.toml
mkdir -p build/holdout/baseline build/holdout/candidate
for version in baseline candidate; do
    cli="../sct-$version/target/release/solidcraft-trivariate"
    "$cli" build "$RECIPE" "build/holdout/$version/model.json" &&
    "$cli" audit "build/holdout/$version/model.json"
done
```

Set `RECIPE` to the path of one private recipe. Use a fresh output directory per recipe when collecting a corpus. Record unsuccessful builds and run the independent geometry checks before comparing block counts. The [public 3D controls](results/propagated-bores.md) include the corresponding STEP readback.

### Additional oblique-separation diagnostic

After the freeze, a separate through-hole diagnostic used a 100 × 100 × 12 mm square with diamond centers at (42, 42) and (58, 58) and half-diagonal 12 mm. The X and Y projections overlap. The two holes are strictly separated by the line x + y = 100. The two historical executables returned 72 blocks and exact volume 113088 mm³. The current construction uses ten blocks for the same exact volume. This example shows why finding a separating line alone does not settle the host shapes, corner correspondence or admissibility of the ring maps. The diagnostic has a separate denominator from the 45-case evaluation. [Reproducible recipe](../tests/regressions/overlapping_projections_diamonds.recipe.json)
