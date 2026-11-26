# Development instructions

From the repository root, read the [architecture contracts](docs/architecture.md) and [Luna roadmap](docs/LUNA.md). The [specification](docs/kernel-requirements.md) defines the scientific guarantees. Historical discussion notes preserve the context for the original choices.

When the code graph is available, prefer its index, symbol search, and call tracing. If those results are insufficient, search text and source files. Project commands run from a standard clone. On systems with RTK, follow its local command conventions as well.

The production kernel, builders, and audits are developed in Rust. Python provides independent oracles and historical references. Each geometric increment must provide Rust functionality exposed by the crate and tested with Cargo.

1. Divide work into reusable capabilities. Before implementation, record the starting revision, affected files, preconditions, thresholds, independent oracle, rejected inputs, and exit criteria.
2. Write documentation, report text, figure labels, gallery copy, docstrings, technical comments and commit messages in English. Use English for code identifiers. Preserve existing recipes and versioned formats. Each family extension must have an identifier and explicit conventions.
3. Keep the library independent of the user interface. The CLI and future adapters use the same builder. Model adoption goes through audit and atomic replacement. Test that rejection preserves the previous state.
4. Distinguish representation fidelity, local regularity, injectivity, trace compatibility, coverage, and analysis suitability. Sampling is diagnostic. A certificate names its domain and method. Incomplete proof yields an indeterminate result.
5. Develop oracles independently of the builder. Keep analytical references in `tests/references`, invalid inputs in `tests/counterexamples`, and each reproduced numerical defect in `tests/regressions` with its recipe and the assertion that fails on the old revision.
6. Run `make check` at integration. Report revision, recipe, environment, commands, measurements, coverage, and limitations. Rejection of a required case counts as a coverage failure. Keep regenerable results under `build/`.
7. `upstream/` contains read-only references. Before reusing code, document its origin, revision, license, required notices, and numerical contract. The original code's license choice before distribution belongs to the project owner (NF-8).
8. If construction conflicts with its oracle, preserve the counterexample and identify the assumption at fault before changing thresholds or scope.

Rebuild tracked geometry from its source recipe whenever the chosen partition changes. Run `make check-published-geometry` before publishing. Keep historical measurements linked to their original revision, and distribute the current geometry and its native audit together.
