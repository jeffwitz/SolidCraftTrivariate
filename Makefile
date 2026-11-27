.PHONY: bootstrap refs test demo demo-native check check-strict check-published-geometry
bootstrap:
	bash scripts/references.sh clone
refs:
	bash scripts/references.sh check
test:
	cargo test --locked --offline
	python3 -m unittest discover -s tests -v
demo-native:
	mkdir -p build
	cargo run --quiet --locked --offline -- plan-sketch tests/references/g2_cross_sketches.recipe.json > build/g2_cross_sketch_plan.json
	cargo run --quiet --locked --offline -- import-solidcraft-sketch tests/references/solidcraft_sketch_bracket.timeline.json build/solidcraft_sketch_bracket.json
	cargo run --quiet --locked --offline -- audit build/solidcraft_sketch_bracket.json
	cargo run --quiet --locked --offline -- export-step build/solidcraft_sketch_bracket.json build/solidcraft_sketch_bracket.step
	cargo run --quiet --locked --offline -- build tests/references/sketch_boolean_blind_pocket.recipe.json build/sketch_boolean_blind_pocket.json
	cargo run --quiet --locked --offline -- audit build/sketch_boolean_blind_pocket.json
	cargo run --quiet --locked --offline -- export-step build/sketch_boolean_blind_pocket.json build/sketch_boolean_blind_pocket.step
demo:
	cargo run --quiet --locked --offline -- build tests/references/affine_pair.json build/affine_pair.json
	python3 -m trivariate audit build/affine_pair.json
	cargo run --quiet --locked --offline -- build tests/references/annular_pair.json build/annular_pair.json
	python3 -m trivariate audit build/annular_pair.json
demo: demo-native
.PHONY: check-axial-knot-reduction
check-axial-knot-reduction:
	cargo build --locked --offline --release
	python3 scripts/check_axial_knot_reduction.py

check-published-geometry: check-axial-knot-reduction
	cargo build --locked --offline --release --example annulus_four_block_chamfer --example annulus_original_edges --example annulus_prepared_volumes --example cylinder_tensor_covers
	python3 scripts/check_operation_geometry.py
	python3 scripts/build_box_cylinder_tall_patches.py --check
	python3 scripts/check_tall_booleans.py
	python3 scripts/check_gallery_geometry.py
	python3 scripts/check_reference_geometry.py --cli target/release/solidcraft-trivariate
	python3 tests/references/replay_geometry_recipes.py --cli target/release/solidcraft-trivariate --work build/reference-recipe-replay --output build/reference-recipe-replay.json
check: test demo check-published-geometry
check-strict: refs check

.PHONY: check-cylinder-iga
check-cylinder-iga: check-axial-knot-reduction
	python3 tests/references/cylinder_torsion.py
	python3 tests/references/tall_cylinder_iga.py --check
	python3 tests/references/axial_knot_torsion.py --check

.PHONY: check-crossed-bending
check-crossed-bending:
	cargo test --locked --offline --release --test crossed_bending --test conforming_bending
	cargo build --locked --offline --release
	python3 scripts/check_native_crossed_bending.py target/release/solidcraft-trivariate
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 docs/results/crossed-interface-iga/check_bending.py target/release/solidcraft-trivariate

.PHONY: check-inherited-space check-inherited-corpus
check-inherited-space:
	cargo build --release --locked --offline
	python3 scripts/check_spline_parent.py --cli target/release/solidcraft-trivariate
	python3 scripts/check_spline_parent.py --cli target/release/solidcraft-trivariate --nonuniform-source
	python3 scripts/check_spline_parent.py --cli target/release/solidcraft-trivariate --reverse-source
	python3 scripts/check_spline_eligibility.py --cli target/release/solidcraft-trivariate --check-native-report
	python3 scripts/check_parent_field.py
	python3 scripts/check_parent_spectrum.py
	python3 scripts/check_inherited_joins.py
	python3 scripts/measure_inherited_fields.py

check-inherited-corpus: check-inherited-space
	python3 scripts/census_support_interfaces.py
	python3 scripts/census_transverse_transitions.py
	python3 scripts/index_support_roles.py
	python3 scripts/check_named_supports.py
	python3 scripts/check_pair_parent_polygons.py
	cmp build/pair-parent-expanded/summary.json docs/results/support-census/pair-parent-expanded.json
	python3 scripts/check_curved_supports.py
	python3 scripts/check_s3_spline_job.py
	python3 scripts/audit_fieldcraft_parent_support.py --cli target/release/solidcraft-trivariate --output build/support-census/fieldcraft-parent-support.json
	python3 scripts/check_pair_parent_corpus.py --cli target/release/solidcraft-trivariate
	cmp build/support-census/role-index.json docs/results/support-census/role-index.json
	python3 scripts/check_polygon_join_mechanics.py

.PHONY: check-local-thb
check-local-thb:
	cargo build --release --locked --offline
	cargo test --release --locked --offline --test thb_strip --test fieldcraft_nurbs_h1
	cargo test --release --locked --offline --test cube_convergence cube_convergence_replay -- --ignored
	python3 scripts/check_thb_mechanics.py
	target/release/solidcraft-trivariate export-thb-extraction build/thb-mechanics/package.json build/thb-mechanics/extraction.json
	python3 scripts/check_thb_function_values.py build/thb-mechanics/extraction.json --model build/thb-mechanics/model.json > build/thb-mechanics/function-values.json
	python3 scripts/check_thb_traction.py

.PHONY: docs docs-inputs
SPHINXBUILD ?= sphinx-build
docs-inputs:
	cargo build --locked --offline --bin solidcraft-trivariate
	RUSTDOCFLAGS="-D rustdoc::broken_intra_doc_links" cargo doc --locked --offline --no-deps --lib --target-dir build/docs-cargo
	python3 scripts/prepare_manual.py
docs: docs-inputs
	$(SPHINXBUILD) -E -a -W --keep-going -b html docs/rtd build/manual
	python3 scripts/check_manual.py

.PHONY: check-continuity-followup
check-continuity-followup:
	cargo test --release --locked --offline --test transverse_transitions
	python3 -m unittest discover -s tests/references -p test_pocket_lateral_c1_sweep.py
	cargo test --release --locked --offline --test cube_convergence matched_dimension_gaussian -- --ignored --exact
	cargo test --release --locked --offline --test pocket_convergence -- --include-ignored
	python3 scripts/check_continuity_followup.py
	cargo test --release --locked --offline --test thb_reentrant_edge -- --include-ignored
	cargo test --release --locked --offline --test reentrant_ritz -- --ignored
	OPENBLAS_NUM_THREADS=1 python3 scripts/check_reentrant_thb.py
