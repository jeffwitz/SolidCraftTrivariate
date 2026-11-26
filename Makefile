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
check-published-geometry:
	cargo build --locked --offline --release
	cargo build --locked --offline --release --example annulus_four_block_chamfer --example annulus_original_edges --example annulus_prepared_volumes
	python3 scripts/check_operation_geometry.py
	python3 scripts/build_box_cylinder_tall_patches.py --check
	python3 scripts/check_tall_booleans.py
	python3 scripts/check_gallery_geometry.py
	python3 scripts/check_reference_geometry.py --cli target/release/solidcraft-trivariate
	python3 tests/references/replay_geometry_recipes.py --cli target/release/solidcraft-trivariate --work build/reference-recipe-replay --output build/reference-recipe-replay.json
check: test demo check-published-geometry
check-strict: refs check
