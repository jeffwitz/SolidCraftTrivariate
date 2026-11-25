.PHONY: bootstrap refs test demo demo-native check
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
	cargo run --quiet --locked --offline -- build tests/references/filleted_bracket.json build/filleted_bracket.json
	python3 -m trivariate audit build/filleted_bracket.json
demo: demo-native
check: refs test demo
