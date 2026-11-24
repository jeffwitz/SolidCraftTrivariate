.PHONY: bootstrap refs test demo check
bootstrap:
	bash scripts/references.sh clone
refs:
	bash scripts/references.sh check
test:
	cargo test --locked --offline
	python3 -m unittest discover -s tests -v
demo:
	cargo run --quiet --locked --offline -- build tests/references/affine_pair.json build/affine_pair.json
	python3 -m trivariate audit build/affine_pair.json
	cargo run --quiet --locked --offline -- build tests/references/annular_pair.json build/annular_pair.json
	python3 -m trivariate audit build/annular_pair.json
	cargo run --quiet --locked --offline -- build tests/references/filleted_bracket.json build/filleted_bracket.json
	python3 -m trivariate audit build/filleted_bracket.json
check: refs test demo
