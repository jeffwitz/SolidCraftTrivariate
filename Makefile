.PHONY: bootstrap refs test demo check
bootstrap:
	bash scripts/references.sh clone
refs:
	bash scripts/references.sh check
test:
	python3 -m unittest discover -s tests -v
demo:
	python3 -m trivariate build tests/references/affine_pair.json build/affine_pair.json
	python3 -m trivariate audit build/affine_pair.json
check: refs test demo
