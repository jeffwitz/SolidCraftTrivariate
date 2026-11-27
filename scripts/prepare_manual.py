"""Generate the Sphinx catalogue from the local Rust executable and its fixtures."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / "target/debug/solidcraft-trivariate"
SOURCES = ROOT / "docs/rtd"
GENERATED = SOURCES / "_generated"


def write_examples(families, replay):
    catalog = json.loads((SOURCES / "recipe-examples.json").read_bytes())
    if set(catalog) != set(families):
        raise ValueError("Every native build family must have one downloadable recipe")
    output = GENERATED / "recipes"
    output.mkdir(parents=True, exist_ok=True)
    example_fields = {}
    with tempfile.TemporaryDirectory() as temporary:
        model = Path(temporary) / "model.json"
        for family, entry in catalog.items():
            if "example" in entry:
                recipe = entry["example"]
            else:
                recipe = json.loads((ROOT / entry["path"]).read_bytes())
                for key in entry["pointer"]:
                    recipe = recipe[key]
            if recipe["family"] != family:
                raise ValueError(f"Fixture no longer describes {family}")
            example_fields[family] = sorted(recipe.keys() - {"family"})
            path = output / f"{family}.json"
            path.write_text(json.dumps(recipe, ensure_ascii=False, indent=2) + "\n")
            if replay:
                subprocess.run([CLI, "build", path, model], stdout=subprocess.DEVNULL, check=True)
                audit = json.loads(subprocess.check_output([CLI, "audit", model]))
                if audit.get("status") != "valid":
                    raise ValueError(f"Example audit failed: {family}")

    return example_fields


def catalogue(replay):
    capabilities = json.loads(subprocess.check_output([CLI, "capabilities"]))
    message = subprocess.check_output([CLI, "--help"], text=True).strip()
    commands = message.removeprefix("usage: solidcraft-trivariate ").split(" | ")
    synopses = {command.split()[0]: command for command in commands}
    if set(synopses) != set(capabilities["commands"]):
        raise ValueError("Native help and capabilities disagree on the command inventory")
    purposes = json.loads((SOURCES / "reference/command-purposes.json").read_bytes())
    if purposes.keys() != synopses.keys():
        raise ValueError("Each command needs a reference entry: "
                         + str(purposes.keys() ^ synopses.keys()))
    lines = ["Command reference", "-----------------", ""]
    for name in sorted(synopses):
        lines.extend([f"``{synopses[name]}``", f"   {purposes[name]}", ""])
    lines += ["", "Recipe families", "---------------", "",
              "Download a recipe and pass it directly to ``build RECIPE MODEL``. "
              "Run ``audit MODEL`` to verify the resulting control nets and incidence. "
              "The inputs listed below are the top-level members supplied by that example. "
              "A built-in fixture can be constructed using only its family name.", "",
              ".. list-table::", "   :header-rows: 1", "", "   * - Family", "     - Example inputs", "     - Download"]
    fields = write_examples(capabilities["families"], replay)
    for family in sorted(capabilities["families"]):
        ref = f":download:`Recipe JSON <../_generated/recipes/{family}.json>`"
        inputs = ", ".join(f"``{key}``" for key in fields[family]) or "Built-in fixture"
        lines.extend([f"   * - ``{family}``", f"     - {inputs}", f"     - {ref}"])
    (GENERATED / "catalog.inc").write_text("\n".join(lines) + "\n")
    print(f"Catalogue: {len(synopses)} commands and {len(capabilities['families'])} families.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replay", action="store_true", help="build and audit every downloadable recipe")
    catalogue(parser.parse_args().replay)
    api = ROOT / "build/docs-extra/api"
    if api.exists():
        shutil.rmtree(api)
    shutil.copytree(ROOT / "build/docs-cargo/doc", api)
