"""Replay the operation gallery against the native CLI and Rust examples."""

import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

from check_gallery_geometry import step_content
from check_tall_booleans import STATES

ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / "docs/site"
RELEASE = ROOT / "target/release"
CLI = RELEASE / "solidcraft-trivariate"


def compare(data, path):
    if json.loads(data) != json.loads(path.read_bytes()):
        raise ValueError(f"Native replay differs from {path.relative_to(ROOT)}")


def check():
    states = ("sharp", "chamfer", "fillet")
    for state in states:
        original = "annulus_original_edges" if state == "fillet" else "annulus_four_block_chamfer"
        for example, prefix in [(original, "edges-four"), ("annulus_prepared_volumes", "edges")]:
            compare(subprocess.check_output([RELEASE / "examples" / example, state]), SITE / f"{prefix}-{state}.json")

    with tempfile.TemporaryDirectory() as directory:
        model = Path(directory) / "model.json"
        for family in ("reference_rectangle", "pierced_square"):
            for state in ("a", "b", "intersection", "union", "difference"):
                base = SITE / f"booleans-{family}-{state}"
                subprocess.run([CLI, "build", base.with_suffix(".recipe.json"), model], check=True, capture_output=True)
                compare(model.read_bytes(), base.with_suffix(".model.json"))
                compare(subprocess.check_output([CLI, "audit", model]), base.with_suffix(".audit.json"))

        box = SITE / "box-cylinder"
        box_audit = json.loads((box / "audit.json").read_bytes())
        for recipe, result in (("boolean", "result"), ("classic", "classic")):
            subprocess.run([CLI, "build", box / f"{recipe}.recipe.json", model],
                           check=True, capture_output=True)
            compare(model.read_bytes(), box / f"{result}.model.json")
            record = box_audit["two_operand_route"] if result == "result" else box_audit
            if (record["model_sha256"] != hashlib.sha256(model.read_bytes()).hexdigest()
                    or record["recipe_sha256"] != hashlib.sha256((box / f"{recipe}.recipe.json").read_bytes()).hexdigest()):
                raise ValueError(f"The box/cylinder audit identifies a different {result} model or recipe")
            if result == "result":
                candidate = Path(directory) / "result.step"
                subprocess.run([CLI, "export-step", model, candidate], check=True, capture_output=True)
                if step_content(candidate.read_bytes()) != step_content((box / "result.step").read_bytes()):
                    raise ValueError("Box minus cylinder STEP differs from its native replay")

        tall = box / "tall"
        for state in STATES:
            subprocess.run([CLI, "build", tall / f"{state}.recipe.json", model], check=True, capture_output=True)
            compare(model.read_bytes(), tall / f"{state}.model.json")
            compare(subprocess.check_output([CLI, "audit", model]), tall / f"{state}.audit.json")

        crossed = SITE / "crossed-tools"
        for rise in (-3, 0, 3):
            for op in ("intersection", "union", "a-minus-b", "b-minus-a"):
                subprocess.run([CLI, "crossed-boolean", crossed / f"rise-{rise}.recipe.json", op, model],
                               check=True, capture_output=True)
                compare(model.read_bytes(), crossed / f"{op}-{rise}.json")
    print("Native replay matched six edge models, ten axial models and audits, twelve crossed-tool models, the box/cylinder comparison with its STEP, and seven tall-cylinder states.")


if __name__ == "__main__":
    check()
