"""Replay the operation gallery against the native CLI and Rust examples."""

import json
from pathlib import Path
import subprocess
import tempfile

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

        crossed = SITE / "crossed-tools"
        for rise in (-3, 0, 3):
            for op in ("intersection", "union", "a-minus-b", "b-minus-a"):
                subprocess.run([CLI, "crossed-boolean", crossed / f"rise-{rise}.recipe.json", op, model],
                               check=True, capture_output=True)
                compare(model.read_bytes(), crossed / f"{op}-{rise}.json")
    print("Native replay matched six edge models, ten axial models and audits, and twelve crossed-tool models.")


if __name__ == "__main__":
    check()
