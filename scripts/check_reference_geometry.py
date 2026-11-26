#!/usr/bin/env python3
"""Rebuild non-gallery geometry deliverables and check their source parity."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import tempfile
import zipfile
from pathlib import Path

from check_gallery_geometry import step_content

ROOT = Path(__file__).resolve().parents[1]


def run(command: list[str], *, cwd: Path = ROOT) -> bytes:
    result = subprocess.run(command, cwd=cwd, check=True, stdout=subprocess.PIPE)
    return result.stdout


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_json(data: bytes) -> object:
    return json.loads(data)


def build(cli: Path, recipe: Path, output: Path) -> None:
    run([str(cli), "build", str(recipe), str(output)])


def audit(cli: Path, model: Path) -> dict[str, object]:
    return load_json(run([str(cli), "audit", str(model)]))


def export_step(cli: Path, model: Path, output: Path) -> None:
    run([str(cli), "export-step", str(model), str(output)])


def compare(label: str, expected: bytes, actual: bytes) -> None:
    if expected != actual:
        raise RuntimeError(
            f"{label} differs: expected sha256={sha256(expected)}, "
            f"actual sha256={sha256(actual)}"
        )


def compare_step(label: str, expected: bytes, actual: bytes) -> None:
    compare(label, step_content(expected), step_content(actual))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cli", type=Path, default=ROOT / "target/release/solidcraft-trivariate")
    parser.add_argument("--cargo", default="cargo")
    args = parser.parse_args()
    cli = args.cli.resolve()
    if not cli.is_file():
        raise SystemExit(f"CLI binary not found: {cli}")

    checks = 0
    with tempfile.TemporaryDirectory(prefix="solidcraft-reference-geometry-") as temp:
        work = Path(temp)

        # Current cross-bore and pocket deliverables are rebuilt from their source recipes.
        direct_cases = [
            (
                "g2-cross",
                ROOT / "tests/references/g2_cross_sketches.recipe.json",
                ROOT / "docs/results/g2-cross-plate.model.json",
                ROOT / "docs/results/g2-cross-trivariate.step",
            ),
            (
                "s3-face-pocket",
                ROOT / "tests/references/s3_face_pocket.recipe.json",
                ROOT / "docs/results/s3-face-pocket.model.json",
                ROOT / "docs/results/s3-face-pocket.step",
            ),
        ]
        for name, recipe, saved_model, saved_step in direct_cases:
            model, step = work / f"{name}.model.json", work / f"{name}.step"
            build(cli, recipe, model)
            result = audit(cli, model)
            if result.get("status") != "valid":
                raise RuntimeError(f"{name} audit did not pass: {result}")
            compare(f"{name} model", saved_model.read_bytes(), model.read_bytes())
            export_step(cli, model, step)
            compare_step(f"{name} STEP", saved_step.read_bytes(), step.read_bytes())
            checks += 2

        # Preserve and distinguish the independently exported SolidCore B-rep input.
        reference_step = ROOT / "docs/results/g2-cross-solidcore.step"
        if not reference_step.is_file() or not direct_cases[0][3].is_file():
            raise RuntimeError("independent SolidCore STEP reference is missing")
        if reference_step.read_bytes() == direct_cases[0][3].read_bytes():
            raise RuntimeError("independent SolidCore and trivariate STEP unexpectedly coincide byte-for-byte")
        checks += 1

        open_recipe = ROOT / "tests/references/s3_face_pocket_open_side.recipe.json"
        open_model = ROOT / "docs/results/s3-face-pocket-open-side.model.json"
        rebuilt = work / "s3-open-side.model.json"
        build(cli, open_recipe, rebuilt)
        if audit(cli, rebuilt).get("status") != "valid":
            raise RuntimeError("S3 open-side model audit did not pass")
        compare("S3 open-side model", open_model.read_bytes(), rebuilt.read_bytes())
        checks += 1

        # Rebuild all five compact FieldCraft variants from gallery source recipes.
        compact_path = ROOT / "docs/results/fieldcraft-replay/fieldcraft-compact.zip"
        cases = {
            "x40-y80": "fieldcraft_x40_y80",
            "x50-y74": "fieldcraft_x50_y74",
            "x51-y74": "fieldcraft_x51_y74",
            "x64-y64": "fieldcraft_x64_y64",
            "x80-y40": "fieldcraft_x80_y40",
        }
        # Every downloadable solid in results must have a replay above or below.
        results = ROOT / "docs/results"
        covered = {p for _, _, model, step in direct_cases for p in (model, step)}
        covered.update((compact_path, open_model, reference_step))
        rim_sources = [
            ("s3-stock", "edit-s3-stock-rim", work / "s3-face-pocket.model.json"),
            ("fieldcraft", "edit-fieldcraft-bore-rims", ROOT / "docs/site/data/fieldcraft_x64_y64/fieldcraft_x64_y64.model.json"),
        ]
        for stem, command, source in rim_sources:
            for state in ("chamfer", "fillet"):
                filename = f"{stem}-{state}.step"
                saved = results / "fixed-rims" / filename
                regenerated = work / filename
                run([str(cli), command, str(source), state, str(regenerated)])
                compare_step(filename, saved.read_bytes(), regenerated.read_bytes())
                covered.add(saved)
                checks += 1
        shipped = {p for p in results.rglob("*") if p.is_file()
                   and p.name.endswith((".model.json", ".step", ".zip"))}
        if covered != shipped:
            raise RuntimeError(f"Uncovered or missing geometry in results: {covered ^ shipped}")

        with zipfile.ZipFile(compact_path) as archive:
            geometries = {f"{case}/{case}.{ext}" for case in cases for ext in ("step", "model.json")}
            geometries.update(("plate/four-sector-plate.model.json", "plate/four-sector-plate.step"))
            members = [name for name in archive.namelist() if name.endswith((".step", ".model.json"))]
            if sorted(members) != sorted(geometries):
                raise RuntimeError("Missing, duplicate or unverified solids inside the FieldCraft bundle")
            for case, site_id in cases.items():
                stem = case
                recipe_bytes = (ROOT / f"docs/site/data/{site_id}/{site_id}.recipe.json").read_bytes()
                compare(f"{case} compact recipe", recipe_bytes, archive.read(f"{case}/{stem}.recipe.json"))
                recipe = work / f"{case}.recipe.json"
                recipe.write_bytes(recipe_bytes)
                model, step = work / f"{case}.model.json", work / f"{case}.step"
                build(cli, recipe, model)
                result = audit(cli, model)
                if result.get("status") != "valid":
                    raise RuntimeError(f"{case} audit did not pass: {result}")
                expected_model = archive.read(f"{case}/{stem}.model.json")
                if load_json(expected_model) != load_json(model.read_bytes()):
                    raise RuntimeError(f"{case} compact model geometry differs from current build")
                expected_audit = load_json(archive.read(f"{case}/{stem}.audit.json"))
                if expected_audit != result:
                    raise RuntimeError(f"{case} compact audit differs from current build")
                export_step(cli, model, step)
                compare_step(f"{case} compact STEP", archive.read(f"{case}/{stem}.step"), step.read_bytes())
                checks += 4

            # The control-prism recipe is also checked against the canonical test recipe.
            plate_recipe = ROOT / "tests/references/four_sector_plate.recipe.json"
            compare("control-prism recipe", plate_recipe.read_bytes(), archive.read("plate/four_sector_plate.recipe.json"))
            model, step = work / "control-prism.model.json", work / "control-prism.step"
            build(cli, plate_recipe, model)
            plate_audit = audit(cli, model)
            if plate_audit.get("status") != "valid":
                raise RuntimeError(f"control-prism audit did not pass: {plate_audit}")
            compare("control-prism model", archive.read("plate/four-sector-plate.model.json"), model.read_bytes())
            control_audit_name = "plate/four-sector-plate.audit.json"
            if control_audit_name not in archive.namelist():
                raise RuntimeError("compact bundle is missing the control-prism audit")
            if load_json(archive.read(control_audit_name)) != plate_audit:
                raise RuntimeError("control-prism audit differs from current build")
            export_step(cli, model, step)
            compare_step("control-prism STEP", archive.read("plate/four-sector-plate.step"), step.read_bytes())
            checks += 3

        # The annulus prototype files are complete model nets produced by Rust examples.
        annulus = ROOT / "docs/results/annulus-fixed-topology"
        commands = [
            ("annulus_four_block_chamfer", "sharp", "four-sharp.json", 4),
            ("annulus_four_block_chamfer", "chamfer", "four-chamfer.json", 4),
            ("annulus_prepared_volumes", "sharp 2", "ten-sharp.json", 10),
            ("annulus_prepared_volumes", "chamfer 2", "ten-chamfer.json", 10),
            ("annulus_prepared_volumes", "fillet 2", "ten-fillet.json", 10),
        ]
        for example, mode, filename, expected_blocks in commands:
            output = run(
                [args.cargo, "run", "--locked", "--offline", "--quiet", "--example", example, "--", *mode.split()]
            )
            saved = (annulus / filename).read_bytes()
            compare(f"{filename} Rust example", saved, output)
            model = load_json(output)
            blocks = model.get("blocks", model.get("charts", []))
            if len(blocks) != expected_blocks or any(not b.get("xyzw") for b in blocks):
                raise RuntimeError(f"{filename} is not a complete {expected_blocks}-map model")
            checks += 1

    print(f"geometry parity passed: {checks} checks; STEP comparison ignores only the header timestamp")
    print("annulus four-/ten-map prototypes are full Rust-example models; they have no CLI audit route")
    print(f"independent SolidCore reference retained: {reference_step.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
