#!/usr/bin/env python3
"""Check every published gallery download against a fresh release-kernel build."""
import sys
# FreeCAD initializes its embedded main module on first import.
if __name__ == "__main__" and "--read-step" in sys.argv:
    import FreeCAD
    import Part

import argparse
import hashlib
import json
import re
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal
from pathlib import Path
from urllib.parse import urlsplit
from gallery_curved_patches import build_catalog as build_patch_catalog
from check_gallery_patch_maps import verify_patch_samples, verify_patch_topology

ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / "docs/site"
POLYGON_COMPARISON = ROOT / "docs/results/sketch-nurbs-comparison.json"
TIMESTAMP = re.compile(rb"(?m)^(FILE_NAME\('SolidCraft boundary',')\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ(')")

def load(data):
    return json.loads(data, parse_float=Decimal)

def first_difference(expected, actual, path="$"):
    if isinstance(expected, dict) and isinstance(actual, dict) and expected.keys() == actual.keys():
        key = next(key for key in expected if expected[key] != actual[key])
        return first_difference(expected[key], actual[key], f"{path}.{key}")
    if isinstance(expected, list) and isinstance(actual, list) and len(expected) == len(actual):
        index = next(i for i, item in enumerate(expected) if item != actual[i])
        return first_difference(expected[index], actual[index], f"{path}[{index}]")
    return f"{path}: expected {str(expected)[:160]}, rebuilt {str(actual)[:160]}"

def step_content(data):
    normalized, count = TIMESTAMP.subn(rb"\1UTC\2", data)
    if count != 1:
        raise ValueError("missing or ambiguous native STEP timestamp")
    return normalized

def asset(url):
    link = urlsplit(url)
    path = (SITE / link.path).resolve()
    if link.netloc or link.scheme or not path.is_relative_to(SITE / "data"):
        raise ValueError(f"invalid gallery asset path: {url}")
    expected = hashlib.sha256(path.read_bytes()).hexdigest()[:12]
    if link.query != f"v={expected}":
        raise ValueError(f"stale download hash: {url}")
    return path

def check(cli, entry, patch_catalog):
    links = {k: asset(entry[k]) for k in ("source", "model", "audit", "step")}
    if "representation" in entry:
        links["representation"] = asset(entry["representation"])
    with tempfile.TemporaryDirectory(prefix="sct-gallery-replay-") as tmp:
        model, step = Path(tmp) / "model.json", Path(tmp) / "solid.step"
        def run(*args):
            proc = subprocess.run([cli, *args], capture_output=True, timeout=180)
            if proc.returncode:
                raise ValueError((proc.stderr + proc.stdout).decode(errors="replace")[-1500:])
            return proc.stdout
        run("build", links["source"], model)
        rebuilt, published = load(model.read_bytes()), load(links["model"].read_bytes())
        if rebuilt != published:
            raise ValueError(f"model mismatch at {first_difference(published, rebuilt)}")
        native = load(run("audit", model))
        evidence = load(links["audit"].read_bytes())
        if native.get("status") != "valid" or native != evidence.get("model_audit", evidence):
            raise ValueError("published audit differs from native replay")
        metrics = {row["label"]: row["value"] for row in entry["metrics"]}
        if (len(metrics) != len(entry["metrics"])
                or int(metrics.get("Bézier cells", metrics.get("Rust model blocks", -1))) != native["blocks"]
                or int(metrics.get("Rust shared interfaces", -1)) != native["shared_faces"]):
            raise ValueError("displayed partition counts differ from the native audit")
        expected_patch = patch_catalog.get(entry["id"])
        if expected_patch is not None:
            if "representation" not in links:
                raise ValueError("exact patch map is absent from the gallery manifest")
            representation = load(links["representation"].read_bytes())
            expected_kind = expected_patch[1]
            row = expected_patch[0]
            if (representation.get("schema") != "gallery_exact_patch_map_v1"
                    or representation.get("case_id") != entry["id"]
                    or representation.get("representation_kind") != expected_kind
                    or representation.get("patches") != row["patches"]
                    or representation.get("interfaces") != row["interfaces"]
                    or representation.get("boundary_faces") != row["boundary_faces"]
                    or representation.get("source_model_sha256") != hashlib.sha256(links["model"].read_bytes()).hexdigest()
                    or representation.get("source_recipe_sha256") != hashlib.sha256(links["source"].read_bytes()).hexdigest()):
                raise ValueError("published exact patch map differs from a fresh extraction")
            verify_patch_samples(representation["patches"], published)
            verify_patch_topology(representation, published)
            recorded = evidence.get("exact_patch_representation", {})
            if (recorded.get("sha256") != hashlib.sha256(links["representation"].read_bytes()).hexdigest()
                    or recorded.get("path", "").split("?", 1)[0] != entry["representation"].split("?", 1)[0]):
                raise ValueError("audit does not identify the published exact patch map")
        elif "representation" in entry:
            raise ValueError("unexpected exact patch representation on an unqualified model")
        if entry["id"].startswith("bench_poly_"):
            comparison = json.loads(POLYGON_COMPARISON.read_bytes())
            row = next((item for item in comparison["cases"] if item["case"] == entry["id"]), None)
            if row is None:
                raise ValueError("polygon case is absent from the frozen NURBS comparison")
            expected_mode = min(row["modes"], key=lambda mode: (
                int(row["modes"][mode]["patches"]), int(row["modes"][mode]["cells"]), mode != "assemble"
            ))
            selection = evidence.get("partition_selection", {})
            if (selection.get("selected_mode") != expected_mode
                    or selection.get("native_model_sha256") != hashlib.sha256(links["model"].read_bytes()).hexdigest()
                    or metrics.get("NURBS partition choice") != expected_mode):
                raise ValueError("published NURBS partition choice does not match the frozen two-mode study")
            recipe = load(links["source"].read_bytes())
            if recipe.get("nurbs_partition") != expected_mode:
                raise ValueError("published polygon recipe does not select the audited NURBS partition")
            expected_patch_count = (len(load(links["representation"].read_bytes())["patches"])
                                    if "representation" in links else native.get("nurbs_patch_count"))
            if (metrics.get("NURBS volumes") != expected_patch_count
                    or ("representation" not in links
                        and len(published["nurbs_patches"]) != metrics["NURBS volumes"])):
                raise ValueError("displayed NURBS count differs from the downloaded model")
            if ("representation" in links
                    and len(published["nurbs_patches"]) != metrics.get("Native NURBS volumes before assembly")):
                raise ValueError("native NURBS count differs from the downloaded model")
            report = evidence["extrusion_quality"]
            for label, key in (("Sampled in-plane κ₂ max", "planar_spectral_kappa"),
                               ("Minimum corner sine", "minimum_corner_sine"),
                               ("Maximum in-plane span / normal thickness", "thickness_direction_aspect")):
                if str(metrics.get(label)) != format(float(report[key]), ".6g"):
                    raise ValueError(f"displayed {label} differs from the quality report")
        for key, field in (("model", "model_sha256"), ("step", "step_sha256"), ("source", "recipe_sha256")):
            if evidence.get("kernel_rebuild", {}).get(field) != hashlib.sha256(links[key].read_bytes()).hexdigest():
                raise ValueError(f"{key} digest differs from the published rebuild record")
        run("export-step", model, step)
        if step_content(step.read_bytes()) != step_content(links["step"].read_bytes()):
            raise ValueError("published STEP differs from current export")
    return links.values()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, default=ROOT / "target/release/solidcraft-trivariate")
    parser.add_argument("--read-step", action="store_true",
                        help="Also read each shipped STEP with FreeCAD/OCCT (requires FreeCAD Python modules)")
    args = parser.parse_args()
    entries = json.loads((SITE / "experiments.json").read_bytes())
    ids = [entry["id"] for entry in entries]
    if not ids or len(set(ids)) != len(ids):
        raise SystemExit("empty gallery or duplicate model IDs")
    referenced, failures = set(), []
    patch_report = build_patch_catalog(SITE / "experiments.json")
    patch_catalog = {row["id"]: (row, "curved") for row in patch_report["entries"]}
    patch_catalog.update({row["id"]: (row, "singleton_linear") for row in patch_report["linear_patch_assemblies"]})
    with ThreadPoolExecutor(max_workers=2) as pool:
        tasks = [(e["id"], pool.submit(check, args.cli.resolve(), e, patch_catalog)) for e in entries]
        for name, future in tasks:
            try:
                referenced.update(future.result())
            except (ValueError, OSError, subprocess.SubprocessError) as exc:
                failures.append(f"{name}: {exc}")
    if failures:
        raise SystemExit("\n".join(failures))
    orphans = sorted(p for p in (SITE / "data").rglob("*") if p.is_file() and p not in referenced)
    if orphans:
        raise SystemExit("Unreferenced gallery files: " + ", ".join(map(str, orphans)))
    if args.read_step:
        import FreeCAD  # Initialize the OCCT-backed Part module.
        import Part
        for entry in entries:
            shape = Part.read(str(asset(entry["step"])))
            if (len(shape.Solids) != 1 or not shape.Shells
                    or not shape.isValid() or not all(shell.isClosed() for shell in shape.Shells)):
                raise SystemExit(f"{entry['id']}: invalid or open STEP solid")
        print(f"{len(entries)} valid closed solids reread by FreeCAD {'.'.join(FreeCAD.Version()[:3])}.")
    print(f"{len(entries)} fresh builds, audits and STEP exports matched; {len(referenced)} downloads checked.")

if __name__ == "__main__":
    main()
