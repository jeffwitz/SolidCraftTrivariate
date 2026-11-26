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

ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / "docs/site"
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

def check(cli, entry):
    links = {k: asset(entry[k]) for k in ("source", "model", "audit", "step")}
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
                or int(metrics.get("Rust model blocks", -1)) != native["blocks"]
                or int(metrics.get("Rust shared interfaces", -1)) != native["shared_faces"]):
            raise ValueError("displayed partition counts differ from the native audit")
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
    with ThreadPoolExecutor(max_workers=2) as pool:
        tasks = [(e["id"], pool.submit(check, args.cli.resolve(), e)) for e in entries]
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
