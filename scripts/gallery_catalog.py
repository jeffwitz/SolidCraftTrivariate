#!/usr/bin/env python3
"""Merge benchmark bundles into the WebGL gallery and version its static assets."""
import argparse
import hashlib
import json
import re
from pathlib import Path

SITE = Path(__file__).resolve().parents[1] / "docs" / "site"
REPORTS = "https://github.com/jeffwitz/SolidCraftTrivariate/blob/main/docs/"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()[:12]


def version_modules(filename, dependencies):
    path = SITE / filename
    source = path.read_text()
    for dependency in dependencies:
        url = "./" + dependency
        source, hits = re.subn(re.escape(url) + r"\?v=[a-f0-9]+",
                               url + "?v=" + digest(SITE / dependency), source)
        if hits != 1:
            raise ValueError(f"Missing or ambiguous dependency {dependency} in {filename}")
    path.write_text(source)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("polygons", type=Path)
    parser.add_argument("brackets", type=Path)
    args = parser.parse_args()
    manifest = SITE / "experiments.json"
    existing = [e for e in json.loads(manifest.read_text())
                if not e["id"].startswith(("bench_poly_", "bench-bracket-"))]
    polygons = json.loads(args.polygons.read_text())
    brackets = json.loads(args.brackets.read_text())
    brackets.sort(key=lambda e: e["id"] != "bench-bracket-east-stack")
    for collection, report in [
        (brackets, "results/propagated-bores.md"),
        (polygons, "results/sketch-planner-polygon-evaluation/README.md"),
    ]:
        for entry in collection:
            entry["report"] = REPORTS + report
            if entry["id"] == "bench_poly_projection_diamonds":
                entry["report"] = REPORTS + "sketch-planner-blind-evaluation.md#additional-oblique-separation-diagnostic"
            for key in ("model", "source", "step", "audit"):
                path = entry[key].split("?")[0]
                entry[key] = path + "?v=" + digest(SITE / path)
    entries = brackets + polygons + existing
    if len({e["id"] for e in entries}) != len(entries):
        raise ValueError("Duplicate gallery identifiers")
    manifest.write_text(json.dumps(entries, indent=2, ensure_ascii=False) + "\n")
    version_modules("app.js", ["renderer.js"])
    version_modules("index.html", ["app.js", "experiments.json", "style.css"])


if __name__ == "__main__":
    main()
