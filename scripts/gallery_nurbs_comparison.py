"""Pack the two audited partitions of each polygon sketch for the WebGL comparison."""
import json
from pathlib import Path

from gallery_catalog import version_modules

ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / "docs/site"
WORK = ROOT / "build/sketch-nurbs-evaluation"


def pack(model):
    groups = {cell: index for index, patch in enumerate(model["nurbs_patches"]) for cell in patch["cells"]}
    points, blocks, indices = {}, [], {}
    for block in model["blocks"]:
        assert block["degrees"] == [1, 1, 1]
        indices[block["id"]] = len(blocks)
        controls = [points.setdefault(tuple(p), len(points)) for p in block["controls"]]
        blocks.append([groups[block["id"]], *controls])
    interfaces = []
    for pair in model["topology"]["interfaces"]:
        (a, side), (b, other) = pair["first"], pair["second"]
        interfaces.append([indices[a], side, indices[b], other])
    return {"points": list(points), "blocks": blocks, "interfaces": interfaces}


def main():
    report = json.loads((ROOT / "docs/results/sketch-nurbs-comparison.json").read_text())
    entries = {e["id"]: e for e in json.loads((SITE / "experiments.json").read_text())}
    cases = []
    for row in report["cases"]:
        key = row["case"]
        geometry = {}
        for mode, metrics in row["modes"].items():
            path = WORK / f"{key}.{mode}.model.json"
            model = json.loads(path.read_bytes())
            assert len(model["nurbs_patches"]) == metrics["patches"]
            assert len(model["blocks"]) == metrics["cells"]
            geometry[mode] = pack(model)
        if geometry["optimize"] == geometry["assemble"]:
            del geometry["optimize"]
        cases.append({"id": key, "title": entries[key]["title"], "modes": row["modes"], "geometry": geometry})
    assert len(cases) == 46
    (SITE / "comparison-data.json").write_text(json.dumps(cases, separators=(",", ":"))+"\n")
    version_modules("compare.js", ["renderer.js", "comparison-data.json"])
    version_modules("compare.html", ["compare.js", "style.css", "compare.css"])


if __name__ == "__main__":
    main()
