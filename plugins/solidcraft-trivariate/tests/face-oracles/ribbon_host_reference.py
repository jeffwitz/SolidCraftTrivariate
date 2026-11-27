#!/usr/bin/env python3
"""Reproduce native SolidCraft 0.48.1 edge/face thickness references."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import tempfile


TESTS = Path(__file__).resolve().parents[1]
STOCK = TESTS / "face-box-add.sc"
BODY = "FaceIntegrationBox"
CASES = {
    "fillet-edge": (
        "solid.fillet",
        [{"body": BODY, "index": 8}],
        "radius",
    ),
    "fillet-face": (
        "solid.fillet",
        [{"body": BODY, "index": i} for i in range(4)],
        "radius",
    ),
    "chamfer-face": (
        "solid.chamfer",
        [{"body": BODY, "index": i} for i in range(4)],
        "distance",
    ),
}
EXPECTED = {
    "fillet-edge": 35999.83904862254,
    "fillet-face": 35997.055216642955,
    "chamfer-face": 35993.14583333333,
}


def main() -> None:
    cli = "solidcraft-cli"
    reports = {}
    with tempfile.TemporaryDirectory(prefix="solidcraft-ribbon-reference-") as temp:
        work = Path(temp)
        for label, (command, edges, parameter) in CASES.items():
            script = work / f"{label}.json"
            output = work / f"{label}.sc"
            measured = work / f"{label}.measurement.json"
            script.write_text(json.dumps({"commands": [{"command": command, "params": {
                "edges": edges,
                parameter: 0.25,
                "tangent_chain": False,
                "body": BODY,
                "name": f"Reference {label}",
            }}]}))
            subprocess.run([cli, "run", str(script), "--in", str(STOCK), "--save", str(output),
                            "--quiet"], check=True)
            subprocess.run([cli, "eval", str(output), "--quiet", "--json-out", str(measured)], check=True)
            report = json.loads(measured.read_text())
            body = report["bodies"][0]
            volume = report["total"]["volume_mm3"]
            if report["body_count"] != 1 or body["closed"] is not True:
                raise AssertionError(f"{label}: expected one closed body, got {report}")
            if abs(volume - EXPECTED[label]) > 1e-8:
                raise AssertionError(f"{label}: {volume} != reference {EXPECTED[label]}")
            reports[label] = {
                "volume_mm3": volume,
                "closed": body["closed"],
                "faces": body["faces"],
                "edges": body["edges"],
                "bbox": body["bbox"],
                "radius_or_distance_mm": 0.25,
                "tangent_chain": False,
            }
    print(json.dumps({"host": "SolidCraft 0.48.1", "stock": STOCK.name,
                      "stock_volume_mm3": 36000, "cases": reports}, indent=2))


if __name__ == "__main__":
    main()
