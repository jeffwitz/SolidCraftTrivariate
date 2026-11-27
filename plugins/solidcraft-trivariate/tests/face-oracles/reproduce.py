#!/usr/bin/env python3
"""Build the four host-recorded face-spline designs and plugin test cases.

Run from any directory with `python3 .../face-oracles/reproduce.py`. The only
non-generated input is the checked-in 0.48.1 FieldCraft bracket design.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile


TESTS = Path(__file__).resolve().parents[1]
BRACKET = TESTS / "ribbon-bracket.sc"
BRACKET_SHA256 = "d357f538f963fce76739f95b1e285ea2e87b00aff608b02f58d5e5267837d968"
GRANTS = {
    "folders": [{"name": "companion", "read": True, "write": True, "from": "companion"}],
    "command_groups": ["file"],
}


def call(cli: str, args: list[str]) -> str:
    proc = subprocess.run([cli, *args], capture_output=True, text=True, check=False)
    if proc.returncode:
        raise RuntimeError(
            f"{cli} {' '.join(args)} failed ({proc.returncode})\n"
            f"{proc.stdout}\n{proc.stderr}"
        )
    return proc.stdout


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value):
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def save_design(cli: str, commands: list[dict], name: str, work: Path, base: Path | None = None):
    script = work / f"{name}.commands.json"
    design = work / f"{name}.sc"
    result = work / f"{name}.results.json"
    write_json(script, {"commands": commands})
    args = ["run", str(script)]
    if base is not None:
        args += ["--in", str(base)]
    args += ["--save", str(design), "--quiet", "--results", "--json-out", str(result)]
    call(cli, args)
    return design, read_json(result)


def model_faces(cli: str, body: str, design: Path, work: Path, label: str):
    _, result = save_design(
        cli,
        [{"command": "model.faces", "params": {"body": body}}],
        label,
        work,
        design,
    )
    assert len(result) == 1 and result[0]["body"] == body, result
    upward = [
        face
        for face in result[0]["faces"]
        if face["exact"]
        and face["plane_normal"] is not None
        and face["plane_normal"][2] > 0.999999
        and face["point"] is not None
    ]
    if not upward:
        raise AssertionError(f"No exact +Z planar face on {body}")
    return max(upward, key=lambda face: face["area"])


def model_y_face(cli: str, body: str, design: Path, work: Path, label: str):
    _, result = save_design(
        cli,
        [{"command": "model.faces", "params": {"body": body}}],
        label,
        work,
        design,
    )
    faces = result[0]["faces"]
    candidates = [
        face
        for face in faces
        if face["exact"]
        and face["plane_normal"] is not None
        and face["plane_normal"][1] > 0.999999
        and face["point"] is not None
    ]
    if not candidates:
        raise AssertionError(f"No exact +Y planar XZ face on {body}")
    return min(candidates, key=lambda face: face["area"])


def evaluate(cli: str, design: Path, work: Path, label: str):
    path = work / f"{label}.eval.json"
    call(cli, ["eval", str(design), "--quiet", "--json-out", str(path)])
    report = read_json(path)
    assert report["body_count"] == 1, report
    body = report["bodies"][0]
    assert body["closed"] is True, body
    return body


def base_box(cli: str, work: Path):
    script = read_json(TESTS / "face-oracles/box-source.json")
    return save_design(cli, script["commands"], "face-box-base", work)[0]


def run_case(cli: str, case: str, operation: str, body_name: str, source: Path, work: Path):
    before = evaluate(cli, source, work, f"{case}-{operation}-before")
    face = model_faces(cli, body_name, source, work, f"{case}-{operation}-faces")
    if case == "bracket":
        assert face["index"] == 8 and math.isclose(face["area"], 368.0, abs_tol=1e-8), face
    else:
        assert face["index"] == 0 and math.isclose(face["area"], 3000.0, abs_tol=1e-8), face

    fixture = read_json(TESTS / f"face-oracles/{case}-{operation}-sketch.json")
    create = fixture["commands"][0]
    create["params"]["plane"] = {
        "face": face["point"],
        "body": body_name,
        "index": face["index"],
    }
    design, _ = save_design(
        cli, fixture["commands"], f"face-{case}-{operation}", work, source
    )
    design_name = f"face-{case}-{operation}.sc"
    target = TESTS / design_name
    target.write_bytes(design.read_bytes())

    # A second fixture contains the host plugin's own operation/depth fields.
    plugin_case = {
        "entry": "extrude",
        "params": {
            "source_feature_id": 3 if case.startswith("box") else 2,
            "height_mm": 2,
            "operation": operation,
        },
        "design": design_name,
        "grants": GRANTS,
        "expect": {
            "files": [
                f"companion/face-{3 if case.startswith('box') else 2}-{operation}-{5 if case.startswith('box') else 3}."
                f"{'result' if case.startswith('box') else 'tool'}.model.json"
            ]
        },
    }
    host_command = {
        "command": "solid.extrude",
        "params": {
            "sketch": "FaceIntegrationSpline",
            "distance": 2,
            "direction": "negative" if operation == "cut" else "positive",
            "operation": "cut" if operation == "cut" else "join",
            "name": "FaceHostReference",
        },
    }
    host_ref, _ = save_design(
        cli, [host_command], f"{case}-{operation}-reference", work, design
    )
    after = evaluate(cli, host_ref, work, f"{case}-{operation}-host-reference")
    if operation == "cut":
        assert after["volume_mm3"] < before["volume_mm3"]
    else:
        assert after["volume_mm3"] > before["volume_mm3"]
        assert after["bbox"]["max"][2] > before["bbox"]["max"][2]
    plugin_case["expect"] = {
        "files": plugin_case["expect"]["files"],
        "model": {
            "body_count": 1,
            "total": {"volume_mm3": after["volume_mm3"]},
        }
    }
    write_json(TESTS / f"face-{case}-{operation}.json", plugin_case)

    history_path = work / f"{case}-{operation}.history.json"
    call(cli, ["script", str(host_ref), "--out", str(history_path)])
    history = read_json(history_path)
    features = [
        item["params"]["feature"]
        for item in history["commands"]
        if item.get("command") == "timeline.add"
    ]
    sketch = next(feature for feature in features if feature.get("type") == "sketch" and feature.get("name") == "FaceIntegrationSpline")
    extrude = next(feature for feature in features if feature.get("type") == "extrude" and feature.get("name") == "FaceHostReference")
    curve = next(curve for curve in sketch["sketch"]["curves"] if curve["type"] == "spline")
    assert sketch["plane"]["type"] == "face"
    assert sketch["plane"]["anchor"]["target"]["body"] == body_name
    assert curve["periodic"] is True
    assert bool(curve.get("control", False)) is fixture["commands"][1]["command"].startswith(
        "sketch.spline.control_point"
    )
    expected_direction = "negative" if operation == "cut" else "positive"
    assert extrude["extent"] == {"convention": "total_distance", "direction": expected_direction, "distance": "2"}
    assert extrude["operation"] == ("cut" if operation == "cut" else "join")

    # Replay the exact recorded feature graph and require byte-level numeric agreement.
    replay = work / f"{case}-{operation}-replay.sc"
    call(cli, ["run", str(history_path), "--save", str(replay), "--quiet"])
    replayed = evaluate(cli, replay, work, f"{case}-{operation}-replay")
    assert replayed["volume_mm3"] == after["volume_mm3"]
    assert replayed["bbox"] == after["bbox"]
    return {
        "target": case,
        "operation": operation,
        "face_index": face["index"],
        "face_area_mm2": face["area"],
        "face_normal": face["plane_normal"],
        "face_point_mm": face["point"],
        "spline": "control" if curve.get("control") else "fit",
        "periodic": curve["periodic"],
        "extrusion_direction": extrude["extent"]["direction"],
        "extrusion_operation": extrude["operation"],
        "volume_before_mm3": before["volume_mm3"],
        "volume_after_mm3": after["volume_mm3"],
        "volume_delta_mm3": after["volume_mm3"] - before["volume_mm3"],
        "bounds_after_mm": after["bbox"],
        "replay_volume_identical": True,
    }


def qualify_xz_and_open_refusal(cli: str, bracket: Path, box: Path, work: Path):
    body_name = "SolidCraft boundary"
    face = model_y_face(cli, body_name, bracket, work, "bracket-xz-faces")
    assert face["index"] == 23 and math.isclose(face["area"], 320.0, abs_tol=1e-8), face
    local_points = [[-6, 20], [-4, 18], [-2, 20], [-4, 22]]
    sketch, _ = save_design(
        cli,
        [
            {
                "command": "sketch.create",
                "params": {
                    "plane": {"face": face["point"], "body": body_name, "index": face["index"]},
                    "name": "FieldCraftXZFaceSketch",
                },
            },
            {"command": "sketch.spline.fit_point", "params": {"points": local_points, "periodic": True}},
            {"command": "sketch.finish", "params": {}},
        ],
        "face-bracket-xz-cut",
        work,
        bracket,
    )
    target = TESTS / "face-bracket-xz-cut.sc"
    target.write_bytes(sketch.read_bytes())
    before = evaluate(cli, bracket, work, "bracket-xz-before")
    extruded, _ = save_design(
        cli,
        [{
            "command": "solid.extrude",
            "params": {
                "sketch": "FieldCraftXZFaceSketch",
                "distance": 2,
                "direction": "negative",
                "operation": "cut",
                "name": "FieldCraft XZ cut",
            },
        }],
        "bracket-xz-cut-result",
        work,
        sketch,
    )
    after = evaluate(cli, extruded, work, "bracket-xz-after")
    assert after["volume_mm3"] < before["volume_mm3"]
    write_json(
        TESTS / "face-bracket-xz-cut.json",
        {
            "entry": "extrude",
            "params": {"source_feature_id": 2, "height_mm": 2, "operation": "cut"},
            "design": "face-bracket-xz-cut.sc",
            "grants": GRANTS,
            "expect": {
                "files": ["companion/face-2-cut-3.tool.model.json"],
                "model": {"body_count": 1, "total": {"volume_mm3": after["volume_mm3"]}},
            },
        },
    )
    text = work / "bracket-xz-cut.txt"
    text.write_text(call(cli, ["textconv", str(sketch)]), encoding="utf-8")
    content = text.read_text(encoding="utf-8")
    history = work / "bracket-xz-cut-history.json"
    call(cli, ["script", str(extruded), "--out", str(history)])
    features = [
        item["params"]["feature"]
        for item in read_json(history)["commands"]
        if item.get("command") == "timeline.add"
    ]
    sketch_feature = next(feature for feature in features if feature.get("name") == "FieldCraftXZFaceSketch")
    plane = sketch_feature["plane"]
    assert plane["plane"]["x_dir"] == [-1.0, 0.0, 0.0]
    assert plane["plane"]["y_dir"] == [0.0, 0.0, 1.0]
    x_dir, y_dir = plane["plane"]["x_dir"], plane["plane"]["y_dir"]
    replay = work / "bracket-xz-cut-replay.sc"
    call(cli, ["run", str(history), "--save", str(replay), "--quiet"])
    replayed = evaluate(cli, replay, work, "bracket-xz-replay")
    assert replayed["volume_mm3"] == after["volume_mm3"]

    # An open spline has no bounded profile for an explicit extrusion selection;
    # the host must refuse before writing a changed design.
    unclosed_script = {
        "commands": [
            {
                "command": "sketch.create",
                "params": {
                    "plane": {"face": [30, 25, 12], "body": "FaceIntegrationBox", "index": 0},
                    "name": "UnclosedFaceSpline",
                },
            },
            {"command": "sketch.spline.fit_point", "params": {"points": [[20, 20], [30, 22], [40, 20]]}},
            {"command": "sketch.finish", "params": {}},
        ]
    }
    unclosed, _ = save_design(cli, unclosed_script["commands"], "face-box-open", work, box)
    target = TESTS / "face-box-open.sc"
    target.write_bytes(unclosed.read_bytes())
    digest_before = hashlib.sha256(unclosed.read_bytes()).hexdigest()
    failure = work / "face-box-open-refusal.txt"
    # Write the attempted command as a separate tiny script for a clear reproducible invocation.
    open_attempt = work / "open-extrusion.json"
    write_json(open_attempt, {"commands": [{
        "command": "solid.extrude",
        "params": {
            "sketch": "UnclosedFaceSpline",
            "profiles": [{"curves": ["s1"]}],
            "distance": 2,
            "direction": "negative",
            "operation": "cut",
            "name": "RejectOpenSpline",
        },
    }]})
    refused = subprocess.run(
        [cli, "run", str(open_attempt), "--in", str(unclosed), "--quiet", "--save", str(work / "should-not-exist.sc")],
        capture_output=True,
        text=True,
        check=False,
    )
    failure.write_text(refused.stderr + refused.stdout, encoding="utf-8")
    assert refused.returncode != 0, "host accepted an open spline as an extrusion profile"
    assert "no profile" in failure.read_text(encoding="utf-8").lower()
    assert not (work / "should-not-exist.sc").exists()
    assert hashlib.sha256(unclosed.read_bytes()).hexdigest() == digest_before
    base_body = evaluate(cli, box, work, "box-stock-for-open-spline")
    open_body = evaluate(cli, unclosed, work, "box-open-spline-stock")
    assert open_body["volume_mm3"] == base_body["volume_mm3"] == 36000.0
    write_json(
        TESTS / "face-box-open.json",
        {
            "entry": "extrude",
            "params": {"source_feature_id": 3, "height_mm": 2, "operation": "cut"},
            "design": "face-box-open.sc",
            "grants": GRANTS,
            "expect": {
                "error": "invalid_input",
                "message": "profile_spans_must_join_at_shared_exact_points",
                "model": {"body_count": 1, "total": {"volume_mm3": 36000.0}},
            },
        },
    )
    return {
        "xz_bracket": {
            "face_index": face["index"],
            "face_area_mm2": face["area"],
            "outward_normal": face["plane_normal"],
            "recorded_sketch_axes": {"x_dir": x_dir, "y_dir": y_dir},
            "local_profile_mm": local_points,
            "extrusion_direction": "negative into the body",
            "volume_before_mm3": before["volume_mm3"],
            "volume_after_mm3": after["volume_mm3"],
            "replay_volume_identical": True,
        },
        "open_spline_refusal": {
            "status": "refused as unbounded profile",
            "source_design_sha256_preserved": digest_before,
            "stock_volume_mm3": open_body["volume_mm3"],
            "extrusion_output_written": False,
        },
    }


def main():
    cli = os.environ.get("SOLIDCRAFT_CLI", "solidcraft-cli")
    version = call(cli, ["--version"]).strip()
    assert version == "solidcraft-cli 0.48.1", f"expected SolidCraft 0.48.1, got {version}"
    sha = hashlib.sha256(BRACKET.read_bytes()).hexdigest()
    assert sha == BRACKET_SHA256, f"ribbon-bracket.sc changed: {sha}"
    with tempfile.TemporaryDirectory(prefix="sc-face-fixtures-") as temporary:
        work = Path(temporary)
        box = base_box(cli, work)
        report = [
            run_case(cli, "box", "cut", "FaceIntegrationBox", box, work),
            run_case(cli, "box", "add", "FaceIntegrationBox", box, work),
            run_case(cli, "box-degree5", "cut", "FaceIntegrationBox", box, work),
            run_case(cli, "bracket", "cut", "SolidCraft boundary", BRACKET, work),
            run_case(cli, "bracket", "add", "SolidCraft boundary", BRACKET, work),
        ]
        extra = qualify_xz_and_open_refusal(cli, BRACKET, box, work)
    print(json.dumps({"host": version, "bracket_sha256": sha, "cases": report, "additional_checks": extra}, indent=2))


if __name__ == "__main__":
    main()
