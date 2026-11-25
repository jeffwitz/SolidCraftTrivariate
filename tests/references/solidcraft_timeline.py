#!/usr/bin/env python3
"""Independent semantic and byte-parity oracle for solidcraft/66 bracket histories.

The host is used only to extract a .sc document's canonical history with
``solidcraft-cli script``. This oracle never asks SolidCraft to execute a user
script. The translated native model is compared with both the direct-command
translation and each SC1 ``native.json`` snapshot.
"""

from __future__ import annotations

import copy
import json
import platform
import shutil
import subprocess
import sys
import tempfile
import time
from decimal import Decimal
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
HOST = shutil.which("solidcraft-cli") or "/usr/bin/solidcraft-cli"
NATIVE = ROOT / "target/release/solidcraft-trivariate"
CASE_ROOT = ROOT / "build/solidcraft-bridge"
REPORT = ROOT / "docs/results/solidcraft-timeline-20261125.json"

POINTS = {
    "bottom": ((0, 0), (60, 0)),
    "right": ((60, 0), (60, 8)),
    "top": ((60, 8), (8, 8)),
    "inner_vertical": ((8, 8), (8, 50)),
    "upper_left": ((8, 50), (0, 50)),
    "left": ((0, 50), (0, 0)),
}
FEATURE_KEYS = {
    "sketch": {"id", "name", "plane", "sketch", "suppressed", "type"},
    "extrude": {
        "body_names", "extent", "id", "name", "operation", "param_names",
        "profiles", "sketch", "suppressed", "type",
    },
    "fillet": {
        "edge_names", "edge_refs", "edges", "id", "name", "param_names",
        "radius", "suppressed", "type",
    },
    "hole": {
        "diameter", "direction", "face_names", "hole", "id", "name",
        "param_names", "position", "suppressed", "type",
    },
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def fraction(value: object, *, strings: bool = False) -> Fraction:
    require(not isinstance(value, bool), "boolean is not a numeric parameter")
    if isinstance(value, Decimal):
        return Fraction(value)
    if isinstance(value, int):
        return Fraction(value)
    if strings and isinstance(value, str):
        return Fraction(Decimal(value))
    raise ValueError("expected literal finite decimal")


def vector(value: object, expected: tuple[int, ...]) -> None:
    require(isinstance(value, list) and len(value) == len(expected), "vector shape mismatch")
    require(tuple(fraction(x) for x in value) == tuple(Fraction(x) for x in expected), "vector mismatch")


def parse_history(history: dict) -> dict:
    require(set(history) == {"about", "commands"}, "unexpected history keys")
    about = history["about"]
    require(set(about) == {"design", "format", "replay", "version"}, "unexpected about keys")
    require(about["format"] == "solidcraft/66", "unsupported history format")
    require(about["design"] in ("editable", "original"), "unsupported design source")
    require(isinstance(about["version"], str) and bool(about["version"]), "missing host version")
    require(isinstance(about["replay"], str) and bool(about["replay"]), "missing replay contract")
    commands = history["commands"]
    require(isinstance(commands, list) and len(commands) == 5, "expected five active features")
    features = []
    for index, entry in enumerate(commands):
        require(set(entry) == {"command", "params"}, "unexpected timeline command fields")
        require(entry["command"] == "timeline.add", "only timeline.add is in the qualified subset")
        require(set(entry["params"]) == {"feature"}, "unexpected timeline parameter fields")
        feature = entry["params"]["feature"]
        kind = feature.get("type")
        require(kind == ("sketch", "extrude", "fillet", "hole", "hole")[index], "feature order mismatch")
        require(set(feature) == FEATURE_KEYS[kind], "unexpected feature fields")
        require(feature["id"] == index + 1 and not isinstance(feature["id"], bool), "feature identity mismatch")
        require(feature["suppressed"] is False, "suppressed feature")
        features.append(feature)

    sketch_feature, extrude, fillet, hole_a, hole_b = features
    require(sketch_feature["name"] == "Profile", "sketch identity mismatch")
    require(sketch_feature["plane"] == {"name": "XY", "type": "origin"}, "sketch plane mismatch")
    sketch = sketch_feature["sketch"]
    require(set(sketch) == {"constraints", "counters", "curves", "points"}, "unexpected sketch fields")
    require(sketch["counters"] == {"k": 1, "l": 6}, "sketch counters mismatch")

    sketch_points = sketch["points"]
    require(isinstance(sketch_points, list) and len(sketch_points) == 7, "expected origin and six vertices")
    for point in sketch_points:
        require(set(point) == {"fixed", "id", "pos"}, "unexpected sketch point fields")
        require(isinstance(point["id"], str) and point["id"], "empty sketch point id")
        require(isinstance(point["fixed"], bool), "invalid point fixed flag")
    require(sketch_points[0]["id"] == "origin" and sketch_points[0]["fixed"], "origin point reference mismatch")
    vector(sketch_points[0]["pos"], (0, 0))
    point_ids = [p["id"] for p in sketch_points]
    require(len(set(point_ids)) == len(point_ids), "duplicate sketch point id")
    positions = [tuple(fraction(x) for x in p["pos"]) for p in sketch_points]
    require(all(len(p) == 2 for p in positions), "sketch point dimension mismatch")
    require(all(not sketch_points[i]["fixed"] for i in range(1, 7)), "only the origin may be fixed")
    origin_vertices = [i for i in range(1, 7) if positions[i] == (0, 0)]
    require(len(origin_vertices) == 1, "origin vertex is ambiguous")

    constraints = sketch["constraints"]
    require(isinstance(constraints, list) and len(constraints) == 1, "unknown or missing sketch constraint")
    constraint = constraints[0]
    require(set(constraint) == {"id", "p", "q", "type"}, "unexpected constraint fields")
    require(constraint == {"id": "k1", "p": 0, "q": origin_vertices[0], "type": "coincident"}, "unresolved or altered constraint")

    curves = sketch["curves"]
    require(isinstance(curves, list) and len(curves) == 6, "expected six contour curves")
    require(all(set(c) == {"a", "b", "construction", "id", "type"} for c in curves), "unexpected curve fields")
    require(all(c["type"] == "line" and c["construction"] is False for c in curves), "unsupported or construction curve")
    curve_ids = [c["id"] for c in curves]
    require(set(curve_ids) == {f"l{i}" for i in range(1, 7)} and len(set(curve_ids)) == 6, "unknown curve identifiers")
    outgoing = [0] * 7
    incoming = [0] * 7
    segment_for_id: dict[str, str] = {}
    for curve in curves:
        a, b = curve["a"], curve["b"]
        require(type(a) is int and type(b) is int and 1 <= a <= 6 and 1 <= b <= 6 and a != b, "invalid curve endpoint id")
        outgoing[a] += 1
        incoming[b] += 1
        edge = (positions[a], positions[b])
        matches = [name for name, ends in POINTS.items() if edge == tuple(tuple(Fraction(x) for x in p) for p in ends)]
        require(len(matches) == 1, "curve graph is not the supported directed L contour")
        require(matches[0] not in segment_for_id.values(), "duplicate contour segment")
        segment_for_id[curve["id"]] = matches[0]
    require(all(outgoing[i] == incoming[i] == 1 for i in range(1, 7)), "curve graph is not a single closed loop")
    require(set(segment_for_id.values()) == set(POINTS), "incomplete contour graph")

    # Stable point names emitted by this solidcraft/66 subset also bind each
    # indexed endpoint to the declared lN object id.
    expected_point_names = ["origin", "l1.start", "l1.end", "l2.end", "l3.end", "l4.end", "l5.end"]
    require(point_ids == expected_point_names, "sketch point IDs are inconsistent with curve endpoints")

    require(extrude["name"] == "Extrude1" and extrude["operation"] == "new_body", "extrude operation mismatch")
    require(extrude["body_names"] == ["Bracket"] and extrude["sketch"] == 1, "extrude dependency mismatch")
    require(extrude["param_names"] == ["d1"] and extrude["profiles"] == {"type": "all"}, "extrude selection mismatch")
    extent = extrude["extent"]
    require(set(extent) == {"convention", "direction", "distance"}, "unexpected extrusion extent fields")
    require(extent["convention"] == "total_distance" and extent["direction"] == "positive", "extrusion direction mismatch")
    require(fraction(extent["distance"], strings=True) == 40, "extrusion distance mismatch")

    require(fillet["name"] == "Fillet1" and fillet["param_names"] == ["d2"], "fillet identity mismatch")
    require(fraction(fillet["radius"], strings=True) == 6, "fillet radius mismatch")
    require(len(fillet["edges"]) == 1, "fillet edge count mismatch")
    vector(fillet["edges"][0], (8, 8, 20))
    require(len(fillet["edge_names"]) == len(fillet["edge_refs"]) == 1, "fillet edge reference count mismatch")
    named_curves = fillet["edge_names"][0].split("|")
    require(len(named_curves) == 2 and all(x.startswith("F2:side:") for x in named_curves), "malformed fillet reference")
    fillet_indices = [curve_ids.index(x.removeprefix("F2:side:")) for x in named_curves]
    require({segment_for_id[curve_ids[i]] for i in fillet_indices} == {"top", "inner_vertical"}, "fillet does not reference the inside corner")
    edge_ref = fillet["edge_refs"][0]
    require(set(edge_ref) == {"dir", "max", "min"}, "unexpected analytic fillet edge fields")
    vector(edge_ref["dir"], (0, 0, 1)); vector(edge_ref["min"], (0, 0, 0)); vector(edge_ref["max"], (60, 50, 40))

    holes_by_segment: dict[str, tuple[int, Fraction]] = {}
    for hole in (hole_a, hole_b):
        require(hole["name"] in ("Hole1", "Hole2"), "hole feature name mismatch")
        require(hole["hole"] == {"type": "simple"}, "unsupported hole type")
        require(len(hole["face_names"]) == 1, "hole face reference count mismatch")
        face_name = hole["face_names"][0]
        require(face_name.startswith("F2:side:"), "hole does not reference an extrusion side face")
        line_id = face_name.removeprefix("F2:side:")
        require(line_id in segment_for_id, "hole face references an unknown sketch curve")
        segment = segment_for_id[line_id]
        if segment == "top":
            vector(hole["position"], (40, 8, 20)); vector(hole["direction"], (0, -1, 0))
        elif segment == "inner_vertical":
            vector(hole["position"], (8, 32, 20)); vector(hole["direction"], (-1, 0, 0))
        else:
            raise ValueError("hole support is not one of the two qualified faces")
        diameter = fraction(hole["diameter"], strings=True)
        ticks = diameter * 8
        require(ticks.denominator == 1 and 40 <= ticks <= 80, "diameter outside the eighth-mm SC1 lattice")
        require(segment not in holes_by_segment, "duplicate hole support")
        holes_by_segment[segment] = (hole["id"], diameter)
        require(hole["param_names"] == [f"d{hole['id'] - 1}"], "hole parameter identity mismatch")
    require(set(holes_by_segment) == {"top", "inner_vertical"}, "both qualified hole supports are required")

    # Reconstruct only after resolving every source id and physical support.
    top_diameter = holes_by_segment["top"][1]
    vertical_diameter = holes_by_segment["inner_vertical"][1]
    return {
        "commands": [
            {"command": "sketch.create", "params": {"plane": "XY", "name": "Profile"}},
            {"command": "sketch.line", "params": {"points": [[0, 0], [60, 0], [60, 8], [8, 8], [8, 50], [0, 50]], "closed": True}},
            {"command": "sketch.finish", "params": {}},
            {"command": "solid.extrude", "params": {"distance": 40, "body_name": "Bracket"}},
            {"command": "solid.fillet", "params": {"edges": [[8, 8, 20]], "radius": 6}},
            {"command": "solid.hole", "params": {"position": [40, 8, 20], "diameter": Decimal(top_diameter.numerator) / Decimal(top_diameter.denominator)}},
            {"command": "solid.hole", "params": {"position": [8, 32, 20], "diameter": Decimal(vertical_diameter.numerator) / Decimal(vertical_diameter.denominator)}},
        ]
    }


def direct_recipe(script: dict) -> tuple[int, int]:
    """Read the direct SolidCraft command form into the independent radius recipe."""
    commands = script["commands"]
    require(len(commands) == 7, "direct-command baseline must contain seven operations")
    names = [entry["command"] for entry in commands]
    require(names == ["sketch.create", "sketch.line", "sketch.finish", "solid.extrude", "solid.fillet", "solid.hole", "solid.hole"], "direct command sequence mismatch")
    require(commands[0]["params"]["plane"] == "XY", "direct sketch plane mismatch")
    contour = [[0, 0], [60, 0], [60, 8], [8, 8], [8, 50], [0, 50]]
    source_points = commands[1]["params"]["points"]
    require(len(source_points) == 6 and any(source_points == contour[k:] + contour[:k] for k in range(6)), "direct contour mismatch")
    require(commands[1]["params"]["closed"] is True, "direct contour must be closed")
    require(fraction(commands[3]["params"]["distance"]) == 40, "direct extrusion mismatch")
    require(fraction(commands[4]["params"]["radius"]) == 6, "direct fillet mismatch")
    radii: dict[str, int] = {}
    for command in commands[5:]:
        params = command["params"]
        position = tuple(fraction(x) for x in params["position"])
        support = "top" if position == (40, 8, 20) else "inner_vertical" if position == (8, 32, 20) else None
        require(support is not None, "direct hole position is unsupported")
        ticks = fraction(params["diameter"]) * 8
        require(ticks.denominator == 1 and 40 <= ticks <= 80, "direct diameter outside qualified lattice")
        require(support not in radii, "duplicate direct hole support")
        radii[support] = ticks.numerator
    require(set(radii) == {"top", "inner_vertical"}, "direct command misses a qualified hole")
    return radii["top"], radii["inner_vertical"]


def load_json(path: Path) -> object:
    return json.loads(path.read_text(), parse_float=Decimal, parse_int=int, parse_constant=lambda _: (_ for _ in ()).throw(ValueError("nonfinite JSON")))


def run_cli(args: list[str], timeout: int = 60) -> subprocess.CompletedProcess[str]:
    return subprocess.run(args, cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout)


def verify_rejection(binary: str, history: dict, name: str, scratch: Path) -> dict:
    changed = copy.deepcopy(history)
    commands = changed["commands"]
    if name == "construction_geometry":
        commands[0]["params"]["feature"]["sketch"]["curves"][0]["construction"] = True
    elif name == "unknown_constraint":
        commands[0]["params"]["feature"]["sketch"]["constraints"].append({"id": "k2", "p": 1, "q": 2, "type": "distance"})
    elif name == "suppressed_feature":
        commands[1]["params"]["feature"]["suppressed"] = True
    elif name == "removed_operation":
        commands.pop()
    elif name == "unresolved_expression":
        commands[1]["params"]["feature"]["extent"]["distance"] = "d1"
    elif name == "inconsistent_selection":
        commands[2]["params"]["feature"]["edge_names"] = ["F2:side:l2|F2:side:l4"]
    elif name == "wrong_hole_direction":
        commands[3]["params"]["feature"]["direction"] = [1, 0, 0]
    else:
        raise AssertionError(name)
    source = scratch / f"invalid-{name}.json"
    target = scratch / f"unchanged-{name}.json"
    source.write_text(json.dumps(
        changed,
        separators=(",", ":"),
        default=lambda value: int(value) if value == value.to_integral_value() else float(value),
    ))
    sentinel = b"previous-model-must-survive\n"
    target.write_bytes(sentinel)
    started = time.perf_counter()
    result = run_cli([binary, "import-solidcraft", str(source), str(target)])
    elapsed = time.perf_counter() - started
    require(result.returncode != 0, f"mutation accepted: {name}")
    require(target.read_bytes() == sentinel, f"failed import changed target: {name}")
    return {"case": name, "rejected": True, "target_preserved": True, "elapsed_s": elapsed}


def main() -> None:
    require(Path(HOST).exists(), f"SolidCraft CLI not found: {HOST}")
    require(NATIVE.exists(), f"build release binary first: {NATIVE}")
    cases = sorted(p.parent for p in CASE_ROOT.glob("x*/editable.sc"))
    require(len(cases) == 5, f"expected five host .sc documents, got {len(cases)}")
    results = []
    errors = []
    total_start = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix="solidcraft-timeline-") as temp:
        scratch = Path(temp)
        for case_dir in cases:
            label = case_dir.name
            timeline_path = scratch / f"{label}.json"
            extract_start = time.perf_counter()
            extracted = run_cli([HOST, "script", str(case_dir / "editable.sc"), "--out", str(timeline_path)])
            extract_s = time.perf_counter() - extract_start
            require(extracted.returncode == 0, f"host timeline extraction failed for {label}: {extracted.stderr}")
            history = load_json(timeline_path)
            direct_expected = parse_history(history)
            direct_source_path = case_dir / "script.json"
            direct_source = load_json(direct_source_path)
            require(direct_recipe(direct_source) == direct_recipe(direct_expected), f"semantic recipe differs from direct commands for {label}")

            candidate_path = scratch / f"candidate-{label}.json"
            direct_model_path = scratch / f"direct-{label}.json"
            candidate_start = time.perf_counter()
            candidate = run_cli([str(NATIVE), "import-solidcraft", str(timeline_path), str(candidate_path)])
            translate_s = time.perf_counter() - candidate_start
            require(candidate.returncode == 0, f"native timeline translation failed for {label}: {candidate.stdout} {candidate.stderr}")
            direct_start = time.perf_counter()
            direct = run_cli([str(NATIVE), "import-solidcraft", str(direct_source_path), str(direct_model_path)])
            direct_s = time.perf_counter() - direct_start
            require(direct.returncode == 0, f"direct native translation failed for {label}: {direct.stdout} {direct.stderr}")
            candidate_bytes = candidate_path.read_bytes()
            direct_bytes = direct_model_path.read_bytes()
            baseline_path = case_dir / "native.json"
            baseline_bytes = baseline_path.read_bytes()
            require(candidate_bytes == direct_bytes, f"timeline/direct model byte mismatch for {label}")
            require(candidate_bytes == baseline_bytes, f"direct model differs from saved SC1 native model for {label}")
            results.append({
                "case": label,
                "host_history_format": history["about"]["format"],
                "host_version": history["about"]["version"],
                "timeline_operations": len(history["commands"]),
                "semantic_recipe_equal_to_direct": True,
                "model_byte_equal_to_direct": True,
                "model_byte_equal_to_saved_sc1_native": True,
                "model_bytes": len(candidate_bytes),
                "extract_s": extract_s,
                "timeline_translation_s": translate_s,
                "direct_translation_s": direct_s,
            })

        first = load_json(scratch / f"{cases[0].name}.json")
        for name in (
            "construction_geometry", "unknown_constraint", "suppressed_feature", "removed_operation",
            "unresolved_expression", "inconsistent_selection", "wrong_hole_direction",
        ):
            errors.append(verify_rejection(str(NATIVE), first, name, scratch))

    report = {
        "protocol": "solidcraft_timeline_v66_bracket_subset_v1",
        "host_binary": str(HOST),
        "host_version": "0.48.1",
        "native_binary": str(NATIVE.relative_to(ROOT)),
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "cases": results,
        "invalid_cases": errors,
        "all_five_model_bytes_equal": len(results) == 5 and all(x["model_byte_equal_to_direct"] and x["model_byte_equal_to_saved_sc1_native"] for x in results),
        "all_invalid_targets_preserved": all(x["target_preserved"] for x in errors),
        "elapsed_s": time.perf_counter() - total_start,
        "scope": "Only the canonical five-operation solidcraft/66 L-bracket timeline is supported; scripts are extracted, not executed in the SolidCraft kernel.",
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"report": str(REPORT.relative_to(ROOT)), "cases": len(results), "invalid_cases": len(errors), "elapsed_s": report["elapsed_s"]}, indent=2))


if __name__ == "__main__":
    main()
