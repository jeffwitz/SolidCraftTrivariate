"""Replay adjacent-parent witnesses on every polygon and axial Boolean census state."""

from collections import Counter, defaultdict
import json
from pathlib import Path
import subprocess

from check_pair_parent_report import check
from census_support_interfaces import sha256

ROOT = Path(__file__).resolve().parents[1]


def main():
    cli = ROOT / "target/release/solidcraft-trivariate"
    census = json.loads((ROOT / "build/support-census/census.json").read_text())
    work = ROOT / "build/pair-parent-expanded"
    work.mkdir(parents=True, exist_ok=True)
    strip = work / "strip.model.json"
    subprocess.run([cli, "build", ROOT / "tests/references/nurbs_strip.recipe.json", strip],
                   check=True, stdout=subprocess.DEVNULL)
    cases = []
    for row in census["polygon_recipes"]["entries"]:
        assert row["status"] == "admitted"
        member, = row["memberships"]
        cases.append((member["corpus"], member["case_id"], ROOT / row["model_path"],
                      row["model_sha256"], row["declared_interface_count"]))
    for row in census["named_states"]["states"]:
        if "fieldcraft" not in row["state_id"].lower():
            cases.append(("axial_booleans", row["state_id"], ROOT / row["source"],
                          row["model_sha256"], row["declared_interface_count"]))
    cases.append(("source_strip", "nurbs_strip", strip, sha256(strip.read_bytes()), 1))
    groups = defaultdict(Counter)
    axes = defaultdict(Counter)
    records = []
    for group, name, model, model_hash, faces in cases:
        assert sha256(model.read_bytes()) == model_hash
        report = work / f"{group}-{name}.pairs.json"
        subprocess.run([cli, "recover-pair-parents", model, report],
                       check=True, stdout=subprocess.DEVNULL)
        evidence = check(model, report)
        assert evidence["interfaces"] == faces
        counts = dict(interfaces=faces, recovered=evidence["recovered_and_exactly_replayed"],
                      refused=evidence["refused"], indeterminate=evidence["indeterminate"])
        assert counts["recovered"] + counts["refused"] + counts["indeterminate"] == faces
        for axis, outcome in evidence["counts_by_first_interface_axis"].items():
            axes[group][axis] += outcome["recovered"]
        groups[group].update(counts)
        groups[group]["models"] += 1
        records.append(dict(group=group, name=name, model_sha256=model_hash, **counts))
    assert sum(item["interfaces"] for group, item in groups.items()
               if group not in ("axial_booleans", "source_strip")) == 3136
    result = dict(scope="Exact pairwise single-span parents, with independent coefficient replay",
                  groups=dict(groups), recovered_by_first_axis=dict(axes), models=records)
    (work / "summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result["groups"], sort_keys=True))


if __name__ == "__main__":
    main()
