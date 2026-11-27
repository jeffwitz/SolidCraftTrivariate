"""Index face roles by their ordinal in the original, fingerprinted model topology."""

import argparse
from collections import defaultdict
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def pair(a, b):
    return tuple(sorted((tuple(a), tuple(b))))


def summarize(census_path, transitions):
    census = json.loads(census_path.read_bytes())
    cases = {}
    for entry in census["polygon_recipes"]["entries"]:
        model = json.loads((ROOT / entry["model_path"]).read_bytes())
        state = entry["native_inherited_space"]["pz2_inherited_c1"]
        fields = json.loads((ROOT / state["package_path"]).read_bytes())["field"]
        fingerprint = entry["recipe_sha256"]
        tests = json.loads(
            (transitions / (fingerprint + ".interfaces.json")).read_bytes()
        )
        certificates = {pair(r["first"], r["second"]): r for r in fields["interfaces"]}
        transitions_here = {pair(*r["faces"]): r for r in tests}
        categories = defaultdict(list)
        for i, link in enumerate(model["topology"]["interfaces"]):
            key = pair(link["first"], link["second"])
            field = certificates.pop(key)
            transition = transitions_here.pop(key)
            if field["role"] == "common_parent_subdivision":
                assert (
                    field["basis_evidence"]
                    == "same_homogeneous_section_affine_physical_z"
                )
                role = "inherited_axial"
            else:
                assert field["field_order"] == 0
                direction = (transition["status"], transition.get("reversed_status"))
                role = next(
                    (
                        t
                        for t in ("constant_transition", "affine_transition")
                        if t in direction
                    ),
                    "outside_affine_test",
                )
            categories[role].append(i)
        assert not certificates and not transitions_here
        labels = {(f["block"], f["side"]): f["feature"] for f in model["design_faces"]}
        exterior = defaultdict(list)
        for i, face in enumerate(model["topology"]["exterior"]):
            exterior[labels[tuple(face)]].append(i)
        cases[fingerprint] = dict(
            names=[v["case_id"] for v in entry["memberships"]],
            model_sha256=entry["model_sha256"],
            interfaces=dict(categories),
            exterior_by_feature=dict(exterior),
            scalar_dofs=[fields["c0_scalar_dofs"], fields["scalar_dofs"]],
        )
    return dict(
        schema="support_role_index_v1",
        index_units="zero-based entries of topology.interfaces and topology.exterior in each source model",
        roles=dict(
            inherited_axial="complete section matches in homogeneous coordinates with common physical translation; native field restores C1",
            constant_transition="separate supports; exact constant first-jet transition in at least one direction; retained C0",
            affine_transition="separate supports; exact affine first-jet transition in at least one direction; retained C0",
            outside_affine_test="separate supports; neither normalized affine test succeeded; retained C0",
        ),
        polygons=cases,
        named_summary=census["named_state_summary"],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--census", type=Path, default=ROOT / "build/support-census/census.json"
    )
    parser.add_argument(
        "--transitions", type=Path, default=ROOT / "build/transverse-transitions"
    )
    parser.add_argument(
        "--output", type=Path, default=ROOT / "build/support-census/role-index.json"
    )
    args = parser.parse_args()
    report = summarize(args.census, args.transitions)
    args.output.write_text(
        json.dumps(report, separators=(",", ":"), sort_keys=True) + "\n"
    )
    print(len(report["polygons"]), "models indexed")
