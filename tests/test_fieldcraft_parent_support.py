"""Tests for exact FieldCraft homogeneous normal-jet evidence."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "fieldcraft_parent_support", ROOT / "scripts/audit_fieldcraft_parent_support.py"
)
audit = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(audit)


def affine_block(identifier: str, x0: int, x1: int) -> dict:
    controls = []
    for u in (x0, x1):
        for v in (0, 1):
            for w in (0, 1):
                controls.append([u, v, w, 1])
    return {"id": identifier, "degrees": [1, 1, 1], "controls": controls}


def quadratic_block(identifier: str, xs: tuple[int, int, int]) -> dict:
    controls = []
    for u in xs:
        for v in (0, 1):
            for w in (0, 1):
                controls.append([u, v, w, 1])
    return {"id": identifier, "degrees": [2, 1, 1], "controls": controls}


class ParentSupportTests(unittest.TestCase):
    def test_exact_split_with_proportional_homogeneous_jets_is_candidate(self):
        a = affine_block("a", 0, 1)
        b = affine_block("b", 1, 3)
        da, db, net_a, net_b = audit.common_degree(a, b, "u1", "u0", {})
        evidence = audit.jet_restriction(net_a, net_b, da, db, "u1", "u0", {})
        self.assertTrue(evidence["full_jet_match"])
        self.assertEqual(
            evidence["candidate_interval_length_ratio_first_to_second"], [1, 2]
        )
        self.assertEqual(evidence["compared_jet_orders"], 1)

    def test_nonmatching_derivative_fails_exact_candidate(self):
        a = quadratic_block("a", (0, 1, 2))
        b = quadratic_block("b", (2, 4, 7))
        da, db, net_a, net_b = audit.common_degree(a, b, "u1", "u0", {})
        evidence = audit.jet_restriction(net_a, net_b, da, db, "u1", "u0", {})
        self.assertFalse(evidence["full_jet_match"])

    def test_all_fieldcraft_model_interfaces_are_enumerated_and_no_full_jet_matches(
        self,
    ):
        paths = sorted(
            (ROOT / "docs/site/data").glob("fieldcraft_*/fieldcraft_*.model.json")
        )
        self.assertEqual(len(paths), 5)
        for path in paths:
            model = audit.audit_model(path)
            self.assertEqual((model["block_count"], model["interface_count"]), (10, 11))
            self.assertEqual(model["candidate_count"], 0)
            self.assertTrue(
                all(
                    "face_gauge_second_to_first" in row["jet_evidence"]
                    for row in model["interfaces"]
                )
            )

    def test_reversed_normal_and_zero_higher_jets_are_compatible(self):
        a = quadratic_block("a", (0, 1, 2))
        b = quadratic_block("b", (6, 4, 2))
        da, db, na, nb = audit.common_degree(a, b, "u1", "u1", {})
        result = audit.jet_restriction(na, nb, da, db, "u1", "u1", {})
        self.assertTrue(result["full_jet_match"])
        self.assertTrue(result["jet_scaled_ratios"][1]["both_jets_zero"])
        b["controls"][-1][0] += 1
        da, db, na, nb = audit.common_degree(a, b, "u1", "u1", {})
        self.assertEqual(
            audit.jet_restriction(na, nb, da, db, "u1", "u1", {})["face_gauge"],
            "mismatch",
        )


if __name__ == "__main__":
    unittest.main()
