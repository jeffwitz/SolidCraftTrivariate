"""Focused tests for the independent support/interface census helpers."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "census_support_interfaces", ROOT / "scripts/census_support_interfaces.py")
assert SPEC and SPEC.loader
CENSUS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CENSUS)


class SupportCensusTests(unittest.TestCase):
    def test_recipe_hash_dedup_retains_every_membership(self):
        recipe = {"family": "fixture", "height_mm": 3}
        entries = CENSUS.grouped_recipes([
            {"corpus": "a", "case_id": "first", "recipe": recipe,
             "expected_admissible": True, "raw_sha256": None},
            {"corpus": "b", "case_id": "second", "recipe": recipe,
             "expected_admissible": True, "raw_sha256": "source-hash"},
        ])
        self.assertEqual(len(entries), 1)
        self.assertEqual([m["corpus"] for m in entries[0]["memberships"]], ["a", "b"])

    def test_exact_face_comparison_handles_orientation(self):
        a = [[(0, 0, 0, 2), (0, 0, 2, 2)],
             [(2, 0, 0, 2), (2, 0, 2, 2)]]
        b = [[a[1][0], a[1][1]], [a[0][0], a[0][1]]]
        exact, detail = CENSUS.exact_homogeneous_face_match(
            ((1, 1), a), ((1, 1), b), {"flip_first": True})
        self.assertTrue(exact)
        self.assertEqual(detail["positive_gauge_second_to_first"], [1, 1])
        b[0][0] = (4, 0, 0, 2)
        self.assertFalse(CENSUS.exact_homogeneous_face_match(
            ((1, 1), a), ((1, 1), b), {"flip_first": True})[0])

    def test_frozen_polygon_corpora_have_68_memberships(self):
        rows = CENSUS.load_polygon_recipes()
        entries = CENSUS.grouped_recipes(rows)
        self.assertEqual(len(rows), 68)
        self.assertEqual(len(entries), 68)
        self.assertEqual(CENSUS.Counter(r["corpus"] for r in rows), {
            "historical_45": 45, "supplemental_12": 12, "heldout_11": 11,
        })

    def test_s2_partial_faces_are_not_counted_as_full_trace_matches(self):
        result = CENSUS.s2_census()
        self.assertEqual(len(result["states"]), 3)
        for state in result["states"]:
            self.assertEqual(state["interface_count"], 11)
            self.assertEqual(state["full_exact_trace_count"], 7)
            self.assertEqual(state["partial_trace_count"], 4)
            self.assertEqual(state["non_face_collapsed_polar_line_constraints"], 2)


if __name__ == "__main__":
    unittest.main()
