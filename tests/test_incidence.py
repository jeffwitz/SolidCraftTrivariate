"""Independent exact checks of whole-face incidence reconstruction."""

from copy import deepcopy
from fractions import Fraction
import importlib.util
import unittest

from trivariate.incidence import enumerate_faces
from trivariate.trace import extract_trace, reorient_trace, traces_match


def box(block_id, x0, x1, *, swap_yz=False, flip_y=False, flip_z=False):
    """Return a box with a quadratic transverse grid and exact controls."""
    controls = []
    for i in range(2):
        for j in range(3):
            for k in range(3):
                y_parameter = Fraction(j, 2)
                z_parameter = Fraction(k, 2)
                if flip_y:
                    y_parameter = 1 - y_parameter
                if flip_z:
                    z_parameter = 1 - z_parameter
                y, z = ((z_parameter, y_parameter) if swap_yz
                        else (y_parameter, z_parameter))
                controls.append([Fraction(x0) + (x1-x0)*i, y, z, 1])
    return {"id": block_id, "degrees": [1, 2, 2], "controls": controls}


class IncidenceTests(unittest.TestCase):
    def test_adjacent_boxes_have_one_interface_and_ten_exterior_faces(self):
        result = enumerate_faces([box("left", 0, 1), box("right", 1, 2)])
        self.assertEqual(len(result["interfaces"]), 1)
        self.assertEqual(len(result["exterior"]), 10)
        self.assertEqual(result["interfaces"][0]["first"], ["left", "u1"])
        self.assertEqual(result["interfaces"][0]["second"], ["right", "u0"])
        orientation = result["interfaces"][0]["orientation"]
        first = extract_trace(box("left", 0, 1), "u1")
        second = extract_trace(box("right", 1, 2), "u0")
        self.assertTrue(traces_match(first, reorient_trace(second, **orientation)))

    def test_all_eight_square_orientations_are_resolved_uniquely(self):
        orientations = set()
        for swap in (False, True):
            for flip_first in (False, True):
                for flip_second in (False, True):
                    left = box("left", 0, 1)
                    right = box("right", 1, 2, swap_yz=swap,
                                flip_y=flip_first, flip_z=flip_second)
                    orientation = enumerate_faces([left, right])["interfaces"][0]["orientation"]
                    orientations.add(tuple(orientation[key] for key in
                                           ("swap", "flip_first", "flip_second")))
        self.assertEqual(len(orientations), 8)

    def test_exact_degree_elevation_accepts_different_whole_face_degrees(self):
        left = box("left", 0, 1)
        right = box("right", 1, 2)
        # Elevate only the right block in the transverse directions by sampling
        # the same affine map at its degree-three Bernstein control locations.
        right["degrees"] = [1, 3, 3]
        right["controls"] = [
            [1 + i, Fraction(j, 3), Fraction(k, 3), 1]
            for i in range(2) for j in range(4) for k in range(4)
        ]
        result = enumerate_faces([left, right])
        self.assertEqual(len(result["interfaces"]), 1)
        self.assertEqual(len(result["exterior"]), 10)

    def test_common_positive_homogeneous_scale_is_accepted(self):
        left, right = box("left", 0, 1), box("right", 1, 2)
        right["controls"] = [
            [7 * value for value in point] for point in right["controls"]
        ]
        result = enumerate_faces([left, right])
        self.assertEqual(len(result["interfaces"]), 1)

    def test_common_corners_with_altered_interior_trace_are_rejected(self):
        left = box("left", 0, 1)
        right = box("right", 1, 2)
        # Center control on the u=0 face; preserve the four physical corners.
        right["controls"][4][2] += Fraction(1, 10)
        with self.assertRaisesRegex(ValueError, "incompatible_face_traces"):
            enumerate_faces([left, right])

    def test_triple_face_incidence_is_rejected(self):
        blocks = [box("left", 0, 1), box("right_a", 1, 2), box("right_b", 1, 3)]
        with self.assertRaisesRegex(ValueError, "triple_face_incidence"):
            enumerate_faces(blocks)

    def test_duplicate_identity_and_degenerate_face_vertices_are_rejected(self):
        left = box("duplicate", 0, 1)
        with self.assertRaisesRegex(ValueError, "duplicate block identity"):
            enumerate_faces([left, deepcopy(left)])
        collapsed = box("collapsed", 0, 1)
        for point in collapsed["controls"]:
            point[1] = Fraction(0)
        with self.assertRaisesRegex(ValueError, "degenerate_face_vertices"):
            enumerate_faces([collapsed])

    @unittest.skipUnless(importlib.util.find_spec("trivariate.bracket"),
                         "bracket constructor has not been integrated")
    def test_bracket_expected_whole_face_counts_when_available(self):
        from trivariate.bracket import bracket_blocks
        blocks = bracket_blocks()
        self.assertEqual(len(blocks), 36)
        result = enumerate_faces(blocks)
        self.assertEqual(len(result["interfaces"]), 55)
        self.assertEqual(len(result["exterior"]), 106)


if __name__ == "__main__":
    unittest.main()
