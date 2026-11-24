"""Orientation and coefficient checks for tensor-product face traces."""

import copy
import unittest

from trivariate.annular import build_annular_pair
from trivariate.trace import extract_trace, reorient_trace, traces_match


class TraceTests(unittest.TestCase):
    def test_annular_interface_traces_match(self):
        model = build_annular_pair({"family": "annular_pair_v1", "heights_mm": [1, 2]})
        lower = extract_trace(model["blocks"][0], "w1")
        upper = extract_trace(model["blocks"][1], "w0")
        self.assertEqual(lower["degrees"], [2, 1])
        self.assertTrue(traces_match(lower, upper))

    def test_extract_trace_independently_checks_all_six_asymmetric_faces(self):
        block = {
            "degrees": [2, 1, 3],
            "controls": [[float(100 * i + 10 * j + k), float(i), float(j), 1.0]
                         for i in range(3) for j in range(2) for k in range(4)],
        }
        degrees = block["degrees"]
        for axis, label in enumerate("uvw"):
            free = [candidate for candidate in range(3) if candidate != axis]
            for side_index in (0, 1):
                expected = []
                fixed = side_index * degrees[axis]
                for a in range(degrees[free[0]] + 1):
                    for b in range(degrees[free[1]] + 1):
                        multi = [None, None, None]
                        multi[axis] = fixed
                        multi[free[0]] = a
                        multi[free[1]] = b
                        i, j, k = multi
                        index = ((i * (degrees[1] + 1)) + j) * (degrees[2] + 1) + k
                        expected.append(block["controls"][index])
                actual = extract_trace(block, f"{label}{side_index}")
                self.assertEqual(actual["degrees"], [degrees[free[0]], degrees[free[1]]])
                self.assertEqual(actual["controls"], expected)

    @staticmethod
    def independent_reorientation(trace, swap, flip_first, flip_second):
        """Build the expected control order using explicit index coordinates."""
        first_degree, second_degree = trace["degrees"]
        rows = []
        out_degrees = ([second_degree, first_degree] if swap
                       else [first_degree, second_degree])
        for i in range(out_degrees[0] + 1):
            row = []
            for j in range(out_degrees[1] + 1):
                if swap:
                    source_i = first_degree - j if flip_second else j
                    source_j = second_degree - i if flip_first else i
                else:
                    source_i = first_degree - i if flip_first else i
                    source_j = second_degree - j if flip_second else j
                row.append(trace["controls"][source_i * (second_degree + 1) + source_j])
            rows.extend(row)
        return {"degrees": out_degrees, "controls": rows}

    def test_reorient_and_compare_all_eight_asymmetric_orientations(self):
        source = {
            "degrees": [2, 3],
            "controls": [[float(10 * i + j), float(i + 1), float(j + 2), 1.0]
                         for i in range(3) for j in range(4)],
        }
        for swap in (False, True):
            for flip_first in (False, True):
                for flip_second in (False, True):
                    target = self.independent_reorientation(
                        source, swap, flip_first, flip_second
                    )
                    self.assertEqual(
                        reorient_trace(source, swap=swap, flip_first=flip_first,
                                       flip_second=flip_second), target
                    )
                    inverse_first = flip_second if swap else flip_first
                    inverse_second = flip_first if swap else flip_second
                    scaled_target = {
                        "degrees": target["degrees"],
                        "controls": [[2.0 * value for value in point]
                                     for point in target["controls"]],
                    }
                    self.assertTrue(traces_match(
                        source, target, swap=swap,
                        flip_first=inverse_first, flip_second=inverse_second,
                    ))
                    self.assertTrue(traces_match(
                        source, scaled_target, swap=swap,
                        flip_first=inverse_first, flip_second=inverse_second,
                    ))

    def test_rejects_one_ulp_and_decimal_proportion_mismatches(self):
        source = {
            "degrees": [1, 0],
            "controls": [[0.0, 0.0, 0.0, 1.0], [1.0, 0.0, 0.0, 1.0]],
        }
        ulp = copy.deepcopy(source)
        ulp["controls"][1][0] = 1.0000000000000002
        self.assertFalse(traces_match(source, ulp))

        # Decimal ratios look equal, but their binary64 values are not proportional.
        first = {"degrees": [0, 0], "controls": [[0.1, 0.0, 0.0, 1.0]]}
        second = {"degrees": [0, 0], "controls": [[0.3, 0.0, 0.0, 3.0]]}
        self.assertFalse(traces_match(first, second))

        isolated_weight = copy.deepcopy(source)
        isolated_weight["controls"][1][3] = 1.0000000000000002
        self.assertFalse(traces_match(source, isolated_weight))

    def test_invalid_blocks_traces_and_orientations_are_rejected(self):
        valid = {"degrees": [1, 1], "controls": [
            [0.0, 0.0, 0.0, 1.0], [0.0, 1.0, 0.0, 1.0],
            [1.0, 0.0, 0.0, 1.0], [1.0, 1.0, 0.0, 1.0],
        ]}
        for malformed in (
            {"degrees": [True, 1], "controls": valid["controls"]},
            {"degrees": [1, 1], "controls": valid["controls"][:-1]},
            {"degrees": [1, 1], "controls": [[0.0, 0.0, 0.0, 0.0]] * 4},
        ):
            with self.assertRaises(ValueError):
                reorient_trace(malformed)
        with self.assertRaises(ValueError):
            extract_trace({"degrees": [1, 1, 1], "controls": []}, "w0")
        with self.assertRaises(ValueError):
            extract_trace({"degrees": [1, 1, 1], "controls": [[0, 0, 0, 1]] * 8}, "x0")
        with self.assertRaises(ValueError):
            reorient_trace(valid, flip_first=1)


if __name__ == "__main__":
    unittest.main()
