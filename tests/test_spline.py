"""Independent checks for the native rational tensor-product map."""

import itertools
import json
import math
from pathlib import Path
import unittest

from trivariate.spline import RationalBezierMap


TOL = 1e-12
ROOT = Path(__file__).resolve().parent


def quarter_annulus_oracle(u, v, w):
    """Evaluate the analytic rational quarter-annulus and its derivatives.

    This oracle uses the closed-form rational quadratic basis, independently
    of the tensor-product control-net evaluator under test.
    """
    q = math.sqrt(0.5)
    one_minus_u = 1.0 - u
    a = one_minus_u * one_minus_u + 2.0 * q * u * one_minus_u
    b = 2.0 * q * u * one_minus_u + u * u
    denominator = one_minus_u * one_minus_u + 2.0 * q * u * one_minus_u + u * u
    da = -2.0 * one_minus_u + 2.0 * q * (1.0 - 2.0 * u)
    db = 2.0 * q * (1.0 - 2.0 * u) + 2.0 * u
    dd = -2.0 * one_minus_u + 2.0 * q * (1.0 - 2.0 * u) + 2.0 * u
    radius = 2.0 - v
    x_ratio = a / denominator
    y_ratio = b / denominator
    dx_du = radius * (da * denominator - a * dd) / denominator**2
    dy_du = radius * (db * denominator - b * dd) / denominator**2
    point = (radius * x_ratio, radius * y_ratio, 3.0 * w)
    derivatives = (
        (dx_du, dy_du, 0.0),
        (-x_ratio, -y_ratio, 0.0),
        (0.0, 0.0, 3.0),
    )
    return point, derivatives


class RationalBezierMapTests(unittest.TestCase):
    def test_affine_map_on_vertices_faces_and_interior(self):
        controls = []
        for i, j, k in itertools.product(range(2), repeat=3):
            controls.append((2.0 * i + j, 3.0 * j, 4.0 * k, 1.0))
        mapping = RationalBezierMap((1, 1, 1), controls)

        for u, v, w in (
            (0, 0, 0), (1, 1, 1), (0, 0.4, 0.7), (1, 0.2, 0.9),
            (0.3, 0, 0.6), (0.8, 1, 0.1), (0.2, 0.5, 0),
            (0.6, 0.3, 1), (0.25, 0.5, 0.75),
        ):
            point, columns = mapping.evaluate(u, v, w)
            expected = (2.0 * u + v, 3.0 * v, 4.0 * w)
            self.assertVectorAlmostEqual(point, expected)
            for actual, expected in zip(
                columns, ((2.0, 0.0, 0.0), (1.0, 3.0, 0.0), (0.0, 0.0, 4.0))
            ):
                self.assertVectorAlmostEqual(actual, expected)

    def test_quarter_annulus_against_closed_form_oracle(self):
        reference = json.loads((ROOT / "references" / "quarter_annulus.json").read_text())
        mapping = RationalBezierMap(tuple(reference["degrees"]), reference["controls"])
        max_position_error = 0.0
        max_derivative_error = 0.0
        grid = (0.0, 0.25, 0.5, 0.75, 1.0)

        for u, v, w in itertools.product(grid, repeat=3):
            actual_point, actual_columns = mapping.evaluate(u, v, w)
            expected_point, expected_columns = quarter_annulus_oracle(u, v, w)
            max_position_error = max(
                max_position_error,
                *(abs(a - e) for a, e in zip(actual_point, expected_point)),
            )
            max_derivative_error = max(
                max_derivative_error,
                *(abs(a - e)
                  for actual_column, expected_column in zip(actual_columns, expected_columns)
                  for a, e in zip(actual_column, expected_column)),
            )
            self.assertVectorAlmostEqual(actual_point, expected_point)
            for actual, expected in zip(actual_columns, expected_columns):
                self.assertVectorAlmostEqual(actual, expected)
            radius = 2.0 - v
            self.assertAlmostEqual(actual_point[0] ** 2 + actual_point[1] ** 2,
                                   radius**2, delta=TOL)
            self.assertAlmostEqual(actual_point[2], 3.0 * w, delta=TOL)

        # Expose actual worst errors in unittest output on failure and keep the
        # requested coordinate/Jacobian maxima available to test runners.
        self.max_position_error = max_position_error
        self.max_derivative_error = max_derivative_error
        self.assertLessEqual(max_position_error, TOL)
        self.assertLessEqual(max_derivative_error, TOL)

    def test_invalid_control_dimensions_and_zero_weights(self):
        valid_controls = [(0.0, 0.0, 0.0, 1.0)] * 8
        malformed = list(valid_controls[:-1])
        zero_weight = list(valid_controls)
        zero_weight[3] = (0.0, 0.0, 0.0, 0.0)
        for degrees, controls in (
            ((1, 1, 1), malformed),
            ((1, 1, 1), [*valid_controls[:-1], (1.0, 2.0, 3.0)]),
            ((1, 1, 1), zero_weight),
        ):
            with self.subTest(degrees=degrees, controls=controls):
                with self.assertRaises(ValueError):
                    RationalBezierMap(degrees, controls)

    def test_unsupported_degree_family_is_distinct_from_invalid_input(self):
        controls = [(0.0, 0.0, 0.0, 1.0)] * 16
        with self.assertRaises(NotImplementedError):
            RationalBezierMap((3, 1, 1), controls)

    def test_parameters_must_be_finite_and_inside_the_cube(self):
        mapping = RationalBezierMap((1, 1, 1), [(0, 0, 0, 1)] * 8)
        for point in ((-0.01, 0, 0), (0, 1.01, 0), (0, 0, math.nan)):
            with self.subTest(point=point), self.assertRaises(ValueError):
                mapping.evaluate(*point)

    def assertVectorAlmostEqual(self, actual, expected):
        self.assertEqual(len(actual), len(expected))
        for actual_value, expected_value in zip(actual, expected):
            self.assertAlmostEqual(actual_value, expected_value, delta=TOL)


if __name__ == "__main__":
    unittest.main()
