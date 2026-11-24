"""Exact restriction and affine-field oracle for rational Bézier subdivision."""

from fractions import Fraction
from itertools import product
from math import comb, prod
import unittest

from trivariate.build import build
from trivariate.refinement import subdivide


def evaluate_coefficients(degrees, controls, point):
    """Evaluate a tensor polynomial by direct expansion in Bernstein form."""
    result = [Fraction(0) for _ in controls[0]]
    ranges = (range(degree + 1) for degree in degrees)
    for index, control in zip(product(*ranges), controls, strict=True):
        weight = prod(
            comb(d, i) * x**i * (1 - x)**(d - i)
            for d, i, x in zip(degrees, index, point)
        )
        for c, value in enumerate(control):
            result[c] += weight * Fraction(value)
    return result


class SubdivisionTests(unittest.TestCase):
    def test_exact_geometry_and_affine_field_on_unisolvent_nodes(self):
        model = build({"family": "annular_pair_v1", "heights_mm": [3, 11]})
        block = model["blocks"][1]
        degrees = block["degrees"]
        field = []
        for control in block["controls"]:
            x, y, z, weight = map(Fraction, control)
            field.append(3 + (2*x - 5*y + 7*z) / weight)
        cut = Fraction(2, 7)
        nodes = tuple(tuple(Fraction(i, d) for i in range(d+1)) for d in degrees)
        for axis in range(3):
            children = subdivide(block, axis, cut, field=field)
            for child, (start, end) in zip(children, ((0, cut), (cut, 1)), strict=True):
                self.assertEqual(list(child["degrees"]), degrees)
                homogeneous_field = [
                    [value * ctrl[3]]
                    for value, ctrl in zip(child["field"], child["controls"], strict=True)
                ]
                for point in product(*nodes):
                    parent_point = list(point)
                    parent_point[axis] = start + (end - start)*point[axis]
                    expected = evaluate_coefficients(degrees, block["controls"], parent_point)
                    actual = evaluate_coefficients(degrees, child["controls"], point)
                    self.assertEqual(actual, expected)
                    scalar, = evaluate_coefficients(degrees, homogeneous_field, point)
                    x, y, z, weight = expected
                    self.assertEqual(scalar, 3*weight + 2*x - 5*y + 7*z)

    def test_refuses_malformed_restrictions(self):
        model = build({"family": "annular_pair_v1", "heights_mm": [3, 11]})
        block = model["blocks"][0]
        for axis, cut in ((-1, 0.5), (True, 0.5), (3, 0.5), (0, 0), (0, 1), (0, True)):
            with self.subTest(axis=axis, cut=cut), self.assertRaises(ValueError):
                subdivide(block, axis, cut)
        with self.assertRaises(ValueError):
            subdivide(block, 0, 0.5, field=[1])
        for controls in (block["controls"][:-1], [[0,0,0,0]]*12):
            with self.assertRaises(ValueError):
                subdivide(dict(block, controls=controls), 0, 0.5)
